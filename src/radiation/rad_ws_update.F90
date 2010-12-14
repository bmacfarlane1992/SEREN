! RAD_WS_UPDATE.F90
! D. A. Hubber - 14/7/2008
! Update radiation transport quantities.
! (Ref : Stamatellos et al. 2007; Forgan et al. 2009)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE rad_ws_update
  use interface_module
  use particle_module
  use hydro_module
  use Eos_module
  use type_module
  use time_module
  implicit none

  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer(kind=ILP) :: dn             ! Int. timestep since start of timestep
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
  real(kind=PR) :: dt                 ! Timestep
  real(kind=PR) :: dt_new             ! Latest timestep
  real(kind=PR) :: dt_old             ! Previous timestep
  real(kind=PR) :: mu_bar_p           ! Mean gas particle mass for p
#if defined(OPENMP)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif

  debug2("Updating radiation transport quantities [rad_ws_update.F90]")
  debug_timing("RAD_WS")

! For block timesteps, first make a list of all hydro SPH particles on 
! an acceleration step, and then parallelize over that list.
  acctot = 0
  allocate(acclist(1:ptot))
  do p=pgasstart,pgasend
     if (accdo(p)) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if
  end do
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)) + 1
#endif


! Calculate flux-limited diffusion terms
! ----------------------------------------------------------------------------
#if defined(DIFFUSION)
 if (acctot > 0) then
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call find_idens(rho(p),idens(p))
        call find_itemp(temp(p),itemp(p))
        call conductivity(p)
     end do
     !$OMP END PARALLEL DO
 end if
#endif
! ----------------------------------------------------------------------------


! Calculate itemp, idens and column density to infinity of all particles
! ----------------------------------------------------------------------------
  if (acctot > 0) then
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
     !$OMP PRIVATE(p,dt,dt_old,dt_new)
     do i=1,acctot
        p = acclist(i)
        dt_old = real(laststep(p),PR)
        dt_new = real(2**(level_step - nlevel(p)),PR)*real(timestep,PR)
#if defined(EULER) || defined(LEAPFROG_KDK) 
        dt = dt_new
#elif defined(RUNGE_KUTTA)
        dt = 0.5_PR*dt_new
#elif defined(LEAPFROG_DKD) || defined(PREDICTOR_CORRECTOR)
        dt = 0.5_PR*(dt_old + dt_new)
#endif
#if defined(DIFFUSION)
        call diffusion(p,dt)
#endif
        call solve_implicit_cooling_ws(p)
     end do
     !$OMP END PARALLEL DO
  end if
! ----------------------------------------------------------------------------


! Perform implicit integration of internal energy
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
!$OMP PRIVATE(dn,dt,mu_bar_p) 
  do p=pgasstart,pgasend

     ! Calculate time since last force computation.
     dn = n - nlast(p)
     dt = real(timestep,PR)*real(dn,PR)
     
#if defined(IONIZING_UV_RADIATION)
     if (temp(p) < temp_min(p)) then
        temp(p) = temp_min(p)
        call find_idens(rho(p),idens(p))
        call find_itemp(temp(p),itemp(p))
        u(p) = eosenergy(rho(p),temp(p),idens(p),itemp(p))
        ueq(p) = u(p)
     end if
#endif
     
     ! Perform implicit integration depending on timestep.
     if (dt_therm(p) <= SMALL_NUMBER) then
        u(p) = u_old(p)
     else if (dt < 40.0_PR*dt_therm(p)) then
        u(p) = u_old(p)*exp(-dt/dt_therm(p)) &
             & + ueq(p)*(1.0_PR - exp(-dt/dt_therm(p)))
     else if (dt >= 40.0_PR*dt_therm(p)) then
        u(p) = ueq(p)
     end if
     
     ! Now update all other thermal properties
     call find_idens(rho(p),idens(p))
     call find_temp_from_energy(idens(p),u(p),itemp(p),temp(p))
     mu_bar_p = eosmu(rho(p),temp(p),idens(p),itemp(p))
     press(p) = Pconst2*temp(p)*rho(p)/mu_bar_p
     sound(p) = sqrt(press(p)/rho(p))
  end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------


  deallocate(acclist)

  return
END SUBROUTINE rad_ws_update
