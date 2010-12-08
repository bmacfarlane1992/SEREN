! DIFFUSION.F90
! Dimitris Stamatellos - 02/09/2009
! Calculates thermal conductivity of p
! N.B. CHECK WITH DIMITRIS
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE diffusion(p,dt)
  use interface_module, only : distance2
  use particle_module
  use neighbour_module
  use kernel_module
  use hydro_module
  use constant_module
  use eos_module
  use scaling_module
  implicit none

  integer, intent(in) :: p               ! particle identifier
  real(kind=PR), intent(in) ::dt         ! Timestep

  integer :: i                           ! neighbour counter
  integer :: kern_p                      ! kernel table element for p
  integer :: kern_pp                     ! kernel table element for p
  integer :: pp                          ! neighbouring particle number
  integer :: pp_numb                     ! number of neighbours
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: dr(1:NDIM)            ! Relative displacement vector
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: dr_unit(1:NDIM)       ! Unit displacement vector
  real(kind=PR) :: dt_diff_pp            ! energy diff. timescale from p to pp
  real(kind=PR) :: du_diff               ! ..
  real(kind=PR) :: dudt_diff             ! energy diff. rate from p
  real(kind=PR) :: dudt_diff_pp          ! energy diff. rate from pp
  real(kind=PR) :: hfactor_p             ! invhp ^ NDIM
  real(kind=PR) :: hfactor_pp            ! invhpp ^ NDIM
  real(kind=PR) :: hp                    ! smoothing length of particle p
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: invhpp                ! ( 1 / hpp )
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: rho_p                 ! density of particle p
  real(kind=PR) :: rho_pp                ! density of particle pp
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: skern                 ! 0.5 * (r/h) * KERNTOT
  real(kind=PR) :: vp(1:NDIM)            ! velocity of particle p
  real(kind=PR) :: wmean                 ! (W(p) + W(pp)) / 2
  real(kind=PR) :: kappaT                ! ..
  real(kind=PR) :: kappaT_pp             ! ..
  real(kind=PR) :: tau_p                 ! ..
  real(kind=PR) :: tau_pp                ! ..


! Create local copies of important properties of particle p
  rp(1:NDIM) = parray(1:NDIM,p)
  hp         = parray(SMOO,p)
  invhp      = 1.0_PR / hp
  hfactor_p  = invhp**(NDIMPR)
  vp(1:NDIM) = v(1:NDIM,p)
  rho_p      = rho(p)
#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_numb))
     do i=1,pp_numb
        pp_templist(i) = pplist(i,p)
     end do
  else
     allocate(pp_templist(1:ptot))
     call get_neib_on_fly(p,hp,pp_numb,ptot,pp_templist)
  end if
#else
  allocate(pp_templist(1:ptot))
  call get_neib_on_fly(p,hp,pp_numb,ptot,pp_templist)
#endif


! Initialise variables
  dudt_diff    = 0.0_PR
  dudt_diff_pp = 0.0_PR
  du_diff      = 0.0_PR
  call getkappa(rho_p,temp(p),idens(p),kappaT)
  tau_p        = sqrt(column2(p))*kappaT


! Loop over all neighbours and sum terms
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     
     ! Create local copies for neighbour pp     
     rho_pp   = rho(pp)
     invhpp   = 1.0_PR / parray(SMOO,pp)
     hfactor_pp = invhpp**(NDIMPR)
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drmag)
     drmag    = sqrt(drmag) + SMALL_NUMBER
     invdrmag = 1.0_PR / drmag
     dr_unit(1:NDIM) = dr(1:NDIM)*invdrmag
     skern    = HALFKERNTOT * drmag * invhp
     kern_p   = int(skern)
     kern_p   = min(kern_p,KERNTOT)
     skern    = HALFKERNTOT * drmag * invhpp 
     kern_pp  = int(skern)
     kern_pp  = min(kern_pp,KERNTOT)
     call getkappa(rho(pp),temp(pp),idens(pp),kappaT_pp)
     tau_pp   = sqrt(column2(pp))*kappaT_pp

     if (0.5_PR*(tau_p + tau_pp) < 1.0_PR) then
        du_diff = 0.0_PR
     else
        wmean = 0.5_PR*(hfactor_p*invhp*w1(kern_p) + &
             &hfactor_pp*invhpp*w1(kern_pp))

        ! Compute actual energy diffused between particles 
        ! (done individualy in order to take into account
        ! the different diffusion rates between particles)
        dudt_diff_pp = 4.0_PR*(parray(MASS,pp)/(rho_pp*rho_p)) * &
             &(k_cond(p)*k_cond(pp)/(k_cond(p) + k_cond(pp))) * &
             &(temp(p) - temp(pp))*wmean*invdrmag
        du_diff = du_diff + dudt_diff_pp*dt
     end if
     
  end do
! ----------------------------------------------------------------------------
  
!  du_dt_diff(p) = dudt_diff
  du_dt_diff(p) = du_diff / dt

  deallocate(pp_templist)

  return
END SUBROUTINE diffusion
