! THERMAL.F90
! C. P. Batty & D. A. Hubber - 19/1/2007
! Calculates temperature, pressure, sound speed of all particles depending
! on the chosen equation of state.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE thermal(p)
  use interface_module, only : distance3,eosenergy,eosmu,&
       &find_idens,find_itemp,find_temp_from_energy
  use particle_module
  use hydro_module
  use scaling_module
  use time_module
#if defined(STELLAR_HEAT)
  use sink_module
#endif
#if defined(RAD_WS)
  use type_module
  use Eos_module
#endif
#if defined(IONIZING_UV_RADIATION)
  use HP_module
#endif
  implicit none

  integer, intent(in) :: p     ! particle id

  real(kind=PR) :: press_p     ! pressure of particle p
  real(kind=PR) :: rho_p       ! density of particle p
  real(kind=PR) :: sound_p     ! sound speed of particle p
  real(kind=PR) :: temp_p      ! temperature of particle p
#if defined(INTERNAL_ENERGY) && defined(RAD_WS)
  real(kind=PR) :: mu_bar_p    ! Mean gas particle mass for particle p
#endif
#if defined(IONIZING_UV_RADIATION) && !defined(RAD_WS)
  real(kind=PR) :: fmu         ! 1/sqrt(mu) when mu = 1 in params.dat
#endif
#if defined(STELLAR_HEAT)
  integer       :: j           ! Aux. loop counter
  real(kind=PR) :: dr(1:NDIM)  ! Relative displacement
  real(kind=PR) :: drsqd       ! Distance squared 
#endif

  rho_p = rho(p)


! Isothermal equation of state
! ----------------------------------------------------------------------------
#if defined(ISOTHERMAL) 
#if defined(DIMENSIONLESS)
  temp_p  = 1.0_PR
#else
  temp_p  = isotemp
#endif
  press_p = Pconst*temp_p*rho_p
  sound_p = sqrt(press_p/rho_p)
#endif

! 'Local' Isothermal equation of state (T~To*R^{-q})
! ----------------------------------------------------------------------------
#if defined(LOCAL_ISOTHERMAL) 
  call ambient_temp(p,temp_p) 
  press_p = Pconst*temp_p*rho_p
  sound_p = sqrt(press_p/rho_p)
#endif

! ----------------------------------------------------------------------------
  
  
! Barotropic equation of state
! ----------------------------------------------------------------------------
#if defined(BAROTROPIC)
  temp_p  = isotemp*(1.0_PR + (rho_p/rhobary)**(gamma - 1.0_PR))
  press_p = Pconst*temp_p*rho_p
  sound_p = sqrt(press_p/rho_p)
#endif
! ----------------------------------------------------------------------------
  
  
! Polytropic equation of state
! ----------------------------------------------------------------------------
#if defined(POLYTROPIC)
  press_p = Kpoly*(rho_p**(gamma))
  temp_p  = Kpoly*(rho_p**(gamma - 1.0_PR))/Pconst
  sound_p = sqrt(gamma*press_p/rho_p)
#endif
! ----------------------------------------------------------------------------
  
  
! Stiff equation of state (Monaghan shear flow)
! ----------------------------------------------------------------------------
#if defined(STIFF)
  if (rho_p > 1.0_PR) then
     press_p = (50.0_PR/7.0_PR)*((rho_p**7) - 1.0_PR)
  else
     press_p = (50.0_PR/7.0_PR) ! *** modified for rho < 1
  end if
  sound_p = sqrt(press_p/rho_p)
  temp_p  = press_p / (Pconst*rho_p)
#endif
! ----------------------------------------------------------------------------
  

! Ideal-gas equation of state when using the internal energy equation
! ----------------------------------------------------------------------------
#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  press_p = (gamma - 1.0_PR)*rho_p*u(p)
  temp_p  = press_p / (Pconst*rho_p)
  sound_p = sqrt(press_p/rho_p)
#endif
! ----------------------------------------------------------------------------
  

! Polytropic-cooling approximation
! ----------------------------------------------------------------------------
#if defined(INTERNAL_ENERGY) && defined(RAD_WS)
#if defined(IONIZING_UV_RADIATION)
  if (temp_p < temp_min(p)) then
     temp(p) = temp_min(p)
     call find_idens(rho(p),idens(p))
     call find_itemp(temp_min(p),itemp(p))
     u(p) = eosenergy(rho(p),temp_min(p),idens(p),itemp(p))
  else 
     call find_idens(rho(p),idens(p))
     call find_temp_from_energy(idens(p),u(p),itemp(p),temp(p))
  end if
#else
  call find_idens(rho(p),idens(p))
  call find_temp_from_energy(idens(p),u(p),itemp(p),temp(p))
#endif
  temp_p   = temp(p)
  mu_bar_p = eosmu(rho_p,temp_p,idens(p),itemp(p))
  press_p  = Pconst2*temp_p*rho_p/mu_bar_p
  sound_p  = sqrt(press_p/rho_p)
#endif
! ----------------------------------------------------------------------------
  
  
! Heating from sinks (stars)
! ----------------------------------------------------------------------------
#if defined(STELLAR_HEAT)
  temp_p = 0.0_PR
  do j=1,STARS
     call distance3(sink(j)%r(1:NDIM),parray(1:NDIM,p),dr(1:NDIM),drsqd)
     temp_p = temp_p + (isotemp**4) / drsqd
  end do
  temp_p  = temp_p + 10000.0_PR ! (10K)**4
  temp_p  = temp_p**0.25_PR
  press_p = Pconst*rho_p*temp_p
  sound_p = sqrt(press_p/rho_p)
#endif
! ----------------------------------------------------------------------------
  
  
! Other thermal quantities due to ionizing radiation
! ----------------------------------------------------------------------------
#if defined(IONIZING_UV_RADIATION) && !defined(RAD_WS)
  temp_p = max(temp_p,temp_min(p))
  fmu = ((temp_p - Tneut)/(Tion - Tneut))*0.414213 + 1.0_PR  !0.562_PR + 0.652_PR 
  if (temp_p .ge. Tion) fmu = 1.414213_PR  !1.213_PR
  if (temp_p .le. Tneut) fmu = 1.0_PR      !0.652_PR
  press_p = fmu*fmu*Pconst*temp_p*rho_p
  sound_p = sqrt(press_p/rho_p)
#endif
! ----------------------------------------------------------------------------
  
  
! Record thermodynamical information in main data arrays
  press(p) = press_p
  sound(p) = sound_p
  temp(p)  = temp_p 
#if defined(VISC_BALSARA)
  balsara(p) = min(1.0_PR,(abs(div_v(p)) / &
       &(balsara(p) + abs(div_v(p)) + BAL_DENOM*sound_p/parray(SMOO,p))))
#endif
  
  
  return
END SUBROUTINE thermal
