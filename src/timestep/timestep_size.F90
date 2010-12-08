! TIMESTEP_SIZE.F90
! C. P. Batty & D. A. Hubber - 29/3/2007
! Calculates ideal timestep for particle p based on multiple possible 
! criterion.  Calculates minimum of:
! 1. Acceleration timestep, dt = sqrt(h / accel), 
! 2. Courant condition,     dt = h / (h*div_v + sound)
!    or (with viscosity)    dt = h / (sound + h*div_v_p + 
!                                     TVISC_FAC*(alpha*sound + beta*h*div_v_p)
! 3. Energy condition,      dt = u / (dudt + SMALL_NUMBER)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE timestep_size(p,dt)
  use definitions
  use type_module
  use particle_module, only : parray,a
  use hydro_module
  use time_module, only : accel_mult,courant_mult
#if defined(HYDRO)
#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  use particle_module, only : u,du_dt
#endif
#if defined(IDEAL_MHD)
  use mhd_module, only : B
#endif
#endif
#if defined(DEBUG_BLOCK_TIMESTEP)
  use time_module, only : nacc,ncour,meantacc,meantcour
#endif
  implicit none

  integer, intent(in)        :: p    ! particle counter
  real(kind=DP), intent(out) :: dt   ! Step size for particle p

  real(kind=DP) :: amag              ! magnitude of acceleration
  real(kind=DP) :: ap(1:NDIM)        ! Local copy of acceleration
  real(kind=DP) :: hp                ! Local copy of smoothing length
  real(kind=DP) :: tacc              ! Acceleration timestep
!#if defined(HYDRO)
  real(kind=DP) :: div_v_p           ! Local copy of velocity divergence
  real(kind=DP) :: tcour             ! Courant time
  real(kind=DP) :: vsignal           ! Signal speed of particle p
#if defined(ARTIFICIAL_VISCOSITY)
  real(kind=DP) :: alpha_p           ! Local copy of alpha visc. value 
  real(kind=DP) :: beta_p            ! Local copy of beta visc. value
#endif
#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  real(kind=DP) :: tenergy           ! Energy time
  real(kind=DP) :: dudt_p            ! Local copy of du_dt for particle p
  real(kind=DP) :: up                ! Local copy of specific internal energy 
#endif
!#endif 

  debug3("Calculating timestep [timestep_size.F90] for particle ", p)

! Make local copies of important particle properties
  ap(1:NDIM) = real(a(1:NDIM,p),DP)
  hp = real(parray(SMOO,p),DP)
  div_v_p = real(abs(div_v(p)),DP)
#if defined(HYDRO)
#if defined(IDEAL_MHD)
  vsignal = sqrt(real(sound(p)**2 + &
       & dot_product(B(1:BDIM,p),B(1:BDIM,p))/rho(p),DP))
#elif defined(SIGNAL_VELOCITY)
  vsignal = real(vsigmax(p),DP)
#else
  vsignal = real(sound(p),DP)
#endif
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD) && defined(VISC_BALSARA)
  alpha_p = real(talpha(p),DP)*real(balsara(p),DP)
  beta_p  = 2.0_DP*alpha_p
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD)
  alpha_p = real(talpha(p),DP)
  beta_p  = 2.0_DP*alpha_p
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_BALSARA)
  alpha_p = real(alpha*balsara(p),DP)
  beta_p  = real(beta*balsara(p),DP)
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_PATTERN_REC)
  alpha_p = real(alpha*pattrec(p),DP)
  beta_p = real(beta*pattrec(p),DP)
#elif defined(ARTIFICIAL_VISCOSITY) 
  alpha_p = real(alpha,DP)
  beta_p  = real(beta,DP)
#endif
#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  up = real(u(p),DP)
  dudt_p = real(abs(du_dt(p)),DP)
#endif
#endif
  amag = sqrt(dot_product(ap(1:NDIM),ap(1:NDIM)))
  if (p > phydroend) vsignal = 0.0_DP

! Acceleration condition on timestep (Always calculated)
  tacc = accel_mult * sqrt(hp / (amag + SMALL_NUMBER_DP))
  dt = tacc

! Courant condition (with or without artificial viscosity), or without 
! gravity (no sound speed, but velocity divergence).
#if defined(HYDRO) && defined(SIGNAL_VELOCITY)
  tcour = courant_mult * hp / vsignal
  dt = min(dt,tcour)
#elif defined(HYDRO) && defined(ARTIFICIAL_VISCOSITY)
  tcour = courant_mult * hp / (vsignal + hp*div_v_p + &
       & TVISC_FAC*(alpha_p*vsignal + beta_p*hp*div_v_p))
  dt = min(dt,tcour)
#elif defined(HYDRO) && !defined(ARTIFICIAL_VISCOSITY)
  tcour = courant_mult * hp / (vsignal + hp*div_v_p)
  dt = min(dt,tcour)
#elif !defined(HYDRO)
  tcour = courant_mult / (div_v_p + SMALL_NUMBER)
  dt = min(dt,tcour)
#endif

! Internal energy timestep
#if defined(HYDRO) && defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  if (p < phydroend) then
     tenergy  = min(accel_mult,courant_mult) * up / (dudt_p + SMALL_NUMBER_DP)
     dt = min(dt,tenergy)
  end if
#endif

#if defined(DEBUG_BLOCK_TIMESTEPS)
  meantcour = meantcour + tcour
  meantacc = meantacc + tacc
  if (tcour > tacc) then
     ncour = ncour + 1
  else
     nacc = nacc + 1
  end if
#endif

#if defined(DEBUG_TIMESTEP_SIZE)
#if defined(HYDRO)
  write(6,*) "Courant timestep,      tcour    : ", tcour
#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  write(6,*) "Energy timestep,       tenergy  : ", tenergy
#endif
#endif
  write(6,*) "Acceleration timestep, tacc     : ", tacc
  write(6,*) "Minimum timestep,      dt       : ", dt
#endif

  return
END SUBROUTINE timestep_size
