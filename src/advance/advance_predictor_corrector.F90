! ADVANCE_PREDICTOR_CORRECTOR.F90
! C. P. Batty - 19/3/2007
! Second order predictor-corrector integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_predictor_corrector(p)
  use interface_module, only : check_boundary_conditions
  use particle_module
  use time_module
  use hydro_module
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  use mhd_module, only : B, B_old, dB_dt
#endif
  implicit none

  integer, intent(in) :: p      ! Particle id

  integer :: dn                 ! Integer timestep since beginning of timestep
  integer :: nfull              ! Full integer timestep
  integer :: nhalf              ! Half the full integer timestep
  real(kind=DP) :: dt           ! Physical time since beginning of timestep
  real(kind=PR) :: rp(1:NDIM)   ! Local copy of position 
  real(kind=PR) :: vp(1:VDIM)   ! Local copy of velocity
#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  real(kind=PR) :: up           ! Local copy of specific internal energy 
#endif
#if defined(VISC_TD)
  real(kind=PR) :: alpha_p      ! Local copy of talpha
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQN)
  real(kind=PR) :: Bp(1:3)      ! Local copy of B-field
#endif

  debug3("Advancing particle ",p)

! Work out integer and real time intervals from beginning of step
  accdo(p) = .false.
  dn       = n - nlast(p)
  nfull    = 2**(level_step - nlevel(p))
  nhalf    = nfull / 2
  dt       = real(timestep,PR)*real(dn,PR)

! Advance particles that are below the full timestep
! ----------------------------------------------------------------------------
  if (dn < nfull) then

     ! First half of the Predictor-Corrector integration step
     rp(1:NDIM) = r_old(1:NDIM,p) + &
          &v_old(1:NDIM,p)*dt + 0.5_PR*a_old(1:NDIM,p)*dt*dt
     vp(1:VDIM) = v_old(1:VDIM,p) + a(1:VDIM,p)*dt
     
     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
     up = u_old(p) + du_dt(p)*dt
#endif
#if defined(VISC_TD)
     alpha_p = talpha_old(p) + dalpha_dt(p)*dt
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
     Bp(1:3) = B_old(1:3,p) + dB_dt(1:3,p)*dt
#endif

     ! Only calculate forces at the half timestep
     if (dn == nhalf) accdo(p) = .true.

! Else advance particles that are at the full timestep
! ----------------------------------------------------------------------------
  else

     ! Second half of the Predictor-Corrector integration step
     rp(1:NDIM) = r_old(1:NDIM,p) + &
          &v_old(1:NDIM,p)*dt + 0.5_PR*a_old(1:NDIM,p)*dt*dt
     vp(1:VDIM) = v_old(1:VDIM,p) + 0.5_PR*(a_old(1:VDIM,p) + a(1:VDIM,p))*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
     up = u_old(p) + du_dt(p)*dt
#endif
#if defined(VISC_TD)
     alpha_p = talpha_old(p) + dalpha_dt(p)*dt
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
     Bp(1:3) = B_old(1:3,p) + dB_dt(1:3,p)*dt
#endif

     ! Record as 'old' values
     nlast(p)        = n
     laststep(p)     = timestep*real(nfull,DP)
     r_old(1:NDIM,p) = rp(1:NDIM)
     v_old(1:VDIM,p) = vp(1:VDIM)
     a_old(1:VDIM,p) = a(1:VDIM,p)
#if defined(INTERNAL_ENERGY)
     u_old(p)        = up
#endif
#if defined(ENTROPIC_FUNCTION)
     Aold(p)         = Aent(p)
#endif
#if defined(VISC_TD)
     talpha_old(p)   = alpha_p
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
     B_old(1:3,p)    = Bp(1:3)
#endif
  end if
! ----------------------------------------------------------------------------

! Record positions, velocities and other quantities in arrays for all cases 
  parray(1:NDIM,p) = rp(1:NDIM)
  v(1:VDIM,p) = vp(1:VDIM)
#if defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  u(p) = up
#endif
#if defined(VISC_TD)
  talpha(p) = alpha_p
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  B(1:3,p) = Bp(1:3)
#endif

  return
END SUBROUTINE advance_predictor_corrector
