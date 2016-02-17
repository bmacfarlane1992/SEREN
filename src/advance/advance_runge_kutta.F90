! ADVANCE_RUNGE_KUTTA.F90
! C. P. Batty & D. A. Hubber - 19/3/2007
! Second order Runge-Kutta integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_runge_kutta(p)
  use particle_module
  use time_module
  use hydro_module
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  use mhd_module, only : B, B_old, dB_dt
#endif
  implicit none

  integer, intent(in) :: p     ! Particle id

  integer :: dn                ! Integer timestep since beginning of timestep
  integer :: nfull             ! Full integer timestep
  integer :: nhalf             ! Half the full integer timestep
  real(kind=DP) :: dt          ! Physical time since beginning of timestep
  real(kind=PR) :: rp(1:NDIM)  ! Local copy of position 
  real(kind=PR) :: vp(1:VDIM)  ! Local copy of velocity

  debug3("Advancing particle ",p)

! Work out integer and real time intervals from beginning of step
  accdo(p) = .false.
  dn       = n - nlast(p)
  nfull    = 2**(level_step - nlevel(p))
  nhalf    = nfull / 2
  dt       = timestep*real(dn,DP)


! Advance particles that are below or at the half timestep
! ----------------------------------------------------------------------------
  if (dn <= nhalf) then

     ! First half of the Runge-Kutta integration step
     rp(1:NDIM) = r_old(1:NDIM,p) + v_old(1:NDIM,p)*dt
     vp(1:VDIM) = v_old(1:VDIM,p) + a(1:VDIM,p)*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

     ! Record velocities if at half timestep
     if (dn == nhalf) then
        v_half(1:VDIM,p) = vp(1:VDIM)
        accdo(p) = .true.
     endif

! Else advance particles that are beyond the half timestep
! ----------------------------------------------------------------------------
  else

     ! Second half of the Runge-Kutta integration step
     rp(1:NDIM) = r_old(1:NDIM,p) + v_half(1:NDIM,p)*dt
     vp(1:VDIM) = v_old(1:VDIM,p) + a(1:VDIM,p)*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

  end if
! ----------------------------------------------------------------------------


! Record positions, velocities and other quantities in arrays for all cases 
  parray(1:NDIM,p) = rp(1:NDIM)
  v(1:VDIM,p) = vp(1:VDIM)
#if defined(ENTROPIC_FUNCTION) && defined(INTERNAL_ENERGY)
  Aent(p) = Aold(p) + dA_dt(p)*dt
  u(p) = (Aent(p)*rho(p)**(gamma - 1.0_PR))/(gamma - 1.0_PR)
#elif defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  u(p) = u_old(p) + du_dt(p)*dt
#endif
#if defined(VISC_TD)
  talpha(p) = talpha_old(p) + dalpha_dt(p)*dt
  if (talpha(p) < alpha_min) talpha(p) = alpha_min
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  B(1:3,p) = B_old(1:3,p) + dB_dt(1:3,p)*dt
#endif

! If end of timestep, record as 'old' values
  if (dn == nfull) then
     accdo(p)        = .true.
     nlast(p)        = n
     laststep(p)     = timestep*real(nfull,DP)
     r_old(1:NDIM,p) = rp(1:NDIM)
     v_old(1:VDIM,p) = vp(1:VDIM)
#if defined(ENTROPIC_FUNCTION) && defined(INTERNAL_ENERGY)
     Aold(p)         = Aent(p)
#elif defined(INTERNAL_ENERGY)
     u_old(p)        = u(p)
#endif
#if defined(VISC_TD)
     talpha_old(p)   = talpha(p)
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
     B_old(1:BDIM,p) = B(p,1:BDIM)
#endif
  end if
 
  return
END SUBROUTINE advance_runge_kutta
