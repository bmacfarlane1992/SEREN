! ADVANCE_LEAPFROG_KDK.F90
! C. P. Batty & D. A. Hubber - 19/3/2007
! Second order leapfrog kick-drift-kick integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_leapfrog_kdk(p)
  use interface_module, only : check_boundary_conditions
  use particle_module
  use time_module
  use hydro_module
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  use mhd_module, only : B, B_old, dB_dt
#endif
  implicit none

  integer, intent(in) :: p      ! Particle id

  integer(kind=ILP) :: dn       ! Integer timestep since beginning of timestep
  integer(kind=ILP) :: nfull    ! Full integer timestep
  integer(kind=ILP) :: nhalf    ! Half the full integer timestep
  real(kind=PR) :: rp(1:NDIM)   ! Local copy of position 
  real(kind=PR) :: vp(1:VDIM)   ! Local copy of velocity
  real(kind=PR) :: dt           ! Physical time since beginning of timestep
  real(kind=PR) :: dt_half      ! Time since half-step

  debug3("Advancing particle ",p)

! Calculate correct value of 'v_old' here according to formal implementation 
! of leapfrog scheme
  if (accdo(p)) then
     v_old(1:VDIM,p) = v_half(1:VDIM,p) + &
          &0.5_PR*a(1:VDIM,p)*real(laststep(p),PR)
     rho_old(p) = rho(p)
  end if

! Work out integer and real time intervals from beginning of step 
  accdo(p) = .false.
  dn       = n - nlast(p)
  nfull    = 2**(level_step - nlevel(p))
  nhalf    = nfull / 2
  dt       = real(timestep,PR)*real(dn,PR)
  dt_half  = real(timestep,PR)*real(dn - nhalf,PR)


! Advance particles that are below or at the half timestep
! ----------------------------------------------------------------------------
  if (dn <= nhalf) then 

     rp(1:NDIM) = r_old(1:NDIM,p) + v_old(1:NDIM,p)*dt
     vp(1:VDIM) = v_old(1:VDIM,p) + a(1:VDIM,p)*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif
 
     ! Record halfstep velocity here
     if (dn == nhalf) v_half(1:VDIM,p) = vp(1:VDIM)
  
! Else advance particles that are later than the halfstep
! ----------------------------------------------------------------------------
  else
 
     rp(1:NDIM) = r_old(1:NDIM,p) + v_half(1:NDIM,p)*dt
     vp(1:VDIM) = v_half(1:VDIM,p) + a(1:VDIM,p)*dt_half

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif
     
  end if
! ----------------------------------------------------------------------------

! Record positions, velocities and other quantities in arrays for all cases 
! (Ignore though for the Spiegel test)
#ifndef SPIEGEL_TEST
  parray(1:NDIM,p) = rp(1:NDIM)
  v(1:VDIM,p) = vp(1:VDIM)
#endif
#if defined(ENTROPIC_FUNCTION) && defined(INTERNAL_ENERGY)
  Aent(p) = Aold(p) + dA_dt(p)*dt
  u(p) = (Aent(p)*rho(p)**(gamma - 1.0_PR))/(gamma - 1.0_PR)
#elif defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  u(p) = u_old(p) + du_dt(p)*dt
  if (u(p) < SMALL_NUMBER) u(p) = u_old(p)*exp(-u_old(p)/du_dt(p))
#endif
#if defined(VISC_TD)
  talpha(p) = talpha_old(p) + dalpha_dt(p)*dt
  if (talpha(p) < alpha_min) talpha(p) = alpha_min
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  B(1:BDIM,p) = B_old(1:BDIM,p) + dB_dt(1:BDIM,p)*dt
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
     B_old(1:BDIM,p) = B(1:BDIM,p)
#endif
  end if

  return
END SUBROUTINE advance_leapfrog_kdk
