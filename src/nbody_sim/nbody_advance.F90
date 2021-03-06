! NBODY_ADVANCE.F90
! D. A. Hubber - 13/9/2010
! Advance positions and velocities of all bodies.  Currently uses the 
! 'predictor' step of the Hermite method.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_advance
  use Nbody_module
  use sink_module
  use time_module

  integer(kind=ILP) :: dn         ! Integer timestep since start of timestep
  integer(kind=ILP) :: nfull      ! Full integer timestep
  integer :: s                    ! Star particle counter
  real(kind=PR) :: dt             ! Physical time since beginning of timestep
  real(kind=DP) :: dt2            ! dt*dt
  real(kind=DP) :: dt3            ! dt*dt*dt

  debug2("[nbody_advance.F90]")
  debug_timing("NBODY_ADVANCE")

! Loop over all stars in simulation
! ----------------------------------------------------------------------------
  do s=1,stot
     star(s)%accdo = .false.
     dn            = n - star(s)%nlast
     nfull         = 2**(level_step - star(s)%nlevel)
     dt            = timestep*real(dn,DP)
     dt2           = dt*dt
     dt3           = dt2*dt

     star(s)%r(1:NDIM) = star(s)%rold(1:NDIM) + star(s)%vold(1:NDIM)*dt + &
          & 0.5_DP*star(s)%a0(1:NDIM)*dt2 + ONESIXTH_DP*star(s)%adot0*dt3
     star(s)%v(1:VDIM) = star(s)%vold(1:VDIM) + star(s)%a0(1:VDIM)*dt + &
          & 0.5_DP*star(s)%adot0(1:VDIM)*dt2

#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(star(s)%r(1:NDIM),star(s)%v(1:VDIM))
#endif

     if (dn == nfull) then
        star(s)%accdo = .true.
        star(s)%nlast = n
        star(s)%rold(1:NDIM) = star(s)%r(1:NDIM)
        star(s)%vold(1:VDIM) = star(s)%v(1:VDIM)
     end if

  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE nbody_advance
