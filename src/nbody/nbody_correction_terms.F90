! NBODY_CORRECTION_TERMS.F90
! D. A. Hubber - 13/9/2010
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_correction_terms
  use definitions
  use time_module
  use Nbody_module
  use sink_module, only : stot
  implicit none

  integer :: s                      ! Star counter
  real(kind=DP) :: dt               ! Physical timestep
  real(kind=DP) :: dt2              ! dt*dt
  real(kind=DP) :: dt3              ! dt*dt*dt

  debug2("Calculating correction terms [nbody_correction_terms.F90]")
  debug_timing("NBODY_CORRECTION")

! 4th-order Hermite correction terms
! ----------------------------------------------------------------------------
#if defined(NBODY_HERMITE4)
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(dt,dt2,dt3) IF (stot > 500)
  do s=1,stot
     if (.not. star(s)%accdo) cycle

     dt  = timestep*real(2**(level_step - star(s)%nlevel),DP)
     dt2 = dt*dt
     dt3 = dt2*dt

     star(s)%a2dot0(1:NDIM) = (-6.0_DP*(star(s)%a0(1:NDIM)-star(s)%a(1:NDIM))&
          & - dt*(4.0_DP*star(s)%adot0(1:NDIM) + 2.0_DP*star(s)%adot(1:NDIM)))&
          & / dt2
     star(s)%a3dot(1:NDIM) = (12.0_DP*(star(s)%a0(1:NDIM) - star(s)%a(1:NDIM))&
          & + 6.0_DP*dt*(star(s)%adot0 + star(s)%adot)) / dt3

     star(s)%r(1:NDIM) = star(s)%r(1:NDIM) &
          & + star(s)%a2dot0(1:NDIM)*dt2*dt2/24.0_DP &
          & + star(s)%a3dot(1:NDIM)*dt3*dt2/120.0_DP
     star(s)%v(1:NDIM) = star(s)%v(1:NDIM) &
          & + star(s)%a2dot0(1:NDIM)*dt3/6.0_DP &
          & + star(s)%a3dot(1:NDIM)*dt2*dt2/24.0_DP
     star(s)%a2dot(1:NDIM) = star(s)%a2dot0(1:NDIM) + star(s)%a3dot(1:NDIM)*dt

     star(s)%rold(1:NDIM)  = star(s)%r(1:NDIM)
     star(s)%vold(1:NDIM)  = star(s)%v(1:NDIM)
     star(s)%a0(1:NDIM)    = star(s)%a(1:NDIM)
     star(s)%adot0(1:NDIM) = star(s)%adot(1:NDIM)
  end do
!$OMP END PARALLEL DO


! 6th-order Hermite correction terms
! ----------------------------------------------------------------------------
#elif defined(NBODY_HERMITE6)


#endif
! ----------------------------------------------------------------------------


  return
END SUBROUTINE nbody_correction_terms
