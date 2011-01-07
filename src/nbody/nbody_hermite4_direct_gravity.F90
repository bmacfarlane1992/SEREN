! NBODY_HERMITE4_DIRECT_GRAVITY.F90
! D. A. Hubber - 23/6/2008
! Calculate the gravitational acceleration of star s due to all other 
! stars in the simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_hermite4_direct_gravity(s,invhs,rs,vs,agravs,adots,potp)
  use interface_module, only : gravity_hermite4
  use definitions
  use Nbody_module
  use sink_module, only : stot
  implicit none

  integer, intent(in) :: s                     ! Star i.d.
  real(kind=DP), intent(in) :: invhs           ! Smoothing length of star s
  real(kind=DP), intent(in) :: rs(1:NDIM)      ! Position of star s
  real(kind=DP), intent(in) :: vs(1:NDIM)      ! Velocity of star s
  real(kind=DP), intent(out) :: agravs(1:NDIM) ! Grav. accel of star s
  real(kind=DP), intent(out) :: adots(1:NDIM)  ! Jerk of star s
  real(kind=DP), intent(out) :: potp           ! Potential of star s
  
  integer :: ss                       ! Secondary star counter
  real(kind=DP) :: adottemp(1:NDIM)   ! Aux. variable for calculating jerk
  real(kind=DP) :: atemp(1:NDIM)      ! Auxilary accel variable 
  real(kind=DP) :: dpotp              ! Aux. star pot variable

  debug3("[nbody_gravity.F90] : ",s)

! Zero arrays
  adots(1:NDIM)    = 0.0_DP
  agravs(1:NDIM)   = 0.0_DP
  potp             = 0.0_DP

! Loop over all other stars and calculate net gravitational acceleration
! ----------------------------------------------------------------------------
  do ss=1,stot
     if (s == ss) cycle
     call gravity_hermite4(invhs,star(ss)%h,star(ss)%m,&
          &rs(1:NDIM),star(ss)%r(1:NDIM),vs(1:NDIM),&
          &star(ss)%v(1:NDIM),atemp(1:NDIM),adottemp(1:NDIM),dpotp)
     agravs(1:NDIM) = agravs(1:NDIM) + atemp(1:NDIM)
     adots(1:NDIM) = adots(1:NDIM) + adottemp(1:NDIM)
     potp = potp + dpotp
     
#if defined(DEBUG_HERMITE4)
     write(6,*) "star :",s,ss,atemp(1:NDIM),adots(1:NDIM),potp
#endif

  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE nbody_hermite4_direct_gravity
