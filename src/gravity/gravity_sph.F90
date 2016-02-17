! GRAVITY_SPH.F90
! C. P. Batty & D. A. Hubber - 24/6/2007
! Calculates gravitational force between masses mp and mpp at positions 
! rp and rpp.  Passes positions, masses and smoothing lengths rather 
! than particle or sink ids. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gravity_sph(invhp,hpp,mpp,rp,rpp,atemp,dpotp)
  use interface_module, only : distance3,ewald_force
  use definitions
  use kernel_module
  implicit none

  real(kind=PR), intent(in)  :: hpp            ! Smoothing length of pp
  real(kind=PR), intent(in)  :: invhp          ! Smoothing length of p
  real(kind=PR), intent(in)  :: mpp            ! Mass of pp
  real(kind=PR), intent(in)  :: rp(1:NDIM)     ! Position of p
  real(kind=PR), intent(in)  :: rpp(1:NDIM)    ! Position of pp
  real(kind=PR), intent(out) :: atemp(1:NDIM)  ! Acceleration vector
  real(kind=PR), intent(out) :: dpotp          ! Gravitational potential

  integer :: kern                    ! Kernel table element
  real(kind=PR) :: dr(1:NDIM)        ! Relative position vector
  real(kind=PR) :: drmag             ! Distance
  real(kind=PR) :: drsqd             ! Distance squared
  real(kind=PR) :: grav              ! Aux. variable
  real(kind=PR) :: invdrmag          ! 1.0 / drmag
  real(kind=PR) :: invhpp            ! 1.0 / hpp
#if defined(EWALD)
  real(kind=PR) :: eaccel(1:NDIM)    ! Ewald grav. acceleration
#endif

! Calculate relative displacement vector between bodies
#if defined(EWALD)
  call distance3(rp(1:NDIM),rpp(1:NDIM),dr(1:NDIM),drsqd)
#else
  dr(1) = rpp(1) - rp(1)
#if NDIM==2 || NDIM==3
  dr(2) = rpp(2) - rp(2)
#endif
#if NDIM==3
  dr(3) = rpp(3) - rp(3)
#endif
#if NDIM==1
  drsqd = dr(1)*dr(1)
#elif NDIM==2
  drsqd = dr(1)*dr(1) + dr(2)*dr(2)
#elif NDIM==3
  drsqd = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
#endif
#endif
  drmag = sqrt(drsqd) + SMALL_NUMBER
  invdrmag = 1.0_PR / drmag

! Calculate contribution of gravity from particle p
  if (invdrmag > INVKERNRANGE*invhp) then
     kern = int(HALFKERNTOT*drmag*invhp)
     kern = min(kern,KERNTOT)
     grav = invdrmag*invhp*invhp*w5(kern)
     dpotp = invhp*w6(kern)
  else
     grav = invdrmag**3
     dpotp = invdrmag
  end if

! Calculate contribution of gravity from particle pp
  if (drmag < KERNRANGE*hpp) then
     invhpp = 1.0_PR/hpp
     kern = int(HALFKERNTOT*drmag*invhpp)
     kern = min(kern,KERNTOT)
     grav = grav + invdrmag*invhpp*invhpp*w5(kern)
     dpotp = dpotp + invhpp*w6(kern)
  else
     grav = grav + invdrmag**3
     dpotp = dpotp + invdrmag
  end if

! Record acceleration in output vector      
  atemp(1:NDIM) = 0.5_PR*mpp*grav*dr(1:NDIM)
  dpotp = 0.5_PR*mpp*dpotp

! Add Ewald periodic correction if required
#if defined(EWALD)
  call ewald_force(dr(1:NDIM),mpp,eaccel(1:NDIM))
  atemp(1:NDIM) = atemp(1:NDIM) + eaccel(1:NDIM)
#endif

  return
END SUBROUTINE gravity_sph
