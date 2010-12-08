! BOUNDING_BOX.F90
! C. P. Batty & D. A. Hubber - 15/11/2010
! Finds bounding box for particles between pstart and pend in main arrays.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE bounding_box(pstart,pend,rmax,rmin)
  use particle_module, only : parray,PR
  use type_module
#if defined(GRAD_H_SPH)
  use particle_module, only : rextent
#endif
  implicit none

  integer, intent(in) :: pstart               ! id of first particle
  integer, intent(out) ::  pend               ! id of final particle
  real(kind=PR), intent(out) :: rmax(1:NDIM)  ! maximum r value
  real(kind=PR), intent(out) :: rmin(1:NDIM)  ! minimum r value

  integer :: k                         ! counter to loop over dimensions
  integer :: p                         ! counter to loop over particles

  debug2("Finding bounding box for selected particles [bounding_box.F90]")

  rmax(1:NDIM) = -BIG_NUMBER
  rmin(1:NDIM) =  BIG_NUMBER

  do p=pstart,pend
     do k=1,NDIM
        rmax(k) = max(rmax(k), parray(k,p))
        rmin(k) = min(rmin(k), parray(k,p))
     end do
  end do

! Calculate domain extent for bisection method in grad-h method
#if defined(GRAD_H_SPH)
  rextent = maxval(rmax(1:NDIM) - rmin(1:NDIM))
#endif

  return
END SUBROUTINE bounding_box
