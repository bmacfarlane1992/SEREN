! BHHYDRO_WALK.F90
! D. A. Hubber - 24/01/2008
! Walks hydro tree to find potential neighbour list
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydro_walk(rp,hrange,pp_pot,pp_max,pp_list)
  use interface_module, only : distance3
  use particle_module
  use tree_module
  implicit none

  integer, intent(in)  :: pp_max            ! Length of pp_list array
  integer, intent(out) :: pp_pot            ! No. of potential neighbours
  integer, intent(out) :: pp_list(1:pp_max) ! Potential neighbour list
  real(kind=PR),intent(in) :: hrange        ! Extent of smoothing kernel
  real(kind=PR), intent(in) :: rp(1:NDIM)   ! Position of particle p

  integer :: i                              ! Auxilary leaf counter
  integer :: c                              ! Cell counter
  real(kind=PR) :: dr(1:NDIM)               ! Relative displacement vector
  real(kind=PR) :: drsqd                    ! Distance squared
  real(kind=PR) :: maxdistsqd1              ! Maximum gather distance squared
  real(kind=PR) :: maxdistsqd2              ! Maximum scatter distance squared

  debug3("Walking tree to find potential neighbour list [BHhydro_walk.F90]","")

! Start on cell 0 and work our way down
  c = 0
  pp_pot = 0

! ============================================================================
  do
     call distance3(rp(1:NDIM),BHhydro(c)%r(1:NDIM),dr(1:NDIM),drsqd)
     maxdistsqd1 = (hrange + BHhydro(c)%rmax) * (hrange + BHhydro(c)%rmax)
     maxdistsqd2 = (BHhydro(c)%hrangemax + BHhydro(c)%rmax) * &
          &(BHhydro(c)%hrangemax + BHhydro(c)%rmax) 

     ! Check if search sphere overlaps cell 'neighbour sphere'
     ! -----------------------------------------------------------------------
     if (drsqd <= maxdistsqd1 .or. drsqd <= maxdistsqd2) then 
       
        ! Check if c is a leaf cell or not
        if (BHhydro(c)%leaf > 0) then

           ! If templist has become too long, flag and return
           if (pp_pot + BHhydro(c)%leaf > pp_max) then
              pp_pot = -1
              return
           end if
              
           ! Add all particles in leaf cell to potential neighbour list
           do i=1,BHhydro(c)%leaf
              pp_pot = pp_pot + 1
              pp_list(pp_pot) = BHhydro(c)%plist(i)
           end do
           c = BHhydro(c)%nextcell
           
        else if (BHhydro(c)%leaf == 0) then
           c = BHhydro(c)%ifopen
        else
           c = BHhydro(c)%nextcell
        end if
     else
        c = BHhydro(c)%nextcell
     end if
     ! -----------------------------------------------------------------------

     if (c > cmax_hydro) exit
  end do
! ============================================================================

#if defined(DEBUG_BHTREEWALK)
  write(6,*) hp,pp_pot,pp_list(1:pp_pot)
#endif

  return
END SUBROUTINE BHhydro_walk
