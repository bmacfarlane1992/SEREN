! BHHYDROWALK_HGATHER.F90
! D. A. Hubber - 24/01/2008
! Walks hydro tree to find potential neighbour list.  Used in h_gather to 
! find all potential neighbours by gather. (in contrast, BHhydro_walk.F90 
! finds potential neighbours by gather and scatter.) 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydrowalk_hgather(rp,hrange,pp_pot,pp_max,pp_list)
  use definitions
  use interface_module, only : distance3
  use tree_module, only : BHhydro,cmax_hydro
  implicit none

  integer, intent(in)  :: pp_max            ! Max. no. of potential neighbours
  integer, intent(out) :: pp_pot            ! No. of potential neighbours
  integer, intent(out) :: pp_list(1:pp_max) ! Potential neighbour list
  real(kind=PR),intent(in) :: hrange        ! Search range of tree walk
  real(kind=PR),intent(in) :: rp(1:NDIM)    ! Position of particle p

  integer :: i                              ! Auxilary leaf counter
  integer :: c                              ! Cell counter
  real(kind=PR) :: dr(1:NDIM)               ! Relative displacement vector
  real(kind=PR) :: drsqd                    ! Distance squared
  real(kind=PR) :: maxdistsqd               ! Maximum gather distance squared

! Initialize variables.  Start on root cell (c=0) and work our way down.
  c = 0
  pp_pot = 0

! ============================================================================
  do
     call distance3(rp(1:NDIM),BHhydro(c)%r(1:NDIM),dr(1:NDIM),drsqd)
     maxdistsqd = (hrange + BHhydro(c)%rmax)*(hrange + BHhydro(c)%rmax)

     ! Check if search sphere overlaps cell 'neighbour sphere'
     ! -----------------------------------------------------------------------
     if (drsqd <= maxdistsqd) then 

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
END SUBROUTINE BHhydrowalk_hgather
