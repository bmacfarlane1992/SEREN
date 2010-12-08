! CHECK_NEIGHBOUR_TIMESTEPS.F90
! D. A. Hubber - 2/6/2009
! Check through all neighbour lists of all active particles and calculate 
! theminimum neighbour timestep of all their neighbours.  Also flag any
! particles which are 'neighbouring' sink particles so they can be set to
! the minimum timestep.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE check_neighbour_timesteps
  use interface_module, only : distance2,get_neib_on_fly
  use neighbour_module
  use time_module, only : accdo,n,nminneib,nlevel,nminneib
  use particle_module, only : parray,ptot
#if defined(SINKS)
  use sink_module
#endif
  implicit none

  integer :: i                             ! Aux. particle counter
  integer :: p                             ! Particle id
  integer :: pp                            ! Neighbour id
  integer :: pp_numb                       ! No. of neighbours
  integer, allocatable :: pp_templist(:)   ! List of neighbours
  integer(kind=ILP) :: levelp              ! Timestep level of particle p
#if defined(SINKS)
  integer :: s                             ! Sink id
  real(kind=PR) :: dr(1:NDIM)              ! Relative displacement vector
  real(kind=PR) :: drsqd                   ! Distance squared
#endif

  debug2("Calculate min. neighbour timesteps [check_neighbour_timesteps.F90]")
  debug_timing("CHECK_NEIBSTEPS")

  allocate(pp_templist(1:ptot))


! Loop over all active particles
! ----------------------------------------------------------------------------
  do p=1,ptot
     if (.not. accdo(p)) cycle
     levelp = nlevel(p)

     ! Copy neighbour lists from array or walk the tree
#if defined(NEIGHBOUR_LISTS)
     pp_numb = pptot(p)
     if (pp_numb <= pp_limit) then
        pp_templist(1:pp_numb) = pplist(1:pp_numb,p)
     else
        call get_neib_on_fly(p,parray(SMOO,p),pp_numb,ptot,pp_templist)
     end if
#else
     call get_neib_on_fly(p,parray(SMOO,p),pp_numb,ptot,pp_templist)
#endif

     ! Loop over all neighbours and record maximum neighbouring timestep level
     do i=1,pp_numb
        pp = pp_templist(i)
        if (levelp > nminneib(pp) .and. nminneib(pp) /= -1) &
             &nminneib(pp) = levelp
     end do

  end do
! ----------------------------------------------------------------------------


! Also check sinks if required.  Mark nminneib as -1 to signify that
! minimum timestep (i.e. same as sinks) must be imposed.
! ----------------------------------------------------------------------------
#if defined(SINKS)
  do p=1,ptot
     do s=1,stot
        call distance2(sink(s)%r(1:NDIM),p,dr(1:NDIM),drsqd)
        if (drsqd <= REXTENTSQD*sink(s)%radius*sink(s)%radius) &
             &nminneib(p) = nlevel_sinks
     end do
  end do
#endif
! ----------------------------------------------------------------------------


  return
END SUBROUTINE check_neighbour_timesteps
