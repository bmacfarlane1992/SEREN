! GATHER_NEIB_ON_FLY.F90
! C. P. Batty & D. A. Hubber - 3/6/2007
! Calculates neighbour list for particle p (using pure gather method). 
! Returns no. of neibs plus list rather than storing in arrays.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gather_neib_on_fly(p,hp,pp_tot,pp_list)
  use interface_module, only : BHhydrowalk_hgather,binary_gather_walk,distance2
  use particle_module
  use neighbour_module
  implicit none

  integer, intent(in) :: p                     ! chosen particle
  real(kind=PR), intent(in) :: hp              ! h of particle p
  integer, intent(out) :: pp_tot               ! number of pot. neibs
  integer, intent(out) :: pp_list(1:LISTSIZE)  ! list of pot. neibs

  integer :: pp                         ! neighbouring particles (p')
  real(kind=PR) :: dr(1:NDIM)           ! vector displacements (p'-p)
  real(kind=PR) :: drsqd                ! p'-p separation squared
  real(kind=PR) :: hrangesqd            ! particle radius (2*h_p) squared
  real(kind=PR) :: rp(1:NDIM)           ! Position of particle p
#if defined(BH_TREE) || defined(BINARY_TREE)
  integer :: i                          ! counter in neighbour search
  integer :: pp_max                     ! length of pp_potlist array
  integer :: pp_pot                     ! no. of potential neighbours
  integer, allocatable :: pp_potlist(:) ! list of potential neighbours
#endif
#if defined(DEBUG_GATHER_NEIB)
  integer :: n_gather                   ! Number of neighbours by gather
#if defined(HMASS)
  real(kind=PR) :: mfrac                ! Mass fraction
  real(kind=PR) :: mgather              ! Mass by gather
  real(kind=PR) :: mp                   ! Mass of particle p
  real(kind=PR) :: mtot                 ! Total mass
#endif
#endif

  debug3("Obtaining neighbour list [gather_neib_on_fly.F90] for particle ", p)

! Make local copies of particle position and initialize other variables
  rp(1:NDIM) = parray(1:NDIM,p)
  hrangesqd = KERNRANGESQD*hp*hp
  pp_tot = 0
#if defined(DEBUG_GATHER_NEIB)
  n_gather = 0
#if defined(HMASS)
  mp = parray(MASS,p)
  mtot = 0.0_PR
#endif
#endif

! Obtain potential neighbour list by tree walk if using tree
! ----------------------------------------------------------------------------
#if defined(BH_TREE) || defined(BINARY_TREE)
  pp_max = min(LISTSIZE,ptot)
  allocate(pp_potlist(1:pp_max))
  do 
     pp_pot = 0
#if defined(BH_TREE)
     call BHhydrowalk_hgather(rp(1:NDIM),KERNRANGE*hp,pp_pot,pp_max,pp_potlist)
#elif defined(BINARY_TREE)
     call binary_gather_walk(rp(1:NDIM),hp,pp_pot,pp_max,pp_potlist)
#endif
     if (pp_pot < 0) then
        pp_max = ptot
        deallocate(pp_potlist)
        allocate(pp_potlist(1:pp_max))
     else
        exit
     end if
  end do
#endif


! Loop over potential list or all particles depending on flags
! ----------------------------------------------------------------------------
#if defined(BH_TREE) || defined(BINARY_TREE)
  do i=1,pp_pot
     pp = pp_potlist(i)
#else
  do pp=1,ptot
#endif
     if (pp == p) cycle

     ! Is potential neighbour within 2hp of particle p or is particle p 
     ! within 2hpp within neighbour?
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     if (drsqd <= hrangesqd) then
        pp_tot = pp_tot + 1
        pp_list(pp_tot) = pp
     end if

     ! Record numbers of neighbours for debugging 
#if defined(DEBUG_GATHER_NEIB)
     if (drsqd <= hrangesqd) n_gather = n_gather + 1
#if defined(HMASS)
     if (drsqd <= hrangesqd) mtot = mtot + parray(MASS,pp)
#endif        
#endif

  end do
! ----------------------------------------------------------------------------


#if defined(DEBUG_GATHER_NEIB)
  write(6,*) "h(",p,") : ",hp,"   n_gather : ", n_gather, pp_pot
#if defined(HMASS)
  mgather = mp * real(pp_gather,PR)
  mfrac = mtot / mp
  write(6,*) "mp :",mp," mgather :",mgather,"  mtot :",mtot,"  mfrac :",mfrac
#endif
#endif

#if defined(BH_TREE) || defined(BINARY_TREE)
  deallocate(pp_potlist)
#endif

  return
END SUBROUTINE gather_neib_on_fly
