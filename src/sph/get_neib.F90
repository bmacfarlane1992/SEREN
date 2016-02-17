! GET_NEIB.F90
! C. P. Batty & D. A. Hubber - 3/6/2007
! Calculates neighbour list for particle p (using gather/scatter method).
! Stores no. of neibs plus list in arrays.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE get_neib(p)
  use interface_module, only : BHhydro_walk,binary_neibfind,distance2
  use neighbour_module
  use particle_module, only : parray,ptot
  implicit none

  integer, intent(in) :: p               ! chosen particle

  integer :: pp                          ! neighbouring particles (p')
  integer :: pp_tot                      ! no. of neighbours
  integer :: pp_templist(1:pp_limit)     ! list of neighbours
  real(kind=PR) :: dr(1:NDIM)            ! vector displacements (p'-p)
  real(kind=PR) :: drsqd                 ! p'-p separation squared
  real(kind=PR) :: hrangesqd_p           ! particle radius (2*h_p) squared
  real(kind=PR) :: hrangesqd_pp          ! particle radius (2*h_pp) squared
  real(kind=PR) :: rp(1:NDIM)            ! Position of particle p
#if defined(BH_TREE) || defined(BINARY_TREE)
  integer :: i                           ! aux. counter in neighbour search
  integer :: pp_max                      ! length of pp_potlist array
  integer :: pp_pot                      ! no. of potential neighbours found
  integer, allocatable :: pp_potlist(:)  ! list of potential neighbours
#endif
#if defined(DEBUG_GET_NEIB)
  integer :: n_gather                    ! No. of neighbours by gather
  integer :: n_scatter                   ! No. of neighbours by scatter
#if defined(HMASS)
  real(kind=PR) :: mfrac                 ! Mass fraction
  real(kind=PR) :: mgather               ! Mass by gather
  real(kind=PR) :: mp                    ! Mass of particle p
  real(kind=PR) :: mtot                  ! Total mass
#endif
#endif

  debug3("Obtaining neighbour list [get_neib.F90] for particle ", p)


! Make local copies of particle position and smoothing length
  rp(1:NDIM) = parray(1:NDIM,p)
  hrangesqd_p = KERNRANGESQD*parray(SMOO,p)*parray(SMOO,p)
  pp_tot = 0
#if defined(DEBUG_GET_NEIB)
  n_gather = 0 ; n_scatter = 0
#if defined(HMASS)
  mp   = parray(MASS,p)
  mtot = 0.0_PR
#endif
#endif

! If using a tree, obtain a smaller potential neighbour list
#if defined(BH_TREE) || defined(BINARY_TREE)
  pp_max = min(LISTSIZE,ptot)
  allocate(pp_potlist(1:pp_max))
  do 
#if defined(BH_TREE)
     call BHhydro_walk(rp(1:NDIM),KERNRANGE*parray(SMOO,p),&
          &pp_pot,pp_max,pp_potlist)
#elif defined(BINARY_TREE)
     call binary_neibfind(rp,parray(SMOO,p),pp_pot,pp_max,pp_potlist)
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
     hrangesqd_pp = KERNRANGESQD*parray(SMOO,pp)*parray(SMOO,pp)

     ! Is potential neighbour within 2*hp of particle p or is particle p 
     ! within 2*hpp within neighbour?  If list is larger than what is allowed, 
     ! don't worry about recording list but record number of neighbours
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     if ((drsqd <= hrangesqd_p) .or. (drsqd <= hrangesqd_pp)) then
        pp_tot = pp_tot + 1
        pp_templist(min(pp_tot,pp_limit)) = pp
     end if

     ! Record numbers of neighbours for debugging 
#if defined(DEBUG_GET_NEIB)
     if (drsqd <= hrangesqd_p) n_gather = n_gather + 1
     if (drsqd <= hrangesqd_pp) n_scatter = n_scatter + 1
#if defined(HMASS)
     if (drsqd <= hrangesqd_p) mtot = mtot + parray(MASS,pp)
#endif        
#endif

  end do
! ----------------------------------------------------------------------------


! Record neighbour info in main arrays.  Only record neighbour ids if 
! less than maximum allowed number of neighbours. 
  pptot(p) = pp_tot
  if (pp_tot <= pp_limit) then
     pplist(1:pp_tot,p) = pp_templist(1:pp_tot)
  end if

#if defined(DEBUG_GET_NEIB)
#if defined(BH_TREE) || defined(BINARY_TREE)
  write(6,*) "h(",p,") : ",parray(SMOO,p),"   n_gather : ", n_gather,&
       &"  n_scatter : ",n_scatter,"  pp_max : ",pp_max,"  pp_tot : ",pp_tot
#elif
  write(6,*) "h(",p,") : ",parray(SMOO,p),"   n_gather : ", n_gather,&
       &"  n_scatter : ",n_scatter,"  pp_tot : ",pp_tot
#endif
#if defined(HMASS)
  mgather = mp*real(pp_gather,PR)
  mfrac = mtot / mp
  write(6,*) "mp :",mp," mgather :",mgather,"  mtot :",mtot,"  mfrac :",mfrac
#endif
#endif

#if defined(BH_TREE) || defined(BINARY_TREE)
  if (allocated(pp_potlist)) deallocate(pp_potlist)
#endif

  return
END SUBROUTINE get_neib
