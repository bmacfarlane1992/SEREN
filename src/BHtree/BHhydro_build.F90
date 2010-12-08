! BHHYDRO_BUILD.F90
! D. A. Hubber - 22/01/2008
! Builds Barnes-Hut tree for a given distribution of particles. Only builds 
! tree for hydro particles (i.e. boundary, icm and 
! self-gravitating gas particles). 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydro_build
  use particle_module
  use type_module
  use tree_module
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none

  logical :: alldone                   ! Flag for when all particles are done
  integer :: c                         ! Cell counter
  integer :: ckid                      ! Auxillary child id
  integer :: clist(1:NCHILD)           ! List of new children cells
  integer :: i                         ! Auxilary loop counter
  integer :: j                         ! Auxilary loop counter
  integer :: k                         ! Dimension counter
  integer :: l                         ! Tree level counter
  integer :: nkids                     ! No. of new children cells
  integer :: nlist                     ! No. of cells to divide on level l
  integer :: nptclin(1:NCHILD)         ! No. of child cell divisions in BH tree
  integer :: p                         ! Particle counter
  integer, allocatable :: templist(:)  ! List of cells to divide
  real(kind=PR) :: cellsize            ! Size of cell
  real(kind=PR) :: rc(1:NDIM)          ! Local copy of position of cell c
  real(kind=PR) :: rmax(1:NDIM)        ! Max. vector of particle positions
  real(kind=PR) :: rmin(1:NDIM)        ! Min. vector of particle positions
  real(kind=PR) :: rp(1:NDIM)          ! Local copy of position of particle p

  debug2("Building Barnes-Hut hydro tree [BHhydro_build.F90]")

  nlist = 0
  allocate(templist(1:cmax_hydro))

! Initialize variables before spatial decomposition
  do p=phtreestart,phtreeend
     cellof(p) = 0
     BHnextptcl(p) = p + 1
  end do

  do c=0,cmax_hydro
     BHtemp(c)%childof(1:NCHILD) = -2
     BHtemp(c)%pfirst = -1
     BHtemp(c)%plast = -1
     BHtemp(c)%leaf = 0
  end do

  BHtemp(0)%pfirst = phtreestart
  BHtemp(0)%plast = phtreeend
  first_cell_hydro(0:LMAX) = 0
  last_cell_hydro(0:LMAX) = 0

! Determine square bounding box of all self-gravitating gas particles and 
! set position of root cell at centre
  call bounding_box(phtreestart,phtreeend,rmax(1:NDIM),rmin(1:NDIM))
  cellsize = 0.5_PR*maxval(rmax(1:NDIM) - rmin(1:NDIM))
  BHtemp(0)%r(1:NDIM) = 0.5_PR*(rmin(1:NDIM) + rmax(1:NDIM))

! Calculate domain extent for bisection method in grad-h method
#if defined(GRAD_H_SPH)
  rextent = maxval(rmax(1:NDIM) - rmin(1:NDIM))
#endif

! Start with root cell and work downwards
  c = 0
  l = 0
  ctot_hydro = 0
  alldone = .false.
  BHtemp(0)%leaf = 0
  BHtemp(0)%nextcell = cmax_hydro + 1

!  BHhydro(0)%r(1:NDIM) = BHtemp(0)%r(1:NDIM)
  BHhydro(0)%leaf = BHtemp(0)%leaf
  BHhydro(0)%nextcell = BHtemp(0)%nextcell
  BHhydro(0)%plist(1:LEAFMAX) = 0

#if defined(DEBUG_BHTREEBUILD)
  write(6,*) "Position of root cell :",BHtemp(0)%r(1:NDIM)
  write(6,*) "Root cell size :",cellsize
#endif

  nlist = 1
  templist(1) = 0


! Main loop
! Loop through deeper and deeper levels until all particles are in leafs
! ============================================================================
  do 
     cellsize = 0.5_PR*cellsize

#if defined(DEBUG_BHTREEBUILD)
     write(6,*) "Building level ",l+1,"   cellsize :",cellsize
     write(6,*) "First cell :",first_cell_hydro(l),&
          &"  Last cell :",last_cell_hydro(l)
#endif

     ! Loop over all cells on parent level
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO IF (l > 0) SCHEDULE(GUIDED,2) DEFAULT(SHARED) &
     !$OMP PRIVATE(c,ckid,k,nkids,nptclin,p,rc,rp)
     do j=1,nlist
        c = templist(j)

#if defined(DEBUG_BHTREEBUILD)
        nptclin(1:NCHILD) = 0
#endif
         
        ! Loop over all particles in current cell, starting with the first
        p = BHtemp(c)%pfirst
        rc(1:NDIM) = BHtemp(c)%r(1:NDIM)
        do 
           rp(1:NDIM) = parray(1:NDIM,p)
             
           ! Determine which quadrant/octant particle p is in 
           ckid = 1
           if (rp(1) > rc(1)) ckid = ckid + 1
           if (rp(2) > rc(2)) ckid = ckid + 2
#if NDIM==3
           if (rp(3) > rc(3)) ckid = ckid + 4
#endif
#if defined(DEBUG_BHTREEBUILD)
           nptclin(ckid) = nptclin(ckid) + 1
#endif
           whichchild(p) = ckid
           BHtemp(c)%childof(ckid) = -1
           if (p == BHtemp(c)%plast) exit
           p = BHnextptcl(p)
        end do
               
#if defined(DEBUG_BHTREEBUILD)
        do k=1,NCHILD
           write(6,*) "Cell ",c,"  nptclin(",k,") :",nptclin(k)
        end do
#endif

     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------


     ! Find all new child cells (Must be done in serial)
     ! -----------------------------------------------------------------------
     do j=1,nlist
        c = templist(j)           
        do k=1,NCHILD
           if (BHtemp(c)%childof(k) == -1) then
              ctot_hydro = ctot_hydro + 1
              BHtemp(c)%childof(k) = ctot_hydro
           end if
        end do
     end do
     ! ----------------------------------------------------------------------- 


     ! Find number of particles in each new child cell
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO IF (l > 0) SCHEDULE(GUIDED,2) DEFAULT(SHARED) &
     !$OMP PRIVATE(c,ckid,clist,i,k,nkids,nptclin,p,rc)
     do j=1,nlist
        c = templist(j)
        nptclin(1:NCHILD) = 0
        p = BHtemp(c)%pfirst
        
        do
           ckid = BHtemp(c)%childof(whichchild(p))
           cellof(p) = ckid
           if (BHtemp(ckid)%pfirst == -1) then
              BHtemp(ckid)%pfirst = p
           else
              BHnextptcl(BHtemp(ckid)%plast) = p
           end if
           BHtemp(ckid)%plast = p
           nptclin(whichchild(p)) = nptclin(whichchild(p)) + 1
           if (p == BHtemp(c)%plast) exit
           p = BHnextptcl(p)
        end do
        

        ! Now create child cells if there are particles in that division
        nkids = 0
        rc(1:NDIM) = BHtemp(c)%r(1:NDIM)
        
        ! --------------------------------------------------------------------
        do k=1,NCHILD
           if (nptclin(k) == 0) cycle
           ckid = BHtemp(c)%childof(k)

           ! If leafcell, flag it and store children ids
           if (nptclin(k) > 0 .and. nptclin(k) <= LEAFMAX) then
              BHtemp(ckid)%leaf = nptclin(k)
              i = 0
              p = BHtemp(ckid)%pfirst   
              do
                 i = i + 1
                 BHtemp(ckid)%plist(i) = p
                 if (p == BHtemp(ckid)%plast) exit
                 p = BHnextptcl(p)
              end do

           ! If not, record particles within for next level search
           else if (nptclin(k) > LEAFMAX) then
              BHtemp(ckid)%leaf = 0
           end if
           
           ! Assign positions of child cells
           if (k==1 .or. k==3 .or. k==5 .or. k==7) then
              BHtemp(ckid)%r(1) = rc(1) - cellsize
           else
              BHtemp(ckid)%r(1) = rc(1) + cellsize
           end if
           
           if (k==1 .or. k==2 .or. k==5 .or. k==6) then
              BHtemp(ckid)%r(2) = rc(2) - cellsize
           else
              BHtemp(ckid)%r(2) = rc(2) + cellsize
           end if          
#if NDIM==3
           if (k < 5) then
              BHtemp(ckid)%r(3) = rc(3) - cellsize
           else
              BHtemp(ckid)%r(3) = rc(3) + cellsize
           end if
#endif

           ! Record id of cells for following linked list construction
           nkids = nkids + 1
           clist(nkids) = ckid

        end do
        ! --------------------------------------------------------------------

        ! Set up linked lists from parent to children cells
        BHtemp(c)%ifopen = clist(1)
        do k=1,nkids-1
           ckid = clist(k)
           BHtemp(ckid)%nextcell = clist(k + 1)
        end do
        BHtemp(clist(nkids))%nextcell = BHtemp(c)%nextcell

        BHhydro(c)%ifopen = BHtemp(c)%ifopen
        do k=1,nkids
           ckid = clist(k)
           BHhydro(ckid)%leaf = BHtemp(ckid)%leaf
           BHhydro(ckid)%nextcell = BHtemp(ckid)%nextcell
           BHhydro(ckid)%ifopen = BHtemp(ckid)%ifopen
           BHhydro(ckid)%plist(1:LEAFMAX) = BHtemp(ckid)%plist(1:LEAFMAX)
        end do
        
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

     ! Record first and last cells in newly created level
     l = l + 1
     first_cell_hydro(l) = last_cell_hydro(l-1) + 1
     last_cell_hydro(l) = ctot_hydro
     
#if defined(DEBUG_BHTREEBUILD)
     write(6,*) "first_cell(",l,") :",first_cell_hydro(l),&
          &"  last_cell(",l,") :",last_cell_hydro(l)
#endif

     ! Check that all new cells are children cells
     nlist = 0
     alldone = .true.
     do c=first_cell_hydro(l),last_cell_hydro(l)
        if (BHtemp(c)%leaf == 0) then
           alldone = .false.
           nlist = nlist + 1
           templist(nlist) = c
        end if
     end do
     
#if defined(DEBUG_BHTREEBUILD)
     write(6,*) "Level ",l,alldone
#endif

     ! If all children are leafs, finish building tree
     if (alldone) exit
     alldone = .false.

     ! Check we have not reached maximum level
     if (l >= LMAX) stop 'Reached maximum level'

  end do
! ============================================================================


! Record total number of levels in the tree
  ltot_hydro  =  l

#if defined(DEBUG_BHTREEBUILD)
  write(6,*) "Total number of levels in tree :",ltot_hydro
  do l=0,ltot_hydro
     write(6,*) "first_cell(",l,") :",first_cell_hydro(l),&
          &"     last_cell(",l,") :",last_cell_hydro(l)
  end do
#endif

  deallocate(templist)

  return 
END SUBROUTINE BHhydro_build
