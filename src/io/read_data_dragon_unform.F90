! READ_DATA_DRAGON_UNFORM.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Reads in initial conditions file in DRAGON binary format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data_dragon_unform(in_file)
  use particle_module
  use hydro_module
  use type_module
  use sink_module
  use time_module
  implicit none

  character(len=*), intent(in) :: in_file    ! unformatted DRAGON snapshot

  integer :: boundaryslot                    ! point in list to boundary ptcl
  integer :: gasslot                         ! point in list to insert ptcl
  integer :: icmslot                         ! point in list to icm particle
  integer :: idata(1:20)                     ! dummy integers
  integer :: k                               ! dimension counter
  integer :: p                               ! counter to loop over particles
  integer :: psplit                          ! counter for particle splitting
  integer :: pdead                           ! counter for accreted particles
  integer :: sinkslot                        ! point in list sink particle
  integer, allocatable :: idummy1(:)         ! dummy array for integers
  integer, allocatable :: ptype(:)           ! particle types
  real(kind=PR) :: rdata(1:50)               ! dummy reals
  real(kind=PR), allocatable :: rdummy1(:)   ! dummy array for reals
  real(kind=PR), allocatable :: rdummy3(:,:) ! dummy array for 3D reals
#if defined(SINKS)
  integer :: s                               ! sink counter
#endif

  debug1("Reading in unformatted data file : "//trim(in_file)//" [read_data_dragon_unform.F90]")

! Initialise counters
  pboundary = 0
  psplit    = 0
  pdead     = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 0

! Open snapshot file
  open(1, file=in_file, status="old", form="unformatted")
 
! First, read in header information 
  read(1) idata
  read(1) rdata

! Assign variables for important information  
  ptot      = idata(1)
  nsteps    = idata(2)
  snapshot  = idata(4)
  pgas_orig = idata(20)
  time      = real(rdata(1),DP)
  lastsnap  = real(rdata(2),DP)
  mgas_orig = real(rdata(50),DP)

  allocate(ptype(1:ptot))
  allocate(idummy1(1:ptot))
  allocate(rdummy1(1:ptot))
  allocate(rdummy3(1:NDIM,1:ptot))


! First pass to get numbers of particles of each type
! ----------------------------------------------------------------------------
  read(1) rdummy3  ! Positions
  read(1) rdummy3  ! Velocities
  read(1) rdummy1  ! Temperatures
  read(1) rdummy1  ! Smoothing lengths
  read(1) rdummy1  ! Densities
  read(1) rdummy1  ! Masses
  read(1) idummy1  ! Particle types

  do p=1,ptot
     ptype(p) = idummy1(p)
     if (ptype(p) == -1) stot = stot + 1
     if (ptype(p) == 0)  pdead = pdead + 1
     if (ptype(p) == 1)  pgas = pgas + 1
     if (ptype(p) == 4)  psplit = psplit + 1
     if (ptype(p) == 6)  pboundary = pboundary + 1
     if (ptype(p) == 9)  picm = picm + 1
  end do

! Reset ptot to only account for gas, boundary and intercloud particles
  ptot = ptot - (stot + pdead)
  
  write(6,*) "Particles   = ", ptot, "   Sinks = ", stot
  write(6,*) "Gas         = ", pgas
  write(6,*) "Boundary    = ", pboundary
  write(6,*) "Intercloud  = ", picm
  write(6,*) "Dark matter = ", pcdm
  write(6,*) "Splitting   = ", psplit
  write(6,*) "Accreted    = ", pdead
  if (psplit /= 0) stop "Fatal error: particle splitting not supported"
  if (ptot /= (pgas + picm + pboundary)) stop "Fatal error: particles do not add up"

! Allocate memory now we know 'correct' value of ptot
  call allocate_memory


! Second pass to assign data
! ----------------------------------------------------------------------------
  rewind(1)

  read(1) idata
  read(1) rdata

! Positions
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  sinkslot = 1
  read(1) rdummy3
  do p=1,(ptot + (stot + pdead))
     if (ptype(p) == 6) then
        parray(1:NDIM,boundaryslot) = rdummy3(1:NDIM,p)
        porig(boundaryslot) = p
        boundaryslot = boundaryslot + 1
     end if
     if (ptype(p) == 9) then
        parray(1:NDIM,icmslot) = rdummy3(1:NDIM,p)
        porig(icmslot) = p
        icmslot = icmslot + 1
     end if
     if (ptype(p) == 1) then
        parray(1:NDIM,gasslot) = rdummy3(1:NDIM,p)
        porig(gasslot) = p
        gasslot = gasslot + 1
     end if
     if (ptype(p) == -1) then
        sink(sinkslot)%r(1:NDIM) = rdummy3(1:NDIM,p)
        sinkslot = sinkslot + 1
     end if
  end do

! Velocities
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  sinkslot = 1
  read(1) rdummy3
  do p=1,(ptot + (stot + pdead))
     if (ptype(p) == 6) then
        v(1:VDIM,boundaryslot) = rdummy3(1:VDIM,p)
        boundaryslot = boundaryslot + 1
     end if
     if (ptype(p) == 9) then
        v(1:VDIM,icmslot) = rdummy3(1:VDIM,p)
        icmslot = icmslot + 1
     end if
     if (ptype(p) == 1) then
        do k=1,NDIM
           v(1:VDIM,gasslot) = rdummy3(1:VDIM,p)
        end do
        gasslot = gasslot + 1
     end if
     if (ptype(p) == -1) then
        sink(sinkslot)%v(1:VDIM) = rdummy3(1:VDIM,p)
        sinkslot = sinkslot + 1
     end if
  end do

! Temperatures
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  sinkslot = 1
  read(1) rdummy1
  do p=1,(ptot + (stot + pdead))
     if (ptype(p) == 6) then
        temp(boundaryslot) = rdummy1(p)
        boundaryslot = boundaryslot + 1
     end if
     if (ptype(p) == 9) then
        temp(icmslot) = rdummy1(p)
        icmslot = icmslot + 1
     end if
     if (ptype(p) == 1) then
        temp(gasslot) = rdummy1(p)
        gasslot = gasslot + 1
     end if
  end do

! Smoothing lengths
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  sinkslot = 1
  read(1) rdummy1
  do p=1,(ptot + (stot + pdead))
     if (ptype(p) == 6) then
        parray(SMOO,boundaryslot) = rdummy1(p)
        boundaryslot = boundaryslot + 1
     end if
     if (ptype(p) == 9) then
        parray(SMOO,icmslot) = rdummy1(p)
        icmslot = icmslot + 1
     end if
     if (ptype(p) == 1) then
        parray(SMOO,gasslot) = rdummy1(p)
        gasslot = gasslot + 1
     end if
     if (ptype(p) == -1) then
        sink(sinkslot)%h = rdummy1(p)
        sinkslot = sinkslot + 1
     end if
  end do

! Densities
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  sinkslot = 1
  read(1) rdummy1
  do p=1,(ptot + (stot + pdead))
     if (ptype(p) == 6) then
        rho(boundaryslot) = rdummy1(p)
        boundaryslot = boundaryslot + 1
     end if
     if (ptype(p) == 9) then
        rho(icmslot) = rdummy1(p)
        icmslot = icmslot + 1
     end if
     if (ptype(p) == 1) then
        rho(gasslot) = rdummy1(p)
        gasslot = gasslot + 1
     end if
  end do

! Masses
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  sinkslot = 1
  read(1) rdummy1
  do p=1,(ptot + (stot + pdead))
     if (ptype(p) == 6) then
        parray(MASS,boundaryslot) = rdummy1(p)
        boundaryslot = boundaryslot + 1
     end if
     if (ptype(p) == 9) then
        parray(MASS,icmslot) = rdummy1(p)
        icmslot = icmslot + 1
     end if
     if (ptype(p) == 1) then
        parray(MASS,gasslot) = rdummy1(p)
        gasslot = gasslot + 1
     end if
     if (ptype(p) == -1) then
        sink(sinkslot)%m = rdummy1(p)
        sinkslot = sinkslot + 1
     end if
  end do

  ! Close file once finished
  close(1)
! ----------------------------------------------------------------------------

  deallocate(rdummy3)
  deallocate(rdummy1)
  deallocate(idummy1)
  deallocate(ptype)

! Convert temperatures to internal energies if required
! ----------------------------------------------------------------------------
#if defined(INTERNAL_ENERGY)
  do p=1,ptot
     u(p) = Pconst*temp(p)/(gamma - 1.0_PR)
  end do
#endif

! Initialise some sink variables not recorded by Dragon format
! ----------------------------------------------------------------------------
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        sink(s)%accrete = .true.
        sink(s)%static  = .false.
        sink(s)%radius  = KERNRANGE*sink(s)%h
     end do
  end if
#endif


! Verification that particles add up correctly
  if (boundaryslot-1 /= pboundary) stop 'boundary particles do not match'
  if (icmslot-1 /= picm + pboundary) stop 'intercloud particles do not match'
  if (gasslot-1 /= pgas + picm + pboundary) stop 'gas particles do not match'
  if (sinkslot-1 /= stot) stop 'sinks do not match'

  return
END SUBROUTINE read_data_dragon_unform
