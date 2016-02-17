
! READ_DATA_DRAGON_FORM.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Reads in initial conditions file in DRAGON ASCII format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data_dragon_form(in_file)
  use particle_module
  use hydro_module
  use type_module
  use sink_module
  use time_module
use filename_module, only: restart_log,run_id,run_dir
  implicit none


  character(len=*), intent(in) :: in_file    ! formatted DRAGON snapshot

  integer :: boundaryslot                    ! point in list to boundary ptcl
  integer :: cdmslot                         ! point in list to cdm particle
  integer :: gasslot                         ! point in list to insert ptcl
  integer :: icmslot                         ! point in list to icm particle
  integer :: idata(1:20)                     ! dummy integers
  integer :: p                               ! counter to loop over particles
  integer :: psplit                          ! counter for particle splitting
  integer :: pdead                           ! counter for accreted particles
  integer :: sinkslot                        ! point in list sink particle
  integer, allocatable :: ptype(:)           ! particle types
  real(kind=SP) :: rdata(1:50)               ! dummy reals
#if defined(SINKS)
  integer :: s                               ! sink counter
#endif

 logical :: findid=.TRUE.             		    ! true when we need to assign particle identifiers

  real(kind=PR) :: raux

  debug1("Reading in formatted data file : "//trim(in_file)//" [read_data_dragon_form.F90]")

! check to see if restart file exists
  restart_log = trim(adjustl(run_dir))//trim(adjustl(run_id))//".restart"

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
  open(1, file=in_file, status="old", form="formatted")

! First, read in header information 
  do p=1,20
     read(1,*) idata(p)
  end do
  do p=1,50
     read(1,*) rdata(p)
  end do

! Assign variables for important information
  ptot      = idata(1)
  nsteps    = idata(2)
  snapshot  = idata(4)
  pgas_orig = idata(20)
  time      = real(rdata(1),DP)
  lastsnap  = real(rdata(2),DP)
  mgas_orig = real(rdata(50),DP)

  allocate(ptype(1:ptot))


! First pass to get numbers of particles of each type
! ----------------------------------------------------------------------------
  do p=1,ptot*2
#if NDIM == 1
     read(1,*) rdata(1)
#elif NDIM == 2
     read(1,*) rdata(1), rdata(2)
#else
     read(1,*) rdata(1), rdata(2), rdata(3)
#endif
  end do
  do p=1,ptot*4
     read(1,*) rdata(1)
  end do
  do p=1,ptot
!     read(1,*) ptype(p)
     read(1,*) raux
     ptype(p) = int(raux)
     if (ptype(p) == -1) stot = stot + 1
     if (ptype(p) == 0)  pdead = pdead + 1
     if (ptype(p) == 1)  pgas = pgas + 1
     if (ptype(p) == 4)  psplit = psplit + 1
     if (ptype(p) == 6)  pboundary = pboundary + 1
     if (ptype(p) == 9)  picm = picm + 1
     if (ptype(p) == 10) pcdm = pcdm + 1
  end do

! Reset ptot to only account for gas, boundary and intercloud particles
  ptot = ptot - (stot + pdead)
  call types

  write(6,*) "Particles   = ", ptot, "   Sinks = ", stot
  write(6,*) "Gas         = ", pgas
  write(6,*) "Boundary    = ", pboundary
  write(6,*) "Intercloud  = ", picm
  write(6,*) "Dark matter = ", pcdm
  write(6,*) "Splitting   = ", psplit
  write(6,*) "Accreted    = ", pdead
  if (psplit /= 0) stop "Fatal error: particle splitting not supported"
  if (ptot /= (pgas + picm + pboundary + pcdm)) &
       &stop "Fatal error: particles do not add up"

! Allocate memory now we know 'correct' value of ptot
  call allocate_memory

! Second pass to assign data
! ----------------------------------------------------------------------------
  rewind(1)

  do p=1,20
     read(1,*) idata(p)
  end do
  do p=1,50
     read(1,*) rdata(p)
  end do




#if defined(PARTICLE_ID)

  inquire(file=restart_log,exist=findid)

#endif



! Positions
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot + pdead))
#if NDIM == 1
     read(1,*) rdata(1)
#elif NDIM == 2
     read(1,*) rdata(1), rdata(2)
#else
     read(1,*) rdata(1), rdata(2), rdata(3)
#endif
     if (ptype(p) == 6) then
        parray(1:NDIM,boundaryslot) = rdata(1:NDIM)
        if (.not. findid) porig(boundaryslot) = p
        boundaryslot = boundaryslot + 1
     else if (ptype(p) == 9) then
        parray(1:NDIM,icmslot) = rdata(1:NDIM)
         if (.not. findid) porig(icmslot) = p
        icmslot = icmslot + 1
     else if (ptype(p) == 1) then
        parray(1:NDIM,gasslot) = rdata(1:NDIM)
         if (.not. findid) porig(gasslot) = p
        gasslot = gasslot + 1
     else if (ptype(p) == 10) then
        parray(1:NDIM,cdmslot) = rdata(1:NDIM)
         if (.not. findid) porig(cdmslot) = p
        cdmslot = cdmslot + 1
     else if (ptype(p) == -1) then
        sink(sinkslot)%r(1:NDIM) = rdata(1:NDIM)
        sinkslot = sinkslot + 1
        ! no porig for sinks
     end if

  end do



! Velocities
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot + pdead))
#if NDIM == 1
     read(1,*) rdata(1)
#elif NDIM == 2
     read(1,*) rdata(1), rdata(2)
#else
     read(1,*) rdata(1), rdata(2), rdata(3)
#endif
     if (ptype(p) == 6) then
        v(1:VDIM,boundaryslot) = rdata(1:VDIM)
        boundaryslot = boundaryslot + 1
     else if (ptype(p) == 9) then
        v(1:VDIM,icmslot) = rdata(1:VDIM)
        icmslot = icmslot + 1
     else if (ptype(p) == 1) then
        v(1:VDIM,gasslot) = rdata(1:VDIM)
        gasslot = gasslot + 1
     else if (ptype(p) == 10) then
        v(1:VDIM,cdmslot) = rdata(1:VDIM)
        cdmslot = cdmslot + 1
     else if (ptype(p) == -1) then
        sink(sinkslot)%v(1:VDIM) = rdata(1:VDIM)
        sinkslot = sinkslot + 1
     end if
  end do

! Temperatures
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot + pdead))
     read(1,*) rdata(1)
     if (ptype(p) == 6) then
        temp(boundaryslot) = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (ptype(p) == 9) then
        temp(icmslot) = rdata(1)
        icmslot = icmslot + 1
     else if (ptype(p) == 1) then
        temp(gasslot) = rdata(1)
        gasslot = gasslot + 1
     else if (ptype(p) == 10) then
        temp(cdmslot) = rdata(1)
        cdmslot = cdmslot + 1
     end if
  end do

! Smoothing lengths
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot + pdead))
     read(1,*) rdata(1)
     if (ptype(p) == 6) then
        parray(SMOO,boundaryslot) = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (ptype(p) == 9) then
        parray(SMOO,icmslot) = rdata(1)
        icmslot = icmslot + 1
     else if (ptype(p) == 1) then
        parray(SMOO,gasslot) = rdata(1)
        gasslot = gasslot + 1
     else if (ptype(p) == 10) then
        parray(SMOO,cdmslot) = rdata(1)
        cdmslot = cdmslot + 1
     else if (ptype(p) == -1) then
        sink(sinkslot)%h = rdata(1)
        sinkslot = sinkslot + 1
     end if
  end do

! Densities
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot + pdead))
     read(1,*) rdata(1)
     if (ptype(p) == 6) then
        rho(boundaryslot) = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (ptype(p) == 9) then
        rho(icmslot) = rdata(1)
        icmslot = icmslot + 1
     else if (ptype(p) == 1) then
        rho(gasslot) = rdata(1)
        gasslot = gasslot + 1
     else if (ptype(p) == 10) then
        rho(cdmslot) = rdata(1)
        cdmslot = cdmslot + 1
     end if
  end do

! Masses
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot + pdead))
     read(1,*) rdata(1)
     if (ptype(p) == 6) then
        parray(MASS,boundaryslot) = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (ptype(p) == 9) then
        parray(MASS,icmslot) = rdata(1)
        icmslot = icmslot + 1
     else if (ptype(p) == 1) then
        parray(MASS,gasslot) = rdata(1)
        gasslot = gasslot + 1
     else if (ptype(p) == 10) then
        parray(MASS,cdmslot) = rdata(1)
        cdmslot = cdmslot + 1
     else if (ptype(p) == -1) then
        sink(sinkslot)%m = rdata(1)
        sinkslot = sinkslot + 1
     end if
  end do


#if defined(PARTICLE_ID)

 if (findid) then 
! Particle Types
! ----------------------------------------------------------------------------
  do p=1,(ptot + (stot + pdead))
     read(1,*) idata(1)    
  end do

! Particle ids
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot + pdead))
     read(1,*) idata(1)
     if (ptype(p) == 6) then
        porig(boundaryslot) = idata(1)
        boundaryslot = boundaryslot + 1
     else if (ptype(p) == 9) then
        porig(icmslot) = idata(1)
        icmslot = icmslot + 1
     else if (ptype(p) == 1) then
        porig(gasslot) = idata(1)
        gasslot = gasslot + 1
     else if (ptype(p) == 10) then
        porig(cdmslot) = idata(1)
        cdmslot = cdmslot + 1
     else if (ptype(p) == -1) then
        sinkslot = sinkslot + 1
        ! no porig for sinks (equal to zero)
     end if
  end do

endif 

#endif

! Close file once finished
  close(1)
! ----------------------------------------------------------------------------

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
  do s=1,stot
     sink(s)%accrete = .true.
     sink(s)%static  = .false.
     sink(s)%radius  = KERNRANGE*sink(s)%h
  end do
#endif

! Verification that particles add up correctly
  if (boundaryslot-1 /= pboundary) stop 'boundary particles do not match'
  if (icmslot-1 /= picm + pboundary) stop 'intercloud particles do not match'
  if (gasslot-1 /= pgas + picm + pboundary) stop 'gas particles do not match'
  if (sinkslot-1 /= stot) stop 'sinks do not match'

  return
END SUBROUTINE read_data_dragon_form
