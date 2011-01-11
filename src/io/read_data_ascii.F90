! READ_DATA_ASCII.F90
! D. A. Hubber - 17/06/2010
! Reads in initial conditions file in simple column ASCII format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data_ascii(in_file)
  use particle_module
  use hydro_module
  use type_module
  use sink_module
  use time_module
  implicit none

  character(len=*), intent(in) :: in_file   ! Name of ascii file

  character(len=20) :: auxdata         ! aux. data id variable
  character(len=20) :: data_id(1:100)  ! list of data identifiers 
  integer :: boundaryslot              ! point in list to boundary particle
  integer :: cdmslot                   ! point in list to cdm particle
  integer :: gasslot                   ! point in list to insert particle
  integer :: i                         ! aux. loop variable
  integer :: icmslot                   ! point in list to icm particle
  integer :: j                         ! aux. loop variable
  integer :: ndata                     ! no. of data columns
  integer :: p                         ! counter to loop over particles
  integer :: psplit                    ! counter for particle splitting
  integer :: pdead                     ! counter for accreted particles
  integer :: ptype                     ! particle type
  integer :: s                         ! sink id
  integer :: sinkslot                  ! point in list sink particle
  real(kind=PR) :: pdata(1:100)        ! aux. particle data array

  debug1("Reading in data file : "//trim(in_file)//" [read_data_ascii.F90]")

! Initialise counters
  ndata     = 0
  pboundary = 0
  psplit    = 0
  pdead     = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 0

! First, read-in ascii column data
  open(1, file="asciicolumns.dat", status="old", form="formatted")
  do
     read(1,*,err=5,end=5) auxdata
     ndata = ndata + 1
     data_id(ndata) = auxdata
  end do
5 close(5)

! Check there are a sufficient no. of data columns, or ptype is not 1st column
  if (ndata <= 1) stop "Invalid no. of data columns"
  if (data_id(1) /= "ptype") stop "First column is not particle type id"
  open(1, file=in_file, status="old", form="formatted")

! First, count the total no. of particles and the no. of each particle type
! ----------------------------------------------------------------------------
  do
     read(1,*,err=10,end=10) ptype,pdata(2:ndata)
     if (ptype == -1) stot = stot + 1
     if (ptype == 0)  pdead = pdead + 1
     if (ptype == 1)  pgas = pgas + 1
     if (ptype == 4)  psplit = psplit + 1
     if (ptype == 6)  pboundary = pboundary + 1
     if (ptype == 9)  picm = picm + 1
     if (ptype == 10) pcdm = pcdm + 1
     if (ptype >= 1)  ptot = ptot + 1
  end do
! ----------------------------------------------------------------------------

10 rewind(1)

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

  call allocate_memory

! Initialise aux. id variables
  boundaryslot = 0
  icmslot      = pboundary 
  gasslot      = pboundary + picm
  cdmslot      = pboundary + picm + pgas
  sinkslot     = 0

! Now read all particle data into main arrays
! ============================================================================
  do j=1,ptot+stot
     read(1,*,err=20,end=20) ptype,pdata(2:ndata)

     ! Determine array location of particle depending on type
     if (ptype == 6) then
        boundaryslot = boundaryslot + 1
        p = boundaryslot
     else if (ptype == 9) then
        icmslot = icmslot + 1
        p = icmslot
     else if (ptype == 1) then
        gasslot = gasslot + 1
        p = gasslot
     else if (ptype == 10) then
        cdmslot = cdmslot + 1
        p = cdmslot
     else if (ptype == -1) then
        sinkslot = sinkslot + 1
        s = sinkslot
        p = -s
     else
        write(6,*) "Unrecognised particle type in read_data_ascii.F90 :",ptype
        stop
     end if

     ! If particle is an SPH particle
     ! -----------------------------------------------------------------------
     if (p > 0) then
        do i=2,ndata
           if (data_id(i) == "x") then
              parray(1,p) = pdata(i)
           else if (data_id(i) == "y" .and. NDIM > 1) then
              parray(2,p) = pdata(i)
           else if (data_id(i) == "z" .and. NDIM == 3) then
              parray(3,p) = pdata(i)
           else if (data_id(i) == "m") then
              parray(MASS,p) = pdata(i)
           else if (data_id(i) == "h") then
              parray(SMOO,p) = pdata(i)
           else if (data_id(i) == "vx") then
              v(1,p) = pdata(i)
           else if (data_id(i) == "vy" .and. VDIM > 1) then
              v(2,p) = pdata(i)
           else if (data_id(i) == "vz" .and. VDIM == 3) then
              v(3,p) = pdata(i)
           else if (data_id(i) == "rho") then
              rho(p) = pdata(i)
#if defined(INTERNAL_ENERGY)
           else if (data_id(i) == "u") then
              u(p) = pdata(i)
#endif
           end if
        end do


     ! Else, if a sink particle
     ! -----------------------------------------------------------------------
#if defined(SINKS)
     else if (p < 0) then
        do i=2,ndata
           if (data_id(i) == "x") then
              sink(s)%r(1) = pdata(i)
           else if (data_id(i) == "y" .and. NDIM > 1) then
              sink(s)%r(2) = pdata(i)
           else if (data_id(i) == "z" .and. NDIM == 3) then
              sink(s)%r(3) = pdata(i)
           else if (data_id(i) == "m") then
              sink(s)%m = pdata(i)
           else if (data_id(i) == "h") then
              sink(s)%h = pdata(i)
           else if (data_id(i) == "vx") then
              sink(s)%v(1) = pdata(i)
           else if (data_id(i) == "vy" .and. VDIM > 1) then
              sink(s)%v(2) = pdata(i)
           else if (data_id(i) == "vz" .and. VDIM == 3) then
              sink(s)%v(3) = pdata(i)
           end if
        end do
        sink(s)%accrete = .true.
        sink(s)%static  = .false.
        sink(s)%radius  = KERNRANGE*sink(s)%h
#endif
     end if
     ! -----------------------------------------------------------------------

  end do
! ============================================================================

20 close(1)

  return
END SUBROUTINE read_data_ascii
