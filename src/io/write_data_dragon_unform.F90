! WRITE_DATA_DRAGON_UNFORM.F90
! C. P. Batty & D. A. Hubber - 12/12/2006
! Writes simulation snapshot to binary (i.e. unformatted) file in 
! DRAGON format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_data_dragon_unform(out_file)
  use particle_module
  use hydro_module
  use neighbour_module
  use scaling_module
  use time_module
  use type_module
  use sink_module
  implicit none

  character(len=*), intent(in) :: out_file    ! formatted DRAGON snapshot

  integer       :: p                          ! particle counter
  integer       :: idata(1:20)                ! dummy integers
  integer, allocatable :: idummy1(:)          ! dummy integer array
  real(kind=PR) :: rdata(1:50)                ! dummy reals
  real(kind=PR), allocatable :: rdummy1(:)    ! dummy real array
  real(kind=PR), allocatable :: rdummy3(:,:)  ! dummy real vector array
!
!
  real(kind=DP), parameter :: trestr = 1000.							    !!! Set restrictive values of density and 
  real(kind=DP), parameter :: rhorestr = 1.E-13							    !!! temperature (in cgs)
  integer	:: i, trho_count								    !!! Counter for particles in T-rho plane
!
!
#if defined(SINKS)
  integer       :: s                          ! Sink particle counter
#endif

  debug2("Writing data to binary DRAGON file [write_data_dragon_unform.F90]")


  open(1, file=out_file, status="unknown", form="unformatted")
#if defined(DEBUG1)
  write(6,*) "Output file : ",trim(out_file),"  (DRAGON binary format)"
#endif

  idata(1:20) = 0
  rdata(1:50) = 0.0_PR
  idata(2)    = int(nsteps)
  idata(4)    = int(snapshot)
  idata(20)   = pgas_orig
  rdata(1)    = real(time*tscale,PR)
  rdata(2)    = real(lastsnap*tscale,PR)
  rdata(50)   = real(mgas_orig*mscale,PR)


! Count particles in T-rho plane
! ----------------------------------------------------------------------------

  trho_count = 0
  do p=1, ptot
	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
		trho_count = trho_count + 1	
	endif
  enddo
!
  idata(1)    = trho_count + stot
  idata(3)    = trho_count
!
  write(1) idata
  write(1) rdata
!
#if defined(SINKS)
  allocate(idummy1(1:(trho_count+stot)))
  allocate(rdummy1(1:(trho_count+stot)))
  allocate(rdummy3(1:NDIM,1:(trho_count+stot)))
#else
  allocate(idummy1(1:trho_count))
  allocate(rdummy1(1:trho_count))
  allocate(rdummy3(1:NDIM,1:trho_count))
#endif

! Positions
! ----------------------------------------------------------------------------

  i = 1
  do p=1,ptot

	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
		rdummy3(1:NDIM,i) = parray(1:NDIM,p)*real(rscale,PR)
		i = i + 1
	endif

  end do

#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy3(1:NDIM,trho_count+s) = sink(s)%r(1:NDIM)*real(rscale,PR)
     end do
  end if
#endif
  write(1) rdummy3

! Velocities
! ----------------------------------------------------------------------------

  i = 1
  do p=1,ptot

	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
		rdummy3(1:VDIM,i) = v(1:VDIM,p)*real(vscale,PR)
		i = i + 1
	endif

  end do

#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy3(1:VDIM,trho_count+s) = sink(s)%v(1:VDIM)*real(vscale,PR)
     end do
  end if
#endif
  write(1) rdummy3

! Temperatures
! ----------------------------------------------------------------------------

  i = 1
  do p=1,ptot

	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
		rdummy1(i) = temp(p)
		i = i + 1
	endif

  end do
  
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(trho_count+s) = 0.0_PR
     end do
  end if
#endif
  write(1) rdummy1

! Smoothing lengths
! ----------------------------------------------------------------------------

  i = 1
  do p=1,ptot

	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
		rdummy1(i) = parray(SMOO,p)*real(rscale,PR)
		i = i + 1
	endif

  end do

#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(trho_count+s) = sink(s)%h*real(rscale,PR)
     end do
  end if
#endif
  write(1) rdummy1

! Densities
! ----------------------------------------------------------------------------

  i = 1
  do p=1,ptot

	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
		rdummy1(i) = rho(p)*real(rhoscale,PR)
		i = i + 1
	endif

  end do

#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(trho_count+s) = rhosink*real(rhoscale,PR)
     end do
  end if
#endif
  write(1) rdummy1

! Masses
! ----------------------------------------------------------------------------

  i = 1
  do p=1,ptot

	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
		rdummy1(i) = parray(MASS,p)*real(mscale,PR)
		i = i + 1
	endif

  end do

#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(trho_count+s) = sink(s)%m*real(mscale,PR)
     end do
  end if
#endif
  write(1) rdummy1

! Particle types
! ----------------------------------------------------------------------------

  i = 1
  if (pboundary > 0) then
     do p=1,pboundary

		if( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
	        	idummy1(i) = 6
			i = i + 1
		endif

     end do
  end if

  i = 1
  if (picm > 0) then
     do p=1,picm

		if( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
			idummy1(i+pboundary) = 9
		i = i + 1
		endif
     end do
  end if

  i = 1
  if (pgas > 0) then
     do p=1,pgas

		if( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
			idummy1(i+pboundary+picm) = 1
			i = i + 1
		endif
     end do
  end if

#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        idummy1(trho_count+s) = -1
     end do
  end if
#endif

  write(1) idummy1


#if defined(PARTICLE_ID)
! Particle ids
! ----------------------------------------------------------------------------

  i = 1
  do p=1,ptot

	if( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then	
		idummy1(i) = porig(p)
		i = i + 1
	endif

! if (porig(p)==4768) write(*,*)  real(time*tscale,PR), parray(1,p)*real(rscale,PR)
  end do
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        idummy1(trho_count+s) = 0
     end do
  end if
#endif
  write(1) idummy1

#endif



 
! Close file once finished and deallocate all temporary memory
  close(1)
  deallocate(rdummy3)
  deallocate(rdummy1)
  deallocate(idummy1)


  return
END SUBROUTINE write_data_dragon_unform
