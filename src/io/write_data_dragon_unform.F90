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
#if defined(SINKS)
  integer       :: s                          ! Sink particle counter
#endif

  debug2("Writing data to binary DRAGON file [write_data_dragon_unform.F90]")

#if defined(SINKS)
  allocate(idummy1(1:(ptot+stot)))
  allocate(rdummy1(1:(ptot+stot)))
  allocate(rdummy3(1:NDIM,1:(ptot+stot)))
#else
  allocate(idummy1(1:ptot))
  allocate(rdummy1(1:ptot))
  allocate(rdummy3(1:NDIM,1:ptot))
#endif

  open(1, file=out_file, status="unknown", form="unformatted")
#if defined(DEBUG1)
  write(6,*) "Output file : ",trim(out_file),"  (DRAGON binary format)"
#endif

  idata(1:20) = 0
  rdata(1:50) = 0.0_PR
  idata(1)    = ptot + stot
  idata(2)    = int(nsteps)
  idata(3)    = ptot
  idata(4)    = int(snapshot)
  idata(20)   = pgas_orig
  rdata(1)    = real(time*tscale,PR)
  rdata(2)    = real(lastsnap*tscale,PR)
  rdata(50)   = real(mgas_orig*mscale,PR)

  write(1) idata
  write(1) rdata

! Positions
! ----------------------------------------------------------------------------
  do p=1,ptot
     rdummy3(1:NDIM,p) = parray(1:NDIM,p)*real(rscale,PR)
  end do
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy3(1:NDIM,ptot+s) = sink(s)%r(1:NDIM)*real(rscale,PR)
     end do
  end if
#endif
  write(1) rdummy3

! Velocities
! ----------------------------------------------------------------------------
  do p=1,ptot
     rdummy3(1:VDIM,p) = v(1:VDIM,p)*real(vscale,PR)
  end do
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy3(1:VDIM,ptot+s) = sink(s)%v(1:VDIM)*real(vscale,PR)
     end do
  end if
#endif
  write(1) rdummy3

! Temperatures
! ----------------------------------------------------------------------------
  do p=1,ptot
     rdummy1(p) = temp(p)
  end do
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(ptot+s) = 0.0_PR
     end do
  end if
#endif
  write(1) rdummy1

! Smoothing lengths
! ----------------------------------------------------------------------------
  do p=1,ptot
     rdummy1(p) = parray(SMOO,p)*real(rscale,PR)
  end do
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(ptot+s) = sink(s)%h*real(rscale,PR)
     end do
  end if
#endif
  write(1) rdummy1

! Densities
! ----------------------------------------------------------------------------
  do p=1,ptot
     rdummy1(p) = rho(p)*real(rhoscale,PR)
  end do
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(ptot+s) = rhosink*real(rhoscale,PR)
     end do
  end if
#endif
  write(1) rdummy1

! Masses
! ----------------------------------------------------------------------------
  do p=1,ptot
     rdummy1(p) = parray(MASS,p)*real(mscale,PR)
  end do
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        rdummy1(ptot+s) = sink(s)%m*real(mscale,PR)
     end do
  end if
#endif
  write(1) rdummy1

! Particle types
! ----------------------------------------------------------------------------
  if (pboundary > 0) then
     do p=1,pboundary
        idummy1(p) = 6
     end do
  end if

  if (picm > 0) then
     do p=1,picm
        idummy1(p+pboundary) = 9
     end do
  end if

  if (pgas > 0) then
     do p=1,pgas
        idummy1(p+pboundary+picm) = 1
     end do
  end if

#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        idummy1(ptot+s) = -1
     end do
  end if
#endif

  write(1) idummy1

! Close file once finished and deallocate all temporary memory
  close(1)
  deallocate(rdummy3)
  deallocate(rdummy1)
  deallocate(idummy1)


  return
END SUBROUTINE write_data_dragon_unform
