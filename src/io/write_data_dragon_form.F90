! WRITE_DATA_DRAGON_FORM.F90
! C. P. Batty & D. A. Hubber - 12/12/2006
! Writes simultions snapshot to ASCII (i.e. formatted) file in DRAGON format.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_data_dragon_form(out_file)
  use particle_module
  use hydro_module
  use neighbour_module
  use scaling_module
  use time_module
  use type_module
  use sink_module
  implicit none

  character(len=*), intent(in) :: out_file  ! formatted DRAGON snapshot

  integer       :: p              ! counter to loop over particles
  integer       :: idata(1:20)    ! dummy integers
  real(kind=PR) :: rdata(1:50)    ! dummy reals
#if defined(SINKS)
  integer       :: s              ! Sink particle counter
#endif

  debug2("Writing data to formatted DRAGON file [write_data_dragon_form.F90]")

  open(1, file=out_file, status="unknown", form="formatted")
#if defined(DEBUG1)
  write(6,*) "Output file : ",trim(out_file),"  (DRAGON ASCII format)"
#endif

! Set header values and write to snapshot file
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

  do p=1,20
     write(1,*) idata(p)
  end do
  do p=1,50
     write(1,*) rdata(p)
  end do

! Positions
! ----------------------------------------------------------------------------
  do p=1,ptot
#if NDIM == 1
     write(1,*) parray(1,p)*real(rscale,PR)
#elif NDIM == 2
     write(1,*) parray(1,p)*real(rscale,PR), parray(2,p)*real(rscale,PR)
#else
     write(1,*) parray(1,p)*real(rscale,PR), parray(2,p)*real(rscale,PR), &
          &parray(3,p)*real(rscale,PR)
#endif
  end do
#if defined(SINKS)
  do s=1,stot
     write(1,*) sink(s)%r(1:NDIM)*real(rscale,PR)
  end do
#endif

! Velocities
! ----------------------------------------------------------------------------
  do p=1,ptot
#if NDIM == 1
     write(1,*) v(1,p)*real(vscale,PR)
#elif NDIM == 2
     write(1,*) v(1,p)*real(vscale,PR), v(2,p)*real(vscale,PR)
#else
     write(1,*) v(1,p)*real(vscale,PR), v(2,p)*real(vscale,PR), &
          &v(3,p)*real(vscale,PR)
#endif
  end do
#if defined(SINKS)
  do s=1,stot
     write(1,*) sink(s)%v(1:VDIM)*real(vscale,PR)
  end do
#endif

! Temperatures
! ----------------------------------------------------------------------------
  do p=1,ptot
     write(1,*) temp(p)
  end do
#if defined(SINKS)
  do s=1,stot
     write(1,*) 0.0_PR
  end do
#endif

! Smoothing lengths
! ----------------------------------------------------------------------------
  do p=1,ptot
     write(1,*) parray(SMOO,p)*real(rscale,PR)
  end do
#if defined(SINKS)
  do s=1,stot
     write(1,*) sink(s)%h*real(rscale,PR)
  end do
#endif

! Densities
! ----------------------------------------------------------------------------
  do p=1,ptot
     write(1,*) rho(p)*real(rhoscale,PR)
  end do
#if defined(SINKS)
  do s=1,stot
     write(1,*) rhosink*real(rhoscale,PR)
  end do
#endif

! Masses
! ----------------------------------------------------------------------------
  do p=1,ptot
     write(1,*) parray(MASS,p)*real(mscale,PR)
  end do
#if defined(SINKS)
  do s=1,stot
     write(1,*) sink(s)%m*real(mscale,PR)
  end do
#endif

! Particle types
! ----------------------------------------------------------------------------
  if (pboundary > 0) then
     do p=pboundarystart,pboundaryend
        write(1,*) 6
     end do
  end if
  if (picm > 0) then
     do p=picmstart,picmend
        write(1,*) 9
     end do
  end if
  if (pgas > 0) then
     do p=pgasstart,pgasend
        write(1,*) 1
     end do
  end if
  if (pcdm > 0) then
     do p=pcdmstart,pcdmend
        write(1,*) 10
     end do
  end if
#if defined(SINKS)
  if (stot > 0) then
     do s=1,stot
        write(1,*) -1
     end do
  end if
#endif

#ifdef PARTICLE_ID

! Particle ids (porig)
! ----------------------------------------------------------------------------
  do p=1,ptot
     write(1,*) porig(p)
  end do
#if defined(SINKS)
  do s=1,stot
     write(1,*) 0
  end do
#endif

#endif 

! Close file once finished
  close(1)



  return
END SUBROUTINE write_data_dragon_form
