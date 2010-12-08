! COM.F90
! D. A. Hubber - 24/2/2008
! Calculate the position and velocity of the centre of mass at the 
! beginning of the simulation.  If required, can transform the particle 
! data to the centre of mass frame (logical flag in params.dat).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE COM
  use particle_module
  use filename_module, only : restart
#if defined(SINKS)
  use sink_module
#endif
  implicit none

  integer :: p             ! Particle counter
#if defined(SINKS) 
  integer :: s             ! Sink counter
#endif

  debug2("Calculating the centre of mass of the system [COM.F90]")

  rcom0(1:NDIM) = 0.0_DP
  vcom0(1:VDIM) = 0.0_DP
  mtot0 = 0.0_DP

! Calculate position and velocity of centre of mass
  do p=1,ptot
     rcom0(1:NDIM) = rcom0(1:NDIM) + real(parray(MASS,p)*parray(1:NDIM,p),DP)
     vcom0(1:VDIM) = vcom0(1:VDIM) + real(parray(MASS,p)*v(1:VDIM,p),DP)
     mtot0 = mtot0 + real(parray(MASS,p),DP)
  end do

! Include sinks
#if defined(SINKS)
  do s=1,stot
     rcom0(1:NDIM) = rcom0(1:NDIM) + &
          &real(sink(s)%m*sink(s)%r(1:NDIM),DP)
     vcom0(1:VDIM) = vcom0(1:VDIM) + &
          &real(sink(s)%m*sink(s)%v(1:VDIM),DP)
     mtot0 = mtot0 + real(sink(s)%m,DP)
  end do
#endif

! Normalise rcom0 and vcom0
  rcom0(1:NDIM) = rcom0(1:NDIM) / mtot0
  vcom0(1:VDIM) = vcom0(1:VDIM) / mtot0
  rcom(1:NDIM) = rcom0(1:NDIM)
  vcom(1:VDIM) = vcom0(1:VDIM)


! Only convert to centre of mass frame if com_frame flag is on, if this is 
! not a restarted run, and periodic boundary conditions are not employed.  
! ----------------------------------------------------------------------------
#ifndef PERIODIC
  if (com_frame  .and. (.not. restart)) then
     write(6,*) "Changing to COM frame"
     write(6,*) "rcom0 :",rcom0(1:NDIM)
     write(6,*) "vcom0 :",vcom0(1:VDIM)
     do p=1,ptot
        parray(1:NDIM,p) = parray(1:NDIM,p) - real(rcom0(1:NDIM),PR)
        v(1:VDIM,p) = v(1:VDIM,p) - real(vcom0(1:VDIM),PR)
     end do

#if defined(SINKS)
     do s=1,stot
        sink(s)%r(1:NDIM) = &
             &sink(s)%r(1:NDIM) - real(rcom0(1:NDIM),PR)
        sink(s)%v(1:VDIM) = &
             &sink(s)%v(1:VDIM) - real(vcom0(1:VDIM),PR)
     end do
#endif

     rcom0(1:NDIM) = 0.0_DP
     vcom0(1:VDIM) = 0.0_DP

  end if
#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE COM
