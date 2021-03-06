! WRITE_SINK_DATA.F90
! D. A. Hubber - 30/5/2008
! Writes out sink information to files.  Writes an individual file for
! each sink formed in the simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_sink_data_sync
  use seren_sim_module
  use filename_module
  use sink_module
  use Nbody_module
  use scaling_module
  use time_module
  implicit none

  character(len=8)  :: file_ext   ! filename extension for data output
  character(len=256) :: out_file  ! filename extension for data output
  logical :: ex                   ! Does file exist already?
  integer :: s                    ! Sink counter
  integer :: idummy               ! Dummy integer
  real(kind=PR) :: rdummy(1:21)   ! Dummy real array

  if (stot == 0) return

  debug2("Writing sink information to files [write_sink_data.F90]")

! Loop over all sinks
! ----------------------------------------------------------------------------
  do s=1,stot

     ! Make filename from runid and sink number
     if (s >= 100) then
        write(file_ext,"(I3)") s
     else if (s>=10) then
        write(file_ext,"(I2)") s
     else
        write(file_ext,"(I1)") s
     end if

     ! Use ".sink" or ".star" depending on if we are in SPH or N-body phase.
     if (nbody_sph_sim .or. nbody_sim) then
        out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
             &".star"//trim(adjustl(file_ext))
     else
        out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
             &".sink"//trim(adjustl(file_ext))
     end if

     ! Check if file exists
     inquire(file=out_file,exist=ex)
     
     ! Read through file to synchronise files with time of simulation.
     if (ex) then
        open(1,file=out_file,status="unknown",&
             &form="formatted",position="rewind")
        do
#ifndef EPISODIC_ACCRETION
           read(1,'(I8,21E18.10)',end=10,err=10) idummy,rdummy(1:21)
#else
           read(1,'(I8,27E18.10,I3)',end=10,err=10) idummy,rdummy(1:21)
#endif
           if (rdummy(1) > time*tscale) exit
        end do
10      backspace (1,err=20)
     else
        open(1,file=out_file,status="unknown",form="formatted")
     end if

     ! Write all useful sink/star information to file
20   if (nbody_sph_sim .or. nbody_sim) then
        write(1,'(I8,21E18.10)') nsteps,&
             &time*tscale,&
             &star(s)%r(1:NDIM)*rscale,&
             &star(s)%v(1:VDIM)*vscale,&
             &star(s)%m*mscale,&
             &star(s)%h*rscale,&
             &star(s)%radius*rscale,&
             &star(s)%a(1:VDIM)*ascale,&
             &star(s)%angmom(1:3)*angmomscale,&
             &star(s)%gpe*Escale,&
             &star(s)%dmdt*dmdtscale,&
             &sink(s)%star_radius*rscale,&
             &sink(s)%luminosity*Lscale,&
             &sink(s)%temperature
     else
#ifndef EPISODIC_ACCRETION
        write(1,'(I8,21E18.10)') nsteps,&
             &time*tscale,&
             &sink(s)%r(1:NDIM)*rscale,&
             &sink(s)%v(1:VDIM)*vscale,&
             &sink(s)%m*mscale,&
             &sink(s)%h*rscale,&
             &sink(s)%radius*rscale,&
             &sink(s)%a(1:VDIM)*ascale,&
             &sink(s)%angmom(1:3)*angmomscale,&
             &sink(s)%gpe*Escale,&
             &sink(s)%dmdt*dmdtscale,&
             &sink(s)%star_radius*rscale,&
             &sink(s)%luminosity*Lscale,&
             &sink(s)%temperature
#else
        write(1,'(I8,27E18.10,I3)') nsteps,&             !1
             &time*tscale,&                             !2
             &sink(s)%r(1:NDIM)*rscale,&                !3/4/5
             &sink(s)%v(1:VDIM)*vscale,&                !6/7/8
             &sink(s)%m*mscale,&                        !9
             &sink(s)%h*rscale,&                        !10
             &sink(s)%radius*rscale,&                   !11
             &sink(s)%a(1:VDIM)*ascale,&                !12/13/14
             &sink(s)%angmom(1:3)*angmomscale,&         !15/16/17
             &sink(s)%gpe*Escale,&                      !18
             &sink(s)%dmdt*dmdtscale,&                  !19
             &sink(s)%star_radius*rscale,&              !20
             &sink(s)%luminosity*Lscale,&               !21
             &sink(s)%temperature,&                     !22
             &sink(s)%Mstar*mscale,&                    !23
             &sink(s)%Mdisc*mscale,&                    !24
             &sink(s)%dmdt_star*dmdtscale,&             !25
             &sink(s)%t_episode_start*tscale,&          !26
             &sink(s)%t_episode_duration*tscale,&       !27
             &sink(s)%dmdt_0*dmdtscale,&                !28
             &sink(s)%accretion_flag                    !29
#endif

     end if
     close(1)
  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE write_sink_data_sync
