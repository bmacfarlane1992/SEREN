! SPH_OUTPUT.F90
! C. P. Batty & D. A. Hubber - 23/8/2007
! Calls routines to generate regular and temporary snapshot files, and 
! to output diagnostic and time information to screen.  Also calculates 
! some diagnostic info for certain debug options.
! ============================================================================ 

#include "macros.h"

! ============================================================================
SUBROUTINE sph_output
  use interface_module, only : track_particle,write_data,write_data_debug
  use time_module
  use filename_module
  use scaling_module
  use particle_module
  implicit none
  
  character(len=11)  :: file_ext    ! filename extension for data output
  character(len=6)   :: file_numb   ! filename extension for data output
  character(len=10)  :: file_numb2  ! filename extension for data output
  character(len=256) :: out_data    ! Data snapshot filename
  character(len=256) :: out_debug   ! Debug data snapshot filename
#if defined(DEBUG_TRACK_PARTICLE)
  integer :: p                      ! Particle counter
  integer :: paux                   ! Aux. particle id
  real(kind=PR) :: rorigin(1:NDIM)  ! Position of origin for debug output
#endif

  debug2("Writing output to screen [output.F90]")
  debug_timing("OUTPUT")


! Output temporary snapshot file (alternating between 2 temp files)
! ----------------------------------------------------------------------------
  if (nsteps == ntempnext) then
     if (ntemp == 1) then
        out_data = out_temp1
        out_debug =trim(adjustl(run_dir))//trim(adjustl(run_id))//".debug.tmp1"
        ntemp = 2
     else if (ntemp == 2) then
        out_data = out_temp2
        out_debug =trim(adjustl(run_dir))//trim(adjustl(run_id))//".debug.tmp2"
        ntemp = 1
     endif
     ntempnext = nsteps + ntempstep
     call write_data(out_data,out_file_form)
#if defined(DEBUG_PLOT_DATA)
     call write_data_debug(out_debug,rzero(1:NDIM))
#endif
#if defined(DEBUG_RSPH_OUTPUT)
     call write_data_RSPH
#endif

     ! Write name of temp file to log
     open(1,file=restart_log,status="unknown",form="formatted")
     write(1,'(a)') trim(adjustl(out_data))
     write(1,'(a)') trim(adjustl(out_file_form))
     close(1)
  end if


! Output snapshot file every nsnapstep steps
! ----------------------------------------------------------------------------
  if (nsteps == nsnapnext .and. nsteps < 100000000) then
     nsnapnext = nsteps + nsnapstep

     write(file_numb2,"(I8.8)") nsteps
     file_ext = ".n"//trim(adjustl(file_numb2))

     out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          &trim(adjustl(fileform_ext))//trim(adjustl(file_ext))
     call write_data(out_data,out_file_form) 
#if defined(DEBUG_PLOT_DATA)
     out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          ".debug"//trim(adjustl(file_ext))
     call write_data_debug(out_debug,rzero(1:NDIM))
#endif
#if defined(DEBUG_RSPH_OUTPUT)
     call write_data_RSPH
#endif

     ! Write name of snapshot file to log for potential restart
     open(1,file=restart_log,status="unknown",form="formatted")
     write(1,'(a)') out_data
     write(1,'(a)') out_file_form
     close(1)
  end if


! Output periodic snapshot file
! ----------------------------------------------------------------------------
  if (time >= nextsnap .and. snapshot < 100000) then
     lastsnap = nextsnap
     nextsnap = nextsnap + snaptime
     snapshot = snapshot + 1

     write(file_numb,"(I5.5)") snapshot
     file_ext = "."//file_numb

     out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          &trim(adjustl(fileform_ext))//trim(adjustl(file_ext))
     call write_data(out_data,out_file_form) 
#if defined(DEBUG_PLOT_DATA)
     out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))&
          &//".debug"//trim(adjustl(file_ext))
     call write_data_debug(out_debug,rzero(1:NDIM))
#endif
#if defined(DEBUG_GRID_RESULTS)
     out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))&
          &//".grid"//trim(adjustl(file_ext))
     call write_data_grid_results(out_debug)
#endif
#if defined(DEBUG_RSPH_OUTPUT)
     call write_data_RSPH
#endif

     ! Write name of snapshot file to log for potential restart
     open(1,file=restart_log,status="unknown",form="formatted")
     write(1,'(a)') out_data
     write(1,'(a)') out_file_form
     close(1)
  end if


#if defined(SPH_SPECIFIC_OUTPUT)
 call sph_specific_output
#endif


! Write to screen diagnostic information 
! ----------------------------------------------------------------------------
#if defined(DEBUG_DIAGNOSTICS)
  if (nsteps == ndiagnext) then
     ndiagnext = ndiagnext + ndiagstep
     call diagnostics
#if defined(TIMING)
     call write_timing_stats
#endif
  end if
#endif


! Write sink information to file
! ----------------------------------------------------------------------------
#if defined(SINKS)
  if (nsteps == nsinknext) then
     nsinknext = nsinknext + nsinkstep
     call write_sink_data
  end if
#endif

! Write specific data output (e.g. at certain densities or temps)
! ----------------------------------------------------------------------------
#if defined(SPH_SPECIFIC_OUTPUT)
 call sph_specific_output
#endif



! Write diagnostic data for various test options
! ----------------------------------------------------------------------------
#if defined(RAD_WS) && defined(DEBUG_RAD)
  call write_rad_ws_test_data
#endif
#if defined(BINARY_TEST)
  call write_binary_data
#endif


! Write data to file for particle tracking
! ----------------------------------------------------------------------------
#if defined(DEBUG_TRACK_PARTICLE)
  paux = -1
  rorigin(1:NDIM) = 0.0_PR
  do p=1,ptot
     if (porig(p) == ptrack) paux = p
  end do
  if (paux /= -1) call track_particle(paux,rorigin(1:NDIM))
#endif


! Write regular time information to screen each step
! ----------------------------------------------------------------------------
#if defined(DEBUG1)
100 format (1X,'Step # : ',I8,4X,'Time : ',G14.8,1X,a,4X,'Snapshots : ',I4)
  write (6,100) nsteps,time*tscale,trim(tunit),snapshot
#endif

  return
END SUBROUTINE sph_output
