! NBODY_OUTPUT.F90
! D. A. Hubber - 28/01/2008
! Write all output from N-body integrator to screen and to files
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_output
  use definitions
  use sink_module
  use scaling_module
  use time_module
  use Nbody_module
  use filename_module
  implicit none

  character(len=11)  :: file_ext     ! filename extension for data output
  character(len=6)   :: file_numb    ! filename extension for data output
  character(len=256) :: out_data     ! Data snapshot filename
  character(len=256) :: out_debug    ! Debug data snapshot filename
  integer :: s                       ! Sink counter
  real(kind=DP) :: mtot              ! Total mass
  real(kind=DP) :: rcom(1:NDIM)      ! Position of centre of mass
  real(kind=DP) :: vcom(1:NDIM)      ! Velocity of centre of mass

  debug2("Writing N-body output to screen and file [nbody_output.F90]")
  debug_timing("NBODY_OUTPUT")

! Regular temp snapshots
! ----------------------------------------------------------------------------
  if (nsteps == ntempnext) then
     ntempnext = nsteps + ntempstep

     ! Calculate COM of system
     rcom(1:NDIM) = 0.0_DP
     vcom(1:NDIM) = 0.0_DP
     mtot = 0.0_DP
     
     do s=1,stot
        rcom(1:NDIM) = rcom(1:NDIM) + &
             &real(star(s)%m*star(s)%r(1:NDIM),DP)
        vcom(1:NDIM) = vcom(1:NDIM) + &
             &real(star(s)%m*star(s)%v(1:NDIM),DP)
        mtot = mtot + real(star(s)%m,DP)
     end do
          
     ! Normalise rcom and vcom
     rcom(1:NDIM) = rcom(1:NDIM) / mtot
     vcom(1:NDIM) = vcom(1:NDIM) / mtot

  end if
! ----------------------------------------------------------------------------


! Output periodic snapshot file
! ----------------------------------------------------------------------------
  if (time >= nextsnap .and. snapshot < 100000) then
     lastsnap = nextsnap
     nextsnap = nextsnap + snaptime
     snapshot = snapshot + 1

     call copy_stars_to_sinks

     write(file_numb,"(I5.5)") snapshot
     file_ext = "."//file_numb

     out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          &trim(adjustl(fileform_ext))//trim(adjustl(file_ext))
     call write_data(out_data,out_file_form) 

     ! Write name of snapshot file to log for potential restart
     open(1,file=restart_log,status="unknown",form="formatted")
     write(1,'(a)') out_data
     write(1,'(a)') out_file_form
     close(1)

!100  format (1X,'Step # : ',I7,4X,'Time : ',G12.6,1X,a)
!     write (6,100) nsteps,time*tscale,tunit

#if defined(DEBUG_OUTPUT_STAR_DATA)
     call write_star_data
#endif

     ! Calculate and write binary properties to file
#if defined(BINARY_STATS)
     call binary_search
#endif

  end if
! ----------------------------------------------------------------------------


! Write to screen diagnostic information 
! ----------------------------------------------------------------------------
#if defined(DEBUG_DAIGNOSTICS)
  if (nsteps == ndiagnext) then
     ndiagnext = ndiagnext + ndiagstep
     call nbody_diagnostics
#if defined(TIMING)
     call write_timing_stats
#endif
  end if
#endif
! ----------------------------------------------------------------------------


! Write regular time information to screen each step
! ----------------------------------------------------------------------------
#if defined(DEBUG1)
100 format (1X,'Step # : ',I8,4X,'Time : ',G14.8,1X,a,4X,'Snapshots : ',I4)
  write (6,100) nsteps,time*tscale,trim(tunit),snapshot
#endif


  return
END SUBROUTINE nbody_output
