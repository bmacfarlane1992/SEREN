! WRITE_TIMING_STATS.F90
! D. A. Hubber - 14/2/2008
! Writes out all timing statistics to file for quick analysis of code profile.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_timing_stats
  use definitions
  use timing_module
  use filename_module, only : run_dir,run_id
  implicit none

  integer :: i                    ! Aux. loop counter
  real(kind=DP) :: ipercentage    ! percentage of total integer time
  real(kind=DP) :: rpercentage    ! percentage of total real time
  character(len=256) :: filename  ! filename extension for data output

  debug2("Write timing information to file [write_timing_stats.F90]")


! Loop over all blocks and print statistics
! ----------------------------------------------------------------------------
  if (mark_tot > 0) then

     filename = trim(adjustl(run_dir))//trim(adjustl(run_id))//".timing"
     open(unit=1,file=filename,status="unknown")

10   format(1X,a20,2X,i10,2X,1F8.4,2X,1G14.6,2X,1F8.4)

     write(1,*) "SEREN timing statistics"
     write(1,*) "----------------------------------------&
          &----------------------------"
     write(1,*) "Marker id                  itime        &
          &i%       rtime            r%"
     write(1,*) "----------------------------------------&
          &----------------------------"
     
     do i=1,mark_tot
        ipercentage = 100.0_DP * real(iblock(i),DP) / real(itime,DP)
        rpercentage = 100.0_DP * rblock(i) / rtime        
        write(1,10) marker_id(i),iblock(i),ipercentage,rblock(i),rpercentage
     end do
     write(1,*) "----------------------------------------&
          &----------------------------"
     write(1,10) 'TOTAL               ',itime,100.0,rtime,100.0
     write(1,*) "----------------------------------------&
          &----------------------------"
     write(1,*)

     write(1,*) "SEREN computational rate statistics"
     write(1,*) "----------------------------------------&
          &----------------------------"
     do i=1,mark_tot
#if defined(GRAVITY)
        if (trim(adjustl(marker_id(i))) == "GRAVITY_FORCES") write(1,*) &
             &"Rate of gravity calcs : ",real(ngravcomp,DP)/rblock(i)
#endif
        if (trim(adjustl(marker_id(i))) == "HYDRO_FORCES") write(1,*) &
             &"Rate of hydro calcs   : ",real(nhydrocomp,DP)/rblock(i)
     end do
     write(1,*) "----------------------------------------&
          &----------------------------"
     
     close(1)

  end if
! ----------------------------------------------------------------------------


  return
END SUBROUTINE write_timing_stats
