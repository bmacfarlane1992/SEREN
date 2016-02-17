! READ_SINK_DATA_SYNC.F90
! D.Stamatellos -- 16/5/2012
! reads sinks file, reduces in size (in needed) and synchronises them

! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_sink_data_sync (s)
  use filename_module
  use particle_module
  use sink_module
  use Nbody_module
  use scaling_module
  use time_module
  implicit none

  integer, intent(in) :: s                      ! sink id

  character(len=8)  :: file_ext   ! filename extension for data output
  character(len=40) :: out_file   ! filename extension for data output
  character(len=8)  :: file_ext2  ! filename extension for data output
  character(len=40) :: out_file2  ! filename extension for data output

  logical :: ex                   ! Does file exist already?
  integer :: i                    ! Auxilary counter
  integer :: p                    ! Particle id
  integer :: idummy(1:2)          ! Dummy integer array
  real(kind=PR) :: rdummy(1:2)    ! Dummy real array
   integer :: idummy0             ! Dummy integer
  integer :: idummy1              ! Dummy integer
  real(kind=PR) :: rdummy0(1:27)  ! Dummy real array
  real(kind=PR) :: macc,tacc
  integer :: dn_sink=100          ! after how many steps to write sink info
  integer :: n_next		  ! next n to output sink info

  macc=0.0_PR
  tacc=0.0_PR

     ! Make filename from runid and sink number
     if (s >= 100) then
        write(file_ext,"(I3)") s
     else if (s>=10) then
        write(file_ext,"(I2)") s
     else
        write(file_ext,"(I1)") s
     end if

     ! open file with sink properties
        out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
             &".sink"//trim(adjustl(file_ext))

       out_file2 = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
             &".sink"//trim(adjustl(file_ext))//".sync"

       open(2,file=out_file2,status="unknown",form="formatted")

        ! Check if file exists
        inquire(file=out_file,exist=ex)
        !write(*,*) out_file


        ! Read through file to synchronise files with time of simulation.
        if (ex) then
           open(1,file=out_file,status="unknown",form="formatted",position="rewind")
! read tcreate, ncreate
#ifdef EPISODIC_ACCRETION
!           read(1,'(I8,27E18.10,I3)',end=30,err=30) idummy0,rdummy0(1:27),idummy1
           read(1,*,end=30,err=30) idummy0,rdummy0(1:27),idummy1
         !  write(2,'(I8,27E18.10,I3)') idummy0,rdummy0(1:27),idummy1
#else
!           read(1,'(I8,21E18.10)',end=30,err=30) idummy0,rdummy0(1:21)
           read(1,*,end=30,err=30) idummy0,rdummy0(1:21)
	!   write(2,'(I8,21E18.10)') idummy0,rdummy0(1:21)
#endif
           sink(s)%ncreate=idummy0
           sink(s)%tcreate=rdummy0(1)/tscale

           backspace (1)

           if (sink(s)%ncreate>sink(1)%ncreate) then
           	rdummy0(1:21)=0.0
	   	write(2,'(I8,21E18.10)') sink(1)%ncreate,rdummy0(1:21)
           	n_next=sink(1)%ncreate+dn_sink           
           	do while (n_next<sink(s)%ncreate)
	     		write(2,'(I8,21E18.10)') n_next,rdummy0(1:21)
             		n_next=n_next+dn_sink
           	enddo
           endif

           do
#ifdef EPISODIC_ACCRETION
!              read(1,'(I8,27E18.10,I3)',end=30,err=30) idummy0,rdummy0(1:27),idummy1
              read(1,*,end=30,err=30) idummy0,rdummy0(1:27),idummy1

!---- write sink info ------
	      if (idummy0>=n_next) then
		write(2,'(I8,27E18.10,I3)') idummy0,rdummy0(1:27),idummy1
                n_next=n_next+dn_sink
	      endif 	      

#else
!              read(1,'(I8,21E18.10)',end=30,err=30) idummy0,rdummy0(1:21)
              read(1,*,end=30,err=30) idummy0,rdummy0(1:21)
!---- write sink info ------
	      if (idummy0>=n_next) then
	        write(2,'(I8,21E18.10)') idummy0,rdummy0(1:21)
  	        n_next=n_next+dn_sink
	      endif 
#endif

 
           end do
30         backspace (1,err=40)
        else
     	write(6,*) "warning:: sink",s, "file is missing", out_file
     	GOTO 50
        end if
40   close(1)
     close(2)



!        sink(s)%r(1:NDIM)=rdummy0(2:4)/rscale
!        sink(s)%v(1:VDIM)=rdummy0(5:7)/vscale
!        sink(s)%m=rdummy0(8)/mscale
!        sink(s)%h=rdummy0(9)/rscale
!        sink(s)%radius=rdummy0(10)/rscale
!        sink(s)%a(1:VDIM)=rdummy0(11:13)/ascale
!        sink(s)%angmom(1:3)=rdummy0(14:16)/angmomscale
!        sink(s)%gpe=rdummy0(17)/Escale
!        sink(s)%dmdt=rdummy0(18)/dmdtscale
!        sink(s)%star_radius=rdummy0(19)/rscale
!        sink(s)%luminosity=rdummy0(20)/Lscale
!        sink(s)%temperature=rdummy0(21)

!#ifndef EPISODIC_ACCRETION
        sink(s)%Mstar=sink(s)%m
!#endif

!#ifdef EPISODIC_ACCRETION
!	sink(s)%Mstar=rdummy0(22)/mscale			
!	sink(s)%Mdisc=rdummy0(23)/mscale			
!	sink(s)%dmdt_star=rdummy0(24)/dmdtscale		
!       sink(s)%t_episode_start=rdummy0(25)/tscale		
!       sink(s)%t_episode_duration=rdummy0(26)/tscale		
!	sink(s)%dmdt_0=rdummy0(27)/dmdtscale			
!	sink(s)%accretion_flag=idummy1		
!#endif

50 continue
  return
END SUBROUTINE read_sink_data_sync
