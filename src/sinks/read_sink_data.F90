! READ_SINK_DATA.F90
! D.Stamatellos -- 20/1/2011
! reconstructs the accretion history from sink files 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_sink_data (s)
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
  character(len=8)  :: file_ext2   ! filename extension for data output
  character(len=40) :: out_file2   ! filename extension for data output

  logical :: ex                   ! Does file exist already?
  integer :: i                    ! Auxilary counter
  integer :: p                    ! Particle id
  integer :: idummy(1:2)          ! Dummy integer array
  real(kind=PR) :: rdummy(1:2)    ! Dummy real array
   integer :: idummy0               ! Dummy integer
  integer :: idummy1          ! Dummy integer
  real(kind=PR) :: rdummy0(1:27)   ! Dummy real array
  real(kind=PR) :: macc,tacc

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

        ! Check if file exists
        inquire(file=out_file,exist=ex)
         write(*,*) out_file

        ! Read through file to synchronise files with time of simulation.
        if (ex) then
           open(1,file=out_file,status="unknown",form="formatted",position="rewind")
! read tcreate, ncreate
#ifdef EPISODIC_ACCRETION
!           read(1,'(I8,27E18.10,I3)',end=30,err=30) idummy0,rdummy0(1:27),idummy1
           read(1,*,end=30,err=30) idummy0,rdummy0(1:27),idummy1
#else
!           read(1,'(I8,21E18.10)',end=30,err=30) idummy0,rdummy0(1:21)
           read(1,*,end=30,err=30) idummy0,rdummy0(1:21)
#endif
           sink(s)%ncreate=idummy0
           sink(s)%tcreate=rdummy0(1)/tscale
           backspace (1)

           do
#ifdef EPISODIC_ACCRETION
!              read(1,'(I8,27E18.10,I3)',end=30,err=30) idummy0,rdummy0(1:27),idummy1
              read(1,*,end=30,err=30) idummy0,rdummy0(1:27),idummy1
#else
!              read(1,'(I8,21E18.10)',end=30,err=30) idummy0,rdummy0(1:21)
              read(1,*,end=30,err=30) idummy0,rdummy0(1:21)
#endif
              if (rdummy0(1) > time*tscale) exit
 
! reconstruct (approximate) accretion array

! make nsinkstep spaces at the end of the accretion array 
                do i=DMDT_RANGE,nsinkstep,-1
!               		sink(s)%macc(i) = sink(s)%macc(i-nsinkstep)
               		sink(s)%tacc(i) = sink(s)%tacc(i-nsinkstep)

           	end do

! fill the last nsinkstep spaces by assuming equal accretion rate per timestep
               do i=nsinkstep,1,-1
!            	sink(s)%macc(i) = (rdummy0(8)-macc)/mscale/nsinkstep
            	sink(s)%tacc(i) = (tacc+ (nsinkstep-i+1)*(rdummy0(1)-tacc)/nsinkstep )/tscale
		enddo

                macc=rdummy0(8)
		tacc=rdummy0(1)
           end do
30         backspace (1,err=40)
        else
     	write(6,*) "warning:: sink",s, "file is missing", out_file
     	GOTO 50
        end if
40   close(1)
        sink(s)%r(1:NDIM)=rdummy0(2:4)/rscale
        sink(s)%v(1:VDIM)=rdummy0(5:7)/vscale
        sink(s)%m=rdummy0(8)/mscale
        sink(s)%h=rdummy0(9)/rscale
        sink(s)%radius=rdummy0(10)/rscale
        sink(s)%a(1:VDIM)=rdummy0(11:13)/ascale
        sink(s)%angmom(1:3)=rdummy0(14:16)/angmomscale
        sink(s)%gpe=rdummy0(17)/Escale
        sink(s)%dmdt=rdummy0(18)/dmdtscale
        sink(s)%star_radius=rdummy0(19)/rscale
        sink(s)%luminosity=rdummy0(20)/Lscale
        sink(s)%temperature=rdummy0(21)

#ifndef EPISODIC_ACCRETION
        sink(s)%Mstar=sink(s)%m
#endif
! reconstruct (approximate) accretion array from the accretion rate at the last step
        do i=DMDT_RANGE-1,1,-1
           sink(s)%macc(i) = sink(s)%dmdt*(sink(s)%tacc(i)-sink(s)%tacc(i+1))
       enddo
           sink(s)%macc(DMDT_RANGE)=sink(s)%macc(DMDT_RANGE-1)

#ifdef EPISODIC_ACCRETION
	sink(s)%Mstar=rdummy0(22)/mscale			
	sink(s)%Mdisc=rdummy0(23)/mscale			
	sink(s)%dmdt_star=rdummy0(24)/dmdtscale		
        sink(s)%t_episode_start=rdummy0(25)/tscale		
        sink(s)%t_episode_duration=rdummy0(26)/tscale		
	sink(s)%dmdt_0=rdummy0(27)/dmdtscale			
	sink(s)%accretion_flag=idummy1		
#endif
50 continue

  return
END SUBROUTINE read_sink_data
