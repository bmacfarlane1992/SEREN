! SPH_SPECIFIC_OUTPUT.F90
! user defined output
! ============================================================================ 

#include "macros.h"

! ============================================================================
SUBROUTINE sph_specific_output
  use interface_module, only : track_particle,write_data,write_data_debug
  use time_module
  use filename_module
  use scaling_module
  use particle_module
  use hydro_module, only:rho,temp
  use type_module, only: pgas
  implicit none
  
  integer	     :: p
  integer	     :: i,j
  logical 	     :: ex          ! Does file exist already?
  character(len=11)  :: file_ext    ! filename extension for data output
  character(len=6)   :: file_numb   ! filename extension for data output
  character(len=10)  :: file_numb2  ! filename extension for data output
  character(len=256) :: out_data    ! Data snapshot filename
  character(len=256) :: out_debug   ! Debug data snapshot filename

  real(kind=DP) 	:: density(100)  
  real(kind=DP) 	:: temperature(100)
  real(kind=DP) 	:: maxdens 
  real(kind=DP) 	:: maxtemp
  integer		:: imaxdens
  integer		:: imaxtemp
  real			:: rcounter
  integer		:: tempflag(100) ! temperatures to output data
  integer 		:: densflag(100) ! densities to output data

! find max dens, max temp particle

  maxdens=0.0
  maxtemp=0.0
  imaxdens=0
  imaxtemp=0

  do p=1, pgas

        if (rho(p)>maxdens) THEN
          maxdens= rho(p)
          imaxdens=p
        endif

        if (temp(p)>maxtemp) THEN
          maxtemp=temp(p)
          imaxtemp=p
        endif

  end do

! fix tempflag, densflaf (checks if file exists already)
  do i=1,100,1
  densflag(i)=0
  tempflag(i)=0
  enddo
!  initdensflag=.TRUE.
!  inittempflag=.TRUE.

#if defined(SPH_OUTPUT_DENS)
 j=1
 rcounter=-17.0
  do while (rcounter<=0) 
   density(j)=10**(1.D0*rcounter)

   write(file_numb,"(F5.1)") (log10(density(j)))
 		write(file_numb2,"(I2.2)") j
     		out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          		&trim(adjustl(fileform_ext))//".dens."//trim(adjustl(file_numb2))//'.'//&
			&trim(adjustl(file_numb))
   ! Check if file exists
   inquire(file=out_data,exist=ex)
   if (ex) densflag(j)=1
   rcounter=rcounter+0.2
   j=j+1
 end do

  j=1
  rcounter=-17.0
  do while (rcounter<=0)

   density(j)=10**(1.D0*rcounter)

        if (maxdens*real(rhoscale,PR)>=density(j) .and. densflag(j)==0) then
                densflag(j)=1
                write(file_numb,"(F5.1)") (log10(density(j)))
                write(file_numb2,"(I2.2)") j
                out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
                        &trim(adjustl(fileform_ext))//".dens."//trim(adjustl(file_numb2))//'.'//&
                        &trim(adjustl(file_numb))
                call write_data(out_data,out_file_form)

                exit
        endif

        if (maxdens*real(rhoscale,PR)<density(j).and. densflag(j)==0) exit
   rcounter=rcounter+0.2
   j=j+1
  end do
#endif

#if defined(SPH_OUTPUT_TEMP)
 j=1
 do i=100, 10000, 100
   temperature(j)=i
   write(file_numb,"(I5.5)") INT(temperature(j))
     out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          &trim(adjustl(fileform_ext))//".temp."//trim(adjustl(file_numb))
   ! Check if file exists
   inquire(file=out_data,exist=ex)
   if (ex) tempflag(j)=1
   j=j+1
 end do


  j=1
  do i=100, 10000, 100
   temperature(j)=i*1.D0
  	if (maxtemp>=temperature(j) .and. tempflag(j)==0) then
		tempflag(j)=1
    		write(file_numb,"(I5.5)") INT(temperature(j))
     		out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          		&trim(adjustl(fileform_ext))//".temp."//trim(adjustl(file_numb))
     		call write_data(out_data,out_file_form) 

		exit
	endif
        if (maxtemp<temperature(j) .and. tempflag(j)==0) exit
   j=j+1
  end do
#endif

  return
END SUBROUTINE sph_specific_output
