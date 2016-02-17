! analyse_disc.F90
! D. Stamatellos 
! Reads in a file and assuming a disc centred around the first sink 
!	(i)   centres around the first sink and rotates disc so that xy is the disc midplane "infile.trans"
!	(ii)  calculates disc properties (Q, T, Sigma, Omega, Mass vs radius)
!	(iii) collective files with what happens in a given radius vs time (appends on file)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE analyse_central_region(in_file)
  use type_module
  use particle_module
  use hydro_module
  use sink_module
  use time_module, only : time,nsteps,snapshot
  use scaling_module
  use filename_module, only :run_id,out_file_form,run_dir, in_file_form

  implicit none
  character(len=40) , intent(in) 	:: in_file
  character(len=40)                	:: out_file  ! snapshot name
  character(len=8)  			:: file_ext   ! filename extension for data output
  character(len=6)   :: file_numb   ! filename extension for data output
 
  integer	:: i				! radial bin counter
  integer	:: p				! particle counter
  integer	:: s				! sink counter
 
  integer  	:: NLIST			! list of NLIST densest particles
  real(kind=PR) :: rdensest(1:NDIM)            ! position of denset particle 
  real(kind=PR) :: vdensest(1:NDIM)            ! vel of denset particle 
  integer, allocatable :: imaxdens(:) 		! Integer id of NLIST densest particles
  real(kind=PR), allocatable :: maxdens(:)      ! Density of NLIST densest particles

!constants  (cgs)
  real(kind=DP),parameter :: mH = 1.660539D-24  ! Hydrogen mass
  real(kind=DP),parameter :: Gconst = 6.67428D-8       ! Grav. constant
  real(kind=DP),parameter :: kboltzmann = 1.38065D-16  ! Boltzmann constant 
  real(kind=PR) :: r1,r2,dr,ytmp,ztmp,vytmp,vztmp,sintheta,costheta! auxiliary variables

  NLIST=50

 allocate(maxdens(1:NLIST))
 allocate(imaxdens(1:NLIST))

 maxdens(1:NLIST)=0.0
 imaxdens(1:NLIST)=0

  do p=1, pgas

 	if (rho(p)>maxdens(1)) THEN  

          maxdens(1)= rho(p)

          imaxdens(1)=p

        do i=1,NLIST-1 
        	if (maxdens(i)>maxdens(i+1)) then
		maxdens(i)=maxdens(i+1)
                imaxdens(i)=imaxdens(i+1)
                maxdens(i+1)=rho(p)
                imaxdens(i+1)=p
                endif
        enddo
	endif
  enddo


  p=imaxdens(NLIST)

  rdensest(1:NDIM)=parray(1:NDIM,p)
  vdensest(1:NDIM)=v(1:NDIM,p)

#ifdef ANALYSE_CENTRE_AROUND_DENSEST
! centre around densest particle
     do p=1,ptot
        parray(1:NDIM,p) = parray(1:NDIM,p) - rdensest(1:NDIM)
        v(1:VDIM,p) = v(1:VDIM,p) - vdensest(1:NDIM)
     end do
#if defined(SINKS)
     do s=1,stot
        sink(s)%r(1:NDIM) = &
             &sink(s)%r(1:NDIM)  - rdensest(1:NDIM)
        sink(s)%v(1:VDIM) = &
             &sink(s)%v(1:VDIM) - vdensest(1:NDIM)
     end do
#endif

     out_file = trim(adjustl(in_file))//"."//trim(adjustl(out_file_form))//".cnt"

!  write file with sink at the origin of the system and disc on xy-plane
write(*,*) out_file, out_file_form
  call write_data(out_file,out_file_form)

#endif

#ifdef ANALYSE_CENTRAL_REGION_PROPERTIES


  open(1,file="TD.dat",status="unknown",&
             &form="formatted",position="append")

 write(1,'(E18.10,1X,I8,7(1X,E15.7),1X,F7.1,2X,18A)') &
		& time*tscale, & 				!1
             	& nsteps,& 	            			!2
             	& rdensest(1:NDIM)*real(rscale,PR),&		!3/4/5
		& vdensest(1:NDIM)*real(vscale,PR), &		!6/7/8
		& rho(imaxdens(NLIST))*real(rhoscale,PR),&	!9
 		& temp(imaxdens(NLIST)),&			!10
		& trim(adjustl(in_file))

 close(1)
#endif


 deallocate(maxdens)
 deallocate(imaxdens)

  return
END SUBROUTINE analyse_central_region
