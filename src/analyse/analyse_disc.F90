! analyse_disc.F90
! D. Stamatellos 
! Reads in a file and assuming a disc centred around the first sink 
!	(i)   centres around the first sink and rotates disc so that xy is the disc midplane "infile.trans"
!	(ii)  calculates disc properties (Q, T, Sigma, Omega, Mass vs radius)
!	(iii) collective files with what happens in a given radius vs time (appends on file)
!
!
!	output files  
!---------------------------------------------------------------------------------------------------------------------
!   .du.cdisc.1		disc file centered around the star (1st sink)
!	.du.cdisc.2		disc file centered around the planet (2nd sink)
!	.du.cdisc.cr    disc file on the planet corotational frame (centered around the star) 
!    these files are in the same format as the output seren files (they could also be in the df format)
!---------------------------------------------------------------------------------------------------------------------
!	.pdisc.1        azimuthaly averaged properties (around the star)
!	.pdisc.2		azimuthaly averaged properties (around the planet-- can be used to find CP disc)
!	.pdisc.cr   	azimuthaly averaged properties (around the star--on planets corotational frame)
!
!    the paramaters are shown on the first line of the file 
!---------------------------------------------------------------------------------------------------------------------
!
!    .rdisc.x.dis   property evolution over time (x:1,2, cr; dis: distance from star/planet)
!
!    the parameters are not shown in the file (but listed here below... search "rdisc.x.dis")
!
!---------------------------------------------------------------------------------------------------------------------


!	 .rdisc.x       file that contains disc size and mass (with different definitions) search "rdisc.x"
!
!
! ============================================================================
		! IMPORT MODULES
! ============================================================================

#include "macros.h"

SUBROUTINE analyse_disc(in_file)
  use type_module
  use particle_module
  use hydro_module
  use sink_module
  use time_module, only : time,snapshot
  use scaling_module
#ifdef ANALYSE_DISC
  use analyse_module
#endif
  use filename_module, only :run_id,out_file_form,run_dir
!
 implicit none
!
! ============================================================================
		! VARIABLE DEFINITION
! ============================================================================
!
 character(len=40) , intent(in) :: in_file
 character(len=40) :: out_file				! snapshot name
 character(len=8) :: file_ext				! filename extension for data output
 character(len=8) :: file_ext2				! filename extension for data output
 character(len=8) :: file_ext3				! filename extension for data output
 character(len=6) :: file_numb				! filename extension for data output
!
 real(kind=PR) :: radius_in				!inner radius in AU
 real(kind=PR) :: radius_out				!outer radius in AU
 integer :: rbins					!number of radial bins	
 real(kind=PR) :: allowance=10 
!
 real(kind=PR) :: sigma_cutoff1=2			! PdBI min surface density
 real(kind=PR) :: sigma_cutoff2=20			! ALMA min surface density
!
 integer :: i,j						! radial bin counter
 integer :: p						! particle counter
 integer :: s,nsin,counter				! sink counter
 integer :: correct
 real(kind=PR),allocatable :: radius(:)			! radius of bin
 real(kind=PR),allocatable :: rsigma(:)			! surface density 	
 real(kind=PR),allocatable :: rtoomre(:)		! toomre parameter
 real(kind=PR),allocatable :: rtemp(:)			! temperature at midplane (<0.05)
 real(kind=PR),allocatable :: rmass(:)			! mass in rbin
 real(kind=PR),allocatable :: romega (:)		! rotational velocity
 real(kind=PR),allocatable :: rsound(:)			! sound speed
 real(kind=PR),allocatable :: rinmass(:)		! mass within radius r (no sink included)
!
 real(kind=PR),allocatable :: rv_r(:)			! average radial velocity for particles at radius r
 real(kind=PR),allocatable :: rv_theta(:)		! average azimuthal velocity ...
 real(kind=PR),allocatable :: rv_z(:)			! average velocity 
 real(kind=PR),allocatable :: v_r(:)			! radial velocity
 real(kind=PR),allocatable :: v_theta(:)		! azimuthal velocity
 real(kind=PR),allocatable :: v_z(:)			! z velocity
 real(kind=PR),allocatable :: v_kepler(:)		! keplerian velocity 
 real(kind=PR),allocatable :: smooth_av(:)		! average smoothing length
 real(kind=PR),allocatable :: smooth_midplane_av(:)	! average smoothing length in midplane
!  
 integer, allocatable :: rparticles(:)
 integer, allocatable :: midplane_rparticles(:)
 real(kind=PR) :: pradius				! distance on the disc midplane
!
 real(kind=PR) :: rdisc_kepler									    
 real(kind=PR) :: mdisc_kepler
 real(kind=PR) :: rdisc_comp1
 real(kind=PR) :: mdisc_comp1
 real(kind=PR) :: rdisc_comp2
 real(kind=PR) :: mdisc_comp2
 real(kind=PR) :: rdisc_sigma1
 real(kind=PR) :: mdisc_sigma1
 real(kind=PR) :: rdisc_sigma2
 real(kind=PR) :: mdisc_sigma2
 real(kind=PR) :: dummy_rsink(1:3)
 real(kind=PR) :: dummy_vsink(1:3)
 real(kind=PR) :: r_hill=0.0				! hill radius
 real(kind=PR) :: mdisc_hill=0.0			! mass within the hill radius (sink/planet not included) 
!
 real(kind=DP), parameter :: mH = 1.660539D-24		! Hydrogen mass
 real(kind=DP), parameter :: Gconst = 6.67428D-8	! Grav. constant
 real(kind=DP), parameter :: kboltzmann = 1.38065D-16	! Boltzmann constant 
 real(kind=DP), parameter :: pi = 3.14159
!
 real(kind=PR) :: r1,r2,dr,ytmp,ztmp,vxtmp,vytmp,vztmp	! auxiliary variables
 real(kind=PR) :: sintheta,costheta,xtmp,costheta_rotation, sintheta_rotation
!
!
! BM defined variables - FUORs project (Q4 2015 -> Q1 2016) 
!
!
 real(kind=PR) :: theta_r, theta_d, theta_incl		!!! Parameters for rotation of disc
 real(kind=PR) :: sum_Lx, sum_Ly, sum_Lz, L_d, L_xy	!!! to x-y plane
 real(kind=PR) :: cosphi, sinphi
 integer :: align_count
 real(kind=PR) :: sini,cosi,cosi_rotation, sini_rotation
 integer :: rot_diagn = 1				!!! Choose whether (1) or not (0) to check rotation angles
 integer :: rotate_ref = 1				!!! Define which feature rotation based on (0 for central sink, 1 for disc) 
!
 real(kind=DP), parameter :: trestr = 1000.		!!! Set restrictive values of density and temperature
 real(kind=DP), parameter :: rhorestr = 1.E-13		!!! (in cgs) to exclude "weird" rho-T plane particles
 integer :: trho_count					!!! Counter for particles within restriction
 integer :: trho_diagn = 0				!!! Choose whether (1) or not (0) to check rho-T plane restricted particles
!
 real(kind=DP), parameter :: delt_r = 4.		!!! Define parameters for PV
 real(kind=DP), parameter :: delt_v = 0.5		!!! diagram
 real(kind=DP), parameter :: r_range = 200.
 real(kind=DP), parameter :: v_range = 60.
 real(kind=DP), parameter :: y_min = -30.
 real(kind=DP), parameter :: y_max = 30. 
!
 real(kind=PR), allocatable :: pv_x(:)			!!! Define arrays for PV diagram  
 real(kind=PR), allocatable :: pv_y(:)			!!! outputs
 real(kind=PR), allocatable :: pv_vz(:)
 integer :: rint, r_iter, vint, v_iter, r_count, v_count, hist_ind, hist_count
 real(kind=PR), allocatable :: r_arr1(:)
 real(kind=PR), allocatable :: r_arr2(:) 
 real(kind=PR), allocatable :: v_arr1(:)
 real(kind=PR), allocatable :: v_arr2(:) 
!
 character(len=40) :: outfile_rawpv			!!! PV diagram snapshot filenames
 character(len=40) :: outfile_histpv
 character(len=40) :: outfile_fitpv1
 character(len=40) :: outfile_fitpv2
!
 integer :: pv_print = 0				!!! Choose whether (1) or not (0) to generate P-V data
!
 integer :: script_params = 1				!!! Choose whether (1) or not (0) parameters fed from file
!
 real(kind=PR) :: inclin, restrkep	!!! Define inclination and Keplerian restriction parameters, and read in
!					!!! either from file or in this script as defined in script_params
 if (script_params .eq. 1) then
	open(unit = 11, file = "/home/ben/Documents/WORK_PLANETS/PROJECTS/FUORS/GENERATE/varfile.dat", &
	   status = "old")
	read(11,*) inclin, restrkep
	close(11)
!
 endif
!
 if (script_params .eq. 0) then
	inclin = 90.
	restrkep = 90.
 endif
!
 restrkep = restrkep/100.				!!! Convert % to fraction as required by code below
!
! ============================================================================
		! MAIN PROGRAM
! ============================================================================
!
 ! initialize bins for sinks and convert to AU

 radius_in=disc_in_radius/206265.            	!inner radius in AU
 radius_out=disc_out_radius/206265.             !outer radius in AU
 rbins=disc_nbins			!number of radial bins	


 !initialize 
 if (stot==0) stop

 ! allocate arrays
 allocate(radius(1:rbins+1))
 allocate(rsigma(1:rbins+1))
 allocate(rtemp(1:rbins+1)) 
 allocate(romega(1:rbins+1))
 allocate(rsound(1:rbins+1))
 allocate(rtoomre(1:rbins+1))
 allocate(rmass(1:rbins+1))
 allocate(rparticles(1:rbins+1))
 allocate(midplane_rparticles(1:rbins+1))
 allocate(rinmass(1:rbins+1))
 allocate(rv_r(1:rbins+1))
 allocate(rv_theta(1:rbins+1))
 allocate(rv_z(1:rbins+1))
 allocate(v_r(1:ptot))
 allocate(v_theta(1:ptot))
 allocate(v_z(1:ptot))
 allocate(v_kepler(1:ptot))
 allocate(smooth_av(1:ptot))
 allocate(smooth_midplane_av(1:ptot))

 counter=1

! create files both for central star and planet ------
#ifdef PLANET_IN_DISC


 do counter=1,3  ! 3 is for converting to corotational frame of the planet


    if (counter==2) then
  ! initialize bins for planet disc

    radius_in=pdisc_in_radius/206265.             	!inner radius in AU
    radius_out=pdisc_out_radius/206265.             !outer radius in AU
    rbins=pdisc_nbins			!number of radial bins	

    endif


    if (counter==3) then
  ! initialize bins for planet disc

    radius_in=disc_in_radius/206265.             	!inner radius in AU
    radius_out=disc_out_radius/206265.             !outer radius in AU
    rbins=disc_nbins			!number of radial bins	

    endif

!calculate planet angle after star is in the origin of coordinates (after first transformation is done)
    if (counter==2) then

       costheta_rotation=sink(2)%r(1)/sqrt(sink(2)%r(1)**2+sink(2)%r(2)**2)! angle to co-rotational frame
       sintheta_rotation=-sink(2)%r(2)/sqrt(sink(2)%r(1)**2+sink(2)%r(2)**2)! angle to co-rotational frame

    endif
    
    nsin=counter
    
    if (counter==3) then
          nsin=2
          dummy_rsink(1:3)=sink(1)%r(1:3)   
          dummy_vsink(1:3)=sink(1)%v(1:3)
       else
          dummy_rsink(1:3)=sink(nsin)%r(1:3)
          dummy_vsink(1:3)=sink(nsin)%v(1:3)
 	endif

#else
do nsin=1,1 ! no planet in disc
	
    dummy_rsink(1:3)=sink(nsin)%r(1:3)
    dummy_vsink(1:3)=sink(nsin)%v(1:3)
	
#endif


   if (counter==3) then
	   file_ext3="cr"
    else   
   	  write(file_ext3,"(I1)") counter
   endif
   

! stop when the second sink forms in the disc

! if (stot>1)  stop 


 ! write(*,*) "Number of sinks",stot
 
! Centre around the sink

! a. particles 
  do p=1,ptot
   parray(1:3,p)=parray(1:3,p)-dummy_rsink(1:3)
   v(1:3,p)=v(1:3,p)-dummy_vsink(1:3)
  enddo

 
! b. sinks  
  do s=1,stot   
    sink(s)%r(1:3)=sink(s)%r(1:3)-dummy_rsink(1:3)
    sink(s)%v(1:3)=sink(s)%v(1:3)-dummy_vsink(1:3)
  end do
!
!
!
!!! Select which reference vector to be used to rotate disc - either central sink or 
!!! disc midplane vector
!
 if (rotate_ref .eq. 0) then				!!! Rotate disc with reference
!							!!! to central sink L vector
	costheta = abs(sink(nsin)%angmom(3))/ &	
	   sqrt(sink(nsin)%angmom(1)**2+sink(nsin)%angmom(2)**2+sink(nsin)%angmom(3)**2)
	cosphi = abs(sink(nsin)%angmom(1))/ &
	   sqrt(sink(nsin)%angmom(1)**2+sink(nsin)%angmom(2)**2+sink(nsin)%angmom(3)**2)
!
	if (rot_diagn .eq. 1) then
		write(*,*) "Rotation taken about central sink L vector"
		write(*,*) "Angle disc rotated about is: ", costheta / (pi / 180.)
		write(*,*) "Angle centering onto z = 0 is: ", cosphi / (pi / 180.)
	endif
!
!							!!! Ensure rotation about y-axis 
!							!!! is in correct direction (clockwise/anticlockwise)
	if ( sink(nsin)%angmom(1) .gt. 0.) then
		cosphi = -cosphi
	endif
 endif
!
 if (rotate_ref .eq. 1) then				!!! Rotate disc with reference
	sum_Lx = 0.					!!! to disc midplane summation
	sum_Ly = 0.					!!! of L vectors
	sum_Lz = 0.
!
	align_count = 0
	do p=1,ptot
		if ( ( sqrt(parray(1,p)**2 + parray(2,p)**2 + parray(3,p)**2)*rscale*206265. .lt. 50 ) ) then		
!		if ( ( sqrt(parray(1,p)**2 + parray(2,p)**2 + parray(3,p)**2)*rscale*206265. .lt. 50 ) .and. &
!		   (rho(p)*rhoscale .gt. 1.e-12) ) then 
!
			align_count = align_count + 1
!
			sum_Lx = sum_Lx + parray(MASS, p) * (v(3,p)*parray(2,p) - v(2,p)*parray(3,p) )
			sum_Ly = sum_Ly + parray(MASS, p) * (v(1,p)*parray(3,p) - v(3,p)*parray(1,p) )
			sum_Lz = sum_Lz + parray(MASS, p) * (v(2,p)*parray(1,p) - v(1,p)*parray(2,p) )
		endif
  	enddo
!
	L_d = sqrt(sum_Lx**2 + sum_Ly**2 + sum_Lz**2)
!
	costheta = abs(sum_Lz) / L_d
	cosphi = abs(sum_Lx) / L_d
!
	if (sum_Lx .gt. 0.) then
		cosphi = -cosphi
	endif
!
	if (rot_diagn .eq. 1) then
		write(*,*) "Rotation taken about disc L vector"
		write(*,*) "Count of particles in density restriction is: ", align_count
		write(*,*) "Angle disc rotated about is: ", costheta / (pi / 180.)
		write(*,*) "Angle centering onto z = 0 is: ", cosphi / (pi / 180.)
	endif
 endif	
!
!							!!! Calculate sin(theta) from trig. identity
 sintheta=sqrt(1-costheta**2)
!
!							!!! Rotate disc to lie on z = 0 plane
							!!! prior to inclination
 cosphi = cos( ( ( acos(cosphi) / (pi / 180.) ) + 90.) * (pi / 180.) )
 sinphi = sqrt( 1 - cosphi**2 )
!
!							!!! Ensure rotation about y-axis 
!							!!! is in correct direction (clockwise/anticlockwise)
 if ( (sum_Lx .gt. 0.) .or. (sink(nsin)%angmom(1) .gt. 0.) ) then
 	sinphi = -sinphi
 endif
!
!!! With rotation angles defined for correct alignment of disc, rotate disc about x and z
!!! axis as per previous D. Stamatellos analyses. Then rotate about y-axis to align
!!! disc to z = 0 plane. 
!!! Once disc is properly aligned, rotate with reference to inclin variable to 
!!! provide LOS inclination for disc.
!
 do p=1,ptot						!!! Rotate about x-axis with theta
	ytmp=parray(2,p)   
	ztmp=parray(3,p)
	vytmp=v(2,p)  
	vztmp=v(3,p)
	parray(2,p)=costheta*ytmp-sintheta*ztmp
	parray(3,p)=sintheta*ytmp+costheta*ztmp   
	v(2,p)=costheta*vytmp-sintheta*vztmp
	v(3,p)=sintheta*vytmp+costheta*vztmp 
 enddo
!
 do s=1,stot
	ytmp=sink(s)%r(2)  
	ztmp=sink(s)%r(3) 
	vytmp=sink(s)%v(2)
	vztmp=sink(s)%v(3)
	sink(s)%r(2)=costheta*ytmp-sintheta*ztmp
	sink(s)%r(3)=sintheta*ytmp+costheta*ztmp 
	sink(s)%v(2)=costheta*vytmp-sintheta*vztmp
	sink(s)%v(3)=sintheta*vytmp+costheta*vztmp     
 end do
!
 if (counter==3) then
!
  	costheta=costheta_rotation
  	sintheta=sintheta_rotation
!
	do p=1,ptot
		xtmp=parray(1,p)			!!! Rotate about z-axis with theta
		ytmp=parray(2,p)   
		vytmp=v(2,p)  
		vxtmp=v(1,p)
		parray(1,p)=costheta*xtmp-sintheta*ytmp
		parray(2,p)=sintheta*xtmp+costheta*ytmp
		v(1,p)=costheta*vxtmp-sintheta*vytmp
		v(2,p)=sintheta*vxtmp+costheta*vytmp 
	enddo
! 
	do s=1,stot
		ytmp=sink(s)%r(2)  
		xtmp=sink(s)%r(1)
		vytmp=sink(s)%v(2)  
		vxtmp=sink(s)%v(1) 
		sink(s)%r(1)=costheta*xtmp-sintheta*ytmp
		sink(s)%r(2)=sintheta*xtmp+costheta*ytmp 
		sink(s)%v(1)=costheta*vxtmp-sintheta*vytmp
		sink(s)%v(2)=sintheta*vxtmp+costheta*vytmp     
	end do
!
	do p=1,ptot
		v(1:3,p)=v(1:3,p)-sink(2)%v(1:3)   
	enddo     
 endif
!
 do p=1,ptot
	ztmp=parray(3,p)			!!! Rotate about y-axis with phi 
	xtmp=parray(1,p)   
	vztmp=v(3,p)  
	vxtmp=v(1,p)
	parray(3,p)=cosphi*ztmp-sinphi*xtmp
	parray(1,p)=sinphi*ztmp+cosphi*xtmp
	v(3,p)=cosphi*vztmp-sinphi*vxtmp
	v(1,p)=sinphi*vztmp+cosphi*vxtmp 
 enddo
! 
 do s=1,stot
	xtmp=sink(s)%r(1)  
	ztmp=sink(s)%r(3) 
	vztmp=sink(s)%v(3)
	vxtmp=sink(s)%v(1)
	sink(s)%r(3)=cosphi*ztmp-sinphi*xtmp
	sink(s)%r(1)=sinphi*ztmp+cosphi*xtmp 
	sink(s)%v(3)=cosphi*vztmp-sinphi*vxtmp
	sink(s)%v(1)=sinphi*vztmp+cosphi*vxtmp     
 end do
!
!!! Now incline disc with user defined angle, only rotating about x and z axes
!
 cosi = cos ( inclin * (pi / 180.) )
 sini=sqrt(1-cosi**2)
!
 do p=1,ptot						!!! Rotate about x-axis with inclin
	ytmp=parray(2,p)   
	ztmp=parray(3,p)
	vytmp=v(2,p)  
	vztmp=v(3,p)
	parray(2,p)=cosi*ytmp-sini*ztmp
	parray(3,p)=sini*ytmp+cosi*ztmp   
	v(2,p)=cosi*vytmp-sini*vztmp
	v(3,p)=sini*vytmp+cosi*vztmp 
 enddo
!
 do s=1,stot
	ytmp=sink(s)%r(2)  
	ztmp=sink(s)%r(3) 
	vytmp=sink(s)%v(2)
	vztmp=sink(s)%v(3)
	sink(s)%r(2)=cosi*ytmp-sini*ztmp
	sink(s)%r(3)=sini*ytmp+cosi*ztmp 
	sink(s)%v(2)=cosi*vytmp-sini*vztmp
	sink(s)%v(3)=sini*vytmp+cosi*vztmp     
 end do
!
if (counter==3) then
!
  	cosi=cosi_rotation
  	sini=sini_rotation
!
	do p=1,ptot
		xtmp=parray(1,p)			!!! Rotate about z-axis with inclin
		ytmp=parray(2,p)   
		vytmp=v(2,p)  
		vxtmp=v(1,p)
		parray(1,p)=cosi*xtmp-sini*ytmp
		parray(2,p)=sini*xtmp+cosi*ytmp
		v(1,p)=cosi*vxtmp-sini*vytmp
		v(2,p)=sini*vxtmp+cosi*vytmp 
	enddo
! 
	do s=1,stot
		ytmp=sink(s)%r(2)  
		xtmp=sink(s)%r(1) 
		vytmp=sink(s)%v(2)  
		vxtmp=sink(s)%v(1) 
		sink(s)%r(1)=cosi*xtmp-sini*ytmp
		sink(s)%r(2)=sini*xtmp+cosi*ytmp 
		sink(s)%v(1)=cosi*vxtmp-sini*vytmp
		sink(s)%v(2)=sini*vxtmp+cosi*vytmp     
	end do
!
	do p=1,ptot
		v(1:3,p)=v(1:3,p)-sink(2)%v(1:3)   
	enddo     
 endif
!
!
!
  ! write file with sink at the origin of the system and disc on xy-plane
! do p=1,ptot
!
!	if ( (rho*rhoscale .gt. rhorestr) .or. (temp .lt. trestr) ) then

		write(file_numb,"(I5.5)") snapshot

!	endif

! enddo

	
     file_ext = "."//file_numb


  out_file = trim(adjustl(in_file))//"."//trim(adjustl(out_file_form))//".cdisc."//trim(adjustl(file_ext3))


  call write_data(out_file,out_file_form)


  
  
!  if (counter==3) EXIT ! in the corotational frame these are not needed

  
! find v_r, v_theta for calculating size of disc
 do p=1,ptot
  vxtmp=v(1,p)
  vytmp=v(2,p)
  vztmp=v(3,p)

  pradius=sqrt(parray(1,p)**2+parray(2,p)**2) 

  costheta=parray(1,p)/pradius
  sintheta=parray(2,p)/pradius

  v_r(p)=vxtmp*costheta+vytmp*sintheta
  v_theta(p)=-vxtmp*sintheta+vytmp*costheta
  v_z(p)=vztmp

 enddo												    

 dr=(radius_out-radius_in)/rbins
  
  do i=1,rbins+1
  radius(i)=radius_in+(i-1)*dr
  enddo

!initialize
  rtemp=0.
  rmass=0.
  rparticles=0.
  midplane_rparticles=0.
  rv_r=0.
  rv_theta=0.
  rv_z=0.
  smooth_av=0.
  smooth_midplane_av=0.

!place into radial bins, restricting rho-T plane dependent on variables entered

 if (trho_diagn .eq. 1) then
	write(*,*) "Reference of simulation (scaled) and density/temperatures as follows"
	write(*,*) "Density [minimum, maximum, restriction]"
	write(*,*) "	", minval(rho*(rhoscale)), maxval(rho*(rhoscale)), rhorestr
	write(*,*) "Temperature [minimum, maximum, restriction]"
	write(*,*) "	", minval(temp), maxval(temp), trestr
 endif
!
 trho_count = 0
 do p=1,ptot
 	pradius=sqrt(parray(1,p)**2+parray(2,p)**2) 
!
	if ( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
!
		if ( trho_diagn .eq. 1 ) then
			write(*,*) "Particle ", p, " with rho = ", rho(p)*rhoscale, " and T = ",temp(p), " restricted"
			CYCLE
		endif
!
		trho_count = trho_count + 1
!
		if (pradius>radius_out+0.5*dr) CYCLE
 		if (parray(3,p)>0.5*pradius) CYCLE
 		do i=1,rbins+1
 			if (pradius>radius(i)-0.5*dr .and.pradius<radius(i)+0.5*dr) then
 		            rparticles(i)=rparticles(i)+1
 		            rmass(i)=rmass(i)+parray(MASS,p)
 		            if (parray(3,p)<0.05*pradius) then
				midplane_rparticles(i)=midplane_rparticles(i)+1
 	 			rtemp(i)=rtemp(i)+temp(p)
				smooth_midplane_av(i)=smooth_midplane_av(i)+parray(SMOO,p)
			
 	 		    endif
			    rv_r(i)=rv_r(i)+v_r(p)
			    rv_theta(i)=rv_theta(i)+v_theta(p)
			    rv_z(i)=rv_z(i)+v_z(p)
				smooth_av(i)=smooth_av(i)+parray(SMOO,p)
		
 	 		endif
 	 	enddo
	endif		    
 enddo
!
 if (trho_diagn .eq. 1) then
	write(*,*) "Number of particles after T-rho plane restriction is: ", trho_count
 endif
!
 if (pv_print .eq. 1) then
!
	allocate(pv_x(1:trho_count))								    !!! Loop over particles, if rho and T	
	allocate(pv_y(1:trho_count))								    !!! restrictions met, output data in array
	allocate(pv_vz(1:trho_count))								    !!! for PV diagram (raw, histogram and fit)
!
	i = 1
	do p=1,ptot
		if( (rho(p)*rhoscale .gt. rhorestr) .or. (temp(p) .lt. trestr) ) then
			pv_x(i) = parray(1,p)*rscale*206265.
			pv_y(i) = parray(2,p)*rscale*206265.
			pv_vz(i) = v(3,p)*vscale  
			i = i + 1
		endif
	enddo
!
	close(2)
!
	outfile_histpv = trim(adjustl(in_file))//"."//trim(adjustl(out_file_form))//".hist_pv."//trim(adjustl(file_ext3))
	open(UNIT=2,file=outfile_histpv,form='formatted',position='append') 
!
	write(2, '(E15.7)') r_range
	write(2, '(E15.7)') v_range
	write(2, '(E15.7)') delt_r
	write(2, '(E15.7)') delt_v
	write(2, '(E15.7)') y_min
	write(2, '(E15.7)') y_max
!
	outfile_fitpv1 = trim(adjustl(in_file))//"."//trim(adjustl(out_file_form))//".fit_Lpv."//trim(adjustl(file_ext3))
	open(UNIT=3,file=outfile_fitpv1,form='formatted',position='append') 

	outfile_fitpv2 = trim(adjustl(in_file))//"."//trim(adjustl(out_file_form))//".fit_Rpv."//trim(adjustl(file_ext3))
	open(UNIT=4,file=outfile_fitpv2,form='formatted',position='append') 

	r_iter = int(2 * (r_range / delt_r) )
	v_iter = int(2 * (v_range / delt_v) )

	allocate(r_arr1(1:r_iter))
	allocate(r_arr2(1:r_iter))
	allocate(v_arr1(1:v_iter))
	allocate(v_arr2(1:v_iter))
!
	r_count = 0
	do rint=1, r_iter
		r_count = r_count + 1	
		r_arr1(rint) = ((r_count-1)*delt_r)-r_range
		r_arr2(rint) = ((r_count)*delt_r)-r_range
!
		v_count = 0	
		do vint=1, v_iter
			v_count = v_count + 1
			v_arr1(vint) = ((v_count-1)*delt_v)-v_range
			v_arr2(vint) = ((v_count)*delt_v)-v_range
!
			hist_count = 0
			do i=1, trho_count
				if ( ( pv_x(i) .gt. r_arr1(rint) ) .and. ( pv_x(i) .lt. r_arr2(rint) ) &
				   .and. (pv_vz(i) .gt. v_arr1(vint) ) .and. (pv_vz(i) .lt. v_arr2(vint) ) & 
				   .and. (pv_y(i) .gt. y_min ) .and. (pv_y(i) .lt. y_max ) ) then
					hist_count = hist_count + 1
				endif
			enddo
!
			write(2, '(2E15.7,2X, I7)') r_arr1(rint), v_arr1(vint), hist_count
!
			if ( (hist_count .ne. 0) .and. (r_arr1(rint) .le. 0.) ) then
				write(3, '(2E15.7, 2X, I7)') r_arr1(rint), v_arr1(vint), hist_count				
			endif
!
			if ( (hist_count .ne. 0) .and. (r_arr1(rint) .ge. 0.) ) then
				write(4, '(2E15.7, 2X, I7)') r_arr1(rint), v_arr1(vint), hist_count				
			endif
!
		enddo
	enddo
!
	close(2)
	close(3)
	close(4)
!
 endif
!
!calculate radial average parameters

 do i=1,rbins+1
        r1=radius(i)-0.5*dr
        r2=radius(i)+0.5*dr        
	rtemp(i)=rtemp(i)/midplane_rparticles(i)
    smooth_av(i)=smooth_av(i)/rparticles(i)
    smooth_midplane_av(i)=smooth_midplane_av(i)/midplane_rparticles(i)
	
        romega(i)=sqrt(Gconst*sink(1)%m*2d33/(radius(i)*206265*1.496d13)**3)
        rsigma(i)=2.17d-4*rmass(i)/(3.14*(r2**2-r1**2))
        rsound(i)=sqrt(kboltzmann*rtemp(i)/(2.3*mH))
     	rtoomre(i)=rsound(i)*romega(i)/(3.14*Gconst*rsigma(i))
	rv_r(i)=rv_r(i)/rparticles(i)
	rv_theta(i)=rv_theta(i)/rparticles(i)
	rv_z(i)=rv_z(i)/rparticles(i)
 enddo    
 
! calculate mass within radius r2
       rinmass(1)=rmass(1)
     
do i=2,rbins+1
        rinmass(i)=rmass(i)+rinmass(i-1)
enddo

do i=1,rbins+1
!calculate keplerian velocity
#ifdef EPISODIC_ACCRETION
  v_kepler(i)=sqrt((rinmass(i)+sink(nsin)%Mstar+sink(nsin)%Mdisc)/(radius(i)))
#else
  v_kepler(i)=sqrt((rinmass(i)+sink(nsin)%Mstar+0.0)/(radius(i)))
#endif
enddo

!output files     

!collective files with what happens in a given radius vs time (appends on file)  rdisc.x.dis  
  do i=1,RBINS+1,5
  write(file_ext2,"(I3)") NINT(radius(i)*206265)


  out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          &".rdisc."//trim(adjustl(file_ext3))//"."//trim(adjustl(file_ext2))

  open(UNIT=1,file=out_file,form='formatted',position='append') 

#ifdef EPISODIC_ACCRETION

#ifdef PLANET_IN_DISC
  write(1,'(21E15.7,I7,2X,18A)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	 !# 1,2,3
  			& rtemp(i),rsigma(i), romega(i), &                                   !# 4,5,6
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), &            !# 7,8
			& sink(1)%Mstar*real(mscale,PR), sink(1)%Mdisc*real(mscale,PR),&     !# 9,10
                        & sink(2)%Mstar*real(mscale,PR), sink(2)%Mdisc*real(mscale,PR),&     !# 11, 12
                        & sqrt((sink(2)%r(1)-sink(1)%r(1))**2+(sink(2)%r(2)-sink(1)%r(2))**2)*real(rscale,PR)*206265,& !# 13
			& rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),&  !# 14,15,16
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+sink(nsin)%Mdisc)/real(mscale,PR))/(radius(i)))*real(vscale,PR), & ! #17
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 18, 19
			& r_hill, mdisc_hill,& !# 20, 21
			& rparticles(i), trim(adjustl(in_file)) !# 22, 23
#else
  write(1,'(16E15.7,I7,2X,18A)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	
  			& rtemp(i),rsigma(i), romega(i), &
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), &
			& sink(1)%Mstar*real(mscale,PR), sink(1)%Mdisc*real(mscale,PR),&
			& rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),&
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+sink(nsin)%Mdisc)/real(mscale,PR))/(radius(i)))*real(vscale,PR),& ! #17
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 18, 19
			& rparticles(i),&
			& trim(adjustl(in_file))

#endif
  
#else

#ifdef PLANET_IN_DISC
  write(1,'(21E15.7,I7,2X,18A)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	             
  			& rtemp(i),rsigma(i), romega(i), &                                 
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), &           
			& sink(1)%Mstar*real(mscale,PR), 0.0,&    
                        & sink(2)%Mstar*real(mscale,PR), 0.0,&    
                        & sqrt((sink(2)%r(1)-sink(1)%r(1))**2+(sink(2)%r(2)-sink(1)%r(2))**2)*real(rscale,PR)*206265,&
			& rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),&  
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+0.0)/real(mscale,PR))/(radius(i)))*real(vscale,PR),& ! #17
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 18, 19
			& r_hill, mdisc_hill,& !# 20, 21
			& rparticles(i),&
			& trim(adjustl(in_file))
#else
  write(1,'(16E15.7,I7,2X,18A)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	
  			& rtemp(i),rsigma(i), romega(i), &
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), &
			& sink(1)%Mstar*real(mscale,PR), 0.0,&
			& rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),&
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+sink(nsin)%Mdisc)/real(mscale,PR))/(radius(i)))*real(vscale,PR),& ! #17
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 18, 19
			&rparticles(i),&
			& trim(adjustl(in_file))

#endif

#endif

  close(1)
enddo

! write file that contains disc mass and size with various definitions
 
! initialize variables
   rdisc_kepler=radius(RBINS+1)*206265
   mdisc_kepler=rinmass(RBINS+1) 
   rdisc_comp1=radius(RBINS+1)*206265
   mdisc_comp1=rinmass(RBINS+1)
   rdisc_comp2=radius(RBINS+1)*206265
   mdisc_comp2=rinmass(RBINS+1)
   rdisc_sigma1=radius(RBINS+1)*206265
   mdisc_sigma1=rinmass(RBINS+1)
   rdisc_sigma2=radius(RBINS+1)*206265
   mdisc_sigma2=rinmass(RBINS+1) 

do i=1,RBINS+1
! disc size using the abs keplerian velocity							    !!! Edit of keplerian criterion using 	
  if (abs(rv_theta(i))<(restrkep*abs(v_kepler(i)))) then					    !!! value prescribed in restrkep (L98) 
! check a few more shells ( allowance AU) just in case there is a gap or something
       correct=0 
       do j=i+1,i+NINT(allowance/((radius(2)-radius(1))*206265))
        if (j>RBINS+1) exit
        if (abs(rv_theta(j))>(restrkep*abs(v_kepler(j))))  correct=correct+1		           
       enddo
       
       if (correct==0) then
	rdisc_kepler=radius(i)*206265
        mdisc_kepler=rinmass(i)       
	exit
       endif
  endif
enddo

do i=1,RBINS+1
! disc size using abs v_theta>2*pi*v_r
  if (abs(rv_theta(i))<2*3.14*abs(rv_r(i))) then

       correct=0 
       do j=i+1,i+NINT(allowance/((radius(2)-radius(1))*206265))
         if (j>RBINS+1) exit
         if (abs(rv_theta(j))>2*3.14*abs(rv_r(j)))  correct=correct+1
       enddo

       if (correct==0) then
	rdisc_comp1=radius(i)*206265
        mdisc_comp1=rinmass(i)
	exit
       endif
  endif
enddo

do i=1,RBINS+1
! disc size using abs v_theta>v_r
  if (abs(rv_theta(i))<abs(rv_r(i))) then
       correct=0
       do j=i+1,i+NINT(allowance/((radius(2)-radius(1))*206265))
         if (j>RBINS+1) exit
         if (abs(rv_theta(j))>abs(rv_r(j)))  correct=correct+1
       enddo

       if (correct==0) then
	rdisc_comp2=radius(i)*206265
        mdisc_comp2=rinmass(i)
        exit
       endif
  endif

enddo

do i=1,RBINS+1
! disc size using min surface density 1st cutoff
  if (rsigma(i)<sigma_cutoff1) then
       correct=0
       do j=i+1,i+NINT(allowance/((radius(2)-radius(1))*206265))
         if (j>RBINS+1) exit
         if (rsigma(j)>sigma_cutoff1) correct=correct+1
       enddo

       if (correct==0) then
	rdisc_sigma1=radius(i)*206265
        mdisc_sigma1=rinmass(i)
	exit
      endif

  endif
enddo

do i=1,RBINS+1
! disc size using min surface density 2nd cutoff 
  if (rsigma(i)<sigma_cutoff2) then

       correct=0
       do j=i+1,i+NINT(allowance/((radius(2)-radius(1))*206265))
         if (j>RBINS+1) exit
         if (rsigma(j)>sigma_cutoff2)  correct=correct+1
       enddo

       if (correct==0) then
	rdisc_sigma2=radius(i)*206265
        mdisc_sigma2=rinmass(i)
	exit
       endif
  endif
enddo

#ifdef PLANET_IN_DISC

if (counter>1) then
mdisc_hill=0.0	
r_hill= sqrt((sink(2)%r(1)-sink(1)%r(1))**2+(sink(2)%r(2)-sink(1)%r(2))**2)*real(rscale,PR)*206265*&
		&(sink(2)%Mstar/(3*sink(1)%Mstar))**(1./3.)
		write (*,* ) "ok", sink(1)%r, sink(2)%r
		do i=1,RBINS+1 
		  if (radius(i)*206265<r_hill) mdisc_hill=rinmass(i)
		enddo
	endif
	
#endif


! .rdisc.x

! file that contains disc size and mass (with different definitions)

   out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          &".rdisc."//trim(adjustl(file_ext3))

   open(UNIT=1,file=out_file,form='formatted',position='append') 

#ifdef PLANET_IN_DISC

#ifdef EPISODIC_ACCRETION

  write(1,'(E15.7,14(1X,E11.5),2X,18A)') time*tscale*1e6,sink(nsin)%Mstar*real(mscale,PR), sink(nsin)%Mdisc*real(mscale,PR),&  !1,2,3
			& rdisc_kepler,mdisc_kepler, & !4,5
			& rdisc_comp1,mdisc_comp1, &   !6,7
			& rdisc_comp2,mdisc_comp2, &   !8,9
			& rdisc_sigma1,mdisc_sigma1, & !10,11
			& rdisc_sigma2,mdisc_sigma2, & !12,13
			& r_hill, mdisc_hill, &		   !14,15
			& trim(adjustl(in_file))       !16
			
#else

  write(1,'(E15.7,14(1X,E11.5),2X,18A)') time*tscale*1e6,sink(nsin)%Mstar*real(mscale,PR), 0.0,& 
			& rdisc_kepler,mdisc_kepler, &
			& rdisc_comp1,mdisc_comp1, &
			& rdisc_comp2,mdisc_comp2, &
			& rdisc_sigma1,mdisc_sigma1, &
			& rdisc_sigma2,mdisc_sigma2, &
			& r_hill, mdisc_hill, &		   
			& trim(adjustl(in_file))

#endif


#else

#ifdef EPISODIC_ACCRETION


  write(1,'(E15.7,12(1X,E11.5),2X,18A)') time*tscale*1e6,sink(nsin)%Mstar*real(mscale,PR), sink(nsin)%Mdisc*real(mscale,PR),&  !1,2,3
			& rdisc_kepler,mdisc_kepler, & !4,5
			& rdisc_comp1,mdisc_comp1, &   !6,7
			& rdisc_comp2,mdisc_comp2, &   !8,9
			& rdisc_sigma1,mdisc_sigma1, & !10,11
			& rdisc_sigma2,mdisc_sigma2, & !12,13
			& trim(adjustl(in_file))       !14
			
#else

  write(1,'(E15.7,12(1X,E11.5),2X,18A)') time*tscale*1e6,sink(nsin)%Mstar*real(mscale,PR), 0.0,& 
			& rdisc_kepler,mdisc_kepler, &
			& rdisc_comp1,mdisc_comp1, &
			& rdisc_comp2,mdisc_comp2, &
			& rdisc_sigma1,mdisc_sigma1, &
			& rdisc_sigma2,mdisc_sigma2, &
			& trim(adjustl(in_file))

#endif

#endif
 
  close(1)




! file with disc properties vs radius   .pdisc.x

     out_file = trim(adjustl(in_file))//".pdisc."//trim(adjustl(file_ext3))

  open(UNIT=1,file=out_file,status='replace',form='formatted')
  
#ifdef PLANET_IN_DISC
  write(1,'(A)')"#t(yr)/r(AU)/Q/T/sig(gcm-2)/Omeg/DiscM(Msun)/rInM/M*/Mdis/Mp/Mpdis/Rp/vr/vthe/vz/vkep/h/h_mid/rhill/mhill/rpart"
#else
  write(1,'(A)')"#t(yr)/r(AU)/Q/T/sig(gcm-2)/Omeg/DiscM(Msun)/rInM/M*/Mdis/vr/vthe/vz/vkep/h/h_mid/rpart"
#endif
  
  
  do i=1,RBINS+1

#ifdef EPISODIC_ACCRETION

#ifdef PLANET_IN_DISC
     write(1,'(21E15.7,I7)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	             ! 1,2,3
  			& rtemp(i),rsigma(i), romega(i), &                                   ! 4,5,6
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), &            ! 7,8
			& sink(1)%Mstar*real(mscale,PR), sink(1)%Mdisc*real(mscale,PR),&     ! 9,10
                        & sink(2)%Mstar*real(mscale,PR), sink(2)%Mdisc*real(mscale,PR),&     ! 11, 12
                        & sqrt((sink(2)%r(1)-sink(1)%r(1))**2+(sink(2)%r(2)-sink(1)%r(2))**2)*real(rscale,PR)*206265,& ! 13
			& rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),&   ! 14,15,16
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+sink(nsin)%Mdisc)/real(mscale,PR))/(radius(i)))*real(vscale,PR),& ! #17
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 18, 19,
			& r_hill, mdisc_hill, & !20,21
			& rparticles(i) ! 22
		
#else
  write(1,'(16E15.7,I7)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	! 1,2,3
  			& rtemp(i),rsigma(i), romega(i), & ! 4,5,6
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), & !7,8
			& sink(1)%Mstar*real(mscale,PR), sink(1)%Mdisc*real(mscale,PR),& !9,10
			& rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),& !11,12,13
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+sink(nsin)%Mdisc)/real(mscale,PR))/(radius(i)))*real(vscale,PR),& ! #14
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 15, 16
			& rparticles(i) ! 17
#endif
  
#else
#ifdef PLANET_IN_DISC
 write(1,'(21E15.7,I7)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	
  			& rtemp(i),rsigma(i), romega(i), &
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), &
			& sink(1)%Mstar*real(mscale,PR), 0.0,&
                        & sink(2)%Mstar*real(mscale,PR), 0.0,&
                        & sqrt((sink(2)%r(1)-sink(1)%r(1))**2+(sink(2)%r(2)-sink(1)%r(2))**2)*real(rscale,PR)*206265,& 
                        & rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),&
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+0.0)/real(mscale,PR))/(radius(i)))*real(vscale,PR),& ! #17
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 18, 19
			& r_hill, mdisc_hill, & !20,21
			& rparticles(i) ! 22
#else
 write(1,'(16E15.7,I7)') time*tscale*1e6,radius(i)*206265,rtoomre(i),&	
  			& rtemp(i),rsigma(i), romega(i), &
			& rmass(i)*real(mscale,PR), rinmass(i)*real(mscale,PR), &
			& sink(1)%Mstar*real(mscale,PR), 0.0,&
			& rv_r(i)*real(vscale,PR),rv_theta(i)*real(vscale,PR),rv_z(i)*real(vscale,PR),&
			& sqrt(((rinmass(i)+sink(nsin)%Mstar+0.0)/real(mscale,PR))/(radius(i)))*real(vscale,PR),& ! #14
			& smooth_av(i)*206265, smooth_midplane_av(i)*206265, &!# 15, 16
			& rparticles(i) ! 17
#endif

#endif


  enddo			
  close(1)


  enddo

  return
END SUBROUTINE analyse_disc
