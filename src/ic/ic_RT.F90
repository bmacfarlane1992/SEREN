! IC_RT.F90
! D. A. Hubber - 26/05/2009
! Generate initial conditions to model Rayleigh-Taylor instability.  
! Creates two gas lattices of different densities in pressure balance 
! at the boundary.  If the dense gas is on the top, then the configuration 
! is unstable to the Rayleigh-Taylor instability.  We seed a perturbation 
! at the interface following one of two possible methods
! 1) Velocity perturbation
! 2) Boundary displacement perturbation
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_RT
  use interface_module, only : read_data,write_data
  use particle_module
  use hydro_module
  use filename_module
  use time_module
  use periodic_module
  use type_module
  use scaling_module
  use neighbour_module
  implicit none

  character(len=256) :: out_file      ! Name of output file
  character(len=256) :: in_file1      ! Name of output file
  character(len=256) :: in_file1_form ! Name of output file
  integer :: i                        ! x-dimension counter
  integer :: j                        ! y-dimension counter
  integer :: k                        ! z-dimension counter
  integer :: nlayers1                 ! No. of layers of medium 1 in y
  integer :: nlayers2                 ! No. of layers of medium 2 in y
  integer :: nwall1                   ! ..
  integer :: nwall2                   ! ..
  integer :: p                        ! Particle counter
  integer :: p1                       ! No. of particles in medium 1
  integer :: p2                       ! No. of particles in medium 2
  integer :: pertmode                 ! Mode of adding perturbation
  integer :: ppd1                     ! Particle per dimension
  integer :: ppd2                     ! Particle per dimension
  integer, allocatable :: ptype(:)    ! New order of particles
  real(kind=PR) :: acc_grav           ! External grav. accel. in y direction
  real(kind=PR) :: amp                ! Amplitude of velocity perturbation
  real(kind=PR) :: lambda             ! Wavelength of perturbation
  real(kind=PR) :: mp1                ! Mass of 'bottom' particles
  real(kind=PR) :: mp2                ! Mass of 'top' particles
  real(kind=PR) :: Press1             ! Constant pressure
  real(kind=PR) :: rho1               ! 'Bottom' density
  real(kind=PR) :: rho2               ! 'Top' density
  real(kind=PR) :: x1                 ! x-size of region 1
  real(kind=PR) :: x2                 ! x-size of region 2
  real(kind=PR) :: xsize              ! Inputted size of domain
  real(kind=PR) :: y1                 ! y-size of region 1
  real(kind=PR) :: y2                 ! y-size of region 2
  real(kind=PR) :: ywall1             ! ..
  real(kind=PR) :: ywall2             ! ..

#ifndef DIMENSIONLESS
  write(6,*) "Compiler flag error : Only works with DIMENSIONLESS flag on"
  stop
#endif
#if NDIM==1 || NDIM==3
  write(6,*) "Compiler flag error : Only works in 2D"
  stop
#endif

! Reading parameter file
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "            ic_RT             "
  write(6,*) "------------------------------"
  open(unit=1,file="RTparams.dat",status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) pertmode
  read(1,*) ppd1,ppd2
  read(1,*) nlayers1,nlayers2
  read(1,*) nwall1,nwall2
  read(1,*) rho1,rho2
  read(1,*) Press1
  read(1,*) acc_grav
  read(1,*) gamma
  read(1,*) xsize
  read(1,*) amp
  read(1,*) lambda
  read(1,*) pp_gather
  read(1,*) hmin
  read(1,*) h_fac

  p1   = ppd1*(nlayers1 + nwall1)
  p2   = ppd2*(nlayers2 + nwall2)
  ptot = p1 + p2
  write(6,*) "Number of particles : ",ptot
  write(6,*) "p1 :",p1,"   p2 :",p2,p1+p2

  x1 = xsize
  x2 = xsize
  y1 = xsize*real(nlayers1,PR)/real(ppd1,PR)
  y2 = xsize*real(nlayers2,PR)/real(ppd2,PR)
  ywall1 = xsize*real(nwall1,PR)/real(ppd1,PR)
  ywall2 = xsize*real(nwall2,PR)/real(ppd2,PR)

! Set periodic variables
  periodic_min(1)    = 0.0_PR
  periodic_max(1)    = xsize
  periodic_min(2)    = 0.0_PR - ywall1
  periodic_max(2)    = y1 + y2 + ywall2
  periodic_min(3)    = 0.0_PR
  periodic_max(3)    = 0.0_PR
  periodic_size(1:3) = periodic_max(1:3) - periodic_min(1:3)
  periodic_half(1:3) = 0.5_PR*periodic_size(1:3)

  write(6,*) periodic_min(1:3)
  write(6,*) periodic_max(1:3)
  write(6,*) periodic_size(1:3)
  write(6,*) x1,x2,y1,y2

! Initialise some variables using parameters
  call initialize_seren_variables_1

! Checking compiler flags and parameter values
  call sanitycheck

  rscale = 1.0_DP
  mscale = 1.0_DP

! Set masses of SPH particles in both regions
  mp1 = rho1*x1*(y1 + ywall1)/real(p1,PR)
  mp2 = rho2*x2*(y2 + ywall2)/real(p2,PR)

  allocate(ptype(1:ptot))
  ptype(1:ptot) = GASID
  pgas = ptot


! Perturbation mode 1 (Velocity perturbation)
! ============================================================================
  if (pertmode == 1) then

     ! Allocate memory for particle data
     call allocate_memory

     ! Bottom layer of particles
     ! -----------------------------------------------------------------------
     p = 0
     do j=1,nlayers1 + nwall1
        do i=1,ppd1
           p = p + 1
           parray(1,p) = xsize*(real(i,PR)-0.5_PR)/real(ppd1) + periodic_min(1)
           parray(2,p) = xsize*(real(j)-0.5_PR)/real(ppd1) + periodic_min(2)
           parray(MASS,p) = mp1
           parray(SMOO,p) = 1.0_PR
           v(1:VDIM,p)    = 0.0_PR
           if (j <= nwall1) ptype(p) = BOUNDARYID
        end do
     end do
     
     ! Top layer of particles
     ! -----------------------------------------------------------------------
     do j=1,nlayers2 + nwall2
        do i=1,ppd2
           p = p + 1
           parray(1,p) = xsize*(real(i,PR)-0.5_PR)/real(ppd2) + periodic_min(1)
           parray(2,p) = xsize*(real(j,PR)-0.5_PR)/real(ppd2) &
                & + periodic_min(2) + y1 + ywall1
           parray(MASS,p) = mp2
           parray(SMOO,p) = 1.0_PR
           v(1:VDIM,p)    = 0.0_PR
           if (j > nlayers2) ptype(p) = BOUNDARYID
        end do
     end do
     
     ! Set velocity perturbation
     ! -----------------------------------------------------------------------
     do p=1,ptot
        
        ! Abel one
        !if (parray(2,p) > 0.3_PR .and. parray(2,p) < 0.7_PR) then
        !   v(2,p) = amp*(1.0_PR + cos(8.0_PR*PI*(parray(1,p) + 0.25)))*&
        !        &(1.0_PR + cos(2.0_PR*PI/0.4_PR*(parray(2,p) - 0.5)))/4.0_PR
        !end if

        ! Springel one
        v(2,p) = amp*(1.0_PR - cos(2.0_PR*PI*(parray(1,p))/lambda))*&
             &(1.0_PR - cos(2.0_PR*PI*(parray(2,p))/(y1 + y2)))

        ! Simple one
        ! if (abs(parray(2,p) - y1) < 0.025) then
        !    v(2,p) = -amp*cos(2.*PI*(parray(1,p))/lambda)
        ! end if
     end do



! Perturbation mode 2 (Boundary perturbation)
! ============================================================================
  else if (pertmode == 2) then

     ! Allocate memory for particle data
     call allocate_memory

     ! Set properties both sides of density boundary
     ! -----------------------------------------------------------------------
     p = 0
     do j=1,(nlayers1 + nlayers2)
        do i=1,ppd1
           p = p + 1
           parray(1,p) = xsize*(real(i) - 0.5_PR)/real(ppd1) - 0.5*xsize
           parray(2,p) = xsize*(real(j) - 0.5_PR)/real(ppd1)
           parray(SMOO,p) = 1.0_PR
           v(1:VDIM,p) = 0.0_PR
           if (parray(2,p) - y1 > -amp*cos(2.0_PR*PI*(parray(1,p))/lambda)) then
              parray(MASS,p) = mp2
           else
              parray(MASS,p) = mp1
           end if
        end do
     end do 


! Load from files
! ============================================================================
  else if (pertmode == 3) then

     ! Reading in formatted data file (DRAGON format)
     call read_data(in_file1,in_file1_form)
     
     ! Set velocity perturbation
     ! -----------------------------------------------------------------------
     do p=1,ptot
        v(2,p) = amp*(1.0_PR - cos(2.0_PR*PI*(parray(1,p))/lambda))*&
             &(1.0_PR - cos(2.0_PR*PI*(parray(2,p))/(y1 + y2)))
     end do
 
  end if
! ============================================================================

! Sort particles in correct type order in memory
  call sort_particle_types(ptype(1:ptot))

  call types


! Setting up scaling units for simulation
  call units

! Converting to dimensionless code units
  call convert_to_code_units

! Initialising kernel tables
  call kernel

! Initialize variables
  call initialize_variables_1

! Build and stock tree
  nbuild = nsteps
  nstock = nsteps
#ifdef BINARY_TREE
  nskeleton = nsteps
#endif
  call tree_update

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.NOT. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.NOT. restart) call h_guess
#endif

! Calculating initial SPH quantities
  call sph_update

! Calculate temperature required to ensure constant pressure
  do p=1,ptot
     press(p) = Press1 + acc_grav*(parray(2,p) - y1)*rho(p)
     temp(p)  = press(p) / rho(p) / Pconst
     u(p)     = Pconst*temp(p)/(gamma - 1.0_PR)
  end do

! Calculating initial forces/accelerations 
  call sph_hydro_forces

! Initialize other variables
  call initialize_variables_2

! Write everything to file 
  call write_data(out_file,out_file_form)
#ifdef DEBUG_PLOT_DATA
  rzero(1:NDIM) = 0.0_PR
  call write_data_debug("ICRT.debug.dat",rzero)
#endif

  deallocate(ptype)

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_RT
