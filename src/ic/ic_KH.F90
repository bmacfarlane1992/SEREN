! IC_KH.F90
! D. A. Hubber - 15/3/2008
! Creates initial conditions for the Kelvin-Helmholtz instability test.
! Assumes we either load in unity sheets/cubes.
! out_file        : Output file name
! out_file_form   : Output file format
! in_file1        : 1st input file name
! in_file1_form   : 1st input file format
! in_file2        : 2nd input file name
! in_file2_form   : 2nd input file format
! p1, p2          : No. of particles in file 1, 2
! n1, n2          : No. of replicas for layers 1, 2
! vx1, vx2        : x-velocity of layer 1, 2
! rho1, rho2      : Density of layer 1, 2
! Press1, Press2  : Pressure of layer 1, 2
! x1, x2          : x-size of input file 1, 2
! y1, y2          : y-size of input file 1, 2
! z1, z2          : z-size of input file 1, 2 (redundant?)
! h_fac           : 'grad-h' h_fac
! pertmode        : Mode of perturbation
! amp             : Amplitude of perturbation
! lambda          : Wavelength of perturbation
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_KH
  use interface_module, only : check_boundary_conditions,write_data
  use particle_module
  use hydro_module
  use filename_module
  use time_module
  use periodic_module
  use type_module
  use neighbour_module, only : h_fac
  implicit none

  character(len=256) :: outfile          ! KH IC file
  character(len=50) :: pertmode          ! Perturbation type
  integer :: i                           ! Aux. counter
  integer :: k                           ! Dimension counter
  integer :: p                           ! Particle counter
  integer :: p1                          ! No. of particles in ??
  integer :: p2                          ! No. of particles in ??
  integer :: nx1                         ! No. of replicated cubes in ??
  integer :: nx2                         ! No. of replicated cubes in ??
  integer :: ny1                         ! No. of replicated cubes in ??
  integer :: ny2                         ! No. of replicated cubes in ??
  real(kind=PR) :: amp                   ! Amplitude of velocity perturbation
  real(kind=PR) :: kpert
  real(kind=PR) :: lambda                ! Wavelength of perturbation
  real(kind=PR) :: offset                ! ..
  real(kind=PR) :: mp1                   ! Mass of 'bottom' particles
  real(kind=PR) :: mp2                   ! Mass of 'top' particles
  real(kind=PR) :: Press1                ! 'Bottom' pressure
  real(kind=PR) :: Press2                ! 'Top' pressure
  real(kind=PR) :: rho1                  ! 'Bottom' density
  real(kind=PR) :: rho2                  ! 'Top' density
  real(kind=PR) :: sigmapert             ! ..
  real(kind=PR) :: T1                    ! 'Bottom' temperature
  real(kind=PR) :: T2                    ! 'Top' temperature
  real(kind=PR) :: volume                ! Volume of ..
  real(kind=PR) :: vx1                   ! 'Bottom' x-velocity 
  real(kind=PR) :: vx2                   ! 'Top' x-velocity
  real(kind=PR), allocatable :: r1(:,:)  ! Positions of 'bottom' particles
  real(kind=PR), allocatable :: r2(:,:)  ! Positions of 'top' particles
#if defined(DEBUG_PLOT_DATA)
  real(kind=PR) :: rcentre(1:NDIM)       ! Origin
#endif

! Set parameters to default values
  call default_parameters

! Read in parameter file 
  write(6,*) "Opening KHparams file"
  open(unit=1,file="KHparams.dat",status='old')
  read(1,*) outfile
  read(1,*) out_file_form
  read(1,*) nx1, nx2
  read(1,*) ny1, ny2
  read(1,*) vx1, vx2
  read(1,*) rho1, rho2
  read(1,*) Press1, Press2
  read(1,*) periodic_min(1), periodic_max(1)
  read(1,*) periodic_min(2), periodic_max(2)
  read(1,*) h_fac
  read(1,*) gamma
  read(1,*) pertmode
  read(1,*) amp
  read(1,*) lambda
  close(1)

! Initialise some variables using parameters
  call initialize_seren_variables_1

! Checking compiler flags and parameter values
  call sanitycheck

  p1 = nx1*ny1
  p2 = nx2*ny2
  allocate(r1(1:NDIM,1:p1))
  allocate(r2(1:NDIM,1:p2))

  pgas = p1 + p2
  ptot = pgas


! Calculate temperature of gas
  T1 = Press1 / rho1
  T2 = Press2 / rho2
  write(6,*) "T1 :",T1
  write(6,*) "T2 :",T2

  volume = 0.5_PR*periodic_size(1)*periodic_size(2)
  mp1 = rho1*volume/real(p1,PR)
  mp2 = rho2*volume/real(p2,PR)
  kpert = 2.0_PR*PI/lambda

! Output parameters to screen for verification
  write(6,*) "outfile :",trim(outfile),"   ",trim(out_file_form)
  write(6,*) "Number of particles :",p1,p2
  write(6,*) "vx1 :",vx1,"    vx2 ",vx2
  write(6,*) "rho1 :",rho1,"    rho2 :",rho2
  write(6,*) "Press1 :",Press1,"    Press2 :",Press2
  write(6,*) "amp :",amp,"  lambda :",lambda
  write(6,*) "Closing KHparams.dat file"
  write(6,*) "mp1 : ",mp1,"    mp2 : ",mp2

  call allocate_memory

  call add_uniform_2d_lattice(periodic_min(1),periodic_max(1),&
       &periodic_min(2),periodic_min(2) + periodic_half(2),r1,p1,nx1,ny1)
  call add_uniform_2d_lattice(periodic_min(1),periodic_max(1),&
       &periodic_min(2) + periodic_half(2),periodic_max(2),r2,p2,nx2,ny2)



! First write LHS of shock tube info 
! ----------------------------------------------------------------------------
  do p=1,p1
     parray(1:NDIM,p) = r1(1:NDIM,p)
     parray(2,p) = parray(2,p) + 0.5_PR*periodic_half(2)
     parray(MASS,p)   = mp1
     parray(SMOO,p)   = 1.0_PR
     v(1:NDIM,p)      = 0.0_PR
     v(1,p)           = vx1
     temp(p)          = T1
#if defined(ENTROPIC_FUNCTION)
     Aent(p)          = Press1/rho1**gamma
#endif
     call check_boundary_conditions(parray(1:NDIM,p),v(1:VDIM,p))
  end do


! Now write RHS of shock tube info 
! ----------------------------------------------------------------------------
  do i=1,p2
     p = p1 + i
     parray(1:NDIM,p) = r2(1:NDIM,i)
     parray(2,p) = parray(2,p) + 0.5_PR*periodic_half(2)
     parray(MASS,p)   = mp2
     parray(SMOO,p)   = 1.0_PR
     v(1:NDIM,p)      = 0.0_PR
     v(1,p)           = vx2
     temp(p)          = T2
#if defined(ENTROPIC_FUNCTION)
     Aent(p)          = Press1/rho2**gamma
#endif
     call check_boundary_conditions(parray(1:NDIM,p),v(1:VDIM,p))
  end do


! Add velocity perturbation here
! ----------------------------------------------------------------------------
  if (pertmode == 'price2008') then
     do p=1,ptot
        if (parray(2,p) > periodic_max(2)) &
             &parray(2,p) = parray(2,p) - periodic_size(2)
        if (abs(parray(2,p) - 0.25_PR) < 0.025_PR) then
           v(2,p) = amp*sin(-2.0_PR*PI*(parray(1,p) + 0.5_PR)/lambda)
        else if (abs(parray(2,p) + 0.25_PR) < 0.025_PR) then
           v(2,p) = amp*sin(2.0_PR*PI*(parray(1,p) + 0.5_PR)/lambda)
        end if
     end do
  else if (pertmode == 'arepo') then
     sigmapert = 0.05_PR/sqrt(2.0_PR)
     do p=1,ptot
        v(2,p) = amp*sin(2.0_PR*PI*parray(1,p)/lambda)*&
             &(exp(-(parray(2,p) + 0.25_PR)**2/2.0_PR/sigmapert**2) +& 
             &exp(-(parray(2,p) - 0.25_PR)**2/2.0_PR/sigmapert**2))
     end do
  else
     stop 'No valid velocity perturbation selected'
  end if

! Setting up for different particle types
  call types
  write(6,*) "pgas :",pgas
  write(6,*) "ptot :",ptot

! Setting up scaling units for simulation
  call units

! Initialising kernel tables
  call kernel

! Initialize variables
  call initialize_seren_variables_2
  call initialize_sph_variables_1

! Build and stock tree
  nbuild = nsteps
  nstock = nsteps
#ifdef BINARY_TREE
  nskeleton = nsteps
#endif
  call tree_update

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.not. restart) call h_guess
#endif

! Calculating initial SPH quantities
  call sph_update

! Calculate temperature required to ensure constant pressure.
  do p=1,ptot
     temp(p)  = Press1 / rho(p) / Pconst
     u(p)     = Pconst*temp(p)/(gamma - 1.0_PR)
     press(p) = (gamma - 1.0_PR)*rho(p)*u(p)
     sound(p) = sqrt(press(p)/rho(p))
  end do

! Calculate hydro forces on all SPH particles
  call sph_hydro_forces

! Initialize other variables
  call initialize_sph_variables_2

! Write everything to file 
  call write_data(outfile,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  rcentre(1:NDIM) = 0.0_PR
  call write_data_debug("ICKH.debug.dat",rcentre(1:NDIM))
#endif

  deallocate(r1)
  deallocate(r2)

  call clean_up

  stop
END PROGRAM ic_KH
