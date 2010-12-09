! IC_SEDOV.F90
! D. A. Hubber & M. Kaplan - 16/06/2008
! Prepare initial conditions for Sedov blastwave test.
! in_file         : Input file name
! in_file_form    : Input file format
! out_file        : Output file name
! out_file_form   : Output file format
! rho0            : Cloud density
! radius          : Cloud radius
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_sedov
  use interface_module, only read_data,write_data
  use particle_module
  use hydro_module
  use constant_module
  use filename_module
  use scaling_module
  use type_module
  use time_module
  use kernel_module
  use Eos_module
  use HIIregion_module
  use neighbour_module, only : pp_gather
  implicit none

  character(len=256) :: out_file          ! Name of output file
  integer :: k                            ! Dimension counter
  integer :: kern                         ! Kernel table variable
  integer :: p                            ! Particle counter
  integer :: phot                         ! No. of particles in the center
  integer :: pcold                        ! rest of the particles
  integer :: pp_out                       ! rest of the particles II
  integer, allocatable :: hotlist(:)      ! potential neighbour list
  real(kind=PR) :: center(1:NDIM)         ! Center
  real(kind=PR) :: drmag                  ! Distance
  real(kind=PR) :: mp                     ! Mass of particle p
  real(kind=PR) :: mtot                   ! Total mass
  real(kind=PR) :: radius                 ! ..
  real(kind=PR) :: rho0                   ! Density of initial sphere
  real(kind=PR) :: r_hot                  ! Radius of hot particle region
  real(kind=PR) :: r_neigh                ! Radius for nearest neighbours in 
                                          ! the center (pp_center*rmax/ptot)
  real(kind=PR) :: rp(1:NDIM)             ! Position of particle p
  real(kind=PR) :: umax                   ! ..
  real(kind=PR) :: utot_temp              ! ..
  real(kind=PR), allocatable :: utemp(:)  ! internal energy

! Reading parameter file 
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "         ic_sedovtest         "
  write(6,*) "------------------------------"
  open(unit=1,file="sedovparams.dat",status="old")
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file
  read(1,*) in_file_form
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) rho0
  read(1,*) radius
  close(1)

! Initialise some variables using parameters
  call initialize_seren_variables_1

! Setting up scaling units for simulation 
  call units

! Reading in formatted data file (DRAGON snapshot).  File should contain 
! unit mass sphere of radius 1, so no need to scale to dimensionless units. 
  call read_data(in_file,in_file_form)

! Allocate memory for list of 'hot' particles
  allocate(hotlist(1:ptot))

! Calculating kernel tables 
  call kernel

! Calculate average value of smoothing length and gas cloud radius to know 
! where to select boundary particles.
  r_hot = radius*(pp_gather/real(ptot,PR))**(ONETHIRD)

! Calculate correct masses by scaling input parameters
  rho0  = rho0 / rhoscale
  mtot  = (4.0_PR*PI*rho0*radius**3)/3.0_PR
  mp    = mtot / real(ptot,PR)
  pcold = 0
  phot  = 0

! Set all particle type variables
  pgas      = ptot
  picm      = 0
  pboundary = 0
  pcdm      = 0
  call types
  write(6,*) "rho0 : ",rho0,"   mtot : ",mtot,"   mp : ",mp

! Set some basic particle properties
  do p=1,ptot
     v(1:VDIM,p)      = 0.0_PR
     parray(MASS,p)   = mp
     parray(1:NDIM,p) = parray(1:NDIM,p)*radius
  end do

! Initialize variables
  call initialize_seren_variables_1
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
  if (.NOT. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.NOT. restart) call h_guess
#endif

! Calculating initial SPH quantities
  call sph_update


! Determine which particles are in the core
! ----------------------------------------------------------------------------
  umax = 0.0_PR
  utot_temp = 0.0_PR
  do p=1,ptot
     rp(1:NDIM) = parray(1:NDIM,p)
     drmag = sqrt(rp(1)*rp(1) + rp(2)*rp(2) + rp(3)*rp(3))
     if (drmag < r_hot) then
        hotlist(p) = 1 
        phot       = phot + 1
        kern       = int(HALFKERNTOT*drmag/(INVKERNRANGE*r_hot))
        kern       = min(kern,KERNTOT)
        u(p)       = w0(kern)*parray(MASS,p)
        utot_temp  = utot_temp + u(p)
        umax       = max(umax,u(p))
     else
        hotlist(p) = 0
        pcold      = pcold + 1
     end if
  end do
  umax = umax / utot_temp


! Correction for the exact values of particle numbers
! ----------------------------------------------------------------------------
  do p=1,ptot
     if (hotlist(p) == 1) then
        u(p)    = u(p) / utot_temp / parray(MASS,p)
        temp(p) = (gamma - 1.0_PR)*u(p)/Pconst
     else
        u(p)    = 1.E-6_PR*umax/parray(MASS,p)
        temp(p) = (gamma - 1.0_PR)*u(p)/Pconst
     end if
  end do

! Calculating initial forces/accelerations 
  call sph_hydro_forces

! Initialize other variables
  call initialize_sph_variables_2
  
! Write data to file
  call write_data(out_file,out_file_form)
#ifdef DEBUG_PLOT_DATA
  rzero(1:3) = 0.0_PR
  call write_data_debug("ICSEDOV.debug.dat",rzero)
#endif

  write(6,*) "ptot        : ",ptot
  write(6,*) "pp_gather   : ",pp_gather
  write(6,*) "mtot        : ",mtot,mp
  write(6,*) "r_hot       : ",r_hot
  write(6,*) "phot, pcold : ",phot, pcold
  write(6,*) "total en    : ", sum(parray(MASS,1:ptot)*u(1:ptot))

! Clean up memory
  deallocate(hotlist)
  call clean_up  

  stop
END PROGRAM ic_sedov
