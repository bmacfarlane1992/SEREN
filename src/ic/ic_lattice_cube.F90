! IC_LATTICE_CUBE.F90
! D. A. Hubber - 21/9/2007
! Generates a regularly spaced lattice of particles in 1, 2 or 3 dimensions.
! Does not calculate any further properties of particles (assuming the 
! outputted lattice will be fed into another initial conditions program).
! ppd            : No. of particles per dimension in cube
! length         : Side length of cube
! out_file       : Name of output file
! out_file_form  : Output file format
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_lattice_cube
  use interface_module, only : paramerror,write_data
  use particle_module
  use hydro_module
  use filename_module, only : out_file_form
  use scaling_module, only : rscale, mscale
  use type_module
  implicit none

  character(len=256) :: out_file  ! Name of output file
  integer :: i                    ! x-dimension counter
  integer :: j                    ! y-dimension counter
  integer :: k                    ! z-dimension counter
  integer :: p                    ! Particle counter
  integer :: ppd                  ! Particle per dimension
  real(kind=PR) :: length         ! Side-length of line/square/cube

#ifndef DIMENSIONLESS
  write(6,*) "Compiler flag error : Only works with DIMENSIONLESS flag on"
  stop
#endif

! Set parameters to default values
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "       ic_lattice_cube        "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Reading in input parameters"
  write(6,*) "Number of particles per dimension "
  read(5,*) ppd
  write(6,*) "Side length : "
  read(5,*) length
  write(6,*) "Name of output file : "
  read(5,*) out_file
  write(6,*) "Format of output file : "
  read(5,*) out_file_form

! Sanity check for parameters
  if (ppd <= 0) call paramerror("Invalid value of ppd")
  if (length <= 0.0_PR) call paramerror("Zero or negative value of length")

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call units

  ptot = ppd**(NDIM)
  write(6,*) "Number of particles : ",ptot

! Allocate memory for particle data
  call allocate_memory
  rscale = 1.0_PR
  mscale = 1.0_PR
  pgas   = ptot

! If NDIM=1, generate line of equally spaced particles
! ----------------------------------------------------------------------------
#if NDIM==1
  do p=1,ptot
     parray(1,p) = length * (real(p,PR) - 0.5_PR) / real(ptot,PR)
     parray(MASS,p) = 1.0_PR / real(ptot,PR)
  end do

! Else if NDIM=2, generate equally-spaced grid of particles
! ----------------------------------------------------------------------------
#elif NDIM==2
  p = 0
  do i=1,ppd
     do j=1,ppd
        p = p + 1
        parray(1,p) = length * (real(i,PR) - 0.5_PR) / real(ppd,PR)
        parray(2,p) = length * (real(j,PR) - 0.5_PR) / real(ppd,PR) 
     end do
  end do

! Else if NDIM=3, generate equally-spaced cubic lattice of particles
! ----------------------------------------------------------------------------
#elif NDIM==3
  p = 0
  do i=1,ppd
     do j=1,ppd
        do k=1,ppd
           p = p + 1
           parray(1,p) = length * (real(i,PR) - 0.5_PR) / real(ppd,PR)
           parray(2,p) = length * (real(j,PR) - 0.5_PR) / real(ppd,PR) 
           parray(3,p) = length * (real(k,PR) - 0.5_PR) / real(ppd,PR) 
        end do
     end do
  end do
! ----------------------------------------------------------------------------
#endif

! Set all other quantities
   do p=1,ptot
      parray(MASS,p) = 1.0_PR / real(ptot,PR)
      parray(SMOO,p) = 0.0_PR
      v(1:NDIM,p)    = 0.0_PR
      rho(p)         = 0.0_PR
      temp(p)        = 1.0_PR
   end do

! Write data to file
  call write_data(out_file,out_file_form)

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_lattice_cube
