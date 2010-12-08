! SEREN_SETUP.F90
! C. P. Batty, D. A. Hubber, A.McLeod & A. P. Whitworth - 8/12/2006
! Set-up all Seren routines.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE seren_setup
  use interface_module, only : paramstore,read_parameters
  use definitions
  use filename_module
  use seren_sim_module
  implicit none

  character(len=256) :: store_file   ! parameters output

! Read in and process command line arguments
  call read_arguments

! Set default parameters
  call default_parameters

! Reading parameter file
  call read_parameters(param_file)

! Initialise some variables using parameters before IC file is read in
  call initialize_seren_variables_1

! Checking compiler flags and parameter values before proceeding further
  call sanitycheck

! Setting up scaling units for simulation
  call units

! Reading in formatted data file
  call read_data(in_file,in_file_form)

! Writing compiler flags and parameters to file
  store_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".params"
  call paramstore(store_file)

! Setting up for different particle types
  call types

! Converting to dimensionless code units
  call convert_to_code_units

! Initialising kernel tables
  call kernel

! Read in opacity tables for radiation transport
#if defined(HYDRO) && defined(GRAVITY) && defined(RAD_WS)
  call read_eos
#endif

! Create Ewald correction table
#if defined(GRAVITY) && defined(EWALD)
  call ewald_init
#endif

! Calculate COM of system, and if required, change to COM
  call COM

! Initialise all other particles arrays
  call initialize_seren_variables_2

  return
END SUBROUTINE seren_setup
