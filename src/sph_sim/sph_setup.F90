! SPH_SETUP.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Performs initialisation routines to prepare for main SPH simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_setup
  use time_module
  use filename_module, only : restart,in_file
  implicit none

  debug1("Setting up the SPH simulation [sph_setup.F90]")

! Setting up array-limits for different particle types
  call types

! Set-up HEALPix rays
#if defined(HEALPIX)
  call initialize_HP_sources
#endif

! Initialize certain variables before first force calculation
  call initialize_sph_variables_1

! after read in and initialization of data do the analysis (doesn't the info below)
#ifdef ANALYSE
#ifdef ANALYSE_PICKUP_FRAMES
   write(*,*) "Picking up specific frames ..."
   call analyse_disc(in_file)
 stop
#endif
#ifdef ANALYSE_CENTRAL_REGION
   write(*,*) "Analysing central region ..."
   call analyse_central_region(in_file)
 stop
#endif
#endif

! Build and stock trees for first time
  call tree_update

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.not. restart) call h_guess
#endif

! Calculating initial SPH quantities (h, neibs, rho, etc.)
  call sph_update

! Initialize all thermal properties depending on options used
#if defined(HYDRO)
  call initialize_thermal_properties
#endif

! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Calculate gravitational forces on all sink particles
#if defined(SINKS) && defined(GRAVITY)
  call sph_sink_forces
#endif

! Calculate initial diagnostic quantities
#if defined(DEBUG_DIAGNOSTICS)
  call diagnostics
#endif

! Initialize other key variables (after initial force calculation)
  call initialize_sph_variables_2

! after read in and initialization of data do the analysis (needes the above info)
#ifdef ANALYSE
#ifdef ANALYSE_DISC
   write(*,*) "Analysing disc around first sink ..."
   call analyse_disc(in_file)
 stop
#endif
#endif

  return
END SUBROUTINE sph_setup
