! SPH_INTEGRATE.F90
! C. P. Batty & D. A. Hubber - 10/1/2007
! Main SPH simulation integration subroutine.  
! Controls calls to all routines inside the main SPH integration loop.  
! 1.  - Builds or updates the tree
! 2.  - Computes new smoothing lengths, neighbour lists and densities 
!       for all particles
! 3.  - Calculates thermal properties of all SPH particles
! 4.  - Calculates all hydro forces on all SPH particles
! 5.  - Calculates all gravitational forces on all SPH particles
! 6.  - Calculates all gravitational forces on all sink particles
! 7.  - Writes output files 
! 8.  - Calculates new timesteps for all SPH particles
! 9.  - Checks all neighbour timesteps are comparable
! 10. - Updates thermal properties due to polytropic cooling (if required)
! 11. - Advances the positions and velocities of all particles
! 12. - Reduce the timestep of particles with low-timestep neighbours
! 13. - Searches for new sinks, and accretes particles to existing sinks
! 14. - Removes any escaping, inconsequential particles
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_integrate
  use particle_module
  use sink_module
  use time_module
  implicit none

  integer :: p
#if defined(SINKS)
  integer :: s
#endif

  debug2("Performing next SPH integration step [sph_integrate.F90]")

! Updating tree properties
  call tree_update

! Calculate new smoothing lengths, densities and other SPH properties
  call sph_update

! Calculating thermal properties for all particles
#if defined(HYDRO)
  call update_thermal_properties
#endif

! Zero acceleration array of all active particles here (for now)
  do p=1,ptot
     if (accdo(p)) a(1:VDIM,p) = 0.0_PR
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     if (accdo(p)) a_grav(1:NDIM,p) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(HYDRO)
     if (accdo(p)) a_hydro(1:NDIM,p) = 0.0_PR
#endif
  end do

! Zero relevant star arrays
#if defined(SINKS)
  do s=1,stot
     if (accdo_sinks) sink(s)%a(1:NDIM) = 0.0_DP
  end do
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
#if defined(SINKS)
  call sph_sink_forces
#endif

! Output diagnostics and/or snapshot files
  call sph_output

! Calculate new timesteps for particles
  call sph_timesteps

! Check and flag if any active particles neighbour particles with 
! relatively long timesteps.  
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
  call check_neighbour_timesteps
#endif

! Update radiative cooling terms
#if defined(HYDRO) && defined(RAD_WS)
  call rad_ws_update
#endif

! Integration scheme to advance particle properties.
  call sph_advance
#if defined(SINKS)
  call sink_advance
#endif

! Reduce the timesteps of any particles with short-timestep neighbours 
! to new 'safe' value when the timesteps are properly synchronised.
#if defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(IMMEDIATE_TIMESTEP_REDUCTION)
  call reduce_timesteps
#endif

! Search for new sinks and accrete particles into existing sinks
#if defined(SINKS) && defined(GRAVITY) && defined(HYDRO) && !defined(NO_ACCRETION)
  call sink_update
#endif

! Remove any outlying particles depending on chosen criteria.
#if defined(REMOVE_OUTLIERS)
  call remove_outlying_particles
#endif

  return
END SUBROUTINE sph_integrate
