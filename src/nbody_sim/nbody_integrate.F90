! NBODY_INTEGRATE.F90
! D. A. Hubber - 13/9/2010
! Integration step
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_integrate
  use Nbody_module
  implicit none

  debug2("Performing next N-body integration step [nbody_integrate.F90]")

! Calculate new timesteps for all N-body particles
  call nbody_timesteps

! Integration scheme to advance particle positions and velocities
  call nbody_advance

! Compute force and force derivatives at beginning of timestep
  call nbody_grav_forces

! Calculate higher-order derivatives and correction forces for previous steps
  call nbody_correction_terms

! Output diagnostics and/or snapshot files
  call nbody_output

  return
END SUBROUTINE nbody_integrate
