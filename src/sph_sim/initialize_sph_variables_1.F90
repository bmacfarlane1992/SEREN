! INITIALIZE_SPH_VARIABLES_1.F90
! D. A. Hubber - 1/10/2007
! Sets values for particular variables that need to be initialized BEFORE the 
! first force calculations in sph_setup.F90.  Other variables (after the first 
! force calculation) are initialized in initialize_sph_variables_2.F90.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_sph_variables_1
  use particle_module
  use hydro_module
  use periodic_module
  use time_module
  use scaling_module
  use filename_module
  use sink_module
  use neighbour_module
  use type_module
  use Nbody_module
  use timing_module
  implicit none

  integer :: p             ! Particle counter

  debug2("Initializing variables [initialize_sph_variables_1.F90]")


! Set accdo variables to ensure all particles are included in 
! sph and force subroutines
! ----------------------------------------------------------------------------
  do p=1,ptot
     accdo(p)    = .true.
     nlevel(p)   = 0
     nlast(p)    = n
     laststep(p) = 0.0_DP
  end do
  timestep = 0.0_DP
#if defined(SINKS)
  accdo_sinks    = .true.
  nlevel_sinks   = 0
  laststep_sinks = 0.0_DP
#endif


! Artificial viscosity and conductivity variables
! ----------------------------------------------------------------------------
#if defined(VISC_TD)
  talpha(1:ptot)     = alpha_min
  talpha_old(1:ptot) = alpha_min
  dalpha_dt(1:ptot)  = 0.0_PR
#endif
#if defined(VISC_BALSARA)
  balsara(1:ptot) = 0.0_PR
#endif
#if defined(VISC_PATTERN_REC)
  pattrec(1:ptot) = 0.0_PR
#endif


! Grad-h SPH correction terms
! ----------------------------------------------------------------------------
#if defined(GRAD_H_SPH)
  do p=1,ptot
     if (rho(p) == 0.0_PR) rho(p) = 1.0_PR
     omega(p) = 1.0_PR
#if defined(GRAVITY)
     parray(ZETA,p) = 0.0_PR
#endif
  end do
#endif


! Zero acceleration arrays
! ----------------------------------------------------------------------------
  a(1:VDIM,1:ptot) = 0.0_PR
#if defined(GRAVITY) && !defined(GEOMETRIC_MAC)
  agravmag(1:ptot) = 0.0_PR
#if defined(SINKS)
  sink(1:stot)%agravmag = 0.0_PR
#endif
#endif
#if defined(DEBUG_FORCES) && defined(GRAVITY)
  a_grav(1:VDIM,1:ptot) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(HYDRO)
  a_hydro(1:VDIM,1:ptot) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(IDEAL_MHD)
  a_mag(1:VDIM,1:ptot) = 0.0_PR
#endif


! Misc variables
! ----------------------------------------------------------------------------
  div_v(1:ptot) = 0.0_PR
#if defined(GRAVITY)
  gpot(1:ptot) = 0.0_PR
#endif
#if defined(DIV_A)
  div_a(1:ptot) = 0.0_PR
#endif
#if defined(SINKS) && defined(GRAVITY)
  ispotmin(1:ptot) = .false.
#endif
#if defined(IONIZING_UV_RADIATION)
  temp_min(1:ptot) = 0.0_PR
#endif
#if defined(TIMING)
  ngravcomp = 0_ILP
  nhydrocomp = 0_ILP
  nsphcomp = 0_ILP
#endif

! Tree-building, stocking and ionization time variables
  nbuild = nsteps
  nstock = nsteps
  nionize = nsteps
  nionall = nsteps
#if defined(BINARY_TREE)
  nskeleton = nsteps
#endif

  return
END SUBROUTINE initialize_sph_variables_1
