! SANITYCHECK.F90
! C. P. Batty & D. A. Hubber - 1/4/2007
! Checks all compiler flags and all input parameters for any conflicts, 
! bad or illegal values etc.  Will cause the program to stop during 
! runtime and print an error message to screen, rather than compile time. 
! ============================================================================ 

#include "macros.h"

! ============================================================================ 
SUBROUTINE sanitycheck
  use interface_module, only : comperror,paramerror
  use particle_module
  use periodic_module
  use scaling_module
  use neighbour_module
  use hydro_module
  use time_module
  use sink_module
  use tree_module
  use Nbody_module
  use filename_module
  use HP_module
  implicit none

  debug2("Checking compiler flags and parameter values [sanitycheck.F90]")

! Check main compiler flags 
! ----------------------------------------------------------------------------
#if NDIM < 1 || NDIM > 3
  call comperror("NDIM out of range")
#endif

#if defined(PERIODIC_X) || defined(PERIODIC_Y) || defined(PERIODIC_Z)
#if !defined(PERIODIC)
  call comperror("PERIODIC flag not on with X/Y/Z")
#endif
#endif

#if !defined(SPH_SIMULATION) && !defined(NBODY_SPH_SIMULATION) && !defined(NBODY_SIMULATION)
  call comperror("No simulation mode activated")
#endif

#if defined(SPH_SIMULATION) || defined(NBODY_SPH_SIMULATION)
#if !defined(EULER) && !defined(RUNGE_KUTTA) && !defined(LEAPFROG_KDK) && !defined(LEAPFROG_DKD) && !defined(PREDICTOR_CORRECTOR)
  call comperror("No valid SPH integration scheme selected")
#endif

#if defined(EULER)
  call comperror("WARNING : The Euler integration scheme is a RUBBISH &
       &integration scheme and should not be used unless you want to find &
       &out how rubbish it is, or you are having a laugh.  If you wish to &
       &use this option, then comment out this error message in &
       &setup/sanity_check.F90.")
#endif

#if !defined(M4_KERNEL) && !defined(QUINTIC_KERNEL) && !defined(GAUSSIAN_3H_KERNEL)
  call comperror("No valid SPH kernel function selected")
#endif

#if defined(IONIZING_UV_RADIATION)
#if !defined(ISOTHERMAL) && !defined(BAROTROPIC)
  call comperror("HII flag on without an allowed EOS")
#endif
#endif

#if defined(RIEMANN_SOLVER) && defined(INTERNAL_ENERGY) && defined(RAD_WS)
!  call comperror("Riemann solver does not work with RAD_WS option")
#endif

#if defined(IDEAL_MHD)
  call comperror("IDEAL_MHD under development; currently disabled")
#endif

#if defined(IDEAL_MHD) && !defined(HYDRO)
  call comperror("IDEAL_MHD flag on and HYDRO flag off")
#endif

#if defined(IDEAL_MHD) && !defined(GRAD_H_SPH)
  call comperror("IDEAL_MHD only functions for GRAD_H_SPH currently")
#endif

#if defined(RAD_WS) && !defined(GRAVITY)
  call comperror("RAD_WS flag on and GRAVITY flag off")
#endif

#if defined(RAD_WS) && !defined(INTERNAL_ENERGY)
  call comperror("RAD_WS flag on and INTERNAL_ENERGY flag off")
#endif

#if !defined(GRAVITY) && !defined(HYDRO)
  call comperror("GRAVITY and HYDRO flags both off")
#endif

#endif


#if defined(GRAVITY) && NDIM == 1
  call comperror("NDIM == 1 for GRAVITY flag")
#endif

#if defined(EWALD) && !defined(PERIODIC)
  call comperror("PERIODIC flag off, EWALD flag on")
#endif

#if defined(BH_TREE) && NDIM == 1
  call comperror("NDIM = 1 for BH tree")
#endif

#if defined(BH_TREE) && defined(GRAVITY)
#if !defined(GEOMETRIC_MAC) && !defined(GADGET_MAC) && !defined(GADGET2_MAC) && !defined(EIGEN_MAC)
  call comperror("No recognised MAC for BH_TREE and GRAVITY")
#endif
#endif

#if defined(BINARY_TREE)
  call comperror("Binary tree not completely implemented; currently disabled")
#endif

#if defined(BINARY_TREE) && LEAFMAX == 1
  call comperror("LEAFMAX == 1 does not work for the binary tree")
#endif

#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
!#if !defined(SINKS)
!  call comperror("SINKS flag off and NBODY flag on")
!#endif
#endif

#if defined(SINK_GRAVITY) && !defined(SINKS)
  call comperror("SINKS flag off and SINK_GRAVITY flag on")
#endif

#if defined(SINK_GRAVITY) && defined(GRAVITY)
  call comperror("SINK_GRAVITY and GRAVITY flags on")
#endif

#if defined(BH_TREE) && defined(REORDER_TREE)
  call comperror("REORDER not working correctly with TREE or ALL options")
#endif

#if !defined(ADAPTIVE_TIMESTEP_LEVELS) && !defined(RESTRICTED_TIMESTEP_LEVELS) && !defined(FIXED_TIMESTEP_LEVELS)
  call comperror("No TIMESTEP option selected")
#endif

#if defined(LEAPFROG_DKD) && defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(IMMEDIATE_TIMESTEP_REDUCTION)
  call comperror("Combination of INTEGRATOR=LFDKD and CHECK_NEIB_TIMESTEP=2 &
       &potentially dangerous due to secular error increase.  If you wish to &
       &use this option, then comment out this error message in &
       &setup/sanity_check.F90.")
#endif

#if defined(SINKS) && defined(SINK_REMOVE_ANGMOM)
!  call comperror("SINK_REMOVE_ANGMOM under development - currently disabled")
#endif


! Now check all input parameters 
! ----------------------------------------------------------------------------
  if (rseed < 1) call paramerror("rseed < 1")
  if (sph_endtime < 0.0_DP) call paramerror("sph_endtime < 0")
  if (firstsnap < 0.0_DP) call paramerror("firstsnap < 0")
  if (snaptime <= 0.0_DP) call paramerror("snaptime <= 0")
  if (ntempstep < 1) call paramerror("ntempstep < 1")
  if (ndiagstep < 1) call paramerror("ndiagstep < 1")
  if (nsinkstep < 1) call paramerror("nsinkstep < 1")
  if (nsnapstep < 1) call paramerror("nsnapstep < 1")
  if (courant_mult < 0.0_DP) call paramerror("courant_mult < 0")
  if (courant_mult > 1.0_DP) call paramerror("courant_mult too big (> 1)")
  if (accel_mult < 0.0_DP) call paramerror("accel_mult < 0")
  if (accel_mult > 1.0_DP) call paramerror("accel_mult too big (> 1)")
  if (sink_mult < 0.0_DP) call paramerror("sink_mult < 0")
  if (sink_mult > 1.0_DP) call paramerror("sink_mult too big (> 1)")
  if (nlevels < 1) call paramerror("nlevels < 1")
  if (nlevels > 20) call paramerror("nlevels too large (> 20)")
  if (dt_fixed < 0.0_DP) call paramerror("dt_fixed < 0")
#if defined(NBODY_SIMULATION)
  if (nbody_endtime < 0.0_DP) call paramerror("nbody_endtime < 0")
  if (nbody_timemult < 0.0_DP) call paramerror("nbody_timemult < 0")
  if (nbody_frac < 0.0_DP .OR. nbody_frac > 1.0_DP) &
       & call paramerror("nbody_frac out of range")
#endif
#if defined(SPH_SIMULATION) && defined(NBODY_SIMULATION)
  if (nbody_endtime < endtime) call paramerror("nbody_endtime < endtime")
#endif

  if (rscale <= 0.0_DP) call paramerror("rscale <= 0")
  if (mscale <= 0.0_DP) call paramerror("mscale <= 0")

#if defined(PERIODIC_X)
  if (periodic_min(1) >= periodic_max(1)) call paramerror("periodic_x")
#endif
#if defined(PERIODIC_Y)
  if (periodic_min(2) >= periodic_max(2)) call paramerror("periodic_y")
#endif
#if defined(PERIODIC_Z)
  if (periodic_min(2) >= periodic_max(2)) call paramerror("periodic_z")
#endif
  if (pp_gather < 1 .or. pp_gather > pp_limit) &
	&call paramerror("pp_gather < 1 or > pp_limit")

#if defined(MINIMUM_H)
  if (hmin < 0.0_PR) call paramerror("hmin < 0")
#endif

#if defined(ARTIFICIAL_VISCOSITY)
  if (alpha < 0.0_PR) call paramerror("alpha < 0")
  if (beta < 0.0_PR) call paramerror("beta < 0")
#if defined(VISC_TD)
  if (alpha_min < 0.0_PR) call paramerror("alpha_min < 0")
  if (alpha_min > alpha) call paramerror("alpha_min > alpha")
#endif
#endif

#if defined(BH_TREE)
  if (nbuildstep <= 0) call paramerror("nbuildstep <= 0")
  if (thetamaxsqd < 0.0_PR) call paramerror("thetamaxsqd < 0.")
  if (thetamaxsqd > 1.0_PR) call paramerror("thetamaxsqd too big (> 1)")
  if (abserror < 0.0_PR) call paramerror("abserror < 0")
  if (abserror > 1.0_PR) call paramerror("abserror > 1")
#endif

#if defined(HYDRO)
#if defined(ISOTHERMAL) || defined(BAROTROPIC) 
  if (isotemp <= 0.0_PR) call paramerror("isotemp <= 0.")
  if (mu_bar <= 0.0_PR) call paramerror("mu_bar <= 0.")
#endif
#if defined(BAROTROPIC) || defined(POLYTROPIC) 
  if (rhobary <= 0.0_PR) call paramerror("rhobary <= 0.")
  if (gamma < 1.0_PR) call paramerror("gamma < 1.")
  if (Kpoly < 0.0_PR) call paramerror("Kpoly < 0.")
#endif
#if defined(EXTERNAL_PRESSURE)
  if (Pext < 0.0_PR) call paramerror("Pext < 0.")
#endif
#endif

#if defined(REMOVE_OUTLIERS)
  if ((.not. rho_remove) .and. (.not. energy_remove) .and. (.not. rad_remove))&
       & call paramerror("No removal criteria selected for REMOVE_OUTLIERS")
#endif

#if defined(SINKS)
  if (rhosink < 0.0_PR) call paramerror("rhosink < 0")
  if (nsearchstep <= 0) call paramerror("nsearchstep <= 0")
#if !defined(FIXED_ABSOLUTE_SINKRAD)
  if (sinkrad < KERNRANGE) call paramerror("sinkrad < KERNRANGE")
#endif
#endif

#if defined(IONIZING_UV_RADIATION)
  if (nionallstep <= 0) call paramerror("nionallstep <= 0")
  if (f1 < 0.0_PR) call paramerror("f1 < 0.0")
  if (f2 < 0.0_PR) call paramerror("f1 < 0.0")
  if (f2 > 2.0_PR) call paramerror("f2 > 2.0")
  if (f3 < 0.0_PR) call paramerror("f1 < 0.0")
  if (f4 < 0.0_PR) call paramerror("f1 < 0.0")
  if (Tneut < 0.0_PR) call paramerror("Tneut < 0.0")
  if (Tneut > Tion) call paramerror("Tneut > Tion")
  if (lmax_hp > HP_LEVELS) call paramerror("lmax_hp too big : Too many HEALPix levels")
#endif

  return
END SUBROUTINE sanitycheck


! ============================================================================ 
SUBROUTINE comperror(errmsg)
  implicit none

  character(len=*), intent(in) :: errmsg

  write(6,*) "Compiler flag error : ", errmsg
  stop

  return
END SUBROUTINE comperror


! ============================================================================ 
SUBROUTINE paramerror(errmsg)
  implicit none

  character(len=*), intent(in) :: errmsg

  write(6,*) "Parameter error : ", errmsg
  stop

  return
END SUBROUTINE paramerror
