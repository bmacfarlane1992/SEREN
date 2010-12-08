! SINK_ADVANCE.F90
! C. P. Batty & D. A. Hubber - 23/8/2007
! Advance positions and velocities of sink particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_advance
  implicit none

  debug2("Advancing positions of all sink particles [advance.F90]")
  debug_timing("ADVANCE")

#if defined(EULER)
  call advance_sink_euler
#elif defined(RUNGE_KUTTA)
  call advance_sink_RK
#elif defined(LEAPFROG_KDK)
  call advance_sink_leapfrog_kdk
#elif defined(LEAPFROG_DKD)
  call advance_sink_leapfrog_dkd
#elif defined(PREDICTOR_CORRECTOR)
  call advance_sink_PC
#endif

  return
END SUBROUTINE sink_advance
