! SINK_UPDATE.F90
! D. A. Hubber - 22/04/2010
! Control subroutine for searching for new sinks and accreting SPH particles 
! to existing sinks.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_update
  use time_module, only : n,nsteps,nsearchnext,nsearchstep
  use sink_module, only : nlast_sinks,stot

! Create sink particles once density reaches required threshold.
! Only check at end of sink (i.e. minimum) timestep.
  if (nsteps == nsearchnext) then
     call sink_search
     nsearchnext = nsearchnext + nsearchstep
  end if

! Remove particles that are now bound to sink particles.
! Ensure accretion only occurs at the end of a full sink timestep.
  if (n == nlast_sinks .and. stot > 0) call accrete_particles

  return
END SUBROUTINE sink_update
