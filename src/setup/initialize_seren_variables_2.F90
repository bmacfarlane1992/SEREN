! INITIALIZE_SEREN_VARIABLES_2.F90
! D. A. Hubber - 1/10/2007
! Initializes various variables in Seren before the main simulation begins.  
! Should be called before SPH/N-body simulation has started.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_seren_variables_2
  use particle_module
  use time_module
  use filename_module
  use type_module
  use neighbour_module
  use sink_module
  implicit none

  integer :: p             ! Particle counter

  debug2("Initializing variables [initialize_particles_arrays.F90]")


! Record original particle identifiers if not a restart
! ----------------------------------------------------------------------------
  if (.not. restart) then
     do p=1,ptot
        porig(p) = p
     end do
  end if


! Initialise time counters
! ----------------------------------------------------------------------------
  if (.not. restart) then
     nsteps    = 0
     time      = 0.0_DP
     nextsnap  = firstsnap
     lastsnap  = 0.0_DP
     snapshot  = 0
     ntempnext = ntempstep
     ndiagnext = ndiagstep
     nsnapnext = nsnapstep
     nsinknext = nsinkstep
  else
     if (time < firstsnap) then
        nextsnap = firstsnap
     else
        nextsnap = lastsnap + snaptime
        if (nextsnap < time) nextsnap = time + snaptime
     end if
     if (in_file_form /= "seren_form" .and. in_file_form /= "sf" .and. &
          &in_file_form /= "seren_unform" .and. in_file_form /= "su") then
        if (ntempnext <= nsteps) ntempnext = nsteps + ntempstep
        if (ndiagnext <= nsteps) ndiagnext = nsteps + ndiagstep
        if (nsnapnext <= nsteps) nsnapnext = nsteps + nsnapstep
        if (nsinknext <= nsteps) nsinknext = nsteps + nsinkstep
     end if
  end if
  n = 0
  nresync = n
  timestep = 0.0_DP


! Misc. variables
! ----------------------------------------------------------------------------
  mgas = 0.0_DP
  do p=pgravitystart,ptot
     mgas = mgas + real(parray(MASS,p),DP)
  end do
  if (.not. restart) mgas_orig = mgas
  if (.not. restart) pgas_orig = pgas
#if defined(SINKS) && defined(SMOOTH_ACCRETION)
  mmean = real(mgas_orig,DP) / real(pgas_orig,DP)
#endif


! Set sink radius for fixed (constant) multiple of h
#if defined(FIXED_HMULT_SINKRAD)
  sinkrad = sinkrad*INVKERNRANGE*&
       &((3.0_PR*pp_gather*mgas)/(4.0_PR*PI*pgas*rhosink))**(ONETHIRD)
#endif

! Set minimum smoothing length to 
#if defined(SINKS) && defined(SMOOTH_ACCRETION) && defined(MINIMUM_H)
  if (stot > 0) then
     hmin = INVKERNRANGE*&
          &((3.0_PR*pp_gather*mgas)/(4.0_PR*PI*pgas*rhosink))**(ONETHIRD)
  else
     hmin = 0.0_PR
  end if
#endif


  return
END SUBROUTINE initialize_seren_variables_2
