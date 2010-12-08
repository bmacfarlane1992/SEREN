! SEREN.F90
! C. P. Batty, D. A. Hubber, A.McLeod & A. P. Whitworth - 8/12/2006 
! Main seren program.
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM seren
  use definitions
  implicit none

#if defined(DEBUG1)
  write(6,*) "**********************************************************************"
  write(6,*) "*            ****     ******    *****     ******   *     *           *"
  write(6,*) "*           *    *    *         *    *    *        **    *           *"
  write(6,*) "*           *         *         *    *    *        * *   *           *"
  write(6,*) "*            ****     *****     *****     ******   *  *  *           *"
  write(6,*) "*                *    *         *    *    *        *   * *           *"
  write(6,*) "*           *    *    *         *    *    *        *    **           *"
  write(6,*) "*            ****     ******    *    *    ******   *     *           *"
  write(6,*) "*                                                                    *"
  write(6,*) "*                            Version 1.0.0                           *"
  write(6,*) "*                              08/12/2010                            *"
  write(6,*) "*                                                                    *"
  write(6,*) "*        Coders : David Hubber, Chris Batty & Andrew McLeod          *"
  write(6,*) "*                 Thomas Bisbas, Krisada Rawiraswattana,             *"
  write(6,*) "*                 Dimitrios Stamatellos, Stefanie Walch,             *"
  write(6,*) "*                 Anthony Whitworth                                  *"
  write(6,*) "*                                                                    *"
  write(6,*) "*               http://www.astro.group.shef.ac.uk/seren              *"
  write(6,*) "**********************************************************************"
#endif


! Call all important routines to set-up up Seren prior to individual sims.
  call seren_setup

! SPH simulation
! ============================================================================
#if defined(SPH_SIMULATION)
  call sph_simulation
#endif

! Hybrid SPH-Nbody simulation
! ============================================================================
#if defined(NBODY_SPH_SIMULATION)
  call nbody_sph_simulation
#endif

! N-body simulation
! ============================================================================
#if defined(NBODY_SIMULATION)
  call nbody_simulation
#endif

! Cleaning up memory
  call clean_up


  stop
END PROGRAM seren
