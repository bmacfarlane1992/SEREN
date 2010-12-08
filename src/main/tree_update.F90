! TREE_UPDATE.F90
! D. A. Hubber - 21/9/2007
! Controlling subroutine to call tree build and tree stock routines.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE tree_update
  use tree_module
  use type_module
  use particle_module, only : ptot
  use time_module, only : nsteps, nbuild, nbuildstep, nstock
  implicit none

! Barnes-Hut tree subroutine calls
! ----------------------------------------------------------------------------
#if defined(BH_TREE)
  debug2("Updating BH tree [tree_update.F90]")
  debug_timing("BH_TREE")

! Build new tree(s)
  if (nbuild == nsteps) then
#if defined(GRAVITY)
     call BHgrav_build
#if defined(REORDER_PARTICLES)
     call BH_reorder_particles
#endif
     ! If all particles are gas particles, copy the gravity tree to the
     ! hydro tree.  Else, build the hydro tree from scratch.
     if (pgas + pcdm == ptot) then
        call copy_BHgrav_to_BHhydro
     else
        call BHhydro_build
     end if
#else
     call BHhydro_build
#if defined(REORDER_PARTICLES)
     call BH_reorder_particles
#endif
#endif
     nstock = nsteps
  end if

! Stock trees on every acceleration step
  if (nstock == nsteps) then
#if defined(GRAVITY)
     call BHgrav_stock
#endif
     call BHhydro_stock
  end if
#endif
! ----------------------------------------------------------------------------


! Binary tree
! ----------------------------------------------------------------------------
#if defined(BINARY_TREE)
  debug2("Updating binary tree [tree_update.F90]")
  debug_timing("BINARY_TREE")

! Build new tree
  if (nbuild == nsteps) then
     call binary_treebuild
  end if

! Stock trees every step
  call binary_treestock

#endif
! ----------------------------------------------------------------------------


! Calculate integer time for next tree build and next tree stock depending 
! on the integration scheme
  if (nbuild == nsteps) then
     nbuild = nbuild + nbuildstep
  endif

  if (nstock == nsteps) then
#if defined(EULER) || defined(RUNGE_KUTTA)
     nstock = nstock + 1
#elif defined(LEAPFROG_KDK) || defined(LEAPFROG_DKD) || defined(PREDICTOR_CORRECTOR)
     nstock = nstock + 2
#endif
  endif

  return
END SUBROUTINE tree_update
