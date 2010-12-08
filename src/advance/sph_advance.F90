! ADVANCE.F90
! C. P. Batty & D. A. Hubber - 23/8/2007
! Calls routines to advance positions and velocities of all particles 
! depending the integration scheme employed.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_advance
  use interface_module, only : advance_boundary_particle,advance_euler,&
       &advance_leapfrog_dkd,advance_leapfrog_kdk,&
       &advance_predictor_corrector,advance_runge_kutta
  use particle_module, only : ptot
  use type_module, only : pboundary,phydrostart
  implicit none

  integer :: p          ! Particle counter

  debug_timing("ADVANCE")

! First, 'advance' or frog-march boundary particles
! ----------------------------------------------------------------------------
  if (pboundary > 0) then
     debug2("Advancing positions of boundary particles [advance.F90]")
     do p=1,pboundary
        call advance_boundary_particle(p)
     end do
  end if
! ----------------------------------------------------------------------------


! Next, advance hydro particles depending on choice of integration scheme
! ----------------------------------------------------------------------------
  debug2("Advancing positions of all SPH particles [advance.F90]")

!$OMP PARALLEL DO DEFAULT(SHARED)
  do p=phydrostart,ptot
#if defined(EULER)
    call advance_euler(p)
#elif defined(RUNGE_KUTTA)
    call advance_runge_kutta(p)
#elif defined(LEAPFROG_KDK)
    call advance_leapfrog_kdk(p)
#elif defined(LEAPFROG_DKD)
    call advance_leapfrog_dkd(p)
#elif defined(PREDICTOR_CORRECTOR)
    call advance_predictor_corrector(p)
#endif
  end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------

  return
END SUBROUTINE sph_advance
