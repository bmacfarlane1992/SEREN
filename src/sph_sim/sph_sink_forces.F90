! SPH_SINK_FORCES.F90
! C. P. Batty & D. A. Hubber - 2/8/2010
! Calculates gravitational accelerations for all sink particles
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_sink_forces
  use interface_module, only : add_external_gravitational_force,&
       &BHgrav_accel,direct_sink_gravity,direct_sph_gravity
  use type_module
  use particle_module
  use sink_module
  use time_module
#if defined(OPENMP)
  use omp_lib
#endif
#if defined(CELL_WALK)
  use tree_module
#endif
  implicit none

  integer :: s                        ! sink counter
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  real(kind=DP) :: agravs(1:NDIM)     ! Grav. acceleration of sink s
  real(kind=DP) :: pots               ! Grav. potential of sink s
#endif

  debug2("Calculating grav. forces for all sinks [sink_grav_forces.F90]")


! Zero acceleration arrays here for now
  if (accdo_sinks .and. stot > 0) then
     do s=1,stot
        sink(s)%a(1:NDIM) = 0.0_PR
#if defined(DEBUG_FORCES)
        sink(s)%ahydro(1:NDIM) = 0.0_PR
        sink(s)%agrav(1:NDIM) = 0.0_PR
#endif
     end do
  end if


! Loop over all sink particles and calculate grav. forces from SPH particles
! ============================================================================
#if defined(GRAVITY) && defined(SPH_SELF_GRAVITY)
  if (accdo_sinks .and. stot > 0) then
     debug2("Calculating forces on all sink particles")
     debug_timing("SINK_FORCES")
     
     ! Calculate gravitational forces on all sink particles
     ! (Note : -s signifies particle is a sink in gravity routine)
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) DEFAULT(SHARED) &
     !$OMP PRIVATE(agravs,pots)
     do s=1,stot
#if defined(BH_TREE)
        call BHgrav_accel(-s,1.0_PR/sink(s)%h,&
             &sink(s)%r(1:NDIM),agravs(1:NDIM),pots)
#elif defined(BINARY_TREE)
        call binary_gravacc(-s,1.0_PR/sink(s)%h,&
             &sink(s)%r(1:NDIM),agravs(1:NDIM),pots)
#else    
        call direct_sph_gravity(-s,1.0_PR/sink(s)%h,&
             &sink(s)%r(1:NDIM),agravs(1:NDIM),pots)
#endif
        sink(s)%a(1:NDIM) = real(agravs(1:NDIM),PR)
        sink(s)%gpot      = real(pots,PR)
        sink(s)%gpe       = sink(s)%m*sink(s)%gpot
#if defined(GRAVITY) && defined(BH_TREE) && !defined(GEOMETRIC_MAC)
        sink(s)%agravmag  = &
             &real(sqrt(dot_product(agravs(1:NDIM),agravs(1:NDIM))),PR)
#endif
#if defined(DEBUG_FORCES)
        sink(s)%agrav(1:NDIM) = sink(s)%agrav(1:NDIM) + real(agravs(1:NDIM),PR)
#endif
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------
     
  end if
#endif
! ============================================================================


! Now calculate contributions from other sink particles
! ============================================================================
  if (accdo_sinks .and. stot > 0) then
     debug_timing("GRAVITY_FORCES")
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(agravs,pots)
     do s=1,stot
        call direct_sink_gravity(-s,sink(s)%h,&
             &sink(s)%r(1:NDIM),agravs(1:NDIM),pots)
        sink(s)%a(1:NDIM) = sink(s)%a(1:NDIM) + real(agravs(1:NDIM),PR)
        sink(s)%gpot      = sink(s)%gpot + real(pots,PR)
        sink(s)%gpe       = sink(s)%gpe + sink(s)%m*sink(s)%gpot
#if defined(GRAVITY) && defined(BH_TREE) && !defined(GEOMETRIC_MAC)
        sink(s)%agravmag  = &
             &real(sqrt(dot_product(agravs(1:NDIM),agravs(1:NDIM))),PR)
#endif
#if defined(DEBUG_FORCES)
        sink(s)%agrav(1:NDIM) = sink(s)%agrav(1:NDIM) + real(agravs(1:NDIM),PR)
#endif
     end do
     !$OMP END PARALLEL DO
  end if
! ============================================================================


! Add external gravitational force to all active gas particles
! ============================================================================
#if defined(EXTERNAL_FORCE)
  if (accdo_sinks .and. stot > 0) then
     do s=1,stot
        call add_external_gravitational_force(sink(s)%r(1:NDIM),&
             &agravs(1:NDIM),pots)
        sink(s)%a(1:NDIM) = sink(s)%a(1:NDIM) + real(agravs(1:NDIM),PR)
#if defined(DEBUG_FORCES)
        sink(s)%agrav(1:NDIM) = sink(s)%agrav(1:NDIM) + real(agravs(1:NDIM),PR)
#endif
     end do
  end if
#endif
! ============================================================================


  return
END SUBROUTINE sph_sink_forces
