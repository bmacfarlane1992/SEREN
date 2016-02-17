! REDUCE_TIMESTEP.F90
! D. A. Hubber - 9/9/2010
! Reduce the timesteps of any particles which have neighbours with 
! relatively low timesteps.  Only reduces timesteps of particles when the 
! new timestep is correctly synchronized, and the particle has not just 
! started a new timestep (which is handled in sph_timestep.F90 and other 
! such routines). Also reduces the timestep of any sink neighbouring 
! particles to the sink timestep (i.e. the lowest occupied timestep level).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE reduce_timesteps
  use interface_module, only : reduce_particle_timestep
  use neighbour_module
  use time_module
  use particle_module, only : ptot,parray
  implicit none

  integer :: p            ! Particle counter

  debug2("Reduce timesteps of all flagged particles [reduce_timestep.F90]")

! Loop over all particles and reduce the timestep of any particles if 
! the timesteps are synchronised.
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p)
  do p=1,ptot
     if (nminneib(p) > level_max) then
        write(6,*) "Problem with nminneib > level_max : ",nminneib(p),level_max
        stop
     end if
     if (nminneib(p) - nlevel(p) > TIMESTEP_LEVEL_DIFF_MAX .and. &
             & nminneib(p) /= -1 .and. nlast(p) /= n .and. &
             & mod(n,2**(level_step - nminneib(p) + TIMESTEP_LEVEL_DIFF_MAX)) &
             & == 0_ILP) then
        call reduce_particle_timestep(p,nminneib(p) - TIMESTEP_LEVEL_DIFF_MAX)
     end if
  end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------

  return
END SUBROUTINE reduce_timesteps
