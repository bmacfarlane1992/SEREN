! UPDATE_THERMAL_PROPERTIES.F90
! D. A. Hubber - 17/9/2009
! Update the thermal properties of all SPH particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE update_thermal_properties
  use interface_module, only : HP_walk_all_rays,thermal
  use time_module
  use hydro_module
  use type_module
  use particle_module, only : ptot
  use HP_module
  implicit none

  integer :: acctot                    ! No. of particles on acc. step
  integer :: i                         ! Aux. loop counter
  integer :: p                         ! Particle counter
  integer, allocatable :: acclist(:)   ! List of active particles
  
  debug2("Calculating thermal properties of all particles [update_thermal_properties.F90]")


! Calculate thermal properties of all (active) particles if on an ionization 
! step due to ionizing radiation from sources.
! ----------------------------------------------------------------------------
#if defined(IONIZING_UV_RADIATION)
!  if (nionall == nsteps .or. nionize == nsteps) then
  if (nionall == nsteps) then

     newtemp(1:ptot) = 0
     temp_min(1:ptot) = 0.0_PR
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     whichHPlevel(1:ptot) = HP_LEVELS + 1
#endif

     ! If there are any sources, call ionization routine
     if (HPtot > 0) then
        do i=1,HPtot
           call HP_walk_all_rays(i)
        end do
        call write_ionization_data
     end if
     
     ! Update all time counters
     if (nionall == nsteps) nionall = nionall + nionallstep
#if defined(EULER) || defined(RUNGE_KUTTA)
     if (nionize == nsteps) nionize = nionize + 1
#else
     if (nionize == nsteps) nionize = nionize + 2
#endif

  end if
#endif


! For multiple particle timesteps, first make a list of all hydro SPH 
! particles on an acceleration step, and then parallelize over that list.
! ----------------------------------------------------------------------------
  debug_timing("THERMAL")
  acctot = 0
  allocate(acclist(1:ptot))
!  do p=phydrostart,phydroend
  do p=1,phydroend
     if (accdo(p)) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if
  end do


! Loop over all active particles and calculate the thermal properties.
! ----------------------------------------------------------------------------
  if (acctot > 0) then     
     !$OMP PARALLEL DO DEFAULT(SHARED) 
!PRIVATE(p)
!     do i=1,acctot
!        p = acclist(i)
     do p=1,ptot
        call thermal(p)
     end do
     !$OMP END PARALLEL DO
  end if

  deallocate(acclist)


  return
END SUBROUTINE update_thermal_properties
