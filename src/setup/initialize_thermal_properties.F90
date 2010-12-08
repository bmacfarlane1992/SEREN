! INITIALIZE_THERMAL_PROPERTIES.F90
! D. A. Hubber - 24/06/2009
! Initialize any thermal properties depending on whether the energy equation
! is solved, if Dragon format is used to read in the ICs, and on the 
! equation of state chosen.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_thermal_properties
  use interface_module
  use filename_module
  use particle_module
  use hydro_module
  use type_module
  use Eos_module
  implicit none

  integer :: p                  ! Particle counter

  debug2("Calculate all initial thermal variables [initialize_thermal_properties.F90]")


! Initialise internal energies based on Radiative cooling algorithm
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
!  if (pgas < ptot) then
!     !$OMP PARALLEL DO DEFAULT(SHARED)
!     do p=1,pgravitystart-1
!        u(p) = Pconst*temp(p)/(gamma - 1.)
!        u_old(p) = u(p)
!        du_dt(p) = real(0.0,PR)
!     end do
!     !$OMP END PARALLEL DO
!  end if

!$OMP PARALLEL DO DEFAULT(SHARED)
  do p=1,ptot
!  do p=pgravitystart,ptot
     call find_idens(rho(p),idens(p))
     call find_itemp(temp(p),itemp(p))
     u(p)     = eosenergy(rho(p),temp(p),idens(p),itemp(p))
     u_old(p) = u(p)
     du_dt(p) = 0.0_PR
#if defined(DIFFUSION)
     du_dt_diff(p) = 0.0_PR
     k_cond(p) = 0.0_PR
#endif
     press(p) = Pconst2*temp(p)*rho(p) / &
          &eosmu(rho(p),temp(p),idens(p),itemp(p))
     sound(p) = sqrt(press(p)/rho(p))
     call find_temp_from_energy(idens(p),u(p),itemp(p),temp(p))
  end do
!$OMP END PARALLEL DO

! Still update all thermal properties in case other routines are needed 
! (e.g. ionization)
  call update_thermal_properties


! Otherwise when using Dragon format (which does not record internal 
! energies), simply convert temperatures to internal energies depending 
! on the ratio of specific heats, gamma, chosen in the params.dat file.
! ----------------------------------------------------------------------------
#elif defined(INTERNAL_ENERGY) && !defined(U_IMPLICIT_SOLVER)
  do p=1,ptot
     u_old(p) = u(p)
     du_dt(p) = 0.0_PR
  end do

#if defined(ENTROPIC_FUNCTION) && defined(INTERNAL_ENERGY)
  do p=1,ptot
     Aent(p) = (gamma - 1.0_PR)*u(p)/rho(p)**(gamma - 1.0_PR)
     Aold(p) = Aent(p)
     dA_dt(p) = 0.0_PR
  end do
#endif

  call update_thermal_properties


! Otherwise, if using some other equation of state (e.g. barotropic), 
! calculate the temperatures now we have the densities.
! ----------------------------------------------------------------------------
#else
  temp(1:ptot) = 0.0_PR
  call update_thermal_properties

#endif
! ----------------------------------------------------------------------------


  return
END SUBROUTINE initialize_thermal_properties
