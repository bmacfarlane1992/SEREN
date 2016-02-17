! ADVANCE_SINK_PC.F90
! D. A. Hubber - 21/6/2007
! Advance positions and velocities of all sink particles once 
! accelerations have been computed. Uses 2nd order Predictor-Corrector
! integration scheme with global or multiple particle timesteps.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_sink_PC
  use particle_module
  use time_module
  use sink_module
  implicit none

  integer :: dn                    ! Integer timestep since start of timestep
  integer :: nfull                 ! Full integer timestep
  integer :: nhalf                 ! Half integer timestep
  integer :: s                     ! Sink particle counter
  real(kind=PR) :: as(1:NDIM)      ! Local copy of acceleration
  real(kind=PR) :: asold(1:NDIM)   ! Local copy of old acceleration
  real(kind=DP) :: dt              ! Physical time since beginning of timestep 
  real(kind=PR) :: rs(1:NDIM)      ! Local copy of position 
  real(kind=PR) :: rsold(1:NDIM)   ! Local copy of old position
  real(kind=PR) :: vs(1:NDIM)      ! Local copy of velocity
  real(kind=PR) :: vsold(1:NDIM)   ! Local copy of old velocity

  debug2("Advancing sink particles [advance_sink_RK.F90]")

! Calculate integer and physical time intervals since beginning of step
  accdo_sinks = .false.
  dn          = n - nlast_sinks
!  nfull       = nstep_sinks
  nfull       = 2**(level_step - nlevel_sinks)
  nhalf       = nfull / 2
  dt          = real(timestep,PR)*real(dn,PR)
  sink_dt     =dt
! Loop over all sink particles 
! ----------------------------------------------------------------------------
  if (stot > 0) then
     do s=1,stot

        ! Advance sink particle if below or at half timestep
        ! --------------------------------------------------------------------
        if (dn < nfull) then

           ! Store local copies of old positions and velocities
           rsold(1:NDIM) = sink(s)%rold(1:NDIM)
           vsold(1:NDIM) = sink(s)%vold(1:NDIM)
           asold(1:NDIM) = sink(s)%aold(1:NDIM)
           as(1:NDIM) = sink(s)%a(1:NDIM)

           ! First half of Runge-Kutta integration step
           rs(1:NDIM) = rsold(1:NDIM) + (vsold(1:NDIM)+0.5*asold(1:NDIM)*dt)*dt
           vs(1:NDIM) = vsold(1:NDIM) + as(1:NDIM)*dt

#ifdef PERIODIC
           call check_periodic_box(rs(1:NDIM))
#endif

        ! Advance if beyond half timestep
        ! --------------------------------------------------------------------
        else

           ! Store local copies of old positions and velocities
           rsold(1:NDIM) = sink(s)%rold(1:NDIM)
           vsold(1:NDIM) = sink(s)%vold(1:NDIM)
           asold(1:NDIM) = sink(s)%aold(1:NDIM)
           as(1:NDIM) = sink(s)%a(1:NDIM)

           ! Second half of Runge-Kutta integration step
           rs(1:NDIM) = rsold(1:NDIM) + (vsold(1:NDIM)+0.5*asold(1:NDIM)*dt)*dt
           vs(1:NDIM) = vsold(1:NDIM) + 0.5*(asold(1:NDIM) + as(1:NDIM))*dt

#ifdef PERIODIC
           call check_periodic_box(rs(1:NDIM))
#endif

           ! If end of timestep, record as 'old' values
           sink(s)%rold(1:NDIM) = rs(1:NDIM)
           sink(s)%vold(1:NDIM) = vs(1:NDIM)
           sink(s)%aold(1:NDIM) = as(1:NDIM)

        end if
        ! --------------------------------------------------------------------

        ! Record positions and velocities in arrays for all cases 
        sink(s)%r(1:NDIM) = rs(1:NDIM)
        sink(s)%v(1:NDIM) = vs(1:NDIM)
     
     end do
  end if
! ----------------------------------------------------------------------------

! Update quantities for next timestep/force calculation
  if (dn == nhalf) accdo_sinks = .TRUE.
  if (dn == nfull) then
     accdo_sinks = .TRUE.
     nlast_sinks = n
  end if

  return
END SUBROUTINE advance_sink_PC
