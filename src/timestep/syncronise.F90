! REDUCE_PARTICLE_TIMESTEP.F90
! D. A. Hubber - 9/9/2010
! Reduce the timestep of particle p due to too large a difference in 
! timesteps between neighbours.  Sets to timestep level, lnew
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE syncronise
  use particle_module
  use hydro_module
  use time_module
  use sink_module 
  use seren_sim_module

  integer :: p          ! Particle counter
  integer :: s          ! Particle counter

  nresync=n

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p)
  do p=1,ptot

  r_old(1:NDIM,p) = parray(1:NDIM,p)
  v_old(1:VDIM,p) = v(1:VDIM,p)

  if (n/= nlast(p)) laststep(p) = timestep*real(n - nlast(p),DP)

#if defined(LEAPFROG_KDK)
  if (n /= nlast(p)) v_half(1:VDIM,p) = &
       &v_old(1:VDIM,p) + 0.5_PR*a(1:VDIM,p)*real(laststep(p),PR)
  if (n /= nlast(p)) parray(1:NDIM,p) = &
       &r_old(1:NDIM,p) + v_half(1:NDIM,p)*real(laststep(p),PR)
#endif
#if defined(PREDICTOR_CORRECTOR)
  a_old(1:VDIM,p) = a(1:VDIM,p)
#endif

#if defined(ENTROPIC_FUNCTION) && defined(INTERNAL_ENERGY)
  Aold(p) = Aent(p)
#elif defined(INTERNAL_ENERGY)
  u_old(p) = u(p)
#endif
  rho_old(p) = rho(p)
#if defined(VISC_TD)
  talpha_old(p) = talpha(p)
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  B_old(1:BDIM,p) = B(1:BDIM,p)
#endif

  if (n /= nlast(p)) nlast(p) = n

! Calculate new accelerations for all particles to ensure we do not keep 
! old and invalid accelerations.
!  accdo(p) = .true.

#if defined(RUNGE_KUTTA) || defined(EULER) || defined(LEAPFROG_KDK)
  accdo(p) = .true.
#elif defined(LEAPFROG_DKD) || defined(PREDICTOR_CORRECTOR)
  accdo(p) = .false.
#endif

 enddo

    ! If using sinks in SPH simulation, reduce sink timesteps here.
     ! -----------------------------------------------------------------------
#if defined(SPH_SIMULATION) && defined(SINKS)
     if (sph_sim) then
        if (n /= nlast_sinks) &
             &laststep_sinks = timestep*real(n - nlast_sinks,DP)
        do s=1,stot
           sink(s)%rold(1:NDIM) = sink(s)%r(1:NDIM)
           sink(s)%vold(1:NDIM) = sink(s)%v(1:NDIM)
#if defined(LEAPFROG_KDK)
           if (n /= nlast_sinks) sink(s)%vhalf = sink(s)%vold(1:VDIM) &
                & + 0.5_PR*sink(s)%a(1:VDIM)*real(laststep_sinks,PR)
           if (n /= nlast_sinks) sink(s)%r(1:NDIM) = sink(s)%rold(1:NDIM) &
                & + sink(s)%vhalf(1:NDIM)*real(laststep_sinks,PR)
#endif
        end do
        nlast_sinks = n
     end if
#endif

     ! If using sinks in SPH simulation, reduce sink timesteps here.
     ! -----------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
     if (nbody_sph_sim .or. nbody_sim) then
        do s=1,stot
           star(s)%nlast = n
        end do
     end if
#endif

!Note for future generations the problem was at reduce particle timestes !defined(RAD_WS) should not be there!!!

  return
END SUBROUTINE syncronise
