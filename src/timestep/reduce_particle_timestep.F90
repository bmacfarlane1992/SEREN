! REDUCE_PARTICLE_TIMESTEP.F90
! D. A. Hubber - 9/9/2010
! Reduce the timestep of particle p due to too large a difference in 
! timesteps between neighbours.  Sets to timestep level, lnew
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE reduce_particle_timestep(p,lnew)
  use particle_module
  use hydro_module
  use time_module

  integer, intent(in) :: p                ! Particle id
  integer(kind=ILP), intent(in) :: lnew   ! New timestep level

#if defined(DEBUG_REDUCE_TIMESTEP)
  write(6,*) "Reducing timestep for particle ",p
  write(6,*) "lold :",nlevel(p),"    lnew :",lnew
  write(6,*) "nlast :",nlast(p),"    n :",n
#endif

  if (n /= nlast(p)) laststep(p) = timestep*real(n - nlast(p),DP)

#if defined(LEAPFROG_KDK)
  if (n /= nlast(p)) v_half(1:VDIM,p) = &
       &v_old(1:VDIM,p) + 0.5_PR*a(1:VDIM,p)*real(laststep(p),PR)
  if (n /= nlast(p)) parray(1:NDIM,p) = &
       &r_old(1:NDIM,p) + v_half(1:NDIM,p)*real(laststep(p),PR)
#endif
#if defined(PREDICTOR_CORRECTOR)
  a_old(1:VDIM,p) = a(1:VDIM,p)
#endif
  r_old(1:NDIM,p) = parray(1:NDIM,p)
  v_old(1:VDIM,p) = v(1:VDIM,p)
#if defined(ENTROPIC_FUNCTION) && defined(INTERNAL_ENERGY)
  Aold(p) = Aent(p)
#elif defined(INTERNAL_ENERGY) && !defined(RAD_WS)
  u_old(p) = u(p)
#endif
  rho_old(p) = rho(p)
#if defined(VISC_TD)
  talpha_old(p) = talpha(p)
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQUATION)
  B_old(1:BDIM,p) = B(1:BDIM,p)
#endif

  if (n /= nlast(p)) nlevel(p) = lnew
  if (n /= nlast(p)) nlast(p) = n
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
  nminneib(p) = lnew
#endif

! Calculate new accelerations for all particles to ensure we do not keep 
! old and invalid accelerations.
!  accdo(p) = .true.

#if defined(RUNGE_KUTTA) || defined(EULER) || defined(LEAPFROG_KDK)
  accdo(p) = .true.
#elif defined(LEAPFROG_DKD) || defined(PREDICTOR_CORRECTOR)
  accdo(p) = .false.
#endif


  return
END SUBROUTINE reduce_particle_timestep
