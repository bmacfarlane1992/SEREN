! SINK_GRAVITY.F90
! D. A. Hubber - 23/02/2010
! Calculates gravitational accelerations for particle p (or sink -p) 
! due to sinks only by direct sum.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE direct_sink_gravity(p,invhp,rp,agravp,potp)
  use interface_module, only : gravity_nbody,gravity_sph
  use definitions
  use sink_module
  use particle_module, only : ptot,parray,a
#if defined(GRAD_H_SPH)
  use hydro_module, only : omega
#endif
#if defined(RAD_WS)
  use particle_module, only : sphgpot
#endif
  implicit none

  integer, intent(in) :: p                      ! Id of current particle
  real(kind=PR), intent(in) :: invhp            ! Smoothing length of p
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(out) :: agravp(1:NDIM)  ! Gravitational accelertation
  real(kind=DP), intent(out) :: potp            ! Gravitational potential

  integer       :: s                 ! sink particle identifier
  integer       :: ss                ! sink particle identifier
  real(kind=PR) :: atemp(1:NDIM)     ! temp. grav acceleration var
  real(kind=PR) :: dpotp             ! grav potential of p due to pp
#if defined(GRAD_H_SPH)
  real(kind=PR) :: zo_p              ! local copy of zeta/omega for p
#endif

  debug3("Calculating gravitational force [sink_gravity.F90] for particle",p)

! Initialise gravitational acceleration and potential to zero
  agravp(1:NDIM) = 0.0_DP
  potp = 0.0_DP

! Set sink id if required
  s = -1
  if (p < 0) s = -p

! Make local copies of zeta/omega if required
#if defined(GRAD_H_SPH)
  if (p <= ptot) zo_p = parray(ZETA,p)
#endif

! Record potential due to only SPH particles for polytropic cooling method
#if defined(RAD_WS)
  if (p > 0) sphgpot(p) = potp
#endif


! Add contributions due to all sinks
! ----------------------------------------------------------------------------
  do ss=1,stot
     if (s == ss) cycle
#if defined(N_BODY)
     call gravity_nbody(real(sink(ss)%m,PR),rp(1:NDIM),&
          &sink(s)%r(1:NDIM),atemp,dpotp)
#else
     call gravity_sph(invhp,sink(ss)%h,real(sink(ss)%m,PR),rp(1:NDIM),&
          &sink(ss)%r(1:NDIM),atemp,dpotp)
#endif
     agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
     potp = potp + real(dpotp,DP)
  end do
! ----------------------------------------------------------------------------


  return
END SUBROUTINE direct_sink_gravity
