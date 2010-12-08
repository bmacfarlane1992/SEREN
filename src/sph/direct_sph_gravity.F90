! DIRECT_GRAVITY.F90
! C. P. Batty & D. A. Hubber - 24/8/2007
! Calculates gravitational accelerations for particle p (or sink s = -p) 
! using direct summation.  Also calculates gravitational potential.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE direct_sph_gravity(p,invhp,rp,agravp,potp)
  use interface_module, only : gravity_gradh,gravity_nbody,gravity_sph
  use definitions
  use particle_module
  use type_module, only : pgravitystart, pgravityend
#if defined(SINKS) 
  use sink_module
#endif
#if defined(GRAD_H_SPH)
  use hydro_module, only : omega
#endif
  implicit none

  integer, intent(in) :: p                      ! Id of current particle
  real(kind=PR), intent(in) :: invhp            ! Smoothing length of p
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(out) :: agravp(1:NDIM)  ! Gravitational accelertation
  real(kind=DP), intent(out) :: potp            ! Gravitational potential

  integer       :: pp                ! second particle identifier
  real(kind=PR) :: atemp(1:NDIM)     ! temp. grav. accel. variable
  real(kind=PR) :: dpotp             ! grav potential of p due to pp
#if defined(GRAD_H_SPH)
  real(kind=PR) :: zo_p              ! local copy of zeta/omega for p
#endif
#if defined(SINKS)
  integer       :: s                 ! sink particle id
  integer       :: ss                ! second sink particle id
#endif

  debug3("Calculating gravitational force [direct_sph_gravity.F90] for particle",p)

! Initialise gravitational acceleration and potential to zero
  agravp(1:NDIM) = 0.0_PR
  potp = 0.0_PR

! Make local copies of zeta/omega if required
#if defined(GRAD_H_SPH)
  if (p < 0) then
     zo_p = 0.0_PR
  else if (p <= ptot) then
     zo_p = parray(ZETA,p)
  end if
#endif

#if defined(SINKS)
  s = -1
  if (p < 0) s = -p
#endif


! Loop over all other gravitating particles for SPH particles
! ----------------------------------------------------------------------------
  if (p > 0) then
     do pp=pgravitystart,pgravityend
        if (pp == p) cycle
#if defined(N_BODY)
        call gravity_nbody(parray(MASS,pp),rp(1:NDIM),&
             &parray(1:NDIM,pp),atemp(1:NDIM),dpotp)
#elif defined(GRAD_H_SPH)
        call gravity_gradh(invhp,parray(SMOO,pp),parray(MASS,pp),rp(1:NDIM),&
             &parray(1:NDIM,pp),zo_p,parray(ZETA,pp),atemp(1:NDIM),dpotp)
#else
        call gravity_sph(invhp,parray(SMOO,pp),parray(MASS,pp),&
             &rp(1:NDIM),parray(1:NDIM,pp),atemp(1:NDIM),dpotp)
#endif
        agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
        potp = potp + real(dpotp,DP)

     end do
! ----------------------------------------------------------------------------


! Else loop over all SPH particles if considering sink
! ----------------------------------------------------------------------------
  else if (p < 0) then
     do pp=pgravitystart,pgravityend
#if defined(N_BODY)
        call gravity_nbody(parray(MASS,pp),rp(1:NDIM),&
             &parray(1:NDIM,pp),atemp(1:NDIM),dpotp)
#else
        call gravity_sph(invhp,parray(SMOO,pp),parray(MASS,pp),&
             &rp(1:NDIM),parray(1:NDIM,pp),atemp(1:NDIM),dpotp)
#endif
        agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
        potp = potp + real(dpotp,DP)
        
     end do
  end if
! ----------------------------------------------------------------------------

! Record potential due to only SPH particles for polytropic cooling method
#if defined(RAD_WS)
  if (p > 0) sphgpot(p) = potp
#endif

  return
END SUBROUTINE direct_sph_gravity
