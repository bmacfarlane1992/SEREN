! RECORD_PARTICLE_DATA.F90
! D. A. Hubber - 06/03/2010
! Records and returns all float variables for particle p in a single array 
! (alldata) relative to some specified origin (rorigin).  Also returns no. 
! of filled elements in array (nelements)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE record_particle_data(p,nelements,alldata,rorigin)
  use particle_module
  use hydro_module
  use mhd_module
  use neighbour_module
  use scaling_module
  use filename_module
  use time_module
  use type_module
  use HP_module
#if defined(RAD_WS)
  use Eos_module
#endif
  implicit none

  integer, intent(in) :: p                      ! I.d. of particle
  integer, intent(out) :: nelements             ! No. of filled elements
  real(kind=PR), intent(out) :: alldata(1:100)  ! Particle data array
  real(kind=PR), intent(in) :: rorigin(1:NDIM)  ! Origin for output

  integer :: k                      ! Dimension counter
  real(kind=PR) :: a_rad            ! magnitude of acceleration
  real(kind=PR) :: dr(1:NDIM)       ! relative displacement vector
  real(kind=PR) :: drmag            ! Magnitude of distance
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: dr_unit(1:NDIM)  ! Unit vector
  real(kind=DP) :: dt               ! timestep
  real(kind=PR) :: dv_unit(1:NDIM)  ! Unit vector in direction of velocity
  real(kind=PR) :: vmag             ! Magnitude of velocity
  real(kind=PR) :: v_rad            ! Radial component of velocity
  real(kind=PR) :: vtemp(1:NDIM)    ! temp vector for scalar product call
#if defined(IONIZING_UV_RADIATION)
  real(kind=PR) :: drhodr           ! Density gradient
#endif
#if defined(DEBUG_FORCES)
  real(kind=PR) :: ag_rad           ! radial gravitational acceleration
  real(kind=PR) :: ah_rad           ! radial hydro acceleration
  real(kind=PR) :: av_rad           ! radial viscous acceleration
  real(kind=PR) :: av_vel           ! viscous accel in direction of velocity
#endif

  nelements      = 0
  alldata(1:100) = 0.0_PR
! Calculate unit vector here
  call distance2(rorigin(1:NDIM),p,dr(1:NDIM),drsqd)
  drmag = sqrt(drsqd) + SMALL_NUMBER
  if (drmag < SMALL_NUMBER) then
     dr_unit(1:NDIM) = 0.0_PR
  else
     dr_unit(1:NDIM) = dr(1:NDIM) / drmag
  end if
  vmag = sqrt(dot_product(v(1:VDIM,p),v(1:VDIM,p))) + SMALL_NUMBER
  if (vmag < SMALL_NUMBER) then
     dv_unit(1:NDIM) = 0.0_PR
  else
     dv_unit(1:NDIM) = v(1:NDIM,p) / vmag
  end if
  
! Positions (always written)
  do k=1,NDIM
     nelements = nelements + 1
     alldata(nelements) = parray(k,p)*real(rscale,PR)
  end do
  
! Velocities (always written)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = v(k,p)*real(vscale,PR)
     vtemp(k) = alldata(nelements)
  end do
  v_rad = dot_product(vtemp,dr_unit)

! Smoothed velocities
#if defined(SMOOTHED_VELOCITY)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = v_smooth(k,p)*real(vscale,PR)
     vtemp(k) = alldata(nelements)
  end do
#endif
  
! Magnetic field
#if defined(IDEAL_MHD) && defined(HYDRO)
  do k=1,BDIM
     nelements = nelements + 1
     alldata(nelements) = B(k,p)*real(Bscale,PR)
  end do
#endif
  
! Induction equation
#if defined(IDEAL_MHD) && defined(INDUCTION_EQN) && defined(HYDRO)
  do k=1,BDIM
     nelements = nelements + 1
     alldata(nelements) = dB_dt(k,p)*real(Bscale/tscale,PR)
  end do
#endif
  
! Induction equation
#if defined(IDEAL_MHD) && defined(DEBUG_MHD) && defined(HYDRO)
  nelements = nelements + 1
  alldata(nelements) = div_B(p)*real(Bscale/rscale,PR)
#endif
  
! Accelerations (always written)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = a(k,p)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  a_rad = dot_product(vtemp,dr_unit)
  
! Hydro accelerations
#if defined(DEBUG_FORCES)
#if defined(HYDRO)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = a_hydro(k,p)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  ah_rad = dot_product(vtemp,dr_unit)
#endif
  
! Gravitational acceleration
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = a_grav(k,p)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  ag_rad = dot_product(vtemp,dr_unit)
#endif
  
! Magnetic acceleration
#if defined(IDEAL_MHD)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = a_mag(k,p)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
#endif

! Viscous acceleration
#if defined(ARTIFICIAL_VISCOSITY)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = a_visc(k,p)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  av_rad = dot_product(vtemp(1:NDIM),dr_unit(1:NDIM))
  av_vel = dot_product(vtemp(1:NDIM) - &
       &av_rad*dr_unit(1:NDIM),dv_unit(1:NDIM))
#endif
#endif
  
! Smoothing length (always written)
  nelements = nelements + 1
  alldata(nelements) = parray(SMOO,p)*real(rscale,PR)
  
! Density (always written)
  nelements = nelements + 1
  alldata(nelements) = rho(p)*real(rhoscale*rhocgs,PR)
  
! div_v (always written)
  nelements = nelements + 1
  alldata(nelements) = div_v(p)/real(tscale,PR)
  
! timestep (always written)
  call timestep_size(p,dt)
  nelements = nelements + 1
  alldata(nelements) = real(dt,PR)*real(tscale,PR)
  
!#if defined(IONIZING_UV_RADIATION)
!  do k=1,NDIM
!     nelements = nelements + 1
!     alldata(nelements) = gradrho(k,p)*real(rhoscale/rscale,PR)
!     vtemp(k) = alldata(nelements)
!  end do
!  drhodr = dot_product(vtemp(1:NDIM),dr_unit(1:NDIM))
!#endif
  
#if defined(HYDRO)
! Temperature
  nelements = nelements + 1
  alldata(nelements) = temp(p)
  
! Sound speed
  nelements = nelements + 1
  alldata(nelements) = sound(p)*real(vscale,PR)
  
! Pressure
  nelements = nelements + 1
  alldata(nelements) = press(p)*real(Pscale,PR)
  
! Specific internal energy
#if defined(INTERNAL_ENERGY)
  nelements = nelements + 1
  alldata(nelements) = u(p)*real(uscale,PR)
  nelements = nelements + 1
  alldata(nelements) = du_dt(p)*real(uscale/tscale,PR)
#endif

! Entropic function
#if defined(ENTROPIC_FUNCTION)
  nelements = nelements + 1
  alldata(nelements) = Aent(p)
  nelements = nelements + 1
  alldata(nelements) = dA_dt(p)
#endif
#endif
  
! Radial distance (always written)
  nelements = nelements + 1
  alldata(nelements) = drmag*real(rscale,PR)
  
! Other radial values
  nelements = nelements + 1
  alldata(nelements) = v_rad
  nelements = nelements + 1
  alldata(nelements) = a_rad
#if defined(DEBUG_FORCES)
#if defined(HYDRO)
  nelements = nelements + 1
  alldata(nelements) = ah_rad
#endif
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  nelements = nelements + 1
  alldata(nelements) = ag_rad
  nelements = nelements + 1
  alldata(nelements) = gpot(p)*real(Escale/mscale,PR)
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  nelements = nelements + 1
  alldata(nelements) = av_rad
  nelements = nelements + 1
  alldata(nelements) = av_vel
#endif
#endif
#if defined(IONIZING_UV_RADIATION)
  nelements = nelements + 1
  alldata(nelements) = drhodr
#endif
  
#if defined(GRAD_H_SPH)
  nelements = nelements + 1
  alldata(nelements) = omega(p)
  nelements = nelements + 1
  alldata(nelements) = h_fac*(parray(MASS,p)/rho(p))**(INVNDIM) / &
       &parray(SMOO,p)
#if defined(GRAVITY)
  nelements = nelements + 1
  alldata(nelements) = parray(ZETA,p)
#endif
#endif
  
! Smoothing length (always written)
  nelements = nelements + 1
  alldata(nelements) = parray(MASS,p)*real(mscale,PR)
  
! Particle type
  nelements = nelements + 1
  if (pboundary > 0 .and. p <= pboundary) then
     alldata(nelements) = BOUNDARYID
  else if (picm > 0 .and. p <= pboundary + picm) then
     alldata(nelements) = ICMID
  else if (pgas > 0 .and. p <= pboundary + picm + pgas) then
     alldata(nelements) = GASID
  end if
  
! Radiative transfer routines
#if defined(RAD_WS) && defined(DEBUG_RAD) 
  nelements = nelements + 1
  alldata(nelements) = rad_info(1,p)*real(kappascale,PR)
  nelements = nelements + 1
  alldata(nelements) = rad_info(2,p)*real(kappascale,PR)
  nelements = nelements + 1
  alldata(nelements) = rad_info(3,p)
  nelements = nelements + 1
  alldata(nelements) = rad_info(4,p)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(5,p)*real(mscale*mcgs/rscale/rcgs/rscale/rcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(6,p)/parray(MASS,p)*real(Escale*Ecgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(7,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(8,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(9,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(10,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
#endif
  
! Flux-limited diffusion info
#if defined(DIFFUSION)
  nelements = nelements + 1
  alldata(nelements) = ueq(p)*real(Escale,PR)
  nelements = nelements + 1
  alldata(nelements) = dt_therm(p)*tscale
  nelements = nelements + 1
  alldata(nelements) = du_dt_diff(p)*Escale/tscale
  nelements = nelements + 1
  alldata(nelements) = k_cond(p)*Escale*Ecgs/rscale/rcgs/tscale
#endif
  
#if defined(DEBUG_SINK_BOUNDARY_PROPERTIES)
  nelements = nelements + 1
  alldata(nelements) = a_rad
  nelements = nelements + 1
  alldata(nelements) = &
       &dot_product(a(1:NDIM,p) - a_rad*dr_unit(1:NDIM),dr_unit(1:NDIM))
#endif
  
#if defined(VISC_TD)
  nelements = nelements + 1
  alldata(nelements) = talpha(p)
#endif
  
#if defined(VISC_BALSARA)
  nelements = nelements + 1
  alldata(nelements) = balsara(p)
#endif

#if defined(VISC_PATTERN_REC)
  nelements = nelements + 1
  alldata(nelements) = pattrec(p)
#endif
  
#if defined(DIV_A)
  nelements = nelements + 1
  alldata(nelements) = div_a(p)
#endif

#if defined(SIGNAL_VELOCITY)
  nelements = nelements + 1
  alldata(nelements) = vsigmax(p)*real(vscale,PR)
#endif

#if defined(DEBUG_HP_WALK_ALL_RAYS)
  nelements = nelements + 1
  alldata(nelements) = whichHPlevel(p)
#endif

! Rate of change of density (always written)
  nelements = nelements + 1
  alldata(nelements) = drhodt(p)*real(rhoscale/tscale,PR)

! Time
  nelements = nelements + 1
  alldata(nelements) = real(time*tscale,PR)
  

  return
END SUBROUTINE record_particle_data
