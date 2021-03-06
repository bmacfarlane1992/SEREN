! AMBIENTTEMP.F90
! D. Stamatellos - 3/1/2008
! Gets ambient (i.e. heating) temperature of particle. 
! SCALING NEEDS SORTING IN THIS SUBROUTINE *****
! CHECK WITH DIMITRIS
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE ambient_temp(p,atemp)
  use particle_module
  use sink_module
  use Tprof_module
  use scaling_module
  implicit none

  integer, intent(in) :: p             ! particle id
  real(kind=PR), intent(out) :: atemp  ! Ambient temperature

#if defined(HDISC_HEATING) || defined(HDISC_HEATING_3D_SINGLE) || defined(STAR_SIMPLE_HEATING) || defined(HDISC_HEATING_PLUS_STAR_SIMPLE_HEATING) || defined(LOCAL_ISOTHERMAL)
  integer :: s                         ! sink id
  real(kind=PR) :: drsqd               ! radius squared
  real(kind=PR) :: dr(1:NDIM)          ! relative displacement vector
#endif


! Assumes one star/disc with disc in the x-y plane 
! ----------------------------------------------------------------------------
#if defined(HDISC_HEATING) || defined(LOCAL_ISOTHERMAL)
  do s=1,stot
     if (sink(s)%m > 0.5) then
        dr(1:NDIM) = sink(s)%r(1:NDIM) - parray(1:NDIM,p)
        drsqd= (dr(1)*dr(1) + dr(2)*dr(2))*(6.684e-14*rscale*rcgs)**2 
        ! x-y distance in AU
        atemp = sqrt(ptemp0*ptemp0*&
             &(drsqd + ptemp_r0*ptemp_r0)**(-ptemp_q) + &
             &temp_inf*temp_inf)
     end if
  end do


! Assumes one star/disc with disc in the x-y plane 
! ----------------------------------------------------------------------------
#elif defined(HDISC_HEATING_3D_SINGLE)
  s = 1
  dr(1:NDIM) = sink(s)%r(1:NDIM) - parray(1:NDIM,p)
  drsqd= (dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))*&
       &(6.684e-14*rscale*rcgs)**2 
  atemp = sqrt(ptemp0*ptemp0*&
       &(drsqd + ptemp_r0*ptemp_r0)**(-ptemp_q) + temp_inf*temp_inf)



! --- central star as HDISC_HEATING; secondary objects include star_simple heating
#elif defined(HDISC_HEATING_PLUS_STAR_SIMPLE_HEATING)

  do s=1,stot

     dr(1:NDIM) = sink(s)%r(1:NDIM) - parray(1:NDIM,p)
 
    if (sink(s)%m > 0.5) then
        drsqd= (dr(1)*dr(1) + dr(2)*dr(2))*(6.684e-14*rscale*rcgs)**2 
        ! x-y distance in AU
        atemp = sqrt(ptemp0*ptemp0*&
             &(drsqd + ptemp_r0*ptemp_r0)**(-ptemp_q) + &
             &temp_inf*temp_inf)
     else
       drsqd= dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
       atemp = (atemp**4 + 0.25_PR*(sink(s)%star_radius**2/drsqd)&
          &*sink(s)%temperature**4)**(0.25_PR)

     endif

  end do




! ..
! ----------------------------------------------------------------------------
#elif defined(STAR_SIMPLE_HEATING)
  atemp = 0.0_PR
  do s=1,stot       
     dr(1:NDIM) = sink(s)%r(1:NDIM) - parray(1:NDIM,p)
     drsqd= dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
     atemp = (atemp**4 + 0.25_PR*(sink(s)%star_radius**2/drsqd)&
          &*sink(s)%temperature**4)**(0.25_PR)
  end do

#if defined(AMBIENT_HEATING)
  atemp = (temp_inf**4 + atemp**4)**0.25_PR
#endif

! ...
! ----------------------------------------------------------------------------
#elif defined(STAR_HEATING) 
  write(*,*) 'Not done yet'
  stop

! Constant background temperature
! ----------------------------------------------------------------------------
#elif defined(AMBIENT_HEATING)
  atemp = temp_inf
#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE ambient_temp

