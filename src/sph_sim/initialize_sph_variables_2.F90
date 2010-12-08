! INITIALIZE_SPH_VARIABLES_2.F90
! D. A. Hubber - 1/10/2007
! Sets values for particular variables that need to be initialized AFTER 
! the first force calculations in setup.  Other variables (before the first 
! force calculation) are initialized in initialize_variables_1.F90
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_sph_variables_2
  use interface_module, only : eos
  use particle_module
  use hydro_module
  use scaling_module
  use type_module
  use time_module
  use neighbour_module
  use timing_module
  use filename_module, only : restart
#if defined(SINKS)
  use sink_module
#endif
#if defined(IDEAL_MHD)
  use mhd_module
#endif
#if defined(RAD_WS)
  use Eos_module, only : column2,fcolumn
#endif
  implicit none

#if defined(RAD_WS)
  integer :: p               ! Particle counter
#endif
#if defined(SINKS)
  integer :: s               ! Sink counter
#endif

  debug2("Initializing variables [initialize_variables_2.F90]")


! Set initial column densities
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
  do p=pgasstart,pgasend
#if defined(RAD_WS_SINK_POT)
     column2(p) = (fcolumn**2)*gpot(p)*rho(p)
#else
     column2(p) = (fcolumn**2)*sphgpot(p)*rho(p)
#endif
     call eos(p)
  end do
#endif


! Record 'old' particle properties for integration scheme
! ----------------------------------------------------------------------------
  r_old(1:NDIM,1:ptot)  = parray(1:NDIM,1:ptot)
  v_old(1:VDIM,1:ptot)  = v(1:VDIM,1:ptot)
  rho_old(1:ptot)       = rho(1:ptot)
#if defined(INTERNAL_ENERGY)
  u_old(1:ptot)         = u(1:ptot)
#endif
#if defined(ENTROPIC_FUNCTION)
  Aold(1:ptot)          = Aent(1:ptot)
#endif
#if defined(LEAPFROG_KDK)
  v_half(1:VDIM,1:ptot) = v(1:VDIM,1:ptot)
#endif
#if defined(PREDICTOR_CORRECTOR)
  a_old(1:VDIM,1:ptot)  = a(1:VDIM,1:ptot)
#endif
#if defined(IDEAL_MHD)
  B_old(1:BDIM,1:ptot)  = B(1:BDIM,1:ptot)
#endif


! Initialise sink particles info not contained in IC file
! ----------------------------------------------------------------------------
#if defined(SINKS)
  do s=1,stot
     sink(s)%rold(1:NDIM) = sink(s)%r(1:NDIM)
     sink(s)%vold(1:VDIM) = sink(s)%v(1:VDIM)
#if defined(LEAPFROG_KDK)
     sink(s)%vhalf(1:VDIM) = sink(s)%v(1:VDIM)
#endif
#if defined(PREDICTOR_CORRECTOR)
     sink(s)%aold(1:VDIM) = sink(s)%a(1:VDIM)
#endif
#if defined(SMOOTH_ACCRETION)
     sink(s)%mmax = 4.0_PR*PI/3.0_PR*rhosink*sink(s)%radius**3
     sink(s)%menc = sink(s)%menc
     sink(s)%trot = TWOPI*sqrt(sink(s)%radius**3/sink(s)%menc)
#endif
  end do
  if (.not. restart) then
     do s=1,stot
        sink(s)%macc(1:DMDT_RANGE) = 0.0_DP
        sink(s)%tacc(1:DMDT_RANGE) = 0.0_DP
        sink(s)%dmdt               = 0.0_DP
        sink(s)%angmom(1:3)        = 0.0_DP
        sink(s)%angmomnet(1:3)     = 0.0_DP
     end do
  end if
#endif

! Set accdo to FALSE since we have already calculated the initial 
! accelerations for integration schemes which require them.
  accdo(1:ptot) = .false.

! Set accdo_sinks depending on integration scheme used
#if defined(SINKS)
  nsearchnext = nsteps + nsearchstep
#if defined(RUNGE_KUTTA) || defined(EULER) || defined(LEAPFROG_KDK)
  accdo_sinks = .true.
#elif defined(LEAPFROG_DKD) || defined(PREDICTOR_CORRECTOR)
  accdo_sinks = .false.
#endif
#endif

! Set time for next tree build and stock
#if defined(LEAPFROG_KDK)
  nbuild = nsteps + 2
  nstock = nsteps + 2
#else
  nbuild = nsteps + 1
  nstock = nsteps + 1
#endif

! Set time for next ioniziing radiation HEALPix walk
#if defined(LEAPFROG_KDK)
  nionall = nsteps + nionallstep
  nionize = nsteps + 2
#else
  nionall = nsteps + nionallstep
  nionize = nsteps + 1
#endif

  return
END SUBROUTINE initialize_sph_variables_2
