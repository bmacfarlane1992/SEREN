! INITIALIZE_SPH_VARIABLES_2.F90
! D. A. Hubber - 1/10/2007
! Sets values for particular variables that need to be initialized AFTER 
! the first force calculations in setup.  Other variables (before the first 
! force calculation) are initialized in initialize_variables_1.F90.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_sph_variables_2
  use interface_module, only : find_equilibrium_temp_ws
  use particle_module
  use hydro_module
  use scaling_module
  use type_module
#ifdef SINK_PROPERTIES_FIX
  use filename_module, only : out_file_form
#endif
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
  use Eos_module, only : column2,fcolumn,idens,itemp
#endif
  implicit none

#if defined(RAD_WS)
  integer :: p               ! Particle counter
 real(kind=PR) :: dt                 ! Timestep
  real(kind=PR) :: dt_new             ! Latest timestep
  real(kind=PR) :: dt_old             ! Previous timestep
  real(kind=PR) :: mu_bar_p           ! Mean gas particle mass for p
#endif
#if defined(SINKS)
  integer :: s               ! Sink counter
#endif

  debug2("Initializing variables [initialize_variables_2.F90]")


! Calculate flux-limited diffusion terms
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
  do p=pgasstart,pgasend
        call find_idens(rho(p),idens(p))
        call find_itemp(temp(p),itemp(p))
#if defined(RAD_WS_SINK_POT)
     column2(p) = (fcolumn**2)*gpot(p)*rho(p)
#else
     column2(p) = (fcolumn**2)*sphgpot(p)*rho(p)
#endif
#if defined(DIFFUSION)
   du_dt_diff(p)=0.0_PR
#endif
 call find_equilibrium_temp_ws(p)
 
 end do
#endif
! ----------------------------------------------------------------------------


! Record 'old' particle properties for integration scheme
!! ----------------------------------------------------------------------------
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
        sink(s)%Mstar              = 0.0_DP
        sink(s)%dmdt_star          = 0.0_DP
#ifdef EPISODIC_ACCRETION
       sink(s)%dmdt_0             = 0.0_DP
       sink(s)%Mdisc              = 0.0_DP
#endif
     end do
  end if

#ifdef SINK_PROPERTIES_FIX
! read in sink properties and calculate accretion history from sink files

#ifndef SINK_PROPERTIES_SYNC
if (restart) then
  if (out_file_form=="dragon_form" .or. out_file_form=="df".or.&
      out_file_form=="dragon_unform" .or. out_file_form=="du") then
     do s=1,stot
        call read_sink_data(s)
     end do
  endif
endif

#endif
#endif
#endif

#ifdef SINK_PROPERTIES_SYNC
 if (out_file_form=="dragon_form" .or. out_file_form=="df".or.&
      out_file_form=="dragon_unform" .or. out_file_form=="du") then
     do s=1,stot
        call read_sink_data_sync(s)
     end do
  endif
stop
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

! Set time for next tree build and stock, and next ionizing HEALPix walk
#if defined(LEAPFROG_KDK)
  nbuild = nsteps + 2
  nstock = nsteps + 2
  nionize = nsteps + 2
#else
  nbuild = nsteps + 1
  nstock = nsteps + 1
  nionize = nsteps + 1
#endif
  nionall = nsteps + nionallstep

  return
END SUBROUTINE initialize_sph_variables_2
