! MODULES.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Contains all global variable modules
! - constant_module
! - Eos_module; Tprof_module; Eos_functions
! - ewald_module
! - filename_module
! - HIIregion_module
! - hydro_module
! - kernel_module
! - mhd_module
! - Nbody_module
! - neighbour_modules
! - particle_module
! - periodic_module
! - scaling_module
! - sink_module
! - sink_correction_module
! - time_module
! - timing_module
! - tree_module
! - type_module
! ============================================================================

#include "macros.h"

! ============================================================================
MODULE constant_module
  use definitions

  ! Physical constants in SI units (unless stated otherwise)
  real(kind=DP),parameter :: r_pc    = 3.08568E16_DP     ! parsec
  real(kind=DP),parameter :: r_au    = 1.49597870E11_DP  ! astronomical unit
  real(kind=DP),parameter :: r_sun   = 6.96E8_DP         ! solar radius
  real(kind=DP),parameter :: r_earth = 6.371E6_DP        ! Earth radius
  real(kind=DP),parameter :: pc_au   = 206265.0_DP       ! parsec (in AU)
  real(kind=DP),parameter :: km_cm   = 1.E5_DP           ! km (in cm)
  real(kind=DP),parameter :: m_sun   = 1.98892E30_DP     ! solar mass
  real(kind=DP),parameter :: m_jup   = 1.8986E27_DP      ! Jupiter mass
  real(kind=DP),parameter :: m_earth = 5.9736E24_DP      ! Earth mass
  real(kind=DP),parameter :: myr = 3.1556952E13_DP       ! megayear
  real(kind=DP),parameter :: yr  = 3.1556952E7_DP        ! year
  real(kind=DP),parameter :: day = 8.64E4_DP             ! day
  real(kind=DP),parameter :: amu = 1.660538782E-27_DP    ! atomic mass unit
  real(kind=DP),parameter :: m_hydrogen = 1.66054E-27_DP ! Hydrogen mass in kg
  real(kind=DP),parameter :: G_const = 6.67428E-11_DP    ! Grav. constant
  real(kind=DP),parameter :: k_boltzmann = 1.3807E-23_DP ! Boltzmann constant 
  real(kind=DP),parameter :: stefboltz = 5.6704E-8_DP    ! stefan-boltzmann
  real(kind=DP),parameter :: e_charge = 1.6021765E-19_DP ! electron charge
  real(kind=DP),parameter :: mu_0 = 1.25663706144E-6_DP  ! permeability of
                                                         ! free space
  real(kind=DP),parameter :: kappa_const = 2.09E-4_DP    ! opacity kappa (cgs)
  real(kind=DP),parameter :: L_sun = 3.839E26_DP         ! solar luminosity

  ! Numerical constants
  real(kind=PR),parameter :: pi = 3.1415926536_PR
  real(kind=PR),parameter :: twopi = 6.283185307_PR
  real(kind=PR),parameter :: invpi = 0.318309886_PR
  real(kind=PR),parameter :: invlogetwo = 1.442695041_PR
  real(kind=PR),parameter :: invlog10two = 3.321928095_PR
  real(kind=PR),parameter :: invsqrttwo = 0.707106781_PR
  real(kind=PR),parameter :: onethird = 0.333333333_PR
  real(kind=PR),parameter :: onesixth = 0.166666666_PR
  real(kind=PR),parameter :: twothirds = 0.666666666_PR
  real(kind=PR),parameter :: big_number = 9.9e20_PR
  real(kind=PR),parameter :: small_number = 1.0e-20_PR

  real(kind=DP),parameter :: pi_dp = 3.1415926536_DP
  real(kind=DP),parameter :: twopi_dp = 6.283185307_DP
  real(kind=DP),parameter :: invpi_dp = 0.318309886_DP
  real(kind=DP),parameter :: invlogetwo_dp = 1.442695041_DP
  real(kind=DP),parameter :: invlog10two_dp = 3.321928095_DP
  real(kind=DP),parameter :: invsqrttwo_dp = 0.707106781_DP
  real(kind=DP),parameter :: onethird_dp = 0.333333333_DP
  real(kind=DP),parameter :: onesixth_dp = 0.166666666_DP
  real(kind=DP),parameter :: twothirds_dp = 0.666666666_DP
  real(kind=DP),parameter :: big_number_dp = 9.9e20_DP
  real(kind=DP),parameter :: small_number_dp = 1.0e-20_DP  

END MODULE constant_module


! ============================================================================
MODULE Eos_module
  use definitions

  integer :: dim_temp                 ! dimension of temp table
  integer :: dim_dens                 ! dimension of dens table
  integer, allocatable :: idens(:)    ! indices for dens from table
  integer, allocatable :: itemp(:)    ! indices for temp from table

  real(kind=PR) :: bdens                        ! log10 density factor
  real(kind=PR) :: btemp                        ! log10 temp. factor
  real(kind=PR) :: densmax                      ! Max. density in eos table
  real(kind=PR) :: densmin                      ! Min. density in eos table
  real(kind=PR) :: fcolumn                      ! column density
                                                ! polytrope correction
  real(kind=PR) :: rad_const                    ! scaling constant
  real(kind=PR) :: tempmin                      ! Max. temp. in eos table
  real(kind=PR) :: tempmax                      ! Min. temp. in eos table

  real(kind=PR) :: z_factor                     ! metallicity

  real(kind=PR), allocatable :: eos_dens(:)     ! density of EOS table
  real(kind=PR), allocatable :: eos_temp(:)     ! temperature of EOS table
  real(kind=PR), allocatable :: eos_energy(:,:) ! energy from EOS table
  real(kind=PR), allocatable :: eos_mu(:,:)     ! mu from EOS table
  real(kind=PR), allocatable :: kappa(:,:)      ! mean opacity EOS table
  real(kind=PR), allocatable :: kappar(:,:)     ! Rosseland opacity EOS table
  real(kind=PR), allocatable :: kappap(:,:)     ! planck opacity EOS table
  real(kind=PR), allocatable :: column2(:)      ! squared column density
  real(kind=PR), allocatable :: ueq(:)          ! Equilibrium u
  real(kind=PR), allocatable :: dt_therm(:)     ! Thermal cooling timescale
#if defined(DEBUG_RAD)
  real(kind=PR), allocatable :: rad_info(:,:)   ! Debug info
#endif


END MODULE Eos_module


!-----------------------------------------------------------------------------
MODULE Tprof_module
  use definitions

! Temperature profile away from a luminosity source
  real(kind=PR) :: ptemp0        ! temperature at r=1AU from the star
  real(kind=PR) :: ptemp_r0      ! temperature softening radius (<<1AU)
  real(kind=PR) :: ptemp_q       ! temperature power law index
  real(kind=PR) :: temp_inf      ! temperature at infinity

END MODULE Tprof_module


! ============================================================================
MODULE ewald_module
  use definitions
  implicit none

! Ewald correction grid size
#if NDIM==2
  integer, parameter :: ewsize(1:NDIM)=(/64,64/)
#elif NDIM==3
  integer, parameter :: ewsize(1:NDIM)=(/64,64,64/)
#endif
  real(kind=PR) :: eforce(1:NDIM)    ! Ewald correction acceleration
  real(kind=PR) :: ewsizeil(1:NDIM)  ! (ewsize-1)/L
  real(kind=PR) :: L(1:NDIM)         ! Periodic box size
  real(kind=PR), allocatable :: fcorr(:,:,:,:)  ! Ewald correction grids

END MODULE ewald_module


! ============================================================================
MODULE filename_module
  use definitions

  logical :: inifile                  ! Has .ini file been written yet?
  logical :: restart                  ! is this simulation a restart?
  integer :: nparams                  ! no. of params available
  integer :: ntemp                    ! next temp file to be written
  integer :: ptrack                   ! (original) id of particle to track
  character(len=256) :: in_file       ! initial conditions file
  character(len=256) :: run_id        ! simulation run identifier (name)
  character(len=256) :: run_dir       ! name of run directory
  character(len=256) :: fileform_ext  ! file format extension for filenames
  character(len=256) :: in_file_form  ! format of initial conditions file
  character(len=256) :: out_init      ! initial snapshot output
  character(len=256) :: out_final     ! final snapshot output
  character(len=256) :: out_temp      ! temporary snapshot output
  character(len=256) :: out_file_form ! format of snapshot files
  character(len=256) :: out_temp1     ! Temporary snapshot 1
  character(len=256) :: out_temp2     ! Temporary snapshot 2
  character(len=256) :: restart_log   ! Restart file that records last snapshot
  character(len=256) :: param_file    ! parameter file
#if defined(DEBUG_PLOT_DATA) 
  real(kind=PR) :: rzero(1:NDIM)      ! Position of origin for debug output
#endif
#if defined(RAD_WS)
 character(len=256) :: eos_opa_file   ! filename of table with eos and opacities
#endif

  type seren_param                             ! Seren parameter structure
     character(len=256)          :: var_name   ! Variable name
     character(len=1)            :: var_type   ! Variable type
     logical                     :: done       ! Has variable been read?

     character(len=256), pointer :: var_c      ! 256-character string
     character(len=20), pointer  :: var_u      ! 20-character unit
     logical, pointer            :: var_l      ! logical
     integer, pointer            :: var_i      ! integer
     integer(kind=ILP), pointer  :: var_j      ! long integer
     real(kind=PR), pointer      :: var_r      ! PR real
     real(kind=DP), pointer      :: var_d      ! DP real
  end type seren_param
  type(seren_param) :: params(1:256)           ! Main parameter array

END MODULE filename_module


! ============================================================================
MODULE HP_module
  use definitions
  use healpix_types

  logical :: HP_ionize                    ! Ionizing radiation?
  integer :: imax                         ! Maximum no. of rays
  integer :: lmax_hp                      ! Max number of HEALPix levels
  integer :: ltot_hp                      ! ..
  integer :: HPtot                        ! No. of UV sources
  integer(kind=i4b) :: x2pix(1:128)       ! ..
  integer(kind=i4b) :: y2pix(1:128)       ! ..
  integer(kind=i4b) :: pix2x(0:1023)      ! ..
  integer(kind=i4b) :: pix2y(0:1023)      ! ..
  real(kind=DP) :: a_star                 ! Recombination coefficient
  real(kind=PR) :: f1                     ! Integration step factor
  real(kind=PR) :: f2                     ! Opening criterion factor
  real(kind=PR) :: f3                     ! Temperature smoothing factor
  real(kind=PR) :: f4                     ! Density interpolation factor
  real(kind=PR) :: HPmaxres               ! ..
  real(kind=DP) :: intmax                 ! Value of integral at IF
  real(kind=DP) :: N_LyC                  ! Flux of UV photons
  real(kind=PR) :: rstatic(1:3)           ! Position of single static source
  real(kind=PR) :: Tion                   ! Temp. of the ionized gas
  real(kind=PR) :: Tneut                  ! Temp. of the neutral gas
  real(kind=PR) :: Xfrac                  ! ..
  real(kind=PR) :: Yfrac                  ! ..

  logical, allocatable :: ionizedo(:)     ! Ionize these particles
  integer, allocatable :: newtemp(:)      ! New temperature calculated
#if defined(DEBUG_HP_WALK_ALL_RAYS)
  integer, allocatable :: whichHPlevel(:) ! HP level of particle p
#endif

#if defined(HEALPIX)
  type HPlevel_node                            ! HEALPix level
     integer :: ifirst                         ! First ray on current level
     integer :: ilast                          ! Last ray on current level
  end type HPlevel_node
  type(HPlevel_node) :: HPlevel(0:HP_LEVELS)   ! HEALPix level info

  type HPray_node                              ! HEALPix ray 
     logical :: done                           ! Is this ray finished?
     integer(kind=I4B) :: ipix                 ! HEALPix pixel id
     integer :: first                          ! 1st particle in linked list
     integer :: last                           ! Last particle in lnked list
     integer :: rayend                         ! Final particle in ray
     real(kind=PR) :: hep                      ! Smoothing length at e.p.
     real(kind=PR) :: rhoep                    ! Density at e.p.
     real(kind=PR) :: gradrhoep(1:NDIM)        ! Density gradient
     real(kind=PR) :: rep(1:NDIM)              ! Position of e.p.
#if defined(IONIZING_UV_RADIATION)
     real(kind=PR) :: integral                 ! Current ionization integral
#endif
  end type HPray_node
  type(HPray_node), allocatable :: HPray(:)    ! HEALPix rays

  type HPsource_node                           ! HEALPix source
     integer :: sinkid                         ! Sink id of source
     integer :: Nlist                          ! Size of list order
     real(kind=DP) :: arot(1:3,1:3)            ! Rotation matrix
     real(kind=DP) :: Aangle                   ! Rotation angle
     real(kind=DP) :: Bangle                   !    "       "
     real(kind=DP) :: Cangle                   !    "       "
     real(kind=PR) :: r(1:NDIM)                ! Position of HP source
     integer, allocatable :: distorder(:)      ! Distance order array
#if defined(IONIZING_UV_RADIATION)
     real(kind=PR) :: intmax                   ! Stromgren integral value
     real(kind=DP) :: N_LyC                    ! Ionizing photons per sec.
#endif
  end type HPsource_node
  type(HPsource_node) :: HPsource(1:SMAX)      ! HP source array


#if defined(IONIZING_UV_RADIATION) && defined(MULTIPLE_SINK_SOURCES)
  integer :: Ntable                            ! No. of table elements
  type table_node                              ! Look-up table node
     real(kind=PR) :: mass                     ! Star mass (M_sun)
     real(kind=PR) :: log_L                    ! log10 (L/L_sun)
     real(kind=PR) :: log_N_LyC                ! log10 (N_LyC/s^{-1})
     real(kind=PR) :: Teff                     ! Effective temperature (K)
  end type table_node
  type(table_node), allocatable :: stellar_table(:)  ! Table of stellar props.
#endif
#endif


END MODULE HP_module


! ============================================================================
MODULE hydro_module
  use definitions

  real(kind=PR) :: alpha                       ! alpha viscosity parameter
  real(kind=PR) :: alpha_min                   ! alpha_min
  real(kind=PR) :: beta                        ! beta viscosity parameter
  real(kind=PR) :: gamma                       ! ratio of specific heats
  real(kind=PR) :: isotemp                     ! isothermal temperature
  real(kind=PR) :: Kpoly                       ! Polytropic constant
  real(kind=PR) :: mu_bar                      ! mean mol. weight
  real(kind=PR) :: newsound_const              ! sound speed scale
  real(kind=PR) :: Pconst                      ! ideal gas law scaling const.
  real(kind=PR) :: Pconst2                     !  "" for RAD_WS EoS.
  real(kind=PR) :: Pext                        ! External pressure
  real(kind=PR) :: rhobary                     ! barotropic density
  real(kind=PR) :: sound_const                 ! sound speed scale
#if defined(ARTIFICIAL_CONDUCTIVITY)
  real(kind=PR) :: alpha_cond                  ! Conductivity alpha
#endif

  real(kind=PR), allocatable :: div_v(:)       ! velocity divergence
  real(kind=PR), allocatable :: press(:)       ! pressure
  real(kind=PR), allocatable :: rho(:)         ! density
  real(kind=PR), allocatable :: sound(:)       ! sound speed
  real(kind=PR), allocatable :: temp(:)        ! temperature
  real(kind=PR), allocatable :: rho_old(:)     ! 'Old' density
  real(kind=PR), allocatable :: drhodt(:)      ! Rate of change of density
#if defined(GRAD_H_SPH)
  real(kind=PR), allocatable :: omega(:)       ! density correction
#endif
#if defined(SIGNAL_VELOCITY)
  real(kind=PR), allocatable :: vsigmax(:)     ! maximum signal velocity
#endif
#if defined(VISC_TD)
  real(kind=PR), allocatable :: talpha(:)      ! time-dependent alpha parameter
  real(kind=PR), allocatable :: talpha_old(:)  ! alpha at start of timestep
  real(kind=PR), allocatable :: dalpha_dt(:)   ! rate of change of alpha
#endif
#if defined(VISC_BALSARA)
  real(kind=PR), allocatable :: balsara(:)     ! Balsara factor
#endif
#if defined(VISC_PATTERN_REC)
  real(kind=PR), allocatable :: pattrec(:)     ! Pattern recognition factor
#endif 
#if defined(IONIZING_UV_RADIATION)
  real(kind=PR), allocatable :: gradrho(:,:)   ! Density gradient vector
  real(kind=PR), allocatable :: temp_min(:)    ! minimum temperature
  real(kind=PR), allocatable :: temp_aux(:)    ! aux. temp array
#endif
#if defined(DIV_A)
  real(kind=PR), allocatable :: div_a(:)       ! acceleration divergence
#endif

END MODULE hydro_module


! ============================================================================
MODULE kernel_module
  use definitions

  real(kind=PR), allocatable :: w0(:)   ! kernel function
  real(kind=PR), allocatable :: w1(:)   ! kernel derivative (orig or TC)
  real(kind=PR), allocatable :: w2(:)   ! kernel W*
  real(kind=PR), allocatable :: w3(:)   ! kernel W** (Combined with w2 now)
  real(kind=PR), allocatable :: w4(:)   ! Original kernel deriviative
  real(kind=PR), allocatable :: w5(:)   ! Modified W* kernel
  real(kind=PR), allocatable :: w6(:)   ! Modified W** kernel
#if defined(GRAD_H_SPH)
  real(kind=PR), allocatable :: wh(:)   ! kernel (w0) derivative w.r.t. h
  real(kind=PR), allocatable :: wg(:)   ! kernel (w3) derivative w.r.t. h
#endif

END MODULE kernel_module


! ============================================================================
MODULE mhd_module
  use definitions

#if defined(IDEAL_MHD)
  real(kind=PR), allocatable :: B(:,:)         ! B-field
  real(kind=PR), allocatable :: B_old(:,:)     ! B-field at start of timestep
  real(kind=PR), allocatable :: dB_dt(:,:)     ! Rate of change of B-field
#if defined(DEBUG_MHD)
  real(kind=PR), allocatable :: div_B(:)       ! Divergence of magnetic field
#endif
#endif

END MODULE mhd_module


! ============================================================================
MODULE Nbody_module
  use definitions

  integer :: nbin                    ! No. of mutiple systems in simulation
  real(kind=DP) :: nbody_endtime     ! Final runtime of N-body simulation
  real(kind=DP) :: nbody_frac        ! Frac. of gas removed before N-body sim
  real(kind=DP) :: nbody_lastsnap    ! Time of last nbody snapshot
  real(kind=DP) :: nbody_timemult    ! N-body timestep multiplier

  type star_node
     logical :: accdo                ! acceleration step?
     integer(kind=ILP) :: ncreate    ! nsteps when sink is created
     integer(kind=ILP) :: nlast      ! n of beginning of current timestep
     integer(kind=ILP) :: nlevel     ! timestep level of star
     real(kind=DP) :: tcreate        ! physical time when sink is created
     real(kind=DP) :: r(1:NDIM)      ! position vectors
     real(kind=DP) :: v(1:VDIM)      ! velocity vectors
     real(kind=DP) :: m              ! particle mass
     real(kind=DP) :: h              ! star softening length
     real(kind=DP) :: radius         ! sink accretion radius
     real(kind=DP) :: rold(1:NDIM)   ! old position vectors
     real(kind=DP) :: vold(1:VDIM)   ! old velocity vectors
     real(kind=DP) :: angmom(1:3)    ! Star internal angular momentum
     real(kind=DP) :: gpe            ! Gravitational potential energy
     real(kind=DP) :: gpot           ! Gravitational potential
     real(kind=DP) :: agravmag       ! Magnitude of grav. acceleration
     real(kind=DP) :: dmdt           ! Accretion rate
     real(kind=DP) :: a(1:VDIM)      ! Acceleration vectors
     real(kind=DP) :: a0(1:VDIM)     ! Acceleration at start of timestep
     real(kind=DP) :: adot(1:VDIM)   ! Jerk vector
     real(kind=DP) :: adot0(1:VDIM)  ! Jerk at start of timestep
     real(kind=DP) :: a2dot(1:VDIM)  ! 2nd deriv accel
     real(kind=DP) :: a2dot0(1:VDIM) ! 2nd deriv accel at start of timestep
     real(kind=DP) :: a3dot(1:VDIM)  ! 3rd deriv accel
     real(kind=DP) :: luminosity     ! luminosity from unresolved star
     real(kind=DP) :: temperature    ! surface temperature of unresolved star
     real(kind=DP) :: star_radius    ! physical radius of unresolved star
     real(kind=DP) :: macc(1:DMDT_RANGE)  ! Masses accreted in previous steps
     real(kind=DP) :: tacc(1:DMDT_RANGE)  ! Times of previous steps
  end type star_node
  type(star_node), allocatable   :: star(:)   ! Star array (N-body)

  type binary_node                   ! Binary star data structure
     integer :: id                   ! System id
     integer :: s1                   ! id of sink/binary 1
     integer :: s2                   ! id of sink/binary 2
     real(kind=DP) :: r(1:NDIM)      ! Position of COM of binary
     real(kind=DP) :: v(1:VDIM)      ! Velocity of COM of binary
     real(kind=DP) :: m              ! Total mass
     real(kind=DP) :: angmom(1:3)    ! Magnitude of angular momentum
     real(kind=DP) :: binen          ! Binding energy
     real(kind=DP) :: ecc            ! Eccentricity
     real(kind=DP) :: period         ! Period
     real(kind=DP) :: q              ! Mass-ratio (=m2/m1)
     real(kind=DP) :: sma            ! Semi-major axis
     real(kind=DP) :: drmag          ! Instantaneous dist. between components
  end type binary_node
  type(binary_node), allocatable :: binary(:)  ! Main binary array

END MODULE Nbody_module


! ============================================================================
MODULE neighbour_module
  use definitions

  integer :: pp_gather                 ! neighbours wanted
  integer :: pp_limit                  ! max. neighbours limit
  integer, allocatable :: pptot(:)     ! number of neighbours
  integer, allocatable :: pplist(:,:)  ! list of neighbours
  real(kind=PR) :: hmin                ! Minimum allowed smoothing length
  real(kind=PR) :: h_fac               ! grad-h density-h factor

END MODULE neighbour_module


! ============================================================================
MODULE particle_module
  use definitions

  logical :: com_frame                        ! flag to change to COM frame
  logical :: rho_remove                       ! Remove particle below min. rho?
  logical :: energy_remove                    ! Remove escaping particles?
  logical :: rad_remove                       ! Remove distant particles?
  integer :: pmax                             ! length of particle arrays
  integer :: ptot                             ! total number of particles
  integer :: rseed                            ! random number seed
  real(kind=PR) :: rholost                    ! Particle removal density
  real(kind=PR) :: rad_lost                   ! Particle removal radius
  real(kind=DP) :: mtot0                      ! Total initial amount of mass
  real(kind=DP) :: rcom(1:NDIM)               ! position of centre of mass
  real(kind=DP) :: vcom(1:VDIM)               ! velocity of centre of mass
  real(kind=DP) :: mtot                       ! total mass
  real(kind=DP) :: rcom0(1:NDIM)              ! Initial position of COM
  real(kind=DP) :: vcom0(1:VDIM)              ! Initial velocity of COM
#if defined(GRAD_H_SPH)
  real(kind=PR) :: rextent                    ! max extent in any 1 dimension
#endif
#if defined(REMOVE_OUTLIERS)
  real(kind=DP) :: mlost                      ! Mass removed from simulation
  real(kind=DP) :: rlost(1:NDIM)              ! COM of removed mass
  real(kind=DP) :: momlost(1:NDIM)            ! Momentum of removed mass
  real(kind=DP) :: angmomlost(1:3)            ! Ang mom. of removed mass
#endif

  integer, allocatable :: porig(:)            ! original particle identifier
  real(kind=PR), allocatable :: parray(:,:)   ! r, m and h particle data
                                              ! grouped into one array
  real(kind=PR), allocatable :: v(:,:)        ! velocity vectors
  real(kind=PR), allocatable :: a(:,:)        ! acceleration vectors
  real(kind=PR), allocatable :: r_old(:,:)    ! old particle positions
  real(kind=PR), allocatable :: v_old(:,:)    ! old particle velocities
#if defined(RUNGE_KUTTA) || defined(LEAPFROG_KDK)
  real(kind=PR), allocatable :: v_half(:,:)   ! half-step velocities
#endif
#if defined(PREDICTOR_CORRECTOR)
  real(kind=PR), allocatable :: a_old(:,:)    ! old particle accelerations
#endif
#if defined(SMOOTHED_VELOCITY)
  real(kind=PR), allocatable :: v_smooth(:,:) ! smoothed velocity
#endif
  real(kind=PR), allocatable :: gpot(:)       ! gravitational potential
#if defined(GRAVITY) && !defined(GEOMETRIC_MAC)
  real(kind=PR), allocatable :: agravmag(:)   ! mag. of grav. acceleration
#endif
#if defined(RAD_WS) && defined(GRAVITY)
  real(kind=PR), allocatable :: sphgpot(:)    ! grav. pot just due to SPH parts.
#endif
#if defined(INTERNAL_ENERGY)
  real(kind=PR), allocatable :: u(:)          ! internal energy
  real(kind=PR), allocatable :: u_old(:)      ! old internal energy
  real(kind=PR), allocatable :: du_dt(:)      ! cooling rate
#if defined(DIFFUSION)
  real(kind=PR), allocatable :: du_dt_diff(:) ! diffused energy rate
  real(kind=PR), allocatable :: k_cond(:)     ! thermal conductivity
  real(kind=PR), allocatable :: lambda_diff(:)    ! flux limiter
#endif
#endif
#if defined(ENTROPIC_FUNCTION)
  real(kind=PR), allocatable :: Aent(:)       ! entropic function
  real(kind=PR), allocatable :: Aold(:)       ! old entropic function
  real(kind=PR), allocatable :: dA_dt(:)      ! rate of change of entropic fn.
#endif
#if defined(SINKS) && defined(GRAVITY)
  logical, allocatable :: ispotmin(:)         ! Does this particle have the min
                                              ! potential of its neibs?
#endif
#if defined(DEBUG_FORCES)
  real(kind=PR), allocatable :: a_hydro(:,:)  ! hydro acceleration
  real(kind=PR), allocatable :: a_grav(:,:)   ! gravitational acceleration
  real(kind=PR), allocatable :: a_mag(:,:)    ! magnetic acceleration
  real(kind=PR), allocatable :: a_visc(:,:)   ! viscous acceleration
#endif
#if defined(DEBUG_DUDTRAD)
  real(kind=PR), allocatable :: dudt_rad(:)   ! radiative cooling rate
#endif


END MODULE particle_module


! ============================================================================
MODULE periodic_module
  use definitions

! Periodic boundary variables
  integer :: psphere                   ! Id of sphere particle
  real(kind=PR) :: periodic_min(1:3)   ! Minimum extent of periodic box
  real(kind=PR) :: periodic_max(1:3)   ! Maximum extent of periodic box
  real(kind=PR) :: periodic_size(1:3)  ! Size of periodic box
  real(kind=PR) :: periodic_half(1:3)  ! Half-size of periodic box
  real(kind=PR) :: rspheremax          ! Radius of spherical mirror

END MODULE periodic_module


! ============================================================================
MODULE scaling_module
  use definitions

! Scaling unit strings
  character(len=20) :: runit         ! length unit
  character(len=20) :: munit         ! mass unit
  character(len=20) :: tunit         ! time unit
  character(len=20) :: vunit         ! velocity unit
  character(len=20) :: aunit         ! acceleration unit
  character(len=20) :: rhounit       ! density unit
  character(len=20) :: sigmaunit     ! column density unit
  character(len=20) :: Punit         ! pressure unit
  character(len=20) :: funit         ! force unit
  character(len=20) :: Eunit         ! energy unit
  character(len=20) :: momunit       ! momentum unit
  character(len=20) :: angmomunit    ! angular momentum unit
  character(len=20) :: angvelunit    ! angular vel unit
  character(len=20) :: dmdtunit      ! accretion rate unit
  character(len=20) :: Lunit         ! luminosity unit
  character(len=20) :: kappaunit     ! opacity unit
  character(len=20) :: Bunit         ! magnetic field unit
  character(len=20) :: Qunit         ! electronic charge unit
  character(len=20) :: Junit         ! current density unit
  character(len=20) :: uunit         ! specific internal energy unit
  character(len=20) :: tempunit      ! temperature unit

! Scaling conversion factors
  real(kind=DP) :: rscale            ! length scaling factor
  real(kind=DP) :: mscale            ! mass scaling factor
  real(kind=DP) :: tscale            ! time scaling factor
  real(kind=DP) :: vscale            ! velocity scaling factor
  real(kind=DP) :: ascale            ! acceleration scaling factor
  real(kind=DP) :: rhoscale          ! density scaling factor
  real(kind=DP) :: sigmascale        ! column density scaling factor
  real(kind=DP) :: Pscale            ! pressure scaling factor
  real(kind=DP) :: fscale            ! force scaling factor
  real(kind=DP) :: Escale            ! energy scaling factor
  real(kind=DP) :: momscale          ! momentum scaling factor
  real(kind=DP) :: angmomscale       ! angular momentum scaling factor
  real(kind=DP) :: angvelscale       ! angular velocity scaling factor
  real(kind=DP) :: dmdtscale         ! accretion rate scaling factor
  real(kind=DP) :: Lscale            ! luminosity scaling factor
  real(kind=DP) :: kappascale        ! opacity scaling factor
  real(kind=DP) :: Bscale            ! magnetic field scaling factor
  real(kind=DP) :: Qscale            ! electronic charge scaling factor
  real(kind=DP) :: Jscale            ! current density scaling factor
  real(kind=DP) :: uscale            ! specific internal energy scaling factor
  real(kind=DP) :: tempscale         ! temperature scaling factor

! Scaling factors for cgs units
  real(kind=DP) :: rcgs              ! length unit in cgs
  real(kind=DP) :: mcgs              ! mass unit in cgs
  real(kind=DP) :: tcgs              ! time unit in cgs
  real(kind=DP) :: vcgs              ! velocity unit in cgs
  real(kind=DP) :: acgs              ! acceleration unit in cgs
  real(kind=DP) :: rhocgs            ! density unit in cgs
  real(kind=DP) :: sigmacgs          ! column density unit in cgs
  real(kind=DP) :: Pcgs              ! pressure unit in cgs
  real(kind=DP) :: fcgs              ! force unit in cgs
  real(kind=DP) :: Ecgs              ! energy unit in cgs
  real(kind=DP) :: momcgs            ! momentum unit in cgs
  real(kind=DP) :: angmomcgs         ! angular momentum unit in cgs
  real(kind=DP) :: angvelcgs         ! angular velocity unit in cgs
  real(kind=DP) :: dmdtcgs           ! accretion rate unit in cgs
  real(kind=DP) :: Lcgs              ! luminosity unit in cgs
  real(kind=DP) :: kappacgs          ! opacity unit in cgs
  real(kind=DP) :: Bcgs              ! magnetic field unit in cgs
  real(kind=DP) :: Qcgs              ! electronic charge unit in cgs
  real(kind=DP) :: Jcgs              ! current density unit in cgs
  real(kind=DP) :: ucgs              ! specific internal energy unit in cgs
  real(kind=DP) :: tempcgs           ! temperature unit in cgs

! Scaling factors for S.I. units
  real(kind=DP) :: r_SI              ! length unit in S.I.
  real(kind=DP) :: m_SI              ! mass unit in S.I.
  real(kind=DP) :: t_SI              ! time unit in S.I.
  real(kind=DP) :: v_SI              ! velocity unit in S.I.
  real(kind=DP) :: a_SI              ! acceleration unit in S.I.
  real(kind=DP) :: rho_SI            ! density unit in S.I.
  real(kind=DP) :: sigma_SI          ! column density unit in S.I.
  real(kind=DP) :: P_SI              ! pressure unit in S.I.
  real(kind=DP) :: f_SI              ! force unit in S.I.
  real(kind=DP) :: E_SI              ! energy unit in S.I.
  real(kind=DP) :: mom_SI            ! momentum unit in S.I.
  real(kind=DP) :: angmom_SI         ! angular momentum unit in S.I.
  real(kind=DP) :: angvel_SI         ! angular velocity unit in S.I.
  real(kind=DP) :: dmdt_SI           ! accretion rate unit in S.I.
  real(kind=DP) :: L_SI              ! luminosity unit in S.I.
  real(kind=DP) :: kappa_SI          ! opacity unit in S.I.
  real(kind=DP) :: B_SI              ! magnetic field unit in S.I.
  real(kind=DP) :: Q_SI              ! electronic charge unit in S.I.
  real(kind=DP) :: J_SI              ! current density unit in S.I.
  real(kind=DP) :: u_SI              ! specific internal energy unit in S.I.
  real(kind=DP) :: temp_SI           ! temperature unit in S.I.

END MODULE scaling_module


! ============================================================================
MODULE seren_sim_module
  use definitions

  logical :: nbody_sim               ! Flag to indicate N-body phase
  logical :: nbody_sph_sim           ! Flag to indicate N-body phase
  logical :: sph_sim                 ! Flag to indicate N-body phase

END MODULE seren_sim_module


! ============================================================================
MODULE sink_module
  use definitions

  logical :: rho_search                 ! Use density to identify sinks?
  logical :: potmin_search              ! Should sink be bottom of pot. well?
  logical :: hill_sphere_search         ! Use hill sphere for selecting sinks
  logical :: energy_search              ! Should new sink be grav. bound?
  logical :: div_v_search               ! Should new sink have div_v < 0?
  logical :: div_a_search               ! Should new sink have div_a < 0?
  logical :: timescale_search           ! Compare timescales for selecting sinks
  logical :: energy_accrete             ! Should accreted particle be bound?
  logical :: first                      ! ???
  logical :: accdo_sinks                ! accdo variable for all sinks
  integer :: stot                       ! Total number of sinks
  integer(kind=ILP) :: nlast_sinks      ! Beginning of current sink timestep
  integer(kind=ILP) :: nlevel_sinks     ! Step size of current sink timestep
  real(kind=PR) :: alpha_ss             ! Sunyaev-Shakura viscosity parameter
  real(kind=PR) :: f_accretion          ! Frac. of accretion lum. radiated away
  real(kind=PR) :: feedback_tdelay      ! Time delay to switch on 
                                        ! stellar radiative feeback
  real(kind=PR) :: feedback_minmass     ! Min. sink mass for switching on 
                                        ! stellar radiative feeback
  real(kind=PR) :: star_radius          ! radius of the protostar
  real(kind=PR) :: BD_radius          ! radius of the brown dwarfs
  real(kind=PR) :: planet_radius          ! radius of the planets
  real(kind=DP) :: laststep_sinks       ! Previous sink timestep
  real(kind=PR) :: rhosink              ! Sink creation density
  real(kind=PR) :: sinkrad              ! Sink creation radius
  real(kind=PR) :: sink_frac            ! Frac. of gas mass in sinks
  real(kind=DP) :: smooth_accrete_frac  ! Critical mass frac. before accretion
  real(kind=DP) :: smooth_accrete_dt    ! Critical timestep before accretion
  real(kind=DP) :: sink_dt              ! sink timestep
#ifdef EPISODIC_ACCRETION
     real(kind=DP) :: alpha_EA           ! episodic accretion effective viscosity alpha
     real(kind=DP) :: dmdt_regular      ! regular accretion rate (NOT EA -- Msun/yr)
     real(kind=DP) :: episodic_time=0.0_DP
#endif

  type sink_node
     logical :: accrete                 ! Is the sink accreting particles?
     logical :: static                  ! Is the sink static or not?
     integer :: id                      ! Sink id
     integer :: ncreate                 ! nsteps when sink is created
     real(kind=DP) :: tcreate           ! Physical time when sink is created
     real(kind=PR) :: r(1:NDIM)         ! Position vectors
     real(kind=PR) :: v(1:VDIM)         ! Velocity vectors
     real(kind=PR) :: m                 ! Particle mass
     real(kind=PR) :: h                 ! Sink softening length
     real(kind=PR) :: radius            ! Sink accretion radius
     real(kind=PR) :: a(1:VDIM)         ! Acceleration vectors
     real(kind=PR) :: rold(1:NDIM)      ! Old position vectors
     real(kind=PR) :: vold(1:VDIM)      ! Old velocity vectors
     real(kind=DP) :: angmom(1:3)       ! Sink internal angular momentum
     real(kind=DP) :: angmomnet(1:3)    ! Total sink ang mom accreted
     real(kind=PR) :: gpot              ! Gravitational potential energy
     real(kind=PR) :: gpe               ! Gravitational potential energy
     real(kind=DP) :: dmdt              ! Accretion rate
     real(kind=DP) :: luminosity        ! Luminosity from unresolved star
     real(kind=DP) :: luminosity_old    ! Luminosity from unresolved star (previous step)
     real(kind=DP) :: temperature       ! Surface temp. of unresolved star
     real(kind=DP) :: star_radius       ! Physical radius of unresolved star
     real(kind=DP) :: macc(1:DMDT_RANGE) ! Masses accreted in previous steps
     real(kind=DP) :: tacc(1:DMDT_RANGE) ! Times of previous steps
     real(kind=DP) :: Mstar           ! mass of the real star
     real(kind=DP) :: dmdt_star       ! accretion onto the real star
#ifdef EPISODIC_ACCRETION
     real(kind=DP) :: dmdt_0          ! initial accretion rate onto the real star
     integer       :: accretion_flag=0  ! 1: episodic accretion is happening
     real(kind=DP) :: Mdisc           ! mass of the disc (for the episodic accretion model)
     real(kind=DP) :: t_episode_start ! starting time of the current episodic accretion event
     real(kind=DP) :: t_episode_duration ! duration of the current episodic accretion event
#endif EPISODIC_ACCRETION

#if defined(SMOOTH_ACCRETION) || defined(SINK_REMOVE_ANGMOM)
     real(kind=PR) :: mmax           ! Accretion mass limit for sink
     real(kind=DP) :: trot           ! Rotational period at edge of sink
     real(kind=DP) :: tvisc          ! ..
     real(kind=DP) :: menc           ! Mass enclosed by sink (sink + gas)
#endif
#if defined(SMOOTH_ACCRETION) || defined(iINK_REMOVE_ANGMOM)
     real(kind=PR) :: cmean          ! Mean sound speed of particles in sink
#endif
#if defined(HEALPIX)
     integer :: HPid                 ! id of corresponding UV source
#endif
#if defined(RUNGE_KUTTA) || defined(LEAPFROG_KDK)
     real(kind=PR) :: vhalf(1:VDIM)  ! Half-step velocity arrays
#endif
#if defined(PREDICTOR_CORRECTOR)
     real(kind=PR) :: aold(1:VDIM)   ! Old accelerations for PC
#endif
#if !defined(GEOMETRIC_MAC)
     real(kind=PR) :: agravmag       ! Mag. of grav. acceleration
#endif
#if defined(DEBUG_FORCES)
     real(kind=PR) :: agrav(1:VDIM)  ! gravitational acceleration
     real(kind=PR) :: ahydro(1:VDIM) ! hydro force (of accreted particles)
#endif
  end type sink_node
  type(sink_node), allocatable :: sink(:)  ! Main sink data structure

END MODULE sink_module


! ============================================================================
MODULE time_module
  use definitions

  integer(kind=ILP) :: level_max     ! Maximum timestep level occupied
  integer(kind=ILP) :: level_step    ! Level of smallest (half) step
  integer(kind=ILP) :: n             ! Current integer time
  integer(kind=ILP) :: nbuild        ! Int. time for next tree build
  integer(kind=ILP) :: nbuildstep    ! Int. time interval between tree builds
  integer(kind=ILP) :: ndiagnext     ! Integer time for next diagnostic
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: ndiagstep     ! Integer steps between diagnostic
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: nionall       ! Int. time for next 'all' ionization
  integer(kind=ILP) :: nionallstep   ! Int. time between 'all' ionization calc.
  integer(kind=ILP) :: nionize       ! Int. time for next HEALPix walk
  integer(kind=ILP) :: nlevels       ! Number of quantized timestep levels
  integer(kind=ILP) :: nresync       ! Integer time to resynchronise
  logical(kind=ILP) :: sync_flag     ! Syncronise flag (sync if 1)
  integer(kind=ILP) :: sync_steps    ! for how many steps should I sync
  integer(kind=ILP) :: nsearchnext   ! Integer time for next sink search
  integer(kind=ILP) :: nsearchstep   ! Int. time interval between sink searches
  integer(kind=ILP) :: nsinknext     ! Integer time for next sink output
  integer(kind=ILP) :: nsinkstep     ! Int. time interval between sink output
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: nsnapnext     ! Integer time for next snapshot
  integer(kind=ILP) :: nsnapstep     ! Int. time interval between snapshots
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: nspare        ! Number of spare levels
  integer(kind=ILP) :: nsteps        ! Number of steps (Euler) or
                                     ! half-steps (RK, LF..)
  integer(kind=ILP) :: nstepsize     ! Current integer step size
  integer(kind=ILP) :: nstock        ! Integer time for next tree stock
  integer(kind=ILP) :: ntempnext     ! Integer time for next temp snapshot
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: ntempstep     ! Integer steps between temp files
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: snapshot      ! Snapshot number

  real(kind=DP) :: accel_mult        ! Acceleration timestep multiplier
  real(kind=DP) :: courant_mult      ! Courant timestep multiplier
  real(kind=DP) :: dt_fixed          ! Reference timestep for creating levels
  real(kind=DP) :: dt_max            ! Maximum stepsize
  real(kind=DP) :: endtime           ! End time
  real(kind=DP) :: firstsnap         ! Time for first snapshot file
  real(kind=DP) :: lastsnap          ! Time of last snapshot
  real(kind=DP) :: nextsnap          ! Time for next snapshot
  real(kind=DP) :: sink_mult         ! Sink timestep multiplier
  real(kind=DP) :: snaptime          ! Snapshot interval (time)
  real(kind=DP) :: time              ! Physical time
  real(kind=DP) :: timestep          ! Real timestep of smallest integer step

  real(kind=DP) :: sph_endtime       ! End time of SPH simulation
  real(kind=DP) :: sph_endmass       ! End gas time of SPH simulation (as a fraction of initial GAS mass)
  real(kind=DP) :: sph_endrho        ! End rho time of SPH simulation (when max dens reaches rho)


  real(kind=DP) :: nbody_sph_endtime ! End time of N-body/SPH simulation

  logical, allocatable :: accdo(:)            ! Update particle properties
  integer(kind=ILP), allocatable :: nlast(:)  ! Integer time of last update
  integer(kind=ILP), allocatable :: nlevel(:) ! Integer time step level
!  integer(kind=ILP), allocatable :: nstep(:)  ! Integer time step size
  real(kind=DP), allocatable :: laststep(:)   ! Previous step size

#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
  integer(kind=ILP), allocatable :: nminneib(:) ! Minimum neighbour stepsize
#endif
#if defined(BINARY_TREE)
  integer(kind=ILP) :: nskeleton    ! Integer time for next skeleton build
#endif

END MODULE time_module


! ============================================================================
MODULE timing_module
  use definitions
  implicit none

  integer :: last_id                         ! Number of last marker
  integer :: last_itime                      ! Last integer time mark
  integer :: itime                           ! Total integer time taken
  integer :: mark_tot                        ! Toal number of markers
  integer :: iblock(1:NBLOCKS)               ! Integer time taken by each block
  integer(kind=ILP) :: ngravcomp             ! ..
  integer(kind=ILP) :: nhydrocomp            ! ..
  integer(kind=ILP) :: nsphcomp              ! ..
  real(kind=DP) :: last_rtime                ! Last real time mark
  real(kind=DP) :: rtime                     ! Total real time taken
  real(kind=DP) :: rblock(1:NBLOCKS)         ! Real time taken by each block
  character(len=40) :: marker_id(1:NBLOCKS)  ! Array of marker strings

END MODULE timing_module

! ============================================================================
MODULE analyse_module
  use definitions
  implicit none
  
  integer :: disc_nbins                  ! number of bins to put particles when calculating azimuthal averages for discs
  real(kind=DP) :: disc_in_radius        ! inner disc radius (AU)
  real(kind=DP) :: disc_out_radius       ! outer disc radius (AU)
  ! for planet discs
  integer :: pdisc_nbins                  ! number of bins to put particles when calculating azimuthal averages for discs
  real(kind=DP) ::pdisc_in_radius        ! inner disc radius (AU)
  real(kind=DP) :: pdisc_out_radius       ! outer disc radius (AU)
END MODULE analyse_module

! ============================================================================
MODULE tree_module
  use definitions

  real(kind=PR) :: thetamaxsqd       ! opening angle criterion squared
  real(kind=PR) :: abserror          ! absolute error parameter (Gadget MAC)


! BH tree variables
! ----------------------------------------------------------------------------
#if defined(BH_TREE)
  integer :: cmax_grav               ! maximum allowed no. of cells
  integer :: ctot_grav               ! total number of cells
  integer :: ltot_grav               ! bottom level of gravity tree
  integer :: cmax_hydro              ! maximum possible number of cells
  integer :: ctot_hydro              ! total number of cells
  integer :: ltot_hydro              ! bottom level of hydro tree

  integer, allocatable :: first_cell_grav(:)   ! id of first cell on level
  integer, allocatable :: last_cell_grav(:)    ! id of last cell on level
  integer, allocatable :: first_cell_hydro(:)  ! id of first cell on level
  integer, allocatable :: last_cell_hydro(:)   ! id of last cell on level

  type BHgrav_node                      ! BH gravity tree node
     integer :: leaf                    ! Number of particles if leaf
     integer :: plist(1:LEAFMAX)        ! List of particle ids
     integer :: nextcell                ! Next cell in list if not opening
     integer :: ifopen                  ! First child cell if opened
     real(kind=PR) :: r(1:NDIM)         ! Centre of mass of cell c
     real(kind=PR) :: m                 ! Total mass of cell c
     real(kind=PR) :: dminsqd           ! Cell-opening distance squared
#ifndef GEOMETRIC_MAC
     real(kind=PR) :: mac               ! Multipole-acceptance criterion
#endif
#if defined(QUADRUPOLE)
     real(kind=PR) :: q(1:NQUAD)        ! Quadrupole moment terms
#endif
#if defined(OCTUPOLE)
     real(kind=PR) :: s(1:NOCT)         ! Octupole moment terms
#endif
  end type BHgrav_node
  type(BHgrav_node), allocatable :: BHgrav(:)  ! BH gravity tree array

  type BHhydro_node                     ! BH hydro tree node
     integer :: leaf                    ! Number of particles if leaf
     integer :: nextcell                ! Next cell in list if not opening
     integer :: ifopen                  ! First child cell if opened
     integer :: plist(1:LEAFMAX)        ! List of particle ids
     real(kind=PR) :: r(1:NDIM)         ! Centre of position of cell c
     real(kind=PR) :: rmax              ! Distance of furthest particle in c
     real(kind=PR) :: hrangemax         ! Maximum extent of kernel
  end type BHhydro_node
  type(BHhydro_node), allocatable :: BHhydro(:) ! BH hydro tree array

  type BH_temp_node                     ! Temp tree node
     integer :: pfirst                  ! First particle in linked list
     integer :: plast                   ! Last particle in linked list
     integer :: leaf                    ! Number of particles if leaf
     integer :: nextcell                ! Next cell in list if not opening
     integer :: ifopen                  ! First child cell if opened
     integer :: plist(1:LEAFMAX)        ! List of particle ids
     integer :: childof(1:NCHILD)       ! List if child cell ids
     real(kind=PR) :: r(1:NDIM)         ! Centre of position of cell c
  end type BH_temp_node
  type(BH_temp_node), allocatable :: BHtemp(:)  ! Temporary gravity array

  integer, allocatable :: cellof(:)     ! Cell particle p is in
  integer, allocatable :: BHnextptcl(:) ! Next particle in linked list
  integer, allocatable :: whichchild(:) ! Child cell particle p is in

  type auxilary_node                     ! Auxilary cell node
     real(kind=PR) :: hmax              ! Max value of h in a cell
     real(kind=PR) :: bbmin(1:NDIM)     ! rmin
     real(kind=PR) :: bbmax(1:NDIM)     ! rmax
  end type auxilary_node
  type(auxilary_node), allocatable :: BHstock(:)  ! bounding box array


! Binary tree variables
! ----------------------------------------------------------------------------
#elif defined(BINARY_TREE)

  integer :: ctot                    ! total no. of cells
  integer :: leaf                    ! max. no. of particles in a leaf-cell
  integer :: link                    ! number of child-cells
  integer :: ltot                    ! total no. of levels
  integer :: pmax                    ! maximum capacity of skeleton
  integer :: pmin                    ! minimum capacity of skeleton
  integer, allocatable :: cp(:)      ! (leaf)cell occupied by particle
  integer, allocatable :: ccp(:)     ! child-cell occupied by particle
  integer, allocatable :: chld(:,:)  ! IDs of cell's child-cells
  integer, allocatable :: lc(:)      ! level of cell
  integer, allocatable :: pnxt(:)    ! Pointer to next particle
  integer, allocatable :: pq(:,:)    ! IDs of sorted values
  integer, allocatable :: prev(:)    ! previous ptcle in leaf-cell-chain

  type binary_node                   ! Binary tree node
     integer :: ifopen               ! First child cell
     integer :: nextcell             ! Next cell if not opening
     integer :: pfirst               ! First particle in cell
     integer :: plast                ! Final particle in cell
     real(kind=PR) :: r(1:NDIM)      ! Position of cell COM
     real(kind=PR) :: rh(1:2*NDIM)   ! Smoothing length bounding box
     real(kind=PR) :: hmax           ! Max. h-value of particles in c
     real(kind=PR) :: rmax           ! Max. particle distance from cell centre
     real(kind=PR) :: m              ! Total mass of cell c
     real(kind=PR) :: dminsqd        ! Cell-opening distance squared
#ifndef GEOMETRIC_MAC
     real(kind=PR) :: mac            ! Multipole-acceptance criterion
#endif
#if defined(QUADRUPOLE)
     real(kind=PR) :: q(1:6)         ! Quadrupole moment terms
#endif
#if defined(OCTUPOLE)
     real(kind=PR) :: s(1:NOCT)      ! Octupole moment terms
#endif
  end type binary_node
  type(binary_node), allocatable :: bin_tree(:)  ! Binary tree array

#if defined(DEBUG_FOLIATE)
  integer, allocatable :: clip(:)    ! ID of next leaf-cell
#endif

#endif


END MODULE tree_module


! ============================================================================
MODULE type_module
  use definitions

  integer :: pboundary           ! No. of boundary particles
  integer :: pcdm                ! No. of cold dark matter particles
  integer :: pdust               ! No. of dust particles
  integer :: pgas                ! No. of gas particles
  integer :: picm                ! No. of ICM particles
  integer :: pion                ! No. of ion particles

  integer :: pboundarystart      ! id of first boundary particle
  integer :: pboundaryend        ! id of last boundary particle
  integer :: pcdmstart           ! id of first cdm particle
  integer :: pcdmend             ! id of last cdm particle
  integer :: pduststart          ! id of first dust particles
  integer :: pdustend            ! id of last dust particles
  integer :: pgas_orig           ! Original no. of gas particles
  integer :: pgasstart           ! id of first gas particle
  integer :: pgasend             ! id of last gas particle
  integer :: picmstart           ! id of first ICM particle
  integer :: picmend             ! id of last ICM particle
  integer :: pionstart           ! id of first ion particle
  integer :: pionend             ! id of last ion particle
  integer :: pionized            ! No. of ionized particles in HII sim
  integer :: pgravityend         ! id of last particle in gravity loop
  integer :: pgravitystart       ! id of first particle in gravity loop
  integer :: phydroend           ! id of last particle in hydro loop
  integer :: phydrostart         ! id of first particle in hydro loop
  integer :: pionizestart        ! id of first particle in ionization routine

  integer :: phtreestart         ! id of first particle in neighbour tree
  integer :: phtreeend           ! id of last particle in neighbour tree
  integer :: pgtreestart         ! id of first particle in gravity tree
  integer :: pgtreeend           ! id of last particle in gravity tree
  integer :: psftreestart        ! id of first particle in 2nd fluid tree
  integer :: psftreeend          ! id of last particle in 2nd fluid tree

  real(kind=DP) :: mgas          ! Total gas mass
  real(kind=DP) :: mgas_orig     ! Total gas mass at start of simulation
  real(kind=DP) :: micm          ! Total mass of icm particles
  real(kind=DP) :: mboundary     ! Total mass of boundary particles
  real(kind=DP) :: mcdm          ! Total mass of cdm particles
  real(kind=DP) :: mmean         ! Mean SPH gas particle mass (N.B. not mu_bar)

END MODULE type_module


! ============================================================================


