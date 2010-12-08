# ----------------------------------------------------------------------------
# Seren Makefile version 1.0.0
# Date : 08/12/2010
# ----------------------------------------------------------------------------
F90                      = gfortran
VERSION_NO               = "1.0.0"
SRCDIR                   = $(PWD)/src
EXEDIR                   = $(PWD)
OPTIMISE                 = 3
OPENMP                   = 1
PROFILE                  = 0
DEBUG                    = 1
NDIM                     = 2
PRECISION                = SINGLE
INFILE_FORMAT            = ALL
OUTFILE_FORMAT           = ALL
PERIODIC                 = 0
X_BOUNDARY               = 0
Y_BOUNDARY               = 0
Z_BOUNDARY               = 0
SPHERICAL_WALL           = 0
CYLINDRICAL_WALL         = 0

# ----------------------------------------------------------------------------
# Simulation selection options
# ----------------------------------------------------------------------------
SPH_SIMULATION           = 1
NBODY_SPH_SIMULATION     = 0
NBODY_SIMULATION         = 0

# ----------------------------------------------------------------------------
# SPH simulation options
# ----------------------------------------------------------------------------
SPH                      = STANDARD
SPH_INTEGRATION          = LFKDK
KERNEL                   = M4TC
HFIND                    = NUMBER
MINIMUM_H                = 0
HYDRO                    = 1
THERMAL                  = ISOTHERMAL
SINK_POTENTIAL_WS        = 0
AMBIENT_HEATING_WS       = 0
SINK_HEATING_WS          = 0
FLUX_LIMITED_DIFFUSION   = 0
IONIZING_RADIATION       = 0
RIEMANN_SOLVER           = 0
ARTIFICIAL_VISCOSITY     = MON97
VISC_TD                  = 0
BALSARA                  = 0
PATTERN_REC              = 0
ARTIFICIAL_CONDUCTIVITY  = 0
EXTERNAL_PRESSURE        = 0
MHD                      = 0
INDUCTION_EQN            = 0
RESISTIVITY              = 0
EXTERNAL_FORCE           = 0
GRAVITY                  = KS
EWALD                    = 0
REMOVE_OUTLIERS          = 0

# ----------------------------------------------------------------------------
# Sink and N-body options
# ----------------------------------------------------------------------------
SINKS                    = 0
SINK_RADIUS              = HMULT
SINK_REMOVE_ANGMOM       = 0
SINK_GRAVITY_ONLY        = 0
NBODY_INTEGRATION        = HERMITE4
BINARY_STATS             = 0

# ----------------------------------------------------------------------------
# Tree options
# ----------------------------------------------------------------------------
TREE                     = BH
MULTIPOLE                = QUADRUPOLE
MAC                      = GADGET
REORDER                  = PARTICLES
CELL_WALK                = 0

# ----------------------------------------------------------------------------
# Misc. options
# ----------------------------------------------------------------------------
SORT                     = INSERTION
TIMESTEP                 = RESTRICTED
CHECK_NEIB_TIMESTEP      = 1
SIGNAL_VELOCITY_DT       = 0
NEIGHBOURLISTS           = 1
TIMING_CODE              = 1
DIMENSIONLESS            = 0
TEST                     = 0

# ----------------------------------------------------------------------------
# Debugging options
# ----------------------------------------------------------------------------
ifneq ($(DEBUG),0)
#DFLAGS += -DDIV_A
#DFLAGS += -DDEBUG_ACCRETE
#DFLAGS += -DDEBUG_ALLOCATE_MEMORY
#DFLAGS += -DDEBUG_BHTREEBUILD
#DFLAGS += -DDEBUG_BHTREESTOCK
#DFLAGS += -DDEBUG_BHTREEWALK
#DFLAGS += -DDEBUG_BHTREEGRAVITY
#DFLAGS += -DDEBUG_BINARY_PROPERTIES
#DFLAGS += -DDEBUG_BINARY_SEARCH
#DFLAGS += -DDEBUG_BLOCK_TIMESTEPS
#DFLAGS += -DDEBUG_COPY_PARTICLE_DATA
#DFLAGS += -DDEBUG_CREATE_SINK
#DFLAGS += -DDEBUG_CREATE_HP_SOURCE
#DFLAGS += -DDEBUG_DENSITY
DFLAGS += -DDEBUG_DIAGNOSTICS
#DFLAGS += -DDEBUG_DIV_V
#DFLAGS += -DDEBUG_DUDTRAD
#DFLAGS += -DDEBUG_ENERGY_EQN
#DFLAGS += -DDEBUG_FOLIATE
#DFLAGS += -DDEBUG_FORCES
#DFLAGS += -DDEBUG_FREEFALL
#DFLAGS += -DDEBUG_GATHER_NEIB
#DFLAGS += -DDEBUG_GET_NEIB
#DFLAGS += -DDEBUG_GRAD_H_SPH
#DFLAGS += -DDEBUG_HEAPSORT
#DFLAGS += -DDEBUG_HERMITE4
#DFLAGS += -DDEBUG_H_GATHER
#DFLAGS += -DDEBUG_H_GATHER_DENSITY
#DFLAGS += -DDEBUG_H_GUESS
#DFLAGS += -DDEBUG_HP_IF
#DFLAGS += -DDEBUG_HP_OUTPUT
#DFLAGS += -DDEBUG_HP_SPLIT_ACTIVE_RAYS
#DFLAGS += -DDEBUG_HP_WALK_ALL_RAYS
#DFLAGS += -DDEBUG_HP_WALK_RAY
#DFLAGS += -DDEBUG_INTEGRATE
#DFLAGS += -DDEBUG_KERNEL
#DFLAGS += -DDEBUG_MHD
#DFLAGS += -DDEBUG_NBODYSETUP
#DFLAGS += -DDEBUG_PARAMETERS
#DFLAGS += -DDEBUG_PLOT_DATA
#DFLAGS += -DDEBUG_RAD
#DFLAGS += -DDEBUG_REDUCE_TIMESTEP
#DFLAGS += -DDEBUG_REMOVE_OUTLIERS
#DFLAGS += -DDEBUG_RSPH_OUTPUT
#DFLAGS += -DDEBUG_SINKCORRECTION
#DFLAGS += -DDEBUG_SINK_REMOVE_ANGMOM
#DFLAGS += -DDEBUG_SINK_SEARCH
#DFLAGS += -DDEBUG_SINK_TIMESTEP
#DFLAGS += -DDEBUG_SKELETON
#DFLAGS += -DSMOOTHED_VELOCITY
#DFLAGS += -DDEBUG_SPH_UPDATE
#DFLAGS += -DDEBUG_SWAP_PARTICLE_DATA
#DFLAGS += -DDEBUG_TIMESTEP_SIZE
#DFLAGS += -DDEBUG_TRACK_PARTICLE
#DFLAGS += -DDEBUG_TREE_BUILD
#DFLAGS += -DDEBUG_TREE_GRAVITY
#DFLAGS += -DDEBUG_TREESTOCK
#DFLAGS += -DDEBUG_TREEWALK
#DFLAGS += -DDEBUG_TYPES
#DFLAGS += -DDEBUG_VISC_BALSARA
#DFLAGS += -DDEBUG_VISC_PATTERN_REC
#DFLAGS += -C=all
#DFLAGS += -C=undefined
#DFLAGS += -Wall -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace
CFLAGS += $(DFLAGS)
endif


# Now include the Makefile tail which contains all the processed options 
# and generates list of object files for compilation
# ----------------------------------------------------------------------------
include makefiletail.mk


# List of possible flags for making the code
# ----------------------------------------------------------------------------
# F90                     : FORTRAN compiler
#		            f95      = NAG f95 compiler
#                           g95      = free f95 compiler (not gnu)
#                           gfortran = gnu f95 compiler
#                           pgf90    = Portland group compiler (coma)
#                           pgf95    = Portland group compiler (iceberg)
#                           ifort    = Intel Fortran compiler

# SRCDIR                  : Absolute path to main Seren directory

# EXEDIR                  : Absolute path to location of Seren executable

# OPTIMISE                : Compiler optimisation level (0, 1, 2 or 3)

# OPENMP                  : OpenMP parallelisation (0 or 1)

# PROFILE                 : 0, 1 (by subroutine) or 2 (by line)

# DEBUG                   : 0 = Debug flags switched off
#                           1 = Debug flags switched on to level 1
#                           2 = Debug flags switched on to level 2
#                           3 = Debug flags switched on to level 3

# NDIM                    : Number of dimensions (1, 2 or 3)

# PRECISION               : Precision of floating point variables
#                           SINGLE    = ..
#                           DOUBLE    = ..
#                           QUADRUPLE = (not functioning yet)

# INFILE_FORMAT           : Available input file format options in SEREN
#                           ALL    = All possible file formats enabled
#                           DRAGON = DRAGON file format
#                           SEREN  = SEREN file format
#                           ASCII  = Simple ASCII file format

# OUTFILE_FORMAT          : Available output file format options in SEREN
#                           ALL    = All possible file formats enabled
#                           DRAGON = DRAGON file format
#                           SEREN  = SEREN file format

# PERIODIC                : Periodic boundary conditions (0 or 1)

# X_BOUNDARY              : Boundary-type in x-direction
#                           0        = No boundaries
#                           PERIODIC = Periodic boundaries
#                           WALL     = LHS and RHS walls in x-dimension
#                           WALL_LHS = LHS wall in x-dimension
#                           WALL_RHS = RHS wall in x-dimension

# Y_BOUNDARY              : as X_BOUNDARY
# Z_BOUNDARY              : as X_BOUNDARY

# SPHERICAL_WALL          : (0 or 1)

# CYLINDRICAL_WALL        : (0 or 1)

# SPH_SIMULATION          : (0 or 1)

# NBODY_SPH_SIMULATION    : (0 or 1)

# NBODY_SIMULATION        : (0 or 1)

# SPH                     : SPH mode
#                           STANDARD   = Traditional SPH (e.g. Monaghan 1992)
#                           GRAD_H_SPH = 'grad-h' conservative SPH 
#                                        (e.g. Price & Monaghan 2005)
#                           RPSPH      = 'Relative pressure' SPH (Abel 2010)

# INTEGRATOR              : Integration scheme used
#                           EULER = 1st order Euler scheme
#                           RK    = 2nd order Runge-Kutta
#                           LFKDK = 2nd order Kick-drift-kick Leap-frog
#                           LFDKD = 2nd order Drift-kick-drift Leap-frog
#                           PC    = 2nd order Predictor-Corrector

# KERNEL                  : Kernel used to calculate SPH quantities
#                           M4          = M4 kernel (Monaghan & Lattanzio 1985)
#                           M4TC        = Modified M4 kernel gradient 
#                                         (Thomas & Couchman 1992)
#                           QUINTIC     = Quintic polynomial kernel 
#                                         (Morris 1996)
#                           QUINTICTC   = As QUINTIC, but with modified grad.
#                           GAUSSIAN_3H = Gaussian kernel truncated at 3h

# HFIND                   : Method to calculate h in standard SPH formulation
#                           NUMBER   = h contains specified no. of neibs
#                           MASS     = h contains specified mass
#                           CONSTANT = Use a constant value of h

# MINIMUM_H               : Enforce a minimum smoothing length (0 or 1)

# HYDRO                   : Hydrodynamical forces (0 or 1)

# THERMAL                 : Equation of state
#                           ISOTHERMAL   = Isothermal EOS
#                           BAROTROPIC   = Barotropic EOS
#                           POLYTROPIC   = Polytropic EOS
#                           ENERGY_EQN   = Solve energy equation
#                           ENTROPY_EQN  = Solve entropy equations
#                           RAD_WS       = WS Radiative approximation scheme
#                           STELLAR_HEAT = Heating from sinks: T ~ R^(-1/2)
#                                          STELLAR_HEAT_1 for 1st sink only
#                                          STELLAR_HEAT_2 for 1st 2 sinks only
#                                          STELLAR_HEAT_ALL for all sinks
#                           STIFF        = Stiff equation of state

# SINK_POTENTIAL_WS       : Use sink potential in RAD_WS method (0 or 1)

# AMBIENT_HEATING_WS      : External heating source (0 or 1)

# SINK_HEATING_WS         : Options for heating gas from sinks
#                           0                   = No heating from sinks
#                           STAR_HEATING        = ..
#                           STAR_SIMPLE_HEATING = ..
#                           HDISC_HEATING = heating by central star in disc 
#                                (assumes one star/disc with disc in x-y plane)
#                                T = sqrt{To^2 [(R^2+Ro^2)/AU^2]^(-q)+Tinf^2}
#                                i.e. T ~ To*(R/Ro)^(-q), 
#                                where R is the distance in disc midplane 
#                                (To, Ro (in AU), Tinf set in params.dat)

# FLUX_LIMITED_DIFFUSION  : Hybrid flux-limited diffusion method (0 or 1)

# IONIZING_RADIAITON      : Include sources of ionizing radiation or not
#                           0                     = No sources of 
#                                                   ionizing radiation
#                           SINGLE_STATIC_SOURCE  = Single static source 
#                                                   located at position given  
#                                                   in params file.
#                           MULTIPLE_SINK_SOURCES = Sink particles become 
#                                                   sources of ionizing rad. 

# RIEMANN_SOLVER          : Use Riemann solver (0 or 1)

# ARTIFICIAL_VISCOSITY    : Artificial viscosity formulation
#                           0      = No artificial viscosity used
#                           AB     = Alpha-Beta (standard)
#                           MON97  = Monaghan Riemann-viscosity

# BALSARA                 : Balsara switch for artificial viscosity (0 or 1)

# VISC_TD                 : Time-dependent artificial viscosity (0 or 1)

# PATTERN_REC             : Keplerian pattern recognition visc. switch (0 or 1)

# ARTIFICIAL_CONDUCTIVITY : Artificial conductivity formulation
#                           PRICE2008   = Price (2008) conductivity
#                           WADSLEY2008 = Wadsley et al. (2008) conductivity

# EXTERNAL_PRESSURE       : Simple external pressure formulation (0 or 1)

# MHD                     : Magneto-hydrodynamics (0)

# INDUCTION_EQN           : (0)

# RESISTIVITY             : Artificial resistivity (0)

# EXTERNAL_FORCE          : External gravitational force options
#                           0        = No external gravitational force
#                           PLUMMER  = Plummer potential
#                           UDS      = Uniform-density sphere
#                           NFW1996  = Navarro, Frenk & White (1996) potential

# GRAVITY                 : Method of computing gravitational force
#		            0     = No gravitational forces computed
#		            KS    = Kernel-softened gravity for 2-body forces
#                           NBODY = Newton's grav. law for all 2-body forces

# SINK_GRAVITY_ONLY       : Sink gravity options
#                           0 = Self-gravity from all SPH gas/sink particles
#                           1 = No self-gravity; sink-gravity only 
#                              (experimental)

# MULTIPOLE               : Gravity tree multipole expansion options
#                           0          = Monopole moments only
#                           QUADRUPOLE = Expand to quadrupole order
#                           OCTUPOLE   = Expand to octupole order

# MAC                     : Multipole-acceptance criterion used
#                           GEOMETRIC = Original Barnes-Hut geometric criterion
#                           GADGET    = GADGET-style MAC
#                           GADGET2   = GADGET2-style MAC
#                           EIGEN     = Eigenvalue MAC

# EWALD                   : Ewald periodic gravity forces (0 or 1)

# SINKS                   : Sink particle options
#                           0           = No sinks used
#                           SIMPLE      = Simple sink particles; accretion only
#                           NO_ACC      = Simple sinks with no accretion
#                           SMOOTH_ACC  = Sinks with smooth accretion 
#                                         (experimental)

# SINK_RADIUS             : Method of selecting the sink radius for new sinks
#                           FIXED_ABSOLUTE = Absolute (constant) sink radius
#                           FIXED_HMULT    = Multiple of h 
#                                            (same sinkrad for all)
#                           HMULT          = Multiple of h 
#                                            (indiv. sinkrad values)

# SINK_REMOVE_ANGMOM      : Redistribute angular momentum from sinks (0 or 1)

# NBODY_INTEGRATION       : Integration scheme for N-body sim
#                           HERMITE4 = 4th order Hermite scheme

# BINARY_STATS            : Automatically calculate binary statistics for 
#                           stars during N-body integrator (0 or 1)

# REMOVE_OUTLIERS         : 0 or 1 (Experimental)

# TREE                    : Tree options
#                           0         = No tree (direct-sum only)
#                           BH        = Barnes-Hut Octal-spatial tree
#                           BINARY    = Binary-number tree (not fully working)

# REORDER                 : Options for re-ordering particles/tree
#                           0         = No reordering of particles
#                           PARTICLES = Only reorder particle array
#                           TREE      = Only reorder tree cells (disabled)
#                           BOTH      = Reorder particle and tree arrays 
#                                       (disabled)

# CELL_WALK               : 0 or 1 (Experimental)

# SORT                    : Sorting algorithm 
#                           INSERTION = Insertion sort
#                           HEAP      = Heapsort

# TIMESTEP                : Options for multiple-particle timestepping
#                           ADAPTIVE   = Block timestep levels adjusted at 
#                                        resync
#                           RESTRICTED = Timestep levels only take certain 
#                                        values (dtmax parameter times 
#                                        integer power of 2)
#                           FIXED      = Fixed block timestep levels

# CHECK_NEIB_TIMESTEP     : Options for reducing particles timesteps
#                           0 = No checking of neighbour timesteps
#                           1 = Only reduce timesteps at end of current 
#                               timestep
#                           2 = Reduce timestep in middle of current timestep 
#                               if required

# NEIGHBOURLISTS          : Store neighbour lists in memory (0 or 1)

# TIMING_CODE             : Use internal timing routines (0 or 1)

# DIMENSIONLESS           : Set all units to dimensionless (0 or 1)

# TEST                    : Specific test options
#                           0        = No test options
#                           SPIEGEL  = Spiegel (ref??) test
#                           FREEFALL = Freefall collapse test
#                           BINARY   = Orbitting binary stars test
#                           PLUMMER  = Plummer sphere stability test
#                           ENTROPY  = Entropy-core test

