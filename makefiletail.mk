# makefiletail.mk
# ============================================================================

# Other directories (for analysis routines)
# ----------------------------------------------------------------------------
#PGPLOTLIBS = -L$(PGPLOT_DIR) -lpgplot -L/sw/lib -lpng -laquaterm -Wl,-framework -Wl,Foundation -lSystemStubs #-lz
PGPLOTLIBS = -L$(PGPLOT_DIR) -lpgplot -lpng
X11LIBS = -L/usr/X11R6/lib -lX11


# Sub-directory linkage
# ----------------------------------------------------------------------------
MODULEDIR = $(SRCDIR)
HEADERS = $(SRCDIR)/headers
VPATH = $(SRCDIR)/advance:$(SRCDIR)/analyse:$(SRCDIR)/BHtree:$(SRCDIR)/binarytree:$(SRCDIR)/dobs:$(SRCDIR)/gravity:$(SRCDIR)/headers:$(SRCDIR)/ic:$(SRCDIR)/io:$(SRCDIR)/healpix:$(SRCDIR)/main:$(SRCDIR)/mhd:$(SRCDIR)/nbody:$(SRCDIR)/nbody_sim:$(SRCDIR)/nbody_sph_sim:$(SRCDIR)/radiation:$(SRCDIR)/setup:$(SRCDIR)/sinks:$(SRCDIR)/sorts:$(SRCDIR)/sph:$(SRCDIR)/sph_sim:$(SRCDIR)/tests:$(SRCDIR)/timestep:$(SRCDIR)/user


# Remove trailing whitespace from user options
# ----------------------------------------------------------------------------
F90                     := $(strip $(F90))
OPTIMISE                := $(strip $(OPTIMISE))
OPENMP                  := $(strip $(OPENMP))
PROFILE                 := $(strip $(PROFILE))
DEBUG                   := $(strip $(DEBUG))
NDIM                    := $(strip $(NDIM))
PRECISION               := $(strip $(PRECISION))
INFILE_FORMAT           := $(strip $(INFILE_FORMAT))
OUTFILE_FORMAT          := $(strip $(OUTFILE_FORMAT))
PERIODIC                := $(strip $(PERIODIC))
X_BOUNDARY              := $(strip $(X_BOUNDARY))
Y_BOUNDARY              := $(strip $(Y_BOUNDARY))
Z_BOUNDARY              := $(strip $(Z_BOUNDARY))
SPHERICAL_WALL          := $(strip $(SPHERICAL_WALL))
CYLINDRICAL_WALL        := $(strip $(CYLINDRICAL_WALL))
SPH_SIMULATION          := $(strip $(SPH_SIMULATION))
NBODY_SPH_SIMULATION    := $(strip $(NBODY_SPH_SIMULATION))
NBODY_SIMULATION        := $(strip $(NBODY_SIMULATION))
SPH                     := $(strip $(SPH))
ANALYSE                 := $(strip $(ANALYSE))
SPH_SPECIFIC_OUTPUT     := $(strip $(SPH_SPECIFIC_OUTPUT))
SPH_OUTPUT_DENS	        := $(strip $(SPH_OUTPUT_DENS))
SPH_OUTPUT_TEMP	        := $(strip $(SPH_OUTPUT_TEMP))
SINK_PROPERTIES_FIX     := $(strip $(SINK_PROPERTIES_FIX))
EPISODIC_ACCRETION      :=$(strip $(EPISODIC_ACCRETION))
SPH_INTEGRATION         := $(strip $(SPH_INTEGRATION))
KERNEL                  := $(strip $(KERNEL))
HFIND                   := $(strip $(HFIND))
MINIMUM_H               := $(strip $(MINIMUM_H))
HYDRO                   := $(strip $(HYDRO))
THERMAL                 := $(strip $(THERMAL))
SINK_POTENTIAL_WS       := $(strip $(SINK_POTENTIAL_WS))
AMBIENT_HEATING_WS      := $(strip $(AMBIENT_HEATING_WS))
SINK_HEATING_WS         := $(strip $(SINK_HEATING_WS))
FLUX_LIMITED_DIFFUSION  := $(strip $(FLUX_LIMITED_DIFFUSION))
IONIZING_RADIATION      := $(strip $(IONIZING_RADIATION))
RIEMANN_SOLVER          := $(strip $(RIEMANN_SOLVER))
ARTIFICIAL_VISCOSITY    := $(strip $(ARTIFICIAL_VISCOSITY))
VISCOSITY_RECEEDING     := $(strip $(VISCOSITY_RECEEDING))
VISC_TD                 := $(strip $(VISC_TD))
BALSARA                 := $(strip $(BALSARA))
PATTERN_REC             := $(strip $(PATTERN_REC))
ARTIFICIAL_CONDUCTIVITY := $(strip $(ARTIFICIAL_CONDUCTIVITY))
MHD                     := $(strip $(MHD))
INDUCTION_EQN           := $(strip $(INDUCTION_EQN))
RESISTIVITY             := $(strip $(RESISTIVITY))
EXTERNAL_FORCE          := $(strip $(EXTERNAL_FORCE))
GRAVITY                 := $(strip $(GRAVITY))
EWALD                   := $(strip $(EWALD))
REMOVE_OUTLIERS         := $(strip $(REMOVE_OUTLIERS))
PARTICLE_ID       	:= $(strip $(PARTICLE_ID))
SINKS                   := $(strip $(SINKS))
KILLING_SINKS		:= $(strip $(KILLING_SINKS))
SINK_RADIUS             := $(strip $(SINK_RADIUS))
SINK_REMOVE_ANGMOM      := $(strip $(SINK_REMOVE_ANGMOM))
NBODY_INTEGRATION       := $(strip $(NBODY_INTEGRATION))
BINARY_STATS            := $(strip $(BINARY_STATS))
SINK_GRAVITY_ONLY       := $(strip $(SINK_GRAVITY_ONLY))
TREE                    := $(strip $(TREE))
MULTIPOLE               := $(strip $(MULTIPOLE))
MAC                     := $(strip $(MAC))
REORDER                 := $(strip $(REORDER))
CELL_WALK               := $(strip $(CELL_WALK))
SORT                    := $(strip $(SORT))
TIMESTEP                := $(strip $(TIMESTEP))
CHECK_NEIB_TIMESTEP     := $(strip $(CHECK_NEIB_TIMESTEP))
SIGNAL_VELOCITY_DT      := $(strip $(SIGNAL_VELOCITY_DT))
NEIGHBOURLISTS          := $(strip $(NEIGHBOURLISTS))
TIMING_CODE             := $(strip $(TIMING_CODE))
DIMENSIONLESS           := $(strip $(DIMENSIONLESS))
TEST                    := $(strip $(TEST))


# Object files always included in compilation list
# ----------------------------------------------------------------------------
MODULE_OBJ += definitions.o HP_types.o modules.o interface.o
SETUP_OBJ += convert_to_code_units.o default_parameters.o
SETUP_OBJ += initialize_seren_variables_1.o initialize_seren_variables_2.o
SETUP_OBJ += paramstore.o read_parameters.o read_arguments.o sanitycheck.o
SETUP_OBJ += seren_setup.o types.o units.o write_makefile_options.o
SETUP_OBJ += write_column_info.o
IO_OBJ += read_data.o write_data.o
IO_OBJ += write_rad_ws_test_data.o
IC_OBJ = ic_subroutines.o velfield.o
GENERIC_OBJ += allocate_memory.o clean_up.o COM.o distance3.o distance3_dp.o
GENERIC_OBJ += heapsort.o insertion_sort.o remove_from_list.o reorder_array.o 


# User object files (Include your own object files here)
# ----------------------------------------------------------------------------
USER_MOD += 
USER_SUB += 


# Compiler flags
# ----------------------------------------------------------------------------
ifeq ($(F90),f95)
CFLAGS += -DCOMPILER_F95
OPT += -I $(MODULEDIR) -I $(HEADERS) -C=all
ifeq ($(OPENMP),1)
OPT += -fopenmp -DOPENMP
endif

else ifeq ($(F90),g95)
CFLAGS += -DCOMPILER_G95
OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -Wall
ifeq ($(OPENMP),1)
OPT += -openmp -DOPENMP
endif

else ifeq ($(F90),pgf90)
CFLAGS += -DCOMPILER_PGF90
OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -Minline=distance,distance2,distance3,gravity_sph,gravity_gradh
#OPT += -fast -fastsse -Mnontemporal
ifeq ($(OPENMP),1)
OPT += -mp=allcores -DOPENMP
endif

else ifeq ($(F90),pgf95)
CFLAGS += -DCOMPILER_PGF95
OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -Minline=distance,distance2,distance3,gravity_sph,gravity_gradh
#OPT += -fast -fastsse -Mnontemporal
ifeq ($(OPENMP),1)
OPT += -mp=allcores -DOPENMP
endif

else ifeq ($(F90),ifort)
CFLAGS += -DCOMPILER_IFORT
OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -finline-functions -inline-forceinline -no-inline-min-size -inline-max-size=100 -ip -ipo #-vec-report=5
ifeq ($(OPENMP),1)
OPT += -openmp -DOPENMP
endif

else ifeq ($(F90),gfortran)
#F90 = /usr/bin/gfortran
CFLAGS += -DCOMPILER_GFORTRAN
OPT += -I $(MODULEDIR) -I $(HEADERS)
#OPT += -mtune=core2 -march=core2 #-mfpmath=sse -msse4.2 -ffast-math 
OPT += -Winline -fexpensive-optimizations -finline-functions -finline-limit=200
ifneq ($(PROFILE),0)
OPT += -Wall -pedantic -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace
endif
ifeq ($(OPENMP),1)
OPT += -fopenmp -DOPENMP
endif
ifeq ($(OPTIMISE),5)
OPT += -O5 -funroll-loops -ftree-vectorize
endif

else ifeq ($(F90),mpif90)
CFLAGS += -DCOMPILER_MPIF90
OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -Minline=distance,distance2,distance3,gravity_sph,gravity_gradh
OPT += -fast -fastsse -Mnontemporal
ifeq ($(OPENMP),1)
OPT += -mp -DOPENMP
endif

else ifeq ($(F90),sunf95)
CFLAGS += -DCOMPILER_SUNF95
OPT += -I$(MODULEDIR) -I$(HEADERS)
OPT += -ansi
OPT += -inline=%auto,distance2,distance3,gravity_sph
ifneq ($(DEBUG),0)
OPT += -C -w4 -u #-xcheck=%all
ifeq ($(OPENMP),1)
OPT += -vpara -loopinfo -DOPENMP
endif
endif
ifeq ($(OPENMP),1)
OPT += -openmp=parallel -autopar -stackvar -DOPENMP
endif

else
ERROR += "Invalid Fortran compiler selected : "$(F90)"\n"
endif


# Version no.
# ----------------------------------------------------------------------------
CFLAGS += -DVERSION_NO="$(VERSION_NO)"


# Optimisation level
# ----------------------------------------------------------------------------
ifeq ($(OPTIMISE),1)
CFLAGS += -DOPTIMISE=1
OPT += -O1
else ifeq ($(OPTIMISE),2)
CFLAGS += -DOPTIMISE=2
OPT += -O2
else ifeq ($(OPTIMISE),3)
CFLAGS += -DOPTIMISE=3
OPT += -O3
else ifeq ($(OPTIMISE),4)
CFLAGS += -DOPTIMISE=4
OPT += -O4
else ifneq ($(OPTIMISE),0)
ERROR += "Invalid value for OPTIMISE : "$(OPTIMISE)"\n"
endif


# Profiling options
# ----------------------------------------------------------------------------
ifeq ($(PROFILE),1)
CFLAGS += -DPROFILE=1
OPT += -pg
else ifeq ($(PROFILE),2)
CFLAGS += -DPROFILE=2
OPT += -pg -g
else ifeq ($(PROFILE),0)
CFLAGS += -DPROFILE=0
else
ERROR += "Invalid value for PROFILE : "$(PROFILE)"\n"
endif


# Dimensionality of the code
# ----------------------------------------------------------------------------
ifeq ($(NDIM),1)
CFLAGS += -DNDIM=1
else ifeq ($(NDIM),2)
CFLAGS += -DNDIM=2
else ifeq ($(NDIM),3)
CFLAGS += -DNDIM=3
else
ERROR += "Invalid value for NDIM : "$(NDIM)"\n"
endif


# Precision of real variables in code
# ----------------------------------------------------------------------------
ifeq ($(PRECISION),DOUBLE)
CFLAGS += -DDOUBLE_PRECISION
else ifeq ($(PRECISION),QUADRUPLE)
CFLAGS += -DQUADRUPLE_PRECISION
else ifneq ($(PRECISION),SINGLE)
ERROR += "Invalid PRECISION option : "$(PRECISION)"\n"
endif


# Input file format
# ----------------------------------------------------------------------------
ifeq ($(INFILE_FORMAT),ALL)
CFLAGS += -DDRAGON_INPUT -DSEREN_INPUT -DASCII_INPUT
IO_OBJ += read_data_ascii.o read_data_dragon_form.o read_data_dragon_unform.o
IO_OBJ += read_data_seren_form.o read_data_seren_unform.o
else ifeq ($(INFILE_FORMAT),DRAGON)
CFLAGS += -DDRAGON_INPUT
IO_OBJ += read_data_dragon_form.o read_data_dragon_unform.o
else ifeq ($(INFILE_FORMAT),SEREN)
CFLAGS += -DSEREN_INPUT
IO_OBJ += read_data_seren_form.o read_data_seren_unform.o
else ifeq ($(INFILE_FORMAT),ASCII)
CFLAGS += -DASCII_INPUT
IO_OBJ += read_data_ascii.o
else
ERROR += "Invalid INFILE_FORMAT option selected : "$(INFILE_FORMAT)"\n"
endif


# Output file format
# ----------------------------------------------------------------------------
ifeq ($(OUTFILE_FORMAT),ALL)
CFLAGS += -DDRAGON_OUTPUT -DSEREN_OUTPUT -DASCII_OUTPUT
IO_OBJ += write_data_ascii.o write_data_dragon_form.o write_data_dragon_unform.o
IO_OBJ += write_data_seren_form.o write_data_seren_unform.o
else ifeq ($(OUTFILE_FORMAT),DRAGON)
CFLAGS += -DDRAGON_OUTPUT
IO_OBJ += write_data_dragon_form.o write_data_dragon_unform.o 
else ifeq ($(OUTFILE_FORMAT),SEREN)
CFLAGS += -DSEREN_OUTPUT
IO_OBJ += write_data_seren_form.o write_data_seren_unform.o
else ifeq ($(OUTFILE_FORMAT),ASCII)
CFLAGS += -DASCII_OUTPUT
IO_OBJ += write_data_ascii.o
else
ERROR += "Invalid OUTFILE_FORMAT option selected : "$(OUTFILE_FORMAT)"\n"
endif


# Periodic boundary conditions
# ----------------------------------------------------------------------------
ifeq ($(PERIODIC),1)
CFLAGS += -DPERIODIC -DBOUNDARY_CONDITIONS
SPH_OBJ += check_boundary_conditions.o
else ifneq ($(PERIODIC),0)
ERROR += "Invalid PERIODIC option selected : "$(PERIODIC)"\n"
endif

ifeq ($(X_BOUNDARY),PERIODIC)
CFLAGS += -DPERIODIC_X
else ifeq ($(X_BOUNDARY),WALL)
CFLAGS += -DWALL_X_LHS -DWALL_X_RHS
else ifeq ($(X_BOUNDARY),WALL_LHS)
CFLAGS += -DWALL_X_LHS
else ifeq ($(X_BOUNDARY),WALL_RHS)
CFLAGS += -DWALL_X_RHS
else ifneq ($(X_BOUNDARY),0)
ERROR += "Invalid X_BOUNDARY option : "$(X_BOUNDARY)"\n"
endif

ifeq ($(Y_BOUNDARY),PERIODIC)
CFLAGS += -DPERIODIC_Y
else ifeq ($(Y_BOUNDARY),WALL)
CFLAGS += -DWALL_Y_LHS -DWALL_Y_RHS
else ifeq ($(Y_BOUNDARY),WALL_LHS)
CFLAGS += -DWALL_Y_LHS
else ifeq ($(Y_BOUNDARY),WALL_RHS)
CFLAGS += -DWALL_Y_RHS
else ifneq ($(Y_BOUNDARY),0)
ERROR += "Invalid Y_BOUNDARY option : "$(Y_BOUNDARY)"\n"
endif

ifeq ($(Z_BOUNDARY),PERIODIC)
CFLAGS += -DPERIODIC_Z
else ifeq ($(Z_BOUNDARY),WALL)
CFLAGS += -DWALL_Z_LHS -DWALL_Z_RHS
else ifeq ($(Z_BOUNDARY),WALL_LHS)
CFLAGS += -DWALL_Z_LHS
else ifeq ($(Z_BOUNDARY),WALL_RHS)
CFLAGS += -DWALL_Z_RHS
else ifneq ($(Z_BOUNDARY),0)
ERROR += "Invalid Z_BOUNDARY option : "$(Z_BOUNDARY)"\n"
endif


# Spherical mirror boundary
# ----------------------------------------------------------------------------
ifeq ($(SPHERICAL_WALL),1)
CFLAGS += -DSPHERICAL_WALL
SPH_OBJ += check_spherical_mirror.o
else ifneq ($(SPHERICAL_WALL),0)
ERROR += "Invalid SPHERICAL_WALL option : "$(SPHERICAL_WALL)"\n"
endif


# Cylindrical mirror boundary
# ----------------------------------------------------------------------------
ifeq ($(CYLINDRICAL_WALL),1)
CFLAGS += -DCYLINDRICAL_WALL
SPH_OBJ += check_cylindrical_mirror.o
else ifneq ($(CYLINDRICAL_WALL),0)
ERROR += "Invalid CYLINDRICAL_WALL option : "$(CYLINDRICAL_WALL)"\n"
endif


# Simulation-mode flags
# ----------------------------------------------------------------------------
ifeq ($(SPH_SIMULATION),1)
CFLAGS += -DSPH_SIMULATION
OBJ += sph_simulation.o sph_setup.o sph_integrate.o
OBJ += sph_output.o sph_timesteps.o syncronise.o
SETUP_OBJ += initialize_sph_variables_1.o initialize_sph_variables_2.o
INCLUDE_SPH_OBJS = 1
else ifneq ($(SPH_SIMULATION),0)
ERROR += "Invalid SPH_SIMULATION option : "$(SPH_SIMULATION)"\n"
endif

ifeq ($(NBODY_SPH_SIMULATION),1)
CFLAGS += -DNBODY_SPH_SIMULATION
OBJ += nbody_sph_simulation.o nbody_sph_setup.o nbody_sph_integrate.o
OBJ += nbody_sph_grav_forces.o
OBJ += nbody_sph_star_forces.o nbody_sph_output.o nbody_sph_timesteps.o
SETUP_OBJ += initialize_nbody_sph_variables_1.o
SETUP_OBJ += initialize_nbody_sph_variables_2.o
INCLUDE_SPH_OBJS = 1
INCLUDE_NBODY_OBJS = 1
else ifneq ($(NBODY_SPH_SIMULATION),0)
ERROR += "Invalid NBODY_SPH_SIMULATION option : "$(NBODY_SPH_SIMULATION)"\n"
endif

ifeq ($(NBODY_SIMULATION),1)
CFLAGS += -DNBODY_SIMULATION
NBODY_OBJ += nbody_simulation.o nbody_setup.o
NBODY_OBJ += nbody_grav_forces.o
NBODY_OBJ += nbody_integrate.o nbody_output.o nbody_timesteps.o
INCLUDE_NBODY_OBJS = 1
else ifneq ($(NBODY_SIMULATION),0)
ERROR += "Invalid NBODY_SIMULATION option : "$(NBODY_SIMULATION)"\n"
endif


# SPH object files
# ============================================================================
ifeq ($(INCLUDE_SPH_OBJS),1)
SPH_OBJ += all_sph.o bounding_box.o diagnostics.o distance2.o \
distance2_dp.o gather_neib_on_fly.o get_neib.o get_neib_on_fly.o \
h_guess.o reduce_particle_timestep.o sph_update.o \
timestep_size.o track_particles.o tree_update.o
IO_OBJ += record_particle_data.o write_data_debug.o write_data_grid_results.o
SETUP_OBJ += initialize_thermal_properties.o


# SPH mode
# ----------------------------------------------------------------------------
ifeq ($(SPH),GRAD_H_SPH)
CFLAGS += -DGRAD_H_SPH
SPH_OBJ += all_sph_gradh.o
else ifeq ($(SPH),RPSPH)
CFLAGS += -DRPSPH
else ifneq ($(SPH),STANDARD)
ERROR += "Invalid SPH option selected : "$(SPH)"\n"
endif


# Integration scheme used in code
# ----------------------------------------------------------------------------
SPH_OBJ += sph_advance.o advance_boundary_particle.o

ifeq ($(SPH_INTEGRATION),EULER)
CFLAGS += -DEULER
SPH_OBJ += advance_euler.o
else ifeq ($(SPH_INTEGRATION),RK)
CFLAGS += -DRUNGE_KUTTA
SPH_OBJ += advance_runge_kutta.o
else ifeq ($(SPH_INTEGRATION),LFKDK)
CFLAGS += -DLEAPFROG_KDK
SPH_OBJ += advance_leapfrog_kdk.o
else ifeq ($(SPH_INTEGRATION),LFDKD)
CFLAGS += -DLEAPFROG_DKD
SPH_OBJ += advance_leapfrog_dkd.o
else ifeq ($(SPH_INTEGRATION),PC)
CFLAGS += -DPREDICTOR_CORRECTOR
SPH_OBJ += advance_predictor_corrector.o
else 
ERROR += "Invalid SPH_INTEGRATION option selected : "$(SPH_INTEGRATION)"\n"
endif


# Method used to calculate h if not using 'grad-h' SPH
# ----------------------------------------------------------------------------
ifeq ($(HFIND),NUMBER)
CFLAGS += -DHGATHER
SPH_OBJ += h_gather.o
else ifeq ($(HFIND),MASS)
CFLAGS += -DHMASS
SPH_OBJ += h_gather_mass.o
else ifeq ($(HFIND),CONSTANT)
CFLAGS += -DCONSTANT_H
else ifeq ($(HFIND),H_RHO)
CFLAGS += -DH_RHO
ifneq ($(SPH),GRAD_H_SPH)
SPH_OBJ += all_sph_gradh.o gather_neib.o 
endif
else 
ERROR += "Invalid HFIND option selected : "$(HFIND)"\n"
endif


# Use a minimum smoothing length
# ----------------------------------------------------------------------------
ifeq ($(MINIMUM_H),1)
CFLAGS += -DMINIMUM_H
else ifneq ($(MINIMUM_H),0)
ERROR += "Invalid MINIMUM_H option selected : "$(MINIMUM_H)"\n"
endif


# Hydrodynamic forces
# ----------------------------------------------------------------------------
ifeq ($(HYDRO),1)
CFLAGS += -DHYDRO
OBJ += sph_hydro_forces.o
ifeq ($(MHD),0)
ifeq ($(SPH),GRAD_H_SPH)
SPH_OBJ += hydro_gradh.o
else
SPH_OBJ += hydro.o
endif
endif
else ifneq ($(HYDRO),0)
ERROR += "Invalid HYDRO option selected : "$(HYDRO)"\n"
endif


# Thermal physics 
# ----------------------------------------------------------------------------
SPH_OBJ += update_thermal_properties.o thermal.o

ifeq ($(THERMAL),ISOTHERMAL)
CFLAGS += -DISOTHERMAL
else ifeq ($(THERMAL),BAROTROPIC)
CFLAGS += -DBAROTROPIC
else ifeq ($(THERMAL),LOCAL_ISOTHERMAL)
CFLAGS += -DLOCAL_ISOTHERMAL
SPH_OBJ += ambienttemp.o 
else ifeq ($(THERMAL),POLYTROPIC)
CFLAGS += -DPOLYTROPIC
else ifeq ($(THERMAL),STIFF)
CFLAGS += -DSTIFF
else ifeq ($(THERMAL),STELLAR_HEAT_1)
CFLAGS += -DSTELLAR_HEAT
CFLAGS += -DSTARS=1
else ifeq ($(THERMAL),STELLAR_HEAT_2)
CFLAGS += -DSTELLAR_HEAT
CFLAGS += -DSTARS=2
else ifeq ($(THERMAL),STELLAR_HEAT_ALL)
CFLAGS += -DSTELLAR_HEAT
CFLAGS += -DSTARS=stot
else ifeq ($(THERMAL),ENERGY_EQN)
CFLAGS += -DINTERNAL_ENERGY
else ifeq ($(THERMAL),ENTROPY_EQN)
CFLAGS += -DENTROPIC_FUNCTION -DINTERNAL_ENERGY
aa
else ifeq ($(THERMAL),RAD_WS)
CFLAGS += -DRAD_WS -DRAD -DINTERNAL_ENERGY -DU_IMPLICIT_SOLVER
SPH_OBJ += rad_ws_update.o find_equilibrium_temp_ws.o 
SPH_OBJ += read_cooling_table_ws.o ambienttemp.o 
ifeq ($(FLUX_LIMITED_DIFFUSION),1)
CFLAGS += -DDIFFUSION
SPH_OBJ += conductivity.o diffusion.o
endif
ifeq ($(AMBIENT_HEATING_WS),1)
CFLAGS += -DAMBIENT_HEATING -DCONST_HEATING
endif
ifeq ($(SINK_HEATING_WS),HDISC_HEATING)
CFLAGS += -DHDISC_HEATING
endif
ifeq ($(SINK_HEATING_WS),HDISC_HEATING_PLUS_STAR_SIMPLE_HEATING)
CFLAGS += -DHDISC_HEATING_PLUS_STAR_SIMPLE_HEATING
endif
ifeq ($(SINK_HEATING_WS),HDISC_HEATING_3D_SINGLE)
CFLAGS += -DHDISC_HEATING_3D_SINGLE
endif
ifeq ($(SINK_HEATING_WS),STAR_HEATING)
CFLAGS += -DSTAR_HEATING
endif
ifeq ($(SINK_HEATING_WS),STAR_SIMPLE_HEATING)
CFLAGS += -DSTAR_SIMPLE_HEATING
endif
ifeq ($(SINK_POTENTIAL_WS),1)
CFLAGS += -DRAD_WS_SINK_POT
endif
else 
ERROR += "Invalid THERMAL option selected : "$(THERMAL)"\n"
endif


# HEALPix routines
# ----------------------------------------------------------------------------
ifneq ($(IONIZING_RADIATION),0)
CFLAGS += -DHEALPIX -DTRAPEZOIDAL_RULE
SPH_OBJ += HP_calculate_basis_vector.o HP_evaluation_point.o
SPH_OBJ += HP_initialize_source.o
SPH_OBJ += HP_inverse_positions.o HP_reorder_lists.o HP_rhoh_ep.o
SPH_OBJ += HP_split_active_rays.o HP_walk_all_rays.o HP_walk_ray.o
SPH_OBJ += create_HP_source.o healpix.o initialize_HP_sources.o
SPH_OBJ += write_ionization_data.o
endif

ifeq ($(IONIZING_RADIATION),SINGLE_STATIC_SOURCE)
CFLAGS += -DIONIZING_UV_RADIATION -DSINGLE_STATIC_SOURCE
SPH_OBJ += HP_ionizing_radiation.o
else ifeq ($(IONIZING_RADIATION),MULTIPLE_SINK_SOURCES)
CFLAGS += -DIONIZING_UV_RADIATION -DMULTIPLE_SINK_SOURCES
SPH_OBJ += HP_ionizing_radiation.o
else ifneq ($(IONIZING_RADIATION),0)
ERROR += "Invalid IONIZING_RADIATION option selected : "$(IONIZING_RADIATION)"\n"
endif


# Riemann solver
# ----------------------------------------------------------------------------
ifeq ($(RIEMANN_SOLVER),1)
CFLAGS += -DRIEMANN_SOLVER
SPH_OBJ += riemann_solver.o effective_gamma.o
ifeq ($(THERMAL),ISOTHERMAL)
CFLAGS += -DISOTHERMAL_SOLVER 
else
CFLAGS += -DITERATIVE_SOLVER 
endif
else ifneq ($(RIEMANN_SOLVER),0)
ERROR += "Invalid RIEMANN_SOLVER option selected : "$(RIEMANN_SOLVER)"\n"
endif


# Artificial viscosity
# ----------------------------------------------------------------------------
ifeq ($(ARTIFICIAL_VISCOSITY),AB)
CFLAGS += -DARTIFICIAL_VISCOSITY -DVISC_AB
else ifeq ($(ARTIFICIAL_VISCOSITY),MON97)
CFLAGS += -DARTIFICIAL_VISCOSITY -DVISC_MON97
else ifneq ($(ARTIFICIAL_VISCOSITY),0)
ERROR += "Invalid ARTIFICIAL_VISCOSITY option selected : "$(ARTIFICIAL_VISCOSITY)"\n"
endif

ifneq ($(ARTIFICIAL_VISCOSITY),0)
ifeq ($(VISCOSITY_RECEEDING),1)
CFLAGS += -DVISCOSITY_RECEEDING
endif
ifeq ($(VISC_TD),1)
CFLAGS += -DVISC_TD
else ifneq ($(VISC_TD),0)
ERROR += "Invalid VISC_TD option selected : "$(VISC_TD)"\n"
endif
ifeq ($(BALSARA),1)
CFLAGS += -DVISC_BALSARA
else ifneq ($(BALSARA),0)
ERROR += "Invalid BALSARA option selected : "$(BALSARA)"\n"
endif
ifeq ($(PATTERN_REC),1)
CFLAGS += -DVISC_PATTERN_REC
else ifneq ($(PATTERN_REC),0)
ERROR += "Invalid PATTERN_REC option selected : "$(PATTERN_REC)"\n"
endif
endif


# Artificial conductivity
# ----------------------------------------------------------------------------
ifeq ($(ARTIFICIAL_CONDUCTIVITY),PRICE2008)
CFLAGS += -DARTIFICIAL_CONDUCTIVITY -DCOND_PRICE2008
else ifeq ($(ARTIFICIAL_CONDUCTIVITY),WADSLEY2008)
CFLAGS += -DARTIFICIAL_CONDUCTIVITY -DCOND_WADSLEY2008
else ifneq ($(ARTIFICIAL_CONDUCTIVITY),0)
ERROR += "Invalid ARTIFICIAL_CONDUCTIVITY option selected\n"
endif


# External pressure flag
# ----------------------------------------------------------------------------
ifeq ($(EXTERNAL_PRESSURE),1)
CFLAGS += -DEXTERNAL_PRESSURE
else ifneq ($(EXTERNAL_PRESSURE),0)
ERROR += "Invalid EXTERNAL_PRESSURE option selected : "$(EXTERNAL_PRESSURE)"\n"
endif


# Magneto-hydrodynamics (MHD)
# ----------------------------------------------------------------------------
ifeq ($(MHD),IDEAL)
CFLAGS += -DIDEAL_MHD
SPH_OBJ += ideal_mhd_gradh.o
ifeq ($(INDUCTION_EQN),STANDARD)
CFLAGS += -DINDUCTION_EQUATION
endif
else ifeq ($(MHD),NONIDEAL)
CFLAGS += -DNON_IDEAL_MHD
else ifneq ($(MHD),0)
ERROR += "Invalid MHD option selected : "$(MHD)"\n"
endif


# Induction equation
# ----------------------------------------------------------------------------
ifeq ($(INDUCTION_EQN),STANDARD)
CFLAGS += -DINDUCTION_EQN
else ifneq ($(INDUCTION_EQN),0)
ERROR += "Invalid value for INDUCTION_EQN : "$(INDUCTION_EQN)"\n"
endif


# Artificial resistivity 
# ----------------------------------------------------------------------------
ifeq ($(RESISTIVITY),1)
CFLAGS += -DARTIFICIAL_RESISTIVITY
SPH_OBJ += artificial_resistivity.o
else ifneq ($(RESISTIVITY),0)
ERROR += "Invalid value for RESISTIVITY : "$(RESISTIVITY)"\n"
endif


# External gravitational force
# ----------------------------------------------------------------------------
ifeq ($(EXTERNAL_FORCE),PLUMMER)
CFLAGS += -DEXTERNAL_FORCE -DPLUMMER_POTENTIAL
OBJ += add_external_gravitational_force.o
else ifeq ($(EXTERNAL_FORCE),UDS)
CFLAGS += -DEXTERNAL_FORCE -DUDS_POTENTIAL
OBJ += add_external_gravitational_force.o
else ifeq ($(EXTERNAL_FORCE),NFW1996)
CFLAGS += -DEXTERNAL_FORCE -DNFW1996_POTENTIAL
OBJ += add_external_gravitational_force.o
else ifneq ($(EXTERNAL_FORCE),0)
ERROR += "Invalid EXTERNAL_FORCE option selected\n"
endif


# Gravity forces
# ----------------------------------------------------------------------------
SPH_OBJ += sph_grav_forces.o

ifeq ($(GRAVITY),KS)
CFLAGS += -DGRAVITY
SPH_OBJ += gravity_sph.o gravity_meanh.o
ifeq ($(SPH),GRAD_H_SPH)
SPH_OBJ += gravity_gradh.o
endif
else ifeq ($(GRAVITY),NBODY)
CFLAGS += -DGRAVITY -DN_BODY
SPH_OBJ += gravity_nbody.o
else ifneq ($(GRAVITY),0)
ERROR += "Invalid GRAVITY option selected\n"
endif

ifneq ($(GRAVITY),0)
ifneq ($(SINK_GRAVITY_ONLY),1)
CFLAGS += -DSPH_SELF_GRAVITY
else ifneq ($(SINK_GRAVITY_ONLY),0)
ERROR += "Invalid value for SINK_GRAVITY_ONLY : "$(SINK_GRAVITY_ONLY)"\n"
endif
endif


# Tree code
# ----------------------------------------------------------------------------
ifeq ($(TREE),BH)
CFLAGS += -DBH_TREE
SPH_OBJ += BHhydro_build.o BHhydro_stock.o BHhydro_update_hmax.o
SPH_OBJ += BHhydro_walk.o BHhydrowalk_hgather.o BHhydro_hguess.o
SPH_OBJ += BH_remove_particles.o
ifneq ($(GRAVITY),0)
SPH_OBJ += BHgrav_build.o BHgrav_stock.o BHgrav_accel.o
SPH_OBJ += copy_BHgrav_to_BHhydro.o
endif
ifeq ($(REORDER),TREE)
CFLAGS += -DREORDER_TREE
endif
ifeq ($(REORDER),PARTICLES)
CFLAGS += -DREORDER_PARTICLES
SPH_OBJ += BH_reorder_particles.o
endif
ifeq ($(REORDER),ALL)
CFLAGS += -DREORDER_TREE -DREORDER_PARTICLES
SPH_OBJ += BH_reorder_particles.o
endif
ifeq ($(CELL_WALK),1)
ifneq ($(GRAVITY),0)
CFLAGS += -DCELL_WALK
SPH_OBJ += BHgrav_grouped_walk.o
endif
else ifneq ($(CELL_WALK),0)
ERROR += "Invalid value for CELL_WALK : "$(CELL_WALK)"\n"
endif

else ifeq ($(TREE),BINARY)
CFLAGS += -DBINARY_TREE
SPH_OBJ += binary_foliate.o binary_neibfind.o
SPH_OBJ += binary_skeleton.o binary_treebuild.o binary_treestock.o
ifneq ($(GRAVITY),0)
SPH_OBJ += binary_gravacc.o
endif
ifeq ($(REORDER),1)
CFLAGS += -DREORDER
SPH_OBJ += swap_particle_data.o
endif

else ifeq ($(TREE),0)
ifneq ($(GRAVITY),0)
SPH_OBJ += direct_sph_gravity.o
endif

else 
ERROR += "Invalid TREE option selected : "$(TREE)"\n"
endif


# Sink-gravity only
# ----------------------------------------------------------------------------
ifneq ($(SINKS),0)
ifeq ($(SINK_GRAVITY_ONLY),1)
CFLAGS += -DSINK_GRAVITY_ONLY
endif
endif


# Multipole moment expansion in gravity tree
# ----------------------------------------------------------------------------
ifeq ($(MULTIPOLE),QUADRUPOLE)
CFLAGS += -DQUADRUPOLE
else ifeq ($(MULTIPOLE),OCTUPOLE)
CFLAGS += -DQUADRUPOLE -DOCTUPOLE
else ifneq ($(MULTIPOLE),0)
ERROR += "Invalid TREE option selected : "$(TREE)"\n"
endif


# Tree gravity MAC
# ----------------------------------------------------------------------------
ifeq ($(MAC),GEOMETRIC)
CFLAGS += -DGEOMETRIC_MAC
else ifeq ($(MAC),GADGET)
CFLAGS += -DGADGET_MAC
else ifeq ($(MAC),GADGET2)
CFLAGS += -DGADGET2_MAC
else ifeq ($(MAC),NEW)
CFLAGS += -DNEW_MAC
else ifeq ($(MAC),EIGEN)
CFLAGS += -DEIGEN_MAC
SPH_OBJ += eigenvalue_mac.o
else
ERROR += "Invalid value for MAC : "$(MAC)"\n"
endif


# SPH kernel-softeneing MAC
# ----------------------------------------------------------------------------
#ifeq ($(SPH_KS_MAC),1)
CFLAGS += -DSPH_KS_MAC
#else
#ifneq ($(GRAVITY),NBODY)
#SPH_OBJ += gravity_nbody.o
#endif
#endif


# Sink particles
# ----------------------------------------------------------------------------
ifneq ($(GRAVITY),0)
ifneq ($(SINKS),0)
CFLAGS += -DSINKS -DACCRETION_RATE -DDEBUG_FORCES
SPH_OBJ += sph_sink_forces.o direct_sink_gravity.o
SPH_OBJ += sink_advance.o sink_search.o sink_timestep.o sink_update.o
SPH_OBJ += write_accreted_particles.o write_sink_data.o

ifeq ($(SPH_INTEGRATION),EULER)
SPH_OBJ += advance_sink_euler.o
else ifeq ($(SPH_INTEGRATION),RK)
SPH_OBJ += advance_sink_RK.o
else ifeq ($(SPH_INTEGRATION),LFKDK)
SPH_OBJ += advance_sink_LFKDK.o
else ifeq ($(SPH_INTEGRATION),LFDKD)
SPH_OBJ += advance_sink_LFDKD.o
else ifeq ($(SPH_INTEGRATION),PC)
SPH_OBJ += advance_sink_PC.o
else ifneq ($(SPH_INTEGRATION),0)
ERROR += "Invalid SPH_INTEGRATION option selected : "$(SPH_INTEGRATION)"\n"
endif

ifeq ($(KILLING_SINKS),1)
CFLAGS +=-DKILLING_SINKS
endif

ifeq ($(SINKS),SMOOTH_ACC)
CFLAGS += -DSMOOTH_ACCRETION -DMINIMUM_H
SPH_OBJ += create_sink.o smooth_accrete_particles.o sink_accretion_properties.o
else
SPH_OBJ += create_sink.o accrete_particles.o sink_accretion_properties.o
endif

ifeq ($(SINK_RADIUS),FIXED_ABSOLUTE)
CFLAGS += -DFIXED_ABSOLUTE_SINKRAD
else ifeq ($(SINK_RADIUS),FIXED_HMULT)
CFLAGS += -DFIXED_HMULT_SINKRAD
else ifeq ($(SINK_RADIUS),HMULT)
CFLAGS += -DHMULT_SINKRAD
else ifneq ($(SINK_RADIUS),0)
ERROR += "Invalid SINK_RADIUS option selected : "$(SINK_RADIUS)"\n"
endif

ifeq ($(SINK_REMOVE_ANGMOM),1)
CFLAGS += -DSINK_REMOVE_ANGMOM #-DCHECK_NEIGHBOUR_TIMESTEPS -DIMMEDIATE_TIMESTEP_REDUCTION
SPH_OBJ += redistribute_sink_angmom.o
endif

ifeq ($(SINKS),NO_ACC)
CFLAGS += -DNO_ACCRETION
endif

endif
endif

# Reconstruct sink properties from dragon format files
# ----------------------------------------------------------------------------
ifneq ($(SINKS),0)
ifeq ($(SINK_PROPERTIES_FIX),1)
CFLAGS += -DSINK_PROPERTIES_FIX
OBJ +=read_sink_data.o
endif
endif

# Sink episodic accretion
# ----------------------------------------------------------------------------
ifneq ($(SINKS),0)
ifneq ($(SINK_HEATING_WS),0)
ifeq ($(EPISODIC_ACCRETION),1)
CFLAGS +=-DEPISODIC_ACCRETION -DSINK_PROPERTIES_FIX
OBJ +=episodic_accretion_model.o
endif
endif
endif

# Check neighbour's timesteps
# ----------------------------------------------------------------------------
ifeq ($(CHECK_NEIB_TIMESTEP),1)
CFLAGS += -DCHECK_NEIGHBOUR_TIMESTEPS
OBJ += check_neighbour_timesteps.o reduce_timesteps.o
else ifeq ($(CHECK_NEIB_TIMESTEP),2)
CFLAGS += -DCHECK_NEIGHBOUR_TIMESTEPS -DIMMEDIATE_TIMESTEP_REDUCTION
OBJ += check_neighbour_timesteps.o reduce_timesteps.o
else ifneq ($(CHECK_NEIB_TIMESTEP),0)
ERROR += "Invalid CHECK_NEIB_TIMESTEP option selected : "$(CHECK_NEIB_TIMESTEP)"\n"
endif


# Neighbour lists
# ----------------------------------------------------------------------------
ifeq ($(HYDRO),1)
ifeq ($(NEIGHBOURLISTS),1)
CFLAGS += -DNEIGHBOUR_LISTS
else ifneq ($(NEIGHBOURLISTS),0)
ERROR += "Invalid NEIGHBOURLISTS option selected : "$(NEIGHBOURLISTS)"\n"
endif
endif


endif
# ============================================================================



# Include N-body source files
# ============================================================================
ifeq ($(INCLUDE_NBODY_OBJS),1)
CFLAGS += -DSINKS
NBODY_OBJ += copy_sinks_to_stars.o copy_stars_to_sinks.o
NBODY_OBJ += nbody_advance.o nbody_correction_terms.o
NBODY_OBJ += nbody_timestep_size.o write_star_data.o
NBODY_OBJ += nbody_diagnostics.o
SETUP_OBJ += nbody_accrete_bound_particles.o nbody_hermite4_extra_terms.o

ifeq ($(NBODY_INTEGRATION),HERMITE4)
CFLAGS += -DNBODY_HERMITE4
NBODY_OBJ += gravity_hermite4_meanh.o
NBODY_OBJ += nbody_hermite4_direct_gravity.o
else
ERROR += "Invalid value for NBODY_INTEGRATION : "$(NBODY_INTEGRATION)"\n"
endif

ifeq ($(NBODY_SPH_SIMULATION),1)
ifeq ($(TREE),BH)
SPH_OBJ += BHgrav_accel_jerk.o
else
SPH_OBJ += sph_hermite4_direct_gravity.o
endif
endif

ifeq ($(NBODY_SIMULATION),1)
ifneq ($(SPH_SIMULATION),1)
SPH_OBJ += gravity_sph.o
endif
endif

ifeq ($(BINARY_STATS),1)
CFLAGS += -DBINARY_STATS
NBODY_OBJ += binary_energy.o binary_properties.o binary_search.o
else ifneq ($(BINARY_STATS),0)
ERROR += "Invalid value for BINARY_STATS : "$(BINARY_STATS)"\n"
endif
endif
# ============================================================================


# Kernel used to calculate SPH quantities
# ----------------------------------------------------------------------------
SPH_OBJ += kernel.o
ifeq ($(KERNEL),M4)
CFLAGS += -DM4_KERNEL
else ifeq ($(KERNEL),M4TC)
CFLAGS += -DM4_KERNEL -DTC_KERNEL
else ifeq ($(KERNEL),QUINTIC)
CFLAGS += -DQUINTIC_KERNEL
else ifeq ($(KERNEL),QUINTICTC)
CFLAGS += -DQUINTIC_KERNEL -DTC_KERNEL
else ifeq ($(KERNEL),GAUSSIAN_3H)
CFLAGS += -DGAUSSIAN_3H_KERNEL
else
ERROR += "Invalid KERNEL option selected : "$(KERNEL)"\n"
endif


# Ewald forces
# ----------------------------------------------------------------------------
ifeq ($(EWALD),1)
CFLAGS += -DEWALD
SPH_OBJ += ewald_init.o ewald_force.o
else ifneq ($(EWALD),0)
ERROR += "Invalid EWALD option selected : "$(EWALD)"\n"
endif


# Remove outliers
# ----------------------------------------------------------------------------
ifeq ($(REMOVE_OUTLIERS),1)
CFLAGS += -DREMOVE_OUTLIERS
SPH_OBJ += remove_outlying_particles.o
else ifneq ($(REMOVE_OUTLIERS),0)
ERROR += "Invalid REMOVE_OUTLIERS option selected : "$(REMOVE_OUTLIERS)"\n"
endif


# Use particle ids
# ----------------------------------------------------------------------------
ifeq ($(PARTICLE_ID),1)
CFLAGS += -DPARTICLE_ID
else ifneq ($(PARTICLE_ID),0)
ERROR += "Invalid PARTICLE_ID option selected : "$(PARTICLE_ID)"\n"
endif

# Sorting algorithm
# ----------------------------------------------------------------------------
ifeq ($(SORT),INSERTION)
CFLAGS += -DINSERTION_SORT
else ifeq ($(SORT),HEAP)
CFLAGS += -DHEAPSORT
else
ERROR += "Invalid SORT option selected : "$(SORT)"\n"
endif


# Multiple particle timestep options
# ----------------------------------------------------------------------------
ifeq ($(TIMESTEP),ADAPTIVE)
CFLAGS += -DADAPTIVE_TIMESTEP_LEVELS
else ifeq ($(TIMESTEP),RESTRICTED)
CFLAGS += -DRESTRICTED_TIMESTEP_LEVELS
else ifeq ($(TIMESTEP),FIXED)
CFLAGS += -DFIXED_TIMESTEP_LEVELS
else ifneq ($(TIMESTEP),0)
ERROR += "Invalid TIMESTEP option selected : "$(TIMESTEP)"\n"
endif


# Signal-velocity timestep
# ----------------------------------------------------------------------------
ifeq ($(SIGNAL_VELOCITY_DT),1)
CFLAGS += -DSIGNAL_VELOCITY
else ifneq ($(SIGNAL_VELOCITY_DT),0)
ERROR += "Invalid SIGNAL_VELOCITY_DT option selected : "$(SIGNAL_VELOCITY_DT)"\n"
endif


# Timing code
# ----------------------------------------------------------------------------
ifeq ($(TIMING_CODE),1)
CFLAGS += -DTIMING
OBJ += timing.o write_timing_stats.o
else ifneq ($(TIMING_CODE),0)
ERROR += "Invalid TIMING_CODE option selected : "$(TIMING_CODE)"\n"
endif


# Dimensionless units
# ----------------------------------------------------------------------------
ifeq ($(DIMENSIONLESS),1)
CFLAGS += -DDIMENSIONLESS
else ifneq ($(DIMENSIONLESS),0)
ERROR += "Invalid value for DIMENSIONLESS\n"
endif


# Test flags and routines
# ----------------------------------------------------------------------------
ifeq ($(TEST),SPIEGEL)
CFLAGS += -DSPIEGEL_TEST -DSPIEGEL_DISPERSION 
else ifeq ($(TEST),FREEFALL)
CFLAGS += -DFREEFALL_TEST
else ifeq ($(TEST),BINARY)
CFLAGS += -DBINARY_TEST
OBJ += write_binary_data.o
else ifeq ($(TEST),PLUMMER)
CFLAGS += -DPLUMMER_TEST
else ifeq ($(TEST),ENTROPY)
CFLAGS += -DENTROPY_CORE_TEST
else ifneq ($(TEST),0)
ERROR += "Invalid TEST option selected : "$(TEST)"\n"
endif


# Debug flags
# ----------------------------------------------------------------------------
ifeq ($(DEBUG),1)
CFLAGS += -DDEBUG1
else ifeq ($(DEBUG),2)
CFLAGS += -DDEBUG1 -DDEBUG2
else ifeq ($(DEBUG),3)
CFLAGS += -DDEBUG1 -DDEBUG2 -DDEBUG3
else ifneq ($(DEBUG),0)
ERROR += "Invalid value for DEBUG\n"
endif

# SPECIFIC_OUTPUT flags
# ----------------------------------------------------------------------------
ifeq ($(SPH_SPECIFIC_OUTPUT),1)
CFLAGS += -DSPH_SPECIFIC_OUTPUT
ifeq ($(SPH_OUTPUT_DENS),1)
CFLAGS += -DSPH_OUTPUT_DENS
endif
ifeq ($(SPH_OUTPUT_TEMP),1)
CFLAGS += -DSPH_OUTPUT_TEMP
endif
USER_SUB += sph_specific_output.o
endif

# ANALYSING flags
# ----------------------------------------------------------------------------
ifneq ($(ANALYSE),0)
CFLAGS += -DANALYSE
endif
ifeq ($(ANALYSE),1)
CFLAGS += -DANALYSE_DISC							# creates
USER_SUB += analyse_disc.o							 
endif
ifeq ($(ANALYSE),2)
CFLAGS += -DANALYSE_CENTRAL_REGION 
CFLAGS += -DANALYSE_CENTRAL_REGION_PROPERTIES 					# creates TD.dat file
USER_SUB += analyse_central_region.o
endif
ifeq ($(ANALYSE),3)
CFLAGS += -DANALYSE_CENTRAL_REGION 	
CFLAGS += -DANALYSE_CENTRE_AROUND_DENSEST 					# creates .cnt file
USER_SUB += analyse_central_region.o
endif
ifeq ($(ANALYSE),4)
CFLAGS += -DANALYSE_CENTRAL_REGION
CFLAGS += -DANALYSE_CENTRAL_REGION_PROPERTIES -DANALYSE_CENTRE_AROUND_DENSEST	# creates TD.dat and cnt. file
USER_SUB += analyse_central_region.o
endif
ifeq ($(ANALYSE),5)
CFLAGS += -DSINK_PROPERTIES_SYNC
USER_SUB +=read_sink_data_sync.o
endif
ifeq ($(ANALYSE),6)
CFLAGS += -DANALYSE_DISC -DPLANET_IN_DISC							# creates
USER_SUB += analyse_disc.o							 
endif


# List of all object files
# ----------------------------------------------------------------------------
OBJ += $(MODULE_OBJ) $(USER_MOD) $(USER_SUB) $(GENERIC_OBJ) $(SPH_OBJ) $(NBODY_OBJ) $(IO_OBJ) $(SETUP_OBJ)


# Code compilation
# ----------------------------------------------------------------------------
%.o: %.F90 definitions.o HP_types.o modules.o interface.o
	$(F90) $(OPT) $(CFLAGS) -c $<

seren :: $(OBJ) seren.o
ifdef ERROR
	@echo -e 'COMPILATION FAILED : Error(s) in Makefile\n' $(ERROR)
	@exit 2
endif
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/seren $(OBJ) seren.o

convert_format :: $(OBJ) convert_format.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/convert_format $(OBJ) convert_format.o

error_norm :: $(OBJ) error_norm.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/error_norm $(OBJ) error_norm.o

ic_BB :: $(OBJ) $(IC_OBJ) ic_BB.o 
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_BB $(OBJ) $(IC_OBJ) ic_BB.o

ic_binary ::  $(OBJ) ic_binary.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_binary $(OBJ) ic_binary.o

ic_binform :: $(OBJ) ic_binform.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_binform $(OBJ) ic_binform.o

ic_blob :: $(OBJ) $(IC_OBJ) ic_blob.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_blob $(OBJ) $(IC_OBJ) ic_blob.o

ic_core :: $(OBJ) $(IC_OBJ) ic_core.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_core $(OBJ) $(IC_OBJ) ic_core.o

ic_entropy_test :: $(OBJ) $(IC_OBJ) ic_entropy_test.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_entropy_test $(OBJ) $(IC_OBJ) ic_entropy_test.o

ic_hcp :: $(OBJ) ic_hcp.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_hcp $(OBJ) ic_hcp.o

ic_KH :: $(OBJ) $(IC_OBJ) ic_KH.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_KH $(OBJ) $(IC_OBJ) ic_KH.o

ic_jeans :: $(OBJ) ic_jeans.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_jeans $(OBJ) ic_jeans.o

ic_lattice_cube :: $(OBJ) ic_lattice_cube.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_lattice_cube $(OBJ) ic_lattice_cube.o

ic_mkdynfric :: $(OBJ) ic_mkdynfric.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_mkdynfric $(OBJ) ic_mkdynfric.o

ic_NTSI :: $(OBJ) ic_NTSI.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_NTSI $(OBJ) ic_NTSI.o

ic_polytrope :: $(OBJ) $(IC_OBJ) ic_polytrope.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_polytrope $(OBJ) $(IC_OBJ) ic_polytrope.o

ic_plummer :: $(OBJ) $(IC_OBJ) ic_plummer.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_plummer $(OBJ) $(IC_OBJ) ic_plummer.o

ic_rad_core :: $(OBJ) ic_radial_core.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_rad_core $(OBJ) ic_radial_core.o

ic_radtest :: $(OBJ) ic_radtest.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_radtest $(OBJ) ic_radtest.o

ic_random_cube :: $(OBJ) ic_random_cube.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_random_cube $(OBJ) ic_random_cube.o

ic_replicate_cubes :: $(OBJ) ic_replicate_cubes.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_replicate_cubes $(OBJ) ic_replicate_cubes.o

ic_RT :: $(OBJ) $(IC_OBJ) ic_RT.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_RT $(OBJ) $(IC_OBJ) ic_RT.o

ic_sedov :: $(OBJ) ic_sedov.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_sedov $(OBJ) ic_sedov.o

ic_shear_flow :: $(OBJ) $(IC_OBJ) ic_shear_flow.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_shear_flow $(OBJ) $(IC_OBJ) ic_shear_flow.o

ic_shocktube :: $(OBJ) smoothed_velocity.o ic_shocktube.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_shocktube $(OBJ) smoothed_velocity.o ic_shocktube.o

ic_SIS :: $(OBJ) $(IC_OBJ) ic_SIS.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_SIS $(OBJ) $(IC_OBJ) ic_SIS.o

ic_sphere :: $(OBJ) ic_sphere.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_sphere $(OBJ) ic_sphere.o

ic_spitzer :: $(OBJ) ic_spitzer.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_spitzer $(OBJ) ic_spitzer.o

ic_vel_pert :: $(OBJ) $(IC_OBJ) ic_vel_pert.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_vel_pert $(OBJ) $(IC_OBJ) ic_vel_pert.o

gravtest :: $(OBJ) direct_sph_gravity.o gravtest.o
	$(F90) $(OPT) $(CFLAGS) -o gravtest $(OBJ) direct_sph_gravity.o gravtest.o

lane_emden :: definitions.o modules.o lane_emden.o
	$(F90) $(OPT) $(CFLAGS) -o lane_emden modules.o lane_emden.o 

nbody_orbits :: definitions.o healpix_types.o modules.o nbody_orbits.o
	$(F90) $(OPT) $(CFLAGS) $(X11LIBS) $(PGPLOTLIBS) -o nbody_orbits definitions.o healpix_types.o modules.o nbody_orbits.o 

radial_average :: $(OBJ) radial_average.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/radial_average $(OBJ) radial_average.o

ic_subdisk :: $(OBJ) ic_subdisk.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_subdisk $(OBJ) ic_subdisk.o

analysedisc :: $(OBJ) seren.o
	$(F90) $(OPT) $(CFLAGS)  -o $(EXEDIR)/analysedisc $(OBJ) seren.o

analysetd :: $(OBJ) seren.o
	$(F90) $(OPT) $(CFLAGS)  -o $(EXEDIR)/analysetd $(OBJ) seren.o

analysecnt ::  $(OBJ) seren.o
	$(F90) $(OPT) $(CFLAGS)  -o $(EXEDIR)/analysecnt $(OBJ)  seren.o

analysetdcnt :: $(OBJ) seren.o
	$(F90) $(OPT) $(CFLAGS)  -o $(EXEDIR)/analysetdcnt $(OBJ)  seren.o

syncsinks :: $(OBJ) seren.o
	$(F90) $(OPT) $(CFLAGS)  -o $(EXEDIR)/syncsinks $(OBJ)  seren.o

# User makefile additions
# ----------------------------------------------------------------------------
#include user/user_makefiletail.mk


clean :: 
	\rm -f *.mod
	\rm -f *.o


# Dependencies
# ----------------------------------------------------------------------------
definitions.o : definitions.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

HP_types.o : HP_types.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

modules.o : modules.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

interface.o : interface.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

interface.o : definitions.o HP_types.o modules.o
modules.o : definitions.o HP_types.o
HP_types.o : definitions.o
