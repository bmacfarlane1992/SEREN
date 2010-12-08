! ALLOCATE_MEMORY.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Allocates main data arrays.  These arrays are deallocated in reverse order 
! when the simulation terminates (in clean_up.F90).  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE allocate_memory
  use particle_module
  use neighbour_module
  use hydro_module
  use kernel_module
  use type_module
  use time_module
  use sink_module
#if defined(BH_TREE)
  use tree_module  
#endif
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  use Nbody_module
#endif
#if defined(RAD_WS)
  use Eos_module
#endif
#if defined(IDEAL_MHD)
  use mhd_module
#endif
#if defined(IONIZING_UV_RADIATION)
  use HP_module
#endif

  debug1("Allocating memory for main simulation arrays [allocate_memory.F90]")


! Only allocate SPH particles when there are some SPH particles
! (in case simulation is N-body only)
! ============================================================================
  if (ptot > 0) then

     ! Set arbitrary maximum no. of particles here
     pmax = min(int(real(PMAXMULT,DP)*real(ptot,DP)),ptot)

     ! Main particle data
     ! -----------------------------------------------------------------------
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating main particle data arrays..."
#endif
     allocate(porig(1:pmax))
     allocate(parray(1:DATATOT,1:pmax))
     allocate(v(1:VDIM,1:pmax))
     allocate(a(1:VDIM,1:pmax))
     allocate(r_old(1:NDIM,1:pmax))
     allocate(v_old(1:VDIM,1:pmax))
#if defined(RUNGE_KUTTA) || defined(LEAPFROG_KDK)
     allocate(v_half(1:VDIM,1:pmax))
#endif
#if defined(PREDICTOR_CORRECTOR)
     allocate(a_old(1:VDIM,1:pmax))
#endif
#if defined(SMOOTHED_VELOCITY)
     allocate(v_smooth(1:VDIM,1:pmax))
#endif
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
     allocate(gpot(1:pmax))
#endif
#if defined(GRAVITY) && !defined(GEOMETRIC_MAC)
     allocate(agravmag(1:pmax)) 
#endif
#if defined(GRAVITY) && defined(RAD_WS)
     allocate(sphgpot(1:pmax))
#endif
#if defined(DEBUG_FORCES)
#if defined(HYDRO)
     allocate(a_hydro(1:VDIM,1:pmax))
#endif
#if defined(EXTERNAL_FORCE) || defined(GRAVITY)
     allocate(a_grav(1:VDIM,1:pmax))
#endif
#if defined(IDEAL_MHD)
     allocate(a_mag(1:VDIM,1:pmax))
#endif
#if defined(ARTIFICIAL_VISCOSITY)
     allocate(a_visc(1:VDIM,1:pmax))
#endif
#endif

     ! Internal energy arrays
     ! -----------------------------------------------------------------------
#if defined(INTERNAL_ENERGY)
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating internal energy arrays..."
#endif
     allocate(u(1:pmax))
     allocate(u_old(1:pmax))
     allocate(du_dt(1:pmax))
#if defined(DIFFUSION)
     allocate(du_dt_diff(1:pmax))
     allocate(k_cond(1:pmax))
#endif
#if defined(DEBUG_DUDTRAD)
     allocate(dudt_rad(1:pmax))
#endif
#endif

     ! Entropy arrays
     ! -----------------------------------------------------------------------
#if defined(ENTROPIC_FUNCTION)
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating neighbour arrays..."
#endif
     allocate(Aent(1:pmax))
     allocate(Aold(1:pmax))
     allocate(dA_dt(1:pmax))
#endif

     ! Neighbour arrays
     ! -----------------------------------------------------------------------
#if defined(NEIGHBOUR_LISTS)
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating neighbour arrays..."
#endif
     allocate(pptot(1:pmax))
     allocate(pplist(1:pp_limit, 1:pmax))
#endif

     ! Hydrodynamical arrays
     ! -----------------------------------------------------------------------
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating hydro arrays..."
#endif
     allocate(rho(1:pmax))
     allocate(temp(1:pmax))
     allocate(press(1:pmax))
     allocate(sound(1:pmax))
     allocate(div_v(1:pmax))
     allocate(rho_old(1:pmax))
     allocate(drhodt(1:pmax))
#if defined(GRAD_H_SPH)
     allocate(omega(1:pmax))
#endif
#if defined(SIGNAL_VELOCITY)
     allocate(vsigmax(1:pmax))
#endif
#if defined(VISC_TD)
     allocate(talpha(1:pmax))
     allocate(talpha_old(1:pmax))
     allocate(dalpha_dt(1:pmax))
#endif
#if defined(VISC_BALSARA)
     allocate(balsara(1:pmax))
#endif
#if defined(VISC_PATTERN_REC)
     allocate(pattrec(1:pmax))
#endif
#if defined(HEALPIX)
     allocate(gradrho(1:NDIM,1:pmax))
#endif
#if defined(DIV_A)
     allocate(div_a(1:pmax))
#endif
#if defined(SINKS) && defined(GRAVITY)
     allocate(ispotmin(1:pmax))
#endif

     ! MHD arrays
     ! -----------------------------------------------------------------------
#if defined(IDEAL_MHD)
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating MHD arrays..."
#endif
     allocate(B(1:BDIM,ptot))
     allocate(B_old(1:BDIM,ptot))
     allocate(dB_dt(1:BDIM,ptot))
#if defined(DEBUG_MHD)
     allocate(div_B(1:pmax))
#endif
#endif

     ! Multiple particle timestep arrays
     ! -----------------------------------------------------------------------
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating timestepping arrays..."
#endif
     allocate(accdo(1:pmax))
     allocate(nlast(1:pmax))
     allocate(nlevel(1:pmax))
     allocate(laststep(1:pmax))
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
     allocate(nminneib(1:pmax))
#endif

     ! BH tree arrays
     ! -----------------------------------------------------------------------
#if defined(BH_TREE)
#if defined(DEBUG_ALLOCATE_MEMORY)
     write(6,*) "Allocating tree arrays..."
#endif
     if (LEAFMAX > 8) then
        cmax_grav  = (2*(pgas + pcdm)) / (LEAFMAX - 4)
        cmax_hydro = (2*pmax) / (LEAFMAX - 4)
     else
        cmax_grav  = 2*(pgas + pcdm)
        cmax_hydro = 2*pmax
     end if
     
#if defined(GRAVITY)
     allocate(first_cell_grav(0:LMAX))
     allocate(last_cell_grav(0:LMAX))
     allocate(BHgrav(0:cmax_grav))
#endif
     allocate(first_cell_hydro(0:LMAX))
     allocate(last_cell_hydro(0:LMAX))
     allocate(BHhydro(0:cmax_hydro))

     allocate(cellof(1:pmax))
     allocate(BHnextptcl(1:pmax))
     allocate(whichchild(1:pmax))
     allocate(BHtemp(0:cmax_hydro))
     allocate(BHstock(0:cmax_hydro))
#endif

     ! Radiative transport approximation arrays
     ! -----------------------------------------------------------------------
#if defined(RAD_WS)
     allocate(idens(1:pmax))
     allocate(itemp(1:pmax))
     allocate(column2(1:pmax))
     allocate(dt_therm(1:pmax))
     allocate(ueq(1:pmax))
#if defined(DEBUG_RAD)
     allocate(rad_info(1:10,ptot))
#endif
#endif

     ! HEALPix arrays
     ! -----------------------------------------------------------------------
#if defined(IONIZING_UV_RADIATION)
     imax = ptot
     write(6,*) "imax : ",imax
     allocate(HPray(1:imax))
     allocate(newtemp(1:pmax))
     allocate(ionizedo(1:pmax))
     allocate(temp_min(1:pmax))
     allocate(temp_aux(1:pmax))
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     allocate(whichHPlevel(1:pmax))
#endif
#endif

  end if
! ============================================================================


! Kernel arrays
! ----------------------------------------------------------------------------
#if defined(DEBUG_ALLOCATE_MEMORY)
  write(6,*) "Allocating kernel arrays..."
#endif
  allocate(w0(0:KERNTOT))
  allocate(w1(0:KERNTOT))
  allocate(w2(0:KERNTOT))
  allocate(w3(0:KERNTOT))
  allocate(w4(0:KERNTOT))
  allocate(w5(0:KERNTOT))
  allocate(w6(0:KERNTOT))
#if defined(GRAD_H_SPH)
  allocate(wh(0:KERNTOT))
  allocate(wg(0:KERNTOT))
#endif

! Sink particle arrays
! ----------------------------------------------------------------------------
#if defined(SINKS) || defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  allocate(sink(1:SMAX))
#endif

! Star particle arrays
! ----------------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  allocate(star(1:SMAX))
#if defined(BINARY_STATS)
  allocate(binary(1:SMAX))
#endif
#endif


  return
END SUBROUTINE allocate_memory
