! CLEAN_UP.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Frees up memory after simulation finishes (i.e. deallocates main arrays). 
! All variables should be deallocated in reverse order to which they are 
! originally allocated in allocate_memory.F90 (or in other routines).   
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE clean_up
  use particle_module
  use neighbour_module
  use hydro_module
  use kernel_module
  use tree_module
  use Nbody_module
  use time_module
  use sink_module
#if defined(RAD_WS)
  use Eos_module
#endif
#if defined(IDEAL_MHD)
  use mhd_module
#endif
#if defined(IONIZING_UV_RADIATION)
  use HP_module
#endif
  implicit none

  debug2("Deallocating all global arrays [clean_up.F90]")
  debug_timing("CLEAN_UP")


! Deallocate all sink arrays only if we are not using the N-body integrator
! (If using N-body integrator, these arrays are deallocated instead in 
! nbody_clean_up.F90)
! ----------------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  if (allocated(star)) deallocate(star)
#if defined(BINARY_STATS)
  if (allocated(binary)) deallocate(binary)
#endif
#endif

! Sink arrays
! ----------------------------------------------------------------------------
#if defined(SINKS) || defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  if (allocated(sink)) deallocate(sink)
#endif

! Kernel arrays
! ----------------------------------------------------------------------------
#if defined(GRAD_H_SPH)
  if (allocated(wg)) deallocate(wg)
  if (allocated(wh)) deallocate(wh)
#endif
  if (allocated(w6)) deallocate(w6)
  if (allocated(w5)) deallocate(w5)
  if (allocated(w4)) deallocate(w4)
  if (allocated(w3)) deallocate(w3)
  if (allocated(w2)) deallocate(w2)
  if (allocated(w1)) deallocate(w1)
  if (allocated(w0)) deallocate(w0)


! Deallocate all particle arrays if allocated (just check main array here)
! ============================================================================
  if (allocated(parray)) then

     ! Ionization routines
     ! -----------------------------------------------------------------------
#if defined(IONIZING_UV_RADIATION)
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     if (allocated(whichHPlevel)) deallocate(whichHPlevel)
#endif
     if (allocated(temp_aux)) deallocate(temp_aux)
     if (allocated(temp_min)) deallocate(temp_min)
     if (allocated(ionizedo)) deallocate(ionizedo)
     if (allocated(newtemp)) deallocate(newtemp)
     if (allocated(HPray)) deallocate(HPray)
#endif

     ! Radiative transport approximation arrays
     ! -----------------------------------------------------------------------
#if defined(RAD_WS)
#if defined(DEBUG_RAD)
     if (allocated(rad_info)) deallocate(rad_info)
#endif
     if (allocated(ueq)) deallocate(ueq)
     if (allocated(dt_therm)) deallocate(dt_therm)
     if (allocated(column2)) deallocate(column2)
     if (allocated(itemp)) deallocate(itemp)
     if (allocated(idens)) deallocate(idens)
#endif

     ! Tree arrays
     ! -----------------------------------------------------------------------
#if defined(BH_TREE)
     if (allocated(BHstock)) deallocate(BHstock)
     if (allocated(BHtemp)) deallocate(BHtemp)
     if (allocated(whichchild)) deallocate(whichchild)
     if (allocated(BHnextptcl)) deallocate(BHnextptcl)
     if (allocated(cellof)) deallocate(cellof)

     if (allocated(BHhydro)) deallocate(BHhydro)
     if (allocated(last_cell_hydro)) deallocate(last_cell_hydro)
     if (allocated(first_cell_hydro)) deallocate(first_cell_hydro)

#if defined(GRAVITY)
     if (allocated(BHgrav)) deallocate(BHgrav)
     if (allocated(last_cell_grav)) deallocate(last_cell_grav)
     if (allocated(first_cell_grav)) deallocate(first_cell_grav)
#endif
#endif

     ! Timestep arrays
     ! -----------------------------------------------------------------------
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
     if (allocated(nminneib)) deallocate(nminneib)
#endif
     if (allocated(laststep)) deallocate(laststep)
     if (allocated(nlevel)) deallocate(nlevel)
     if (allocated(nlast)) deallocate(nlast)
     if (allocated(accdo)) deallocate(accdo)

     ! MHD arrays
     ! -----------------------------------------------------------------------
#if defined(IDEAL_MHD)
#if defined(DEBUG_MHD)
     if (allocated(div_B)) deallocate(div_B)
#endif
     if (allocated(dB_dt)) deallocate(dB_dt)
     if (allocated(B_old)) deallocate(B_old)
     if (allocated(B)) deallocate(B)
#endif

     ! SPH hydro arrays
     ! -----------------------------------------------------------------------
#if defined(SINKS) && defined(GRAVITY)
     if (allocated(ispotmin)) deallocate(ispotmin)
#endif
#if defined(DIV_A)
     if (allocated(div_a)) deallocate(div_a)
#endif
#if defined(HEALPIX)
     if (allocated(gradrho)) deallocate(gradrho)
#endif
#if defined(VISC_PATTERN_REC)
     if (allocated(pattrec)) deallocate(pattrec)
#endif
#if defined(VISC_BALSARA)
     if (allocated(balsara)) deallocate(balsara)
#endif
#if defined(VISC_TD)
     if (allocated(dalpha_dt)) deallocate(dalpha_dt)
     if (allocated(talpha_old)) deallocate(talpha_old)
     if (allocated(talpha)) deallocate(talpha)
#endif
#if defined(SIGNAL_VELOCITY)
     if (allocated(vsigmax)) deallocate(vsigmax)
#endif
#if defined(GRAD_H_SPH)
     if (allocated(omega)) deallocate(omega)
#endif
     if (allocated(drhodt)) deallocate(drhodt)
     if (allocated(rho_old)) deallocate(rho_old)
     if (allocated(div_v)) deallocate(div_v)
     if (allocated(sound)) deallocate(sound)
     if (allocated(press)) deallocate(press)
     if (allocated(temp)) deallocate(temp)
     if (allocated(rho)) deallocate(rho)     
#if defined(NEIGHBOUR_LISTS)
     if (allocated(pplist)) deallocate(pplist)
     if (allocated(pptot)) deallocate(pptot)
#endif

     ! Main particle data arrays
     ! -----------------------------------------------------------------------
#if defined(ENTROPIC_FUNCTION)
     if (allocated(dA_dt)) deallocate(dA_dt)
     if (allocated(Aold)) deallocate(Aold)
     if (allocated(Aent)) deallocate(Aent)
#endif
#if defined(INTERNAL_ENERGY)
#if defined(DEBUG_DUDTRAD)
     if (allocated(dudt_rad)) deallocate(dudt_rad)
#endif
#if defined(DIFFUSION)
     if (allocated(du_dt_diff)) deallocate(du_dt_diff)
     if (allocated(k_cond)) deallocate(k_cond)
     if (allocated(lambda_diff)) deallocate(lambda_diff)
#endif
     if (allocated(du_dt)) deallocate(du_dt)
     if (allocated(u_old)) deallocate(u_old)
     if (allocated(u)) deallocate(u)
#endif
#if defined(DEBUG_FORCES)
#if defined(ARTIFICIAL_VISCOSITY)
     if (allocated(a_visc)) deallocate(a_visc)
#endif
#if defined(IDEAL_MHD)
     if (allocated(a_mag)) deallocate(a_mag)
#endif
#if defined(HYDRO)
     if (allocated(a_hydro)) deallocate(a_hydro)
#endif 
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
     if (allocated(a_grav)) deallocate(a_grav)
#endif
#endif
#if defined(GRAVITY) && defined(RAD_WS)
     if (allocated(sphgpot)) deallocate(sphgpot)
#endif
#if defined(GRAVITY) && !defined(GEOMETRIC_MAC)
     if (allocated(agravmag)) deallocate(agravmag)
#endif
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
     if (allocated(gpot)) deallocate(gpot)
#endif
#if defined(SMOOTHED_VELOCITY)
     if (allocated(v_smooth)) deallocate(v_smooth)
#endif
#if defined(PREDICTOR_CORRECTOR)
     if (allocated(a_old)) deallocate(a_old)
#endif
#if defined(RUNGE_KUTTA) || defined(LEAPFROG_KDK)
     if (allocated(v_half)) deallocate(v_half)
#endif
     if (allocated(v_old)) deallocate(v_old)
     if (allocated(r_old)) deallocate(r_old)
     if (allocated(a)) deallocate(a)
     if (allocated(v)) deallocate(v)
     if (allocated(parray)) deallocate(parray)
     if (allocated(porig)) deallocate(porig)

  end if
! ============================================================================


  return
END SUBROUTINE clean_up
