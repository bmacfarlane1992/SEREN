! SPH_HYDRO_FORCES.F90
! C. P. Batty & D. A. Hubber - 11/1/2007
! Calculates hydrodynamical accelerations for all particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_hydro_forces
  use interface_module, only : hydro,hydro_gradh
  use type_module
  use particle_module
  use time_module
  use timing_module
#if defined(OPENMP)
  use omp_lib
#endif
#if defined(CELL_WALK)
  use tree_module
#endif
  implicit none

  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
#if defined(OPENMP)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif
#if defined(DEBUG_VISC_BALSARA) || defined(DEBUG_VISC_PATTERN_REC)
  real(kind=DP) :: mean               ! Mean Balsara/pattern recog. factor
#endif

  debug2("Calculating hydro forces for all particles [sph_hydro_forces.F90]")

  allocate(acclist(1:ptot))


! Zero acceleration array of all active particles here (for now)
  do p=1,ptot
     if (accdo(p)) a(1:VDIM,p) = 0.0_PR
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     if (accdo(p)) a_grav(1:NDIM,p) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     if (accdo(p)) a_hydro(1:NDIM,p) = 0.0_PR
#endif
  end do


! Loop over all hydro particles (excluding boundary particles) and 
! calculate hydro forces
! ============================================================================
#if defined(HYDRO)

! For multiple particle timesteps, first make a list of all hydro SPH 
! particles on an acceleration step, and then parallelize over that list.
  acctot = 0
  do p=phydrostart,phydroend
     if (accdo(p)) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if
  end do
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
#if defined(TIMING)
  nhydrocomp = nhydrocomp + int(acctot,ILP)
#endif

  ! --------------------------------------------------------------------------
  if (acctot > 0) then
     debug_timing("HYDRO_FORCES")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     !!$OMP PARALLEL DO SCHEDULE(GUIDED,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
#if defined(IDEAL_MHD) && defined(GRAD_H_SPH)
        call ideal_mhd_gradh(p)
#elif defined(GRAD_H_SPH)
        call hydro_gradh(p)
#else
        call hydro(p)
#endif
     end do
     !$OMP END PARALLEL DO
  end if
  ! --------------------------------------------------------------------------

#endif
! ============================================================================


! Calculate and write out mean pattern recognition factor
#if defined(VISC_PATTERN_REC) && defined(DEBUG_VISC_PATTERN_REC)
  mean = 0.0_PR
  do p=1,ptot
     mean = mean + pattrec(p)
  end do
  open(1, file="vg.dat", position="append")
  write(1,*) nsteps,mean/real(ptot,PR)
  close(1)
#endif

! Calculate and write out mean Balsara factor
#if defined(VISC_BALSARA) && defined(DEBUG_VISC_BALSARA)
  mean = 0.0_PR
  do p=1,ptot
     mean = balsara(p) + mean
  end do
  open(1, file="bal.dat", position="append")
  write(1,*) nsteps,mean/real(ptot,PR)
  close(1)
#endif

  deallocate(acclist)

  return
END SUBROUTINE sph_hydro_forces
