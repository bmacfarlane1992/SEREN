! SPH_UPDATE.F90
! C. P. Batty & D. A. Hubber - 18/1/2007
! Control subroutine for calculating all SPH properties, 
! i.e. smoothing lengths, neighbour lists and densities, of all particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_update
  use interface_module, only : all_sph,all_sph_gradh,get_neib,h_gather
  use particle_module
  use neighbour_module
  use type_module
  use timing_module
  use time_module, only : accdo
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none

  integer :: acctot                   ! No. of particles on acc. step
  integer :: acchydrotot              ! No. of active hydro particles
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
#if defined(OPENMP)
  integer :: chunksize                ! Loop chunksize for OMP
#endif

  debug2("Calculating SPH quantities [sph_update.F90]")


! Prepare list of particles which are on an acceleration step
! ----------------------------------------------------------------------------
  allocate(acclist(1:ptot))
  acctot = 0
  do p=1,ptot
     if (accdo(p)) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if
  end do
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
#if defined(TIMING)
  nsphcomp = nsphcomp + int(acctot,ILP)
#endif


! Call all subroutines only when there are active paricles
! ============================================================================
  if (acctot > 0) then
     
     ! 'grad-h' SPH method
     ! -----------------------------------------------------------------------
#if defined(GRAD_H_SPH)
     debug2("Calculating densities and smoothing lengths for all particles")
     debug_timing("ALL_SPH_GRADH")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call all_sph_gradh(p)
     end do
     !$OMP END PARALLEL DO
     
     ! Update hmax in hydro cells now smoothing lengths have been calculated
#if defined(BH_TREE)
     debug_timing("BHHYDRO_UPDATE")
     if (acctot > 0) call BHhydro_update_hmax
#endif
     
     ! Getting neighbour lists for all particles 
#if defined(NEIGHBOUR_LISTS) && defined(HYDRO)
     debug2("Calculating neighbour lists for all particles")
     debug_timing("GET_NEIB")

     ! Reconstruct active particle list for hydro particles only
     acchydrotot = 0
     do p=phydrostart,phydroend
        if (accdo(p)) then
           acchydrotot = acchydrotot + 1
           acclist(acchydrotot) = p
        end if
     end do
#if defined(OPENMP)
     chunksize = int(CHUNKFRAC*real(acchydrotot,PR)/&
          &real(omp_get_max_threads(),PR)) + 1
#endif

     if (acchydrotot > 0) then
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) &
        !$OMP DEFAULT(SHARED) PRIVATE(p)
        do i=1,acchydrotot
           p = acclist(i)
           call get_neib(p)
        end do
        !$OMP END PARALLEL DO
     end if
#endif
     
     
     ! Standard SPH
     ! -----------------------------------------------------------------------
#else
     ! Calculating new smoothing lengths for all particles
#if defined(CONSTANT_H)
     debug2("Setting constant smoothing lengths for all particles")
     parray(SMOO,1:ptot) = hmin

#elif defined(HGATHER) || defined(HMASS)
     debug2("Calculating smoothing lengths for all particles")
     debug_timing("H_GATHER")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call h_gather(p,parray(SMOO,p),parray(1:NDIM,p))
     end do
     !$OMP END PARALLEL DO
#endif
     
     ! Update hmax in hydro cells now smoothing lengths have been calculated
#if defined(BH_TREE)
     debug_timing("BHHYDRO UPDATE")
     call BHhydro_update_hmax
#endif

     ! Getting neighbour lists for all particles
#if defined(NEIGHBOUR_LISTS)
     debug2("Calculating neighbour lists for all particles")
     debug_timing("GET_NEIB")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call get_neib(p)
     end do
     !$OMP END PARALLEL DO
#endif

     ! Calculating SPH densities for all particles
     debug2("Calculating SPH quantities for all particles")
     debug_timing("ALL_SPH")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call all_sph(p)
     end do
     !$OMP END PARALLEL DO
#endif
     ! -----------------------------------------------------------------------

  end if
! ============================================================================

  if (allocated(acclist)) deallocate(acclist)

  return
END SUBROUTINE sph_update
