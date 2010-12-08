! SPH_GRAV_FORCES.F90
! C. P. Batty & D. A. Hubber - 11/1/2007
! Calculates gravitational accelerations for all particles 
! (including sink particles) - controls calls to relevant subroutines.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_grav_forces
  use interface_module, only : add_external_gravitational_force,&
       &BHgrav_accel,direct_sph_gravity,direct_sink_gravity
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

#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
  real(kind=DP) :: agravp(1:NDIM)     ! Grav. acceleration of particle p
  real(kind=DP) :: potp               ! Grav. potential of particle p
#if defined(OPENMP)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif
#if defined(CELL_WALK)
  logical :: activecell               ! Does cell contain active particles?
  integer :: c                        ! Cell counter
  integer :: pcounter                 ! Aux. particle counter
  integer :: pstart                   ! First particle id
#endif

  debug2("Calculating grav. forces for all particles [sph_grav_forces.F90]")

  allocate(acclist(1:ptot))

! Calculate gravitational forces of groups using tree
! ============================================================================
#if defined(GRAVITY) && defined(SPH_SELF_GRAVITY) && defined(BH_TREE) && defined(CELL_WALK)
  debug_timing("GROUPED_GRAVITY")

  acctot = 0
  pstart = pgravitystart
  pcounter = pstart - 1
  c = 0

  ! Walk gravity tree and find leaf cells containing active particles
  do
     if (BHgrav(c)%leaf > 0) then
        activecell = .false.
        do i=1,BHgrav(c)%leaf
           p = BHgrav(c)%plist(i)
           if (accdo(p)) activecell = .true.
        end do
        if (activecell) then
           acctot = acctot + 1
           acclist(acctot) = c
        end if
        c = BHgrav(c)%nextcell
     else if (BHgrav(c)%leaf == 0) then
        c = BHgrav(c)%ifopen
     else
        c = BHgrav(c)%nextcell
     end if
     if (c > ctot_grav) exit
  end do

  if (acctot > 0) then
     !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) PRIVATE(c)
     do i=1,acctot
        c = acclist(i)
        call BHgrav_grouped_walk(c)
     end do
     !$OMP END PARALLEL DO
  end if


! Loop over all self-gravitating particles and calculate forces
! ============================================================================
#elif defined(GRAVITY) && defined(SPH_SELF_GRAVITY)

! For multiple particle timesteps, first make a list of all self-gravitating
! SPH particles on an acceleration step, and then parallelize over that list.
  acctot = 0
  do p=pgravitystart,pgravityend
     if (accdo(p)) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if
  end do
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
#if defined(TIMING)
  ngravcomp = ngravcomp + int(acctot,ILP)
#endif


#if defined(SPH_SELF_GRAVITY)
  ! --------------------------------------------------------------------------
  if (acctot > 0) then
     debug_timing("GRAVITY_FORCES")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
     !$OMP PRIVATE(agravp,p,potp)
     do i=1,acctot
        p = acclist(i)
#if defined(BH_TREE)
        call BHgrav_accel(p,1.0_PR/parray(SMOO,p),&
             &parray(1:NDIM,p),agravp(1:NDIM),potp)
#elif defined(BINARY_TREE)
        call binary_gravacc(p,1.0_PR/parray(SMOO,p),&
             &parray(1:NDIM,p),agravp(1:NDIM),potp)
#else    
        call direct_sph_gravity(p,1.0_PR/parray(SMOO,p),&
             &parray(1:NDIM,p),agravp(1:NDIM),potp)
#endif
        a(1:NDIM,p) = a(1:NDIM,p) + real(agravp(1:NDIM),PR)
        gpot(p) = real(potp,PR)
#if defined(DEBUG_FORCES)
        a_grav(1:NDIM,p) = a_grav(1:NDIM,p) + real(agravp(1:NDIM),PR)
#endif
#if defined(GRAVITY) && defined(BH_TREE) & !defined(GEOMETRIC_MAC)
        agravmag(p) = real(sqrt(dot_product(agravp,agravp)),PR)
#endif
     end do
     !$OMP END PARALLEL DO
  end if
  ! --------------------------------------------------------------------------
#endif


#if defined(SINKS)
  ! --------------------------------------------------------------------------
  if (acctot > 0) then
     debug_timing("GRAVITY_FORCES")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
     !$OMP PRIVATE(agravp,p,potp)
     do i=1,acctot
        p = acclist(i)
        call direct_sink_gravity(p,1.0_PR/parray(SMOO,p),&
             &parray(1:NDIM,p),agravp(1:NDIM),potp)
        a(1:NDIM,p) = a(1:NDIM,p) + real(agravp(1:NDIM),PR)
        gpot(p) = gpot(p) + real(potp,PR)
#if defined(DEBUG_FORCES)
        a_grav(1:NDIM,p) = real(agravp(1:NDIM),PR)
#endif
#if defined(GRAVITY) && defined(BH_TREE) & !defined(GEOMETRIC_MAC)
        agravmag(p) = real(sqrt(dot_product(agravp,agravp)),PR)
#endif
     end do
     !$OMP END PARALLEL DO
  end if
  ! --------------------------------------------------------------------------
#endif

#endif
! ============================================================================


! Add external gravitational force to all active gas particles
! ============================================================================
#if defined(EXTERNAL_FORCE)
  do p=pgravitystart,pgravityend
     if (.not. accdo(p)) cycle
     call add_external_gravitational_force(parray(1:NDIM,p),&
          &agravp(1:NDIM),potp)
#if defined(GRAVITY) || defined(HYDRO)
     a(1:NDIM,p) = a(1:NDIM,p) + real(agravp(1:NDIM),PR)
#else
     a(1:NDIM,p) = real(agravp(1:NDIM),PR)
#endif
#if defined(DEBUG_FORCES)
     a_grav(1:NDIM,p) = real(agravp(1:NDIM),PR)
#endif
  end do
#endif
! ============================================================================


  deallocate(acclist)

! Special routine for ionization simulations
#if defined(HII_NO_INFALL)
  call hii_no_infall
#endif
#endif

  return
END SUBROUTINE sph_grav_forces
