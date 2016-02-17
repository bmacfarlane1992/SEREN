! SPH_SIM_TIMESTEPS.F90
! D. A. Hubber - 6/9/2010
! Computes ideal timesteps for all particles, then quantises to lower 
! block timestep level.  Allows the timestep of a particle, and the global 
! minimum, to increase by one level if synchronised correctly.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_timesteps
  use interface_module, only : distance2,sink_timestep,timestep_size
  use time_module
  use particle_module
  use hydro_module
  use scaling_module
  use sink_module
  use filename_module
  implicit none

  integer :: p                           ! Particle counter
  integer :: p_update                    ! No. of timesteps updated
  integer(kind=ILP) :: level             ! Timestep level
  integer(kind=ILP) :: level_old         ! Old timestep level
  integer(kind=ILP) :: nlastlevel        ! Timestep level
  real(kind=DP) :: dt_min                ! Minimum step size
  real(kind=DP) :: dt                    ! Step size of particle p
  real(kind=DP), allocatable :: step(:)  ! Step size storage
#if defined(SINKS)
  integer :: s                           ! Sink particle counter
#endif
#if defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(SINKS)
  real(kind=PR) :: dr(1:NDIM)            ! Relative displacement vector
  real(kind=PR) :: drsqd                 ! Distance squared
#endif
#if defined(DEBUG_BLOCK_TIMESTEPS)
  integer :: pmin                               ! Particle with dt_min
  integer(kind=ILP), allocatable :: ninlevel(:) ! Number of particles in level
#endif
!  character(len=256) :: out_file         ! Data snapshot filename

  debug_timing("SPH_TIMESTEPS")
  p_update = 0


! Check if we are at a resynchronisation step
! ============================================================================
  if (n == nresync) then

      debug2("Resynchronising block timesteps [sph_timesteps.F90]")

      ! Wrap n back around to 0 for resync step
      n = 0_ILP
      allocate(step(1:ptot))
      dt_min = BIG_NUMBER_DP
      
      ! Loop over all hydro and gas particles and find ideal time steps,
      ! and also the minimum of all particle timesteps
      ! ----------------------------------------------------------------------
      !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(MIN:dt_min) PRIVATE(dt)
      do p=1,ptot
         nlast(p) = 0_ILP
         call timestep_size(p,dt)
         step(p) = dt
         dt_min = min(dt,dt_min)
      end do
      !$OMP END PARALLEL DO
      
      ! If using sink particles, find their timestep including minimum
#if defined(SINKS)
      do s=1,stot
         call sink_timestep(s,dt)
         dt_min = min(dt,dt_min)
      end do
#endif

#if defined(DEBUG_BLOCK_TIMESTEPS)
      write(6,*) "dt_min : ",dt_min*tscale
#endif

      ! Compute min/max timestep values for resync 
      level_max = nlevels - 1_ILP
#if defined(RESTRICTED_TIMESTEP_LEVELS)
      level     = int(INVLOGETWO*log(dt_fixed / dt_min),ILP) + 1_ILP
      dt_min    = dt_fixed / 2.0_DP**(level)
#elif defined(FIXED_TIMESTEP_LEVELS)
      dt_min    = dt_fixed / 2.0_DP**(level_max)
#endif
      dt_max    = dt_min * 2.0_DP**(level_max)

#if defined(DEBUG_BLOCK_TIMESTEPS)
      write(6,*) "dt_min : ",dt_min*tscale,"    dt_max : ",dt_max*tscale
      write(6,*) "level_max : ",&
           &int(INVLOGETWO*log(dt_max / minval(step)),ILP) + 1_ILP
#endif
      
      ! Loop over hydro and gas particles and assign step levels
      ! ----------------------------------------------------------------------
      do p=1,ptot
         dt        = step(p)
         level     = min(int(INVLOGETWO*log(dt_max/dt),ILP) + 1_ILP, level_max)
         level     = max(level,0_ILP)
         nlevel(p) = level
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
         nminneib(p) = level
#endif
#if defined(DEBUG_BLOCK_TIMESTEPS)
         if (level > level_max) stop 'timestep level > level_max'
#endif
      end do

      
      ! If particles are sink neighbours, set to same timestep as sinks.
      ! ----------------------------------------------------------------------
#if defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(SINKS)
      do p=1,ptot
         do s=1,stot
            call distance2(sink(s)%r(1:NDIM),p,dr(1:NDIM),drsqd)
            if (drsqd <= REXTENTSQD*sink(s)%radius*sink(s)%radius) &
                 &nminneib(p) = level_max
         end do
         if (nminneib(p) - nlevel(p) > TIMESTEP_LEVEL_DIFF_MAX) then
            nlevel(p)   = level_max - TIMESTEP_LEVEL_DIFF_MAX
            nminneib(p) = level_max - TIMESTEP_LEVEL_DIFF_MAX
         end if
      end do
#endif

      ! 2 integer steps required per real time step for Runge-Kutta,
      ! leapfrog or predictor-corrector, else 1 integer timestep for Euler.
#ifndef EULER
      level_step = level_max + 1_ILP
#else
      level_step = level_max
#endif
      nresync = 2**(level_step)
      timestep = dt_max / real(nresync,DP)
      nstepsize = 1

      ! Set all sink variables for next step
#if defined(SINKS)
      nlast_sinks = n
      nlevel_sinks = level_max
#endif

      ! Free up memory for temp arrays
      deallocate(step)


! ============================================================================
! If not resync time, check if any timesteps need to be recomputed
! ============================================================================
  else

     debug2("Recalculating particle stepsizes [sph_timesteps.F90]")
     level_old = level_max

     ! Loop over all hydro and gas particles and find which particles 
     ! need their timesteps to be updated.
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dt,level,nlastlevel) &
     !$OMP REDUCTION(+ : p_update)
     do p=1,ptot
        
        ! Skip if timestep not finished
        if (n == nlast(p)) then
           p_update = p_update + 1
           nlastlevel = nlevel(p)
           
           ! Find ideal timestep for particle p
           call timestep_size(p,dt)
           level = max(int(INVLOGETWO*log(dt_max/dt),ILP) + 1_ILP, 0_ILP)
           
           ! Allow timestep to go up one level if synchronised correctly.
           ! Else, allow it to fall to any lower timetep.
           if (level < nlastlevel .and. nlastlevel > 1_ILP .and. &
                mod(n,2**(level_step - (nlastlevel - 1_ILP))) == 0_ILP) then
              nlevel(p) = nlastlevel - 1_ILP
           else if (level > nlastlevel) then
              nlevel(p) = level
           else
              nlevel(p) = nlastlevel
           end if
        end if
        
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------


     ! Now loop over all sinks and determine smallest timestep.
     ! ----------------------------------------------------------------------
#if defined(SINKS)
     if (n == nlast_sinks) then
        if (stot > 0) then
           dt_min = BIG_NUMBER_DP
           p_update = p_update + 1
           do s=1,stot
              call sink_timestep(s,dt)
              dt_min = min(dt_min,dt)
           end do
           nlevel_sinks = max(int(INVLOGETWO * log(dt_max/dt_min),ILP) &
                &+ 1_ILP, 0_ILP)
        end if
     end if
#endif

      ! If no particles or sinks are updated, do not update timesteps
      ! ----------------------------------------------------------------------
      if (p_update > 0) then
         
         ! Calculate maximum level
         level_max = 0
#if defined(SINKS)
         if (stot > 0) level_max = nlevel_sinks
#endif
         do p=1,ptot
            level_max = max(level_max,nlevel(p))
         end do

         ! If minimum stepsize is synchronized, then increase if possible
         if (level_max <= level_old - 1_ILP .and. level_old > 1_ILP .and. &
              & mod(n,2**(level_step - (level_old - 1_ILP))) == 0_ILP) then
#if defined(DEBUG_BLOCK_TIMESTEPS)
            write(6,*) "Timestep level removed : ",level_old,level_max
#endif
            level_max = level_old - 1_ILP
         else if (level_max < level_old) then
            level_max = level_old
         end if


         ! If neighbour's timesteps are more than TIMESTEP_DIFF_MAX times 
         ! times smaller, immediately reduce the timestep if possible. 
         ! Also, skip update if particle is neighbouring a sink.
         ! -------------------------------------------------------------------
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
         do p=1,ptot
#if defined(SINKS)
            do s=1,stot
               call distance2(sink(s)%r(1:NDIM),p,dr(1:NDIM),drsqd)
               if (drsqd <= REXTENTSQD*sink(s)%radius*sink(s)%radius) &
                    &nminneib(p) = level_max
            end do
#endif
            nminneib(p) = min(nminneib(p),level_max)
            if (nminneib(p) - nlevel(p) > TIMESTEP_LEVEL_DIFF_MAX .and. &
                 & nlast(p) == n .and. &
                 & mod(n,2**(level_step - nminneib(p) + &
                 & TIMESTEP_LEVEL_DIFF_MAX)) == 0_ILP) then
               nlevel(p) = min(nminneib(p) - TIMESTEP_LEVEL_DIFF_MAX,level_max)
               nminneib(p) = nlevel(p)
            else if (nlast(p) == n) then
               nminneib(p) = nlevel(p)
            end if
         end do
#endif

         ! Rescale all integer step variables if level has changed
         ! -------------------------------------------------------------------
         if (level_max > level_old) then
            nlast(1:ptot) = nlast(1:ptot) * (2**(level_max - level_old))
            n = n*(2_ILP**(level_max - level_old))
            nresync = nresync*(2_ILP**(level_max - level_old))
#if defined(SINKS)
            nlast_sinks = nlast_sinks*(2**(level_max - level_old))
#endif
         else if (level_max < level_old) then
            nlast(1:ptot) = nlast(1:ptot) / (2**(level_old - level_max))
            n = n / (2_ILP**(level_old - level_max))
            nresync = nresync / (2_ILP**(level_old - level_max))
#if defined(SINKS)
            nlast_sinks = nlast_sinks / (2**(level_old - level_max))
#endif
         end if


         ! 2 integer steps required per real time step for Runge-Kutta,
         ! leapfrog or predictor-corrector, else 1 integer timestep for Euler.
#ifndef EULER
         level_step = level_max + 1_ILP
#else
         level_step = level_max 
#endif
         timestep = dt_max / real(nresync,DP)

         ! Set all sink variables for next step
#if defined(SINKS)
         nlevel_sinks = level_max
#endif

      end if
      ! ----------------------------------------------------------------------

  end if
! ============================================================================


#if defined(DEBUG_BLOCK_TIMESTEPS)
  write(6,*) "n      : ",n,"    nresync : ",nresync
  write(6,*) "nsteps : ",nsteps,"   time : ",time*tscale
  if (n == 0 .or. p_update > 0) then
     allocate(ninlevel(0:level_max))
     ninlevel(0:level_max) = 0
     do p=1,ptot
        level = nlevel(p)
        ninlevel(level) = ninlevel(level) + 1_ILP
     end do
#if defined(SINKS)     
     if (stot > 0) ninlevel(nlevel_sinks) = min(-ninlevel(nlevel_sinks),-1_ILP)
     write(6,*) "nlevel_sinks : ",nlevel_sinks,level_max,nlast_sinks
#endif
     do level=0,level_max
        write(6,*) "No. of particles in level ",level," : ",ninlevel(level)
     end do
     if (ninlevel(level_max) == 0 .and. stot == 0) then
        write(6,*) "No particles in highest level : ",ninlevel(level_max)
        do p=1,ptot
           if (nlevel(p) == level_max) &
                &write(6,*) "But there's one here ... ",p,nlevel(p)
        end do
        stop
     end if
     deallocate(ninlevel)
  end if
#endif

!  out_file =trim(adjustl(run_dir))//trim(adjustl(run_id))//".minstep"
!  open(unit=1,file=out_file,position='append')
!  write(1,*) nsteps,time,timestep*tscale
!  close(1)


! Advance time and integer times
  n      = n + 1
  nsteps = nsteps + 1
  time   = time + timestep


  return
END SUBROUTINE sph_timesteps
