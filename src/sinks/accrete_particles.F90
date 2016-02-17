! ACCRETE_PARTICLES.F90
! D. A. Hubber - 17/6/2007
! Removes particles from computational domain that have been accreted by 
! sink particles.  First, we identify all particles that are within any sink 
! boundaries.  Next, we determine if any of these particles are bound to any 
! of the sinks.  If accreted, the mass, momentum and angular momentum of the 
! particle is added to the sink.  The particle is then removed from the 
! simulation.  Finally, the particles are reordered so all live particles 
! are contiguous in the main data arrays (i.e. no holes in memory with dead 
! particles). 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE accrete_particles
  use particle_module
  use sink_module
  use type_module
  use kernel_module
  use scaling_module
  use time_module, only : nsteps,nstock
  implicit none

  integer :: i                           ! Auxilary counter
  integer :: ndead                       ! Total number of dead particles
  integer :: p                           ! Particle counter
  integer :: pp_tot                      ! Number of particles accreted by s
  integer :: s                           ! Sink counter
  integer :: ss                          ! Aux sink i
  integer, allocatable :: deadlist(:)    ! List of dead (accreted) particles
  integer, allocatable :: pp_templist(:) ! Temp. list for particle ids
  integer, allocatable :: sinkflag(:)    ! Flag sink-particle overlap
  real(kind=DP) :: acom_temp(1:VDIM)     ! Accel. of particle-sink system
  real(kind=DP) :: angmomsink(1:3)       ! Angular momentum of sink s
  real(kind=DP) :: ap(1:VDIM)            ! Acceleration of particle p
  real(kind=DP) :: as(1:VDIM)            ! Acceleration of sink s
  real(kind=PR) :: atemp(1:NDIM)         ! Aux. acceleration vector
  real(kind=PR) :: dpotp                 ! Aux. grav potential variable
  real(kind=DP) :: dr(1:NDIM)            ! Relative position vector
  real(kind=DP) :: dv(1:VDIM)            ! Relative velocity vector
  real(kind=DP) :: drsqd                 ! Distance squared
  real(kind=DP) :: drsqd2                ! Second distance squared variable
  real(kind=DP) :: invhp                 ! 1.0 / hp
  real(kind=DP) :: invhs                 ! 1.0 / hsink
  real(kind=DP) :: mp                    ! Mass of particle p
  real(kind=DP) :: ms                    ! Mass of sink s
  real(kind=DP) :: mtot_temp             ! Mass of particle-sink system
  real(kind=DP) :: rads                  ! Radius of sink particle s
  real(kind=DP) :: rcom_temp(1:NDIM)     ! COM of particle-sink system
  real(kind=PR) :: reduced_mass          ! Reduced mass of 2-body system
  real(kind=DP) :: rp(1:NDIM)            ! Position of particle p
  real(kind=DP) :: rs(1:NDIM)            ! Position of sink s
  real(kind=DP) :: total_energy          ! Total energy of sink-particle system
  real(kind=DP) :: vcom_temp(1:VDIM)     ! COM velocity of particle-sink system
  real(kind=DP) :: vp(1:VDIM)            ! Velocity of particle p
  real(kind=DP) :: vs(1:VDIM)            ! Velocity of sink s
#if defined(LEAPFROG_KDK)
  real(kind=DP) :: vhalf_temp(1:VDIM)    ! Half-step velocity of p-s system
#endif
#if defined(ACCRETION_RATE)
  real(kind=DP) :: maccreted             ! Total mass of accreted particles
#endif
#if defined(DEBUG_ACCRETE)
  logical :: accreteflag                 ! Flag if any particles are accreted
#endif

  debug2("Accreting SPH particles to sinks [accrete_particles.F90]") 
  debug_timing("ACCRETE_PARTICLES")

! Initialize sink accretion arrays
  allocate(deadlist(1:ptot))
  allocate(sinkflag(1:ptot))
  allocate(pp_templist(1:ptot))
  do p=1,ptot
     sinkflag(p) = -1
  end do
  ndead = 0

! Find id of sink to which particle is closest and within the sink radius
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dr,drsqd,drsqd2,rads,rs,s,ss)
  do p=pgasstart,pgasend
     do s=1,stot
        if (.NOT. sink(s)%accrete) cycle
        rs(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
        rads = real(sink(s)%radius,DP)
        call distance2_dp(rs(1:NDIM),p,dr(1:NDIM),drsqd)

        if (drsqd < rads*rads) then
           if (sinkflag(p) == -1) then
              sinkflag(p) = s
           else
              ss = sinkflag(p)
              rs(1:NDIM) = real(sink(ss)%r(1:NDIM),DP)
              call distance2_dp(rs(1:NDIM),p,dr(1:NDIM),drsqd2)
              if (drsqd < drsqd2) sinkflag(p) = s
           end if
#if defined(DEBUG_ACCRETE)
           write(6,*) "Check accretion : ",p,sinkflag(p),&
                &sqrt(drsqd)*rscale,rads*rscale
#endif
        end if

     end do
  end do
!$OMP END PARALLEL DO

! Write some information to screen for debugging purposes
#if defined(DEBUG_ACCRETE)
  do p=1,ptot
     if (sinkflag(p) /= -1) ndead = ndead + 1
  end do
  if (ndead > 0) then
     write(6,*) "ptot (before accretion) : ",ptot
     write(6,*) "sink candidates : ",ndead
     accreteflag = .true.
     call diagnostics 
  else
     accreteflag = .false.
  end if
  ndead = 0
#endif


! If a particle is within the sink radius of and bound to a sink,
! accrete particle to the sink.
! ============================================================================
  do s=1,stot

#if defined(DEBUG_ACCRETE)
     write(6,*) "Accreting sink? : ",s,sink(s)%accrete
#endif

     ! Skip sink if it's not flagged as accreting
     if (.NOT. sink(s)%accrete) cycle

     ! Create local copies of old sink properties
     rs(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
     vs(1:VDIM) = real(sink(s)%v(1:VDIM),DP)
     as(1:VDIM) = real(sink(s)%a(1:VDIM),DP)
     ms    = real(sink(s)%m,DP)
     rads  = real(sink(s)%radius,DP)
     invhs = 1.0_DP / real(sink(s)%h,DP)

     ! Temporary array to store new sink properties
     rcom_temp(1:NDIM) = ms*rs(1:NDIM)
     vcom_temp(1:VDIM) = ms*vs(1:VDIM)
     acom_temp(1:VDIM) = ms*as(1:VDIM)
     mtot_temp = ms
#if defined(LEAPFROG_KDK)
     vhalf_temp(1:VDIM) = ms*real(sink(s)%vhalf(1:VDIM),DP)
#endif
#if defined(DEBUG_FORCES)
     sink(s)%agrav(1:NDIM) = real(ms,PR)*sink(s)%agrav(1:NDIM)
     sink(s)%ahydro(1:NDIM) = real(ms,PR)*sink(s)%ahydro(1:NDIM)
#endif
#if defined(ACCRETION_RATE)
     maccreted = 0.0_DP
#endif
     pp_tot = 0


     ! Search through all gravitating particles and find those bound to sink s
     ! -----------------------------------------------------------------------
     do p=pgasstart,pgasend

        ! If the particle is within radius of sink s (and only sink s), then
        ! check if particle p is bound to sink s
        if (sinkflag(p) /= s) cycle

        ! Make local copy of particle quantities
        invhp      = 1.0_DP / real(parray(SMOO,p),DP)
        mp         = real(parray(MASS,p),DP)
        rp(1:NDIM) = real(parray(1:NDIM,p),DP)
        ap(1:VDIM) = real(a(1:VDIM,p),DP)
        vp(1:VDIM) = real(v(1:VDIM,p),DP)
        dv(1:VDIM) = vs(1:VDIM) - vp(1:VDIM)
        reduced_mass = ms*mp/(ms + mp)

        ! Only perform energy check if option chosen in params file
        if (energy_accrete) then
#if defined(N_BODY)
           call gravity_nbody(real(ms,PR),real(rp(1:NDIM),PR),&
                &real(rs(1:NDIM),PR),atemp(1:NDIM),dpotp)
#elif defined(GRAD_H_SPH)
           call gravity_gradh(real(invhp,PR),real(invhs,PR),real(ms,PR),&
                &real(rp(1:NDIM),PR),real(rs(1:NDIM),PR),&
                &0.0_PR,0.0_PR,atemp(1:NDIM),dpotp)
#else
           call gravity_sph(real(invhp,PR),real(invhs,PR),real(ms,PR),&
                &real(rp(1:NDIM),PR),real(rs(1:NDIM),PR),atemp(1:NDIM),dpotp)
#endif
           total_energy = 0.5_DP*reduced_mass*&
                &dot_product(dv(1:NDIM),dv(1:NDIM)) - real(dpotp,DP)*mp
        else
           total_energy = -BIG_NUMBER_DP
        end if

#if defined(DEBUG_ACCRETE)
        write(6,*) "candidate : ",p,total_energy
#endif

        ! If bound, add particle properties to the sink and then add 
        ! particle id to list of 'dead' (i.e. accreted) particles.
        if (total_energy < 0.0_DP) then
#if !defined(KILLING_SINKS)
           mtot_temp = mtot_temp + mp
#endif           
           rcom_temp(1:NDIM) = rcom_temp(1:NDIM) + mp*rp(1:NDIM)
           vcom_temp(1:VDIM) = vcom_temp(1:VDIM) + mp*vp(1:VDIM)
           acom_temp(1:VDIM) = acom_temp(1:VDIM) + mp*ap(1:VDIM)
#if defined(LEAPFROG_KDK)
           vhalf_temp(1:VDIM) = vhalf_temp(1:VDIM) + &
                &mp*real(v_half(1:VDIM,p),DP)
#endif
#if defined(DEBUG_FORCES)
           sink(s)%agrav(1:NDIM) = sink(s)%agrav(1:NDIM) + &
                &real(mp,PR)*a_grav(1:NDIM,p)
           sink(s)%ahydro(1:NDIM) = sink(s)%ahydro(1:NDIM) + &
                &real(mp,PR)*a_hydro(1:NDIM,p)
#endif
           ndead  = ndead + 1
           pp_tot = pp_tot + 1
           deadlist(ndead) = p
#if defined(ACCRETION_RATE)
           maccreted = maccreted + mp
#endif
        end if
     end do

     ! Store new sink properties if any particles have been accreted
     ! -----------------------------------------------------------------------
     if (pp_tot > 0) then

        ! First calculate centres of mass, velocity and acceleration
        rcom_temp(1:NDIM) = rcom_temp(1:NDIM) / mtot_temp
        vcom_temp(1:VDIM) = vcom_temp(1:VDIM) / mtot_temp
        acom_temp(1:VDIM) = acom_temp(1:VDIM) / mtot_temp
#if defined(LEAPFROG_KDK)
        vhalf_temp(1:VDIM) = vhalf_temp(1:VDIM) / mtot_temp
#endif

        ! Angular momentum of old COM around new COM
        angmomsink(1:3) = 0.0_DP
        call distance3_dp(rcom_temp(1:NDIM),rs(1:NDIM),dr(1:NDIM),drsqd)
        dv(1:VDIM) = vs(1:VDIM) - vcom_temp(1:VDIM)
#if NDIM==3
        angmomsink(1) = angmomsink(1) + ms*dr(2)*dv(3) - ms*dr(3)*dv(2)
        angmomsink(2) = angmomsink(2) + ms*dr(3)*dv(1) - ms*dr(1)*dv(3)
#endif
        angmomsink(3) = angmomsink(3) + ms*dr(1)*dv(2) - ms*dr(2)*dv(1)

        ! Now add angular momentum of accreted particles to new COM
        if (ndead > 0) then
           do i=1,ndead
              p = deadlist(i)
              
              ! Ensure we only include particles accreted by this sink
              if (sinkflag(p) /= s) cycle
              mp         = real(parray(MASS,p),DP)
              rp(1:NDIM) = real(parray(1:NDIM,p),DP)
              dv(1:VDIM) = vcom_temp(1:VDIM) - real(v(1:VDIM,p),DP)
              call distance3_dp(rp(1:NDIM),rcom_temp(1:NDIM),dr(1:NDIM),drsqd)
#if NDIM==3
              angmomsink(1) = angmomsink(1) + mp*dr(2)*dv(3) - mp*dr(3)*dv(2)
              angmomsink(2) = angmomsink(2) + mp*dr(3)*dv(1) - mp*dr(1)*dv(3)
#endif
              angmomsink(3) = angmomsink(3) + mp*dr(1)*dv(2) - mp*dr(2)*dv(1)
           end do
        end if


        ! Now store all sink properties in main data structures
        sink(s)%m              = real(mtot_temp,PR)
        sink(s)%r(1:NDIM)      = real(rcom_temp(1:NDIM),PR)
        sink(s)%v(1:VDIM)      = real(vcom_temp(1:VDIM),PR)
        sink(s)%a(1:VDIM)      = real(acom_temp(1:VDIM),PR)
        sink(s)%rold(1:NDIM)   = real(rcom_temp(1:NDIM),PR)
        sink(s)%vold(1:VDIM)   = real(vcom_temp(1:VDIM),PR)
        sink(s)%angmom(1:3)    = sink(s)%angmom(1:3) + real(angmomsink(1:3),PR)
        sink(s)%angmomnet(1:3) = sink(s)%angmomnet(1:3) + real(angmomsink(1:3),PR)

#ifdef EPISODIC_ACCRETION
! add the mass of the accreted particles onto the sink's "disc" 
sink(s)%Mdisc=sink(s)%Mdisc+(mtot_temp-ms)
#endif


#if defined(LEAPFROG_KDK)
        sink(s)%vhalf(1:VDIM)  = real(vhalf_temp(1:VDIM),PR)
#endif
#if defined(PREDICTOR_CORRECTOR)
        sink(s)%aold(1:VDIM)   = real(acom_temp(1:VDIM),PR)
#endif
#if defined(DEBUG_FORCES)
        sink(s)%agrav(1:NDIM)  = sink(s)%agrav(1:NDIM) / real(mtot_temp,PR)
        sink(s)%ahydro(1:NDIM) = sink(s)%ahydro(1:NDIM) / real(mtot_temp,PR)
#endif

     end if
     ! -----------------------------------------------------------------------

     ! Calculate accretion rate and protostar properties if required
#if defined(ACCRETION_RATE)
     call sink_accretion_properties(s,maccreted)
#endif

  end do
! ============================================================================


! Record accretion history for each sink
  do s=1,stot
     pp_tot = 0
     if (ndead > 0) then
        do i=1,ndead
           p = deadlist(i)
           if (sinkflag(p) /= s) cycle
           pp_tot = pp_tot + 1
           pp_templist(pp_tot) = p
        end do
        if (pp_tot > 0) then
           call write_accreted_particles(s,pp_tot,pp_templist(1:pp_tot))
        end if
     end if
  end do

! If particles have been accreted, need to reorder arrays 
  if (ndead > 0) then
     call insertion_sort_int(ndead,deadlist(1:ndead))
     call remove_from_list(ndead,deadlist(1:ndead))
  end if

! Redistribute angular momentum contained in sinks to neighbouring particles
#if defined(SINK_REMOVE_ANGMOM)
  call redistribute_sink_angmom
#endif

! Free memory
  deallocate(pp_templist)
  deallocate(sinkflag)
  deallocate(deadlist)

! Re-stock tree immediately in next step
  nstock = nsteps

#if defined(DEBUG_ACCRETE)
  if (accreteflag) then
     write(6,*) "ptot (after accretion) : ",ptot
     call diagnostics 
  end if
#endif

! Calculate fraction of gas accreted by sinks in case of switching 
! to the N-body integrator
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  mgas = 0.0_DP
  do p=pgasstart,pgasend
     mgas = mgas + real(parray(MASS,p),DP)
  end do
  sink_frac = (mgas_orig - mgas) / mgas_orig
#endif

  return
END SUBROUTINE accrete_particles
