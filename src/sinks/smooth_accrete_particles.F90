! SMOOTH_ACCRETE_PARTICLES.F90
! D. A. Hubber & S. K. Walch - 03/01/2010
! Accretes mass from SPH gas particles that have entered a sink particle.
! First, identifies which particles are inside sinks and records the id of 
! the closest sink.   Second, calculates the amount of mass to be accreted 
! by all sinks based on the amount of matter contained, and on whether the 
! matter is mainly collapsing or rotating.  Third, transfers the mass, 
! momentum and angular momentum from the SPH particles to the sink particles.
! Particles are completely accreted by sinks if their mass drops below some 
! minimum threshold or their density exceeds several times the sink density.
! ============================================================================

#include "macros.h"
#define ACCRETE_FRAC 0.005_DP
#define ACCRETE_RHO 10.0_DP
#define ACCRETE_DT 0.02_DP

! ============================================================================
SUBROUTINE accrete_particles
  use interface_module, only : distance2_dp,distance3_dp,insertion_sort_int,&
       &insertion_sort_real,remove_from_list,sink_accretion_properties,&
       &timestep_size,write_accreted_particles
  use particle_module
  use sink_module
  use hydro_module
  use type_module
  use kernel_module
  use time_module
  use scaling_module
  implicit none

  integer :: i                           ! Auxilary counter
  integer :: kern                        ! ..
  integer :: ndead                       ! Total number of dead particles
  integer :: p                           ! Particle counter
  integer :: pp_tot                      ! Number of particles accreted by s
  integer :: s                           ! Sink counter
  integer :: ss                          ! Aux sink i
  integer, allocatable :: deadlist(:)    ! List of dead (accreted) particles
  integer, allocatable :: plist(:)       ! List of particle ids
  integer, allocatable :: pp_templist(:) ! Temp. list for particle ids
  integer, allocatable :: sinkflag(:)    ! Flag sink-particle overlap
  real(kind=DP) :: acom_temp(1:VDIM)     ! Accel. of particle-sink system
  real(kind=DP) :: angmomsink(1:3)       ! Angular momentum of sink s
  real(kind=DP) :: ap(1:VDIM)            ! Acceleration of particle p
  real(kind=DP) :: as(1:VDIM)            ! Acceleration of sink s
  real(kind=DP) :: cmean                 ! ..
  real(kind=DP) :: div_v_s               ! div_v of particles inside sink
  real(kind=DP) :: dr(1:NDIM)            ! Relative position vector
  real(kind=DP) :: drmag                 ! ..
  real(kind=DP) :: drsqd                 ! Distance squared
  real(kind=DP) :: drsqd2                ! Second distance squared variable
  real(kind=DP) :: dt                    ! Timestep
  real(kind=DP) :: dv(1:VDIM)            ! Relative velocity vector
  real(kind=DP) :: gpe                   ! ..
  real(kind=DP) :: hfactor_s             ! hfactor for sink
  real(kind=DP) :: invhs                 ! 1 / hsink
  real(kind=DP) :: ketot                 ! Kinetic energy of particles in sink
  real(kind=DP) :: mcontained(1:SMAX)    ! Mass contained within sinks
  real(kind=DP) :: mp                    ! Mass of particle p
  real(kind=DP) :: ms                    ! Mass of sink s
  real(kind=DP) :: mtemp                 ! Temporary mass variable
  real(kind=DP) :: mtot_temp             ! Mass of particle-sink system
  real(kind=DP) :: rads                  ! Radius of sink particle s
  real(kind=DP) :: rcom_temp(1:NDIM)     ! COM of particle-sink system
  real(kind=DP) :: rhomax                ! ..
  real(kind=DP) :: rhomin                ! ..
  real(kind=DP) :: rotke                 ! Rotational k.e. of particles in sink
  real(kind=DP) :: rp(1:NDIM)            ! Position of particle p
  real(kind=DP) :: rs(1:NDIM)            ! Position of sink s
  real(kind=DP) :: taccrete(1:SMAX)      ! Accretion timescale
  real(kind=DP) :: tbh                   ! ..
  real(kind=DP) :: tdiv                  ! Velocity-divergence timescale
  real(kind=DP) :: tff                   ! Freefall timescale
  real(kind=DP) :: tkh                   ! Kelvin-Helmholtz timescale
  real(kind=DP) :: tvisc                 ! Viscous-disk timescale
  real(kind=DP) :: vcom_temp(1:VDIM)     ! COM velocity of particle-sink system
  real(kind=DP) :: vp(1:VDIM)            ! Velocity of particle p
  real(kind=DP) :: vs(1:VDIM)            ! Velocity of sink s
  real(kind=DP) :: weight(1:SMAX)        ! Weighting normalisation factor
  real(kind=DP) :: weight1
  real(kind=DP) :: weight2
  real(kind=DP), allocatable :: macc(:) ! Mass lost from each particle
  real(kind=PR), allocatable :: rlist(:) ! Distance squared list
#if defined(LEAPFROG_KDK)
  real(kind=DP) :: vhalf_temp(1:VDIM)    ! Half-step velocity of p-s system
#endif
#if defined(ACCRETION_RATE)
  real(kind=DP) :: maccreted             ! Total mass of accreted particles
#endif

  debug2("Accreting mass to sinks [smooth_accrete_particles.F90]")
  debug_timing("ACCRETE_PARTICLES")

! Initialize sink accretion arrays
  allocate(deadlist(1:ptot))
  allocate(plist(1:ptot))
  allocate(sinkflag(1:ptot))
  allocate(pp_templist(1:ptot))
  allocate(macc(1:ptot))
  allocate(rlist(1:ptot))
  ndead              = 0
  sinkflag(1:ptot)   = -1
  mcontained(1:SMAX) = 0.0_DP
  taccrete(1:SMAX)   = 0.0_DP
  weight(1:SMAX)     = 0.0_DP
  macc(1:ptot)      = 0.0_DP


! Find id of sink to which particle is closest and within the sink radius.
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dr,drsqd,drsqd2,mtemp,rads,rs,s,ss)
  do p=pgasstart,pgasend
     do s=1,stot
        if (.not. sink(s)%accrete) cycle

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
        end if

     end do
  end do
!$OMP END PARALLEL DO


! Calculate accretion timescale and total mass accreted from all particles
! ============================================================================
  do s=1,stot
     if (.not. sink(s)%accrete) cycle

     rs(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
     rads       = real(sink(s)%radius,DP)
     vs(1:VDIM) = real(sink(s)%v(1:VDIM),DP)
     ketot      = 0.0_DP
     rotke      = 0.0_DP
     gpe        = 0.0_DP
     rhomax     = 0.0_DP
     rhomin     = BIG_NUMBER_DP
     div_v_s    = 0.0_DP
     invhs      = real(1.0_PR / sink(s)%h,DP)
     hfactor_s  = invhs**real(NDIMPR,DP)
     cmean      = 0.0_DP
     tvisc      = 0.0_DP
     tbh        = 0.0_DP
     pp_tot     = 0
     mcontained(s) = 0.0_DP


     ! Now calculate the actual mass lost by all particles onto sinks.
     ! -----------------------------------------------------------------------
     do p=1,ptot
        if (sinkflag(p) == s) then
           call distance2_dp(rs(1:NDIM),p,dr(1:NDIM),drsqd)
           if (drsqd < rads*rads) then
              pp_tot = pp_tot + 1
              plist(pp_tot) = p
              rlist(pp_tot) = real(drsqd,PR)
           end if
        end if
     end do

     ! Cycle to next sink if there are no particles within current sink
     if (pp_tot == 0) cycle
           
     ! Sort particles into list of increasing distance from sink
     call insertion_sort_real(pp_tot,plist(1:pp_tot),rlist(1:pp_tot))


     ! Calculate the weightings for mass-removal onto sinks
     ! -----------------------------------------------------------------------
     do i=1,pp_tot
        p = plist(i)
        mp = real(parray(MASS,p),DP)
        mcontained(s) = mcontained(s) + mp
        rhomax = max(rhomax,real(rho(p),DP))
        rhomin = min(rhomin,real(rho(p),DP))
        cmean = cmean + mp*real(sound(p),DP)

        call distance2_dp(rs(1:NDIM),p,dr(1:NDIM),drsqd)
        drmag = sqrt(drsqd) + SMALL_NUMBER_DP
        dr(1:NDIM) = dr(1:NDIM) / drmag
        dv(1:VDIM) = v(1:VDIM,p) - vs(1:VDIM)
        ketot = ketot + 0.5_DP*mp*dot_product(dv(1:VDIM),dv(1:VDIM))
        
        dv(1:VDIM) = dv(1:VDIM) - &
             &dot_product(dv(1:NDIM),dr(1:NDIM))*dr(1:NDIM)
        rotke = rotke + 0.5_DP*mp*dot_product(dv(1:VDIM),dv(1:VDIM))
        
        kern  = int(HALFKERNTOT*drsqd*invhs)
        if (kern >= 0 .and. kern < KERNTOT) then
           div_v_s = div_v_s - mp*&
                &dot_product(dv(1:NDIM),dr(1:NDIM))*w4(kern)*hfactor_s*invhs
           gpe = gpe + 0.5_DP*mp*sink(s)%menc*invhs*real(w6(kern),DP)
        end if
        
        tvisc = tvisc + mp*real(sqrt(mcontained(s)*drmag)&
             &/sound(p)**2,DP)/alpha_ss
        tbh = tbh + mp*real(sound(p)**3/&
             &(4.0_DP*PI*rho(p)*(sink(s)%m + mcontained(s))),DP)
     end do



     ! Calculate the accretion timescale, and hence the total amount of mass 
     ! to be accreted by sink s this step.  If contained mass is greater than 
     ! maximum allowed, accrete excess mass.  Otherwise, accrete a small 
     ! amount based on freefall/viscous timescale.
     ! -----------------------------------------------------------------------
     tff = sqrt((3.0_DP*PI)/(32.0_DP*rhomax))
     tvisc = tvisc / mcontained(s)
     tbh = tbh / mcontained(s)

     sink(s)%menc = mcontained(s) + real(sink(s)%m,DP)
     sink(s)%trot = real(TWOPI*sqrt(sink(s)%radius**3/(sink(s)%menc)),DP)
     sink(s)%cmean = cmean / mcontained(s)
     sink(s)%tvisc = tvisc


     weight1 = 2.0_DP*rotke/gpe
     weight2 = 2.0_DP*(0.5_DP*gpe - rotke)/gpe
     taccrete(s) = 10.0_DP**(weight1*log10(tvisc) + weight2*log10(tbh))
     if (taccrete(s) > max(tff,tvisc)) taccrete(s) = max(tff,tvisc)
     if (taccrete(s) < min(tff,tvisc)) taccrete(s) = min(tff,tvisc)

     mtemp = mcontained(s)*max(1.0_DP - exp(-timestep*&
          &real(2**(level_step - nlevel_sinks),DP)/taccrete(s)),0.0_DP)

#if defined(DEBUG_SMOOTH_ACCRETE_PARTICLES)
     write(6,*) "smooth sink : ",s,pp_tot,sink(s)%radius*rscale,rotke/gpe
     write(6,*) "mcont       : ",sink(s)%menc*mscale,&
          &mcontained(s)/sink(s)%mmax,mcontained(s)/mmean
     write(6,*) "tff/tvisc/tbh : ",tff*tscale,tvisc*tscale,tbh*tscale
     write(6,*) "taccrete/trot : ",taccrete(s)*tscale,sink(s)%trot*tscale
#endif

     if (mcontained(s) > real(sink(s)%mmax,DP)) then
!        write(6,*) "TOO MUCH MASS : ",mcontained(s)/mmean,&
!             &(mcontained(s) - sink(s)%mmax)/mmean,mtemp/mmean
!        write(6,*) "NEW TIMESCALES : ",taccrete(s)*tscale, mmean/&
!             &(1.0_DP + (mcontained(s) - real(sink(s)%mmax))/mmean)
        taccrete(s) = taccrete(s) / &
             &(1.0_DP + (mcontained(s) - real(sink(s)%mmax))/mmean)
        mtemp = mcontained(s)*max(1.0_DP - exp(-timestep*&
             &real(2**(level_step - nlevel_sinks),DP)/taccrete(s)),0.0_DP)
!        write(6,*) "NEW MASS : ",mtemp/mmean,taccrete(s)*tscale
     end if
     mcontained(s) = mtemp
     mtot_temp = mtemp


     ! Now loop through list and accrete mass from nearest particles 
     ! until we reach the accreted mass limit
     ! -----------------------------------------------------------------------
     do i=1,pp_tot
        p = plist(i)
        mp = parray(MASS,p)
        mtemp = min(mtot_temp,mp)
        call timestep_size(p,dt)

!        if (mp - mtemp < smooth_accrete_frac*mmean) then
!           macc(p) = mp
!#if defined(DEBUG_ACCRETE)
!           write(6,*) "Instant accrete : particle mass too small ",&
!                &mp - mtemp,smooth_accrete_frac*mmean
!#endif
        if (dt < smooth_accrete_dt*sink(s)%trot) then
           macc(p) = mp
!#if defined(DEBUG_ACCRETE)
!           write(6,*) "Instant accrete : particle timestep too small ",&
!                &dt,smooth_accrete_dt*sink(s)%trot
!#endif
        else
           macc(p) = mtemp
        end if
        mtot_temp = mtot_temp - macc(p)
        if (mtot_temp < SMALL_NUMBER) exit

     end do

  end do
! ============================================================================


! Loop through all sinks and accrete some fraction of the mass of any 
! particle that lies within the sink radius.
! ============================================================================
  do s=1,stot

     ! Skip sink if it's not flagged as accreting
     if (.not. sink(s)%accrete) cycle

     ! Create local copies of old sink properties
     rs(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
     vs(1:VDIM) = real(sink(s)%v(1:VDIM),DP)
     as(1:VDIM) = real(sink(s)%a(1:VDIM),DP)
     ms         = real(sink(s)%m,DP)
     rads       = real(sink(s)%radius,DP)
     invhs      = 1.0_DP / real(sink(s)%h,DP)

     ! Temporary array to store new sink properties
     rcom_temp(1:NDIM) = ms*rs(1:NDIM)
     vcom_temp(1:VDIM) = ms*vs(1:VDIM)
     acom_temp(1:VDIM) = ms*as(1:VDIM)
     mtot_temp = 0.0_DP
#if defined(LEAPFROG_KDK)
     vhalf_temp(1:VDIM) = real(ms,PR)*real(sink(s)%vhalf(1:VDIM),DP)
#endif
#if defined(DEBUG_FORCES)
     sink(s)%agrav(1:NDIM) = real(ms,PR)*sink(s)%agrav(1:NDIM)
     sink(s)%ahydro(1:NDIM) = real(ms,PR)*sink(s)%ahydro(1:NDIM)
#endif
#if defined(ACCRETION_RATE)
     maccreted = 0.0_DP
#endif
     pp_tot = 0


     ! Search through all gravitating particles
     ! -----------------------------------------------------------------------
     do p=pgasstart,pgasend
        if (sinkflag(p) /= s) cycle

        ! Make local copy of particle quantities
        mp         = real(parray(MASS,p),DP)
        rp(1:NDIM) = real(parray(1:NDIM,p),DP)
        vp(1:VDIM) = real(v(1:VDIM,p),DP)
        ap(1:VDIM) = real(a(1:VDIM,p),DP)

        pp_tot    = pp_tot + 1
        mtemp     = macc(p)
        mtot_temp = mtot_temp + mtemp
        rcom_temp(1:NDIM) = rcom_temp(1:NDIM) + mtemp*rp(1:NDIM)
        vcom_temp(1:VDIM) = vcom_temp(1:VDIM) + mtemp*vp(1:VDIM)
        acom_temp(1:VDIM) = acom_temp(1:VDIM) + mtemp*ap(1:VDIM)
#if defined(LEAPFROG_KDK)
        vhalf_temp(1:VDIM) = vhalf_temp(1:VDIM) + mtemp*v_half(1:VDIM,p)
#endif
#if defined(DEBUG_FORCES)
        sink(s)%agrav(1:NDIM) = sink(s)%agrav(1:NDIM) &
             & + real(mtemp,PR)*a_grav(1:NDIM,p)
        sink(s)%ahydro(1:NDIM) = sink(s)%ahydro(1:NDIM) &
             & + real(mtemp,PR)*a_hydro(1:NDIM,p)
#endif
#if defined(ACCRETION_RATE)
        maccreted = maccreted + mtemp
#endif           

     end do
     ! -----------------------------------------------------------------------

     ! Add old sink mass
     mtot_temp = mtot_temp + ms


     ! Store new sink properties if any mass has been accreted
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
        ! --------------------------------------------------------------------
        do p=pgasstart,pgasend
           if (sinkflag(p) /= s) cycle

           rp(1:NDIM) = real(parray(1:NDIM,p),DP)
           mp = real(parray(MASS,p),DP)
        
           ! Calculate mass to be accreted again
           if (mp - macc(p) < smooth_accrete_frac*mmean) then
              mtemp = macc(p)
              pp_tot = pp_tot + 1
              ndead = ndead + 1
              deadlist(ndead) = p
#if defined(DEBUG_ACCRETE)
              write(6,*) "dead particle : ",p,ndead
#endif
           else
              mtemp = macc(p)
              pp_tot = pp_tot + 1
           end if

           dv(1:VDIM) = vcom_temp(1:VDIM) - v(1:VDIM,p)
           call distance3_dp(rp(1:NDIM),rcom_temp(1:NDIM),dr(1:NDIM),drsqd)
#if NDIM==3
           angmomsink(1) = angmomsink(1) + mtemp*(dr(2)*dv(3) - dr(3)*dv(2))
           angmomsink(2) = angmomsink(2) + mtemp*(dr(3)*dv(1) - dr(1)*dv(3))
#endif
           angmomsink(3) = angmomsink(3) + mtemp*(dr(1)*dv(2) - dr(2)*dv(1))

        end do
        ! --------------------------------------------------------------------
        

        ! Now store all sink properties in main data structures
        sink(s)%m              = real(mtot_temp,PR)
        sink(s)%r(1:NDIM)      = real(rcom_temp(1:NDIM),PR)
        sink(s)%v(1:VDIM)      = real(vcom_temp(1:VDIM),PR)
        sink(s)%a(1:VDIM)      = real(acom_temp(1:VDIM),PR)
        sink(s)%rold(1:NDIM)   = real(rcom_temp(1:NDIM),PR)
        sink(s)%vold(1:VDIM)   = real(vcom_temp(1:VDIM),PR)
        sink(s)%angmom(1:3)    = sink(s)%angmom(1:3) + angmomsink(1:3)
        sink(s)%angmomnet(1:3) = sink(s)%angmomnet(1:3) + angmomsink(1:3)
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


! Now remove particles with mass=0
  do p=1,ptot
     if (sinkflag(p) /= -1) parray(MASS,p) = parray(MASS,p) - real(macc(p),PR)
     	if    (parray(MASS,p)==0.D0) THEN
          	ndead = ndead + 1
              	deadlist(ndead) = p
	endif
  end do

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
  deallocate(rlist)
  deallocate(macc)
  deallocate(pp_templist)
  deallocate(sinkflag)
  deallocate(plist)
  deallocate(deadlist)

! Re-stock tree immediately in next step
  nstock = nsteps

#if defined(DEBUG_ACCRETE)
  write(6,*) "ptot (after accretion) : ",ptot
  call diagnostics 
  do p=1,ptot
     if (parray(MASS,p) < smooth_accrete_frac*mmean .or. &
          &parray(MASS,p) > (1.0_DP + smooth_accrete_frac)*mmean) then
        write(6,*) "What the hell is going on?? : ",p,parray(MASS,p),mmean
        stop
     end if
  end do
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

