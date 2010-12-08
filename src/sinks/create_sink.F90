! CREATE_SINK.F90
! D. A. Hubber - 12/6/2007
! Creates a new sink particle at 
! i) centre of mass of particle psink and its neighbours (simple sinks), 
! ii) position of particle psink (smoothly accreting sinks).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE create_sink(psink)
  use interface_module, only : distance3_dp,gather_neib_on_fly,&
       &insertion_sort_int,remove_from_list,write_accreted_particles
  use particle_module
  use neighbour_module
  use sink_module
  use type_module
  use time_module, only : nbuild,nsteps,nstock,time
  implicit none

  integer, intent(in) :: psink            ! Particle identifier

  integer :: ndead                        ! Number of dead particles
  integer :: pp_tot                       ! Number of particles in sink radius
  integer :: s                            ! Sink particle id
  integer :: pp_templist(1:LISTSIZE)      ! Temp. list of neighbours
  integer, allocatable :: deadlist(:)     ! List of dead particles
  real(kind=DP) :: asink(1:NDIM)          ! Acceleration of new sink 
  real(kind=DP) :: angmomsink(1:3)        ! Angular momentum of sink 
  real(kind=DP) :: drsqd                  ! Distance squared
  real(kind=DP) :: dr(1:NDIM)             ! Vector displacement
  real(kind=DP) :: dv(1:NDIM)             ! Velocity difference
  real(kind=DP) :: hp                     ! Smoothing length of particle p
  real(kind=DP) :: mp                     ! Mass of particle p
  real(kind=DP) :: msink                  ! Mass of new sink 
  real(kind=DP) :: rads                   ! Radius of new sink
  real(kind=DP) :: rp(1:NDIM)             ! Position of particle p
  real(kind=DP) :: rsink(1:NDIM)          ! Position of new sink 
  real(kind=DP) :: vp(1:NDIM)             ! Velocity of particle p
  real(kind=DP) :: vsink(1:NDIM)          ! Velocity of new sink 
#if defined(LEAPFROG_KDK)
  real(kind=DP) :: vshalf(1:VDIM)         ! Half-step velocity of p-s system
#endif
#if !defined(SMOOTH_ACCRETION)
  integer :: i                            ! counter over pp for copying
  integer :: pp                           ! neighbouring particles (p')
  real(kind=DP) :: mpp                    ! Mass of particle pp
  real(kind=DP) :: rpp(1:NDIM)            ! Position of particle pp
#endif
#if defined(DEBUG_FORCES)
  real(kind=DP) :: as_grav(1:NDIM)        ! New sink grav. force
  real(kind=DP) :: as_hydro(1:NDIM)       ! New sink hydro force
#endif

#if defined(DEBUG2) || defined(DEBUG_CREATE_SINK)
  write(6,*) "Creating sink [create_sink.F90] at position of particle ",psink
  call diagnostics
#endif

! Initialize sink accretion arrays
  allocate(deadlist(1:ptot))


! Include particle psink in sink properties
! ----------------------------------------------------------------------------
  rp(1:NDIM)    = real(parray(1:NDIM,psink),DP)
  mp            = real(parray(MASS,psink),DP)
  hp            = real(parray(SMOO,psink),DP)
  vp(1:NDIM)    = real(v(1:NDIM,psink),DP)
  msink         = mp
  rsink(1:NDIM) = mp*rp(1:NDIM)
  vsink(1:NDIM) = mp*vp(1:NDIM)
  asink(1:NDIM) = mp*real(a(1:NDIM,psink),DP)
#if defined(LEAPFROG_KDK)
  vshalf(1:NDIM) = mp*real(v_half(1:NDIM,psink),DP)
#endif
#if defined(DEBUG_FORCES)
  as_grav(1:NDIM)  = mp*real(a_grav(1:NDIM,psink),DP)
  as_hydro(1:NDIM) = mp*real(a_hydro(1:NDIM,psink),DP)
#endif
#if defined(HMULT_SINKRAD)
  rads = real(sinkrad*hp,DP)
#else
  rads = real(sinkrad,DP)
#endif

! Gather all neighbours within sink radius
  call gather_neib_on_fly(psink,INVKERNRANGE*real(rads,PR),pp_tot,pp_templist)


! Now loop over neighbours and add contributions to sink properties 
! (i.e. total mass, centre of mass, and velocity of centre of mass)
! ----------------------------------------------------------------------------
#if !defined(SMOOTH_ACCRETION)
  do i=1,pp_tot
     pp = pp_templist(i)
     mpp = real(parray(MASS,pp),DP)
     msink = msink + mpp
     rsink(1:NDIM) = rsink(1:NDIM) + mpp*real(parray(1:NDIM,pp),DP)
     vsink(1:NDIM) = vsink(1:NDIM) + mpp*real(v(1:NDIM,pp),DP)
     asink(1:NDIM) = asink(1:NDIM) + mpp*real(a(1:NDIM,pp),DP)
#if defined(LEAPFROG_KDK)
     vshalf(1:NDIM) = vshalf(1:NDIM) + mpp*real(v_half(1:NDIM,pp),DP)
#endif
#if defined(DEBUG_FORCES)
     as_grav(1:NDIM)  = as_grav(1:NDIM) + mpp*real(a_grav(1:NDIM,pp),DP)
     as_hydro(1:NDIM) = as_hydro(1:NDIM) + mpp*real(a_hydro(1:NDIM,pp),DP)
#endif
  end do
#endif

! Position, velocity and acceleration of new sink particle
  rsink(1:NDIM) = rsink(1:NDIM) / msink
  vsink(1:NDIM) = vsink(1:NDIM) / msink
  asink(1:NDIM) = asink(1:NDIM) / msink
#if defined(LEAPFROG_KDK)
  vshalf(1:NDIM) = vshalf(1:NDIM) / msink
#endif
#if defined(DEBUG_FORCES)
  as_grav(1:NDIM)  = as_grav(1:NDIM) / msink
  as_hydro(1:NDIM) = as_hydro(1:NDIM) / msink
#endif


! Compute angular momentum of particles in sink about centre of mass 
! ----------------------------------------------------------------------------
! First include angular momentum of particle psink
  call distance3_dp(rp(1:NDIM),rsink(1:NDIM),dr(1:NDIM),drsqd)
  dv(1:NDIM) = vsink(1:NDIM) - vp(1:NDIM)
#if NDIM==3
  angmomsink(1) = mp*dr(2)*dv(3) - mp*dr(3)*dv(2)
  angmomsink(2) = mp*dr(3)*dv(1) - mp*dr(1)*dv(3)
#else
  angmomsink(1) = 0.0_DP ; angmomsink(2) = 0.0_DP
#endif
  angmomsink(3) = mp*dr(1)*dv(2) - mp*dr(2)*dv(1)


! Loop over neighbours and add contribution to angular momentum
! ----------------------------------------------------------------------------
#if !defined(SMOOTH_ACCRETION)
  do i=1,pp_tot
     pp = pp_templist(i)
     mpp = parray(MASS,pp)
     rpp(1:NDIM) = parray(1:NDIM,pp)
     dv(1:NDIM) = vsink(1:NDIM) - v(1:NDIM,pp)
     call distance3_dp(rpp(1:NDIM),rsink(1:NDIM),dr(1:NDIM),drsqd)
#if NDIM==3
     angmomsink(1) = angmomsink(1) + mpp*dr(2)*dv(3) - mpp*dr(3)*dv(2)
     angmomsink(2) = angmomsink(2) + mpp*dr(3)*dv(1) - mpp*dr(1)*dv(3)
#endif
     angmomsink(3) = angmomsink(3) + mpp*dr(1)*dv(2) - mpp*dr(2)*dv(1)
  end do
#endif


! Set sink particle properties
! ----------------------------------------------------------------------------
! Increase total sink counter
  stot = stot + 1
  s = stot

! Record sink properties
  sink(s)%accrete        = .true.
  sink(s)%static         = .false.
  sink(s)%id             = s
  sink(s)%ncreate        = int(nsteps)
  sink(s)%tcreate        = time
  sink(s)%m              = real(msink,PR)
  sink(s)%r(1:NDIM)      = real(rsink(1:NDIM),PR)
  sink(s)%v(1:NDIM)      = real(vsink(1:NDIM),PR)
  sink(s)%a(1:NDIM)      = real(asink(1:NDIM),PR)
  sink(s)%rold(1:NDIM)   = real(rsink(1:NDIM),PR)
  sink(s)%vold(1:NDIM)   = real(vsink(1:NDIM),PR)
#if defined(LEAPFROG_KDK)
  sink(s)%vhalf(1:NDIM)  = real(vshalf(1:NDIM),PR)
#endif
#if defined(PREDICTOR_CORRECTOR)
  sink(s)%aold(1:NDIM)   = real(asink(1:NDIM),PR)
#endif
#if defined(DEBUG_FORCES)
  sink(s)%agrav(1:NDIM)  = real(as_grav(1:NDIM),PR)
  sink(s)%ahydro(1:NDIM) = real(as_hydro(1:NDIM),PR)
#endif

! Radius and sink h are multiples of value of h of psink
  sink(s)%h              = real(INVKERNRANGE*rads,PR)
  sink(s)%radius         = real(rads,PR)
  sink(s)%angmom(1:3)    = angmomsink(1:3)
  sink(s)%angmomnet(1:3) = angmomsink(1:3)

! Initialize other properties
  sink(s)%dmdt           = 0.0_DP
  sink(s)%star_radius    = 0.0_DP
  sink(s)%luminosity     = 0.0_DP
  sink(s)%temperature    = 0.0_DP
  sink(s)%macc(1:DMDT_RANGE) = 0.0_DP
  sink(s)%tacc(1:DMDT_RANGE) = 0.0_DP
#if defined(BH_TREE) && defined(GADGET_MAC)
  sink(s)%agravmag       = 0.0_DP
#endif
#if defined(SMOOTH_ACCRETION)
  sink(s)%mmax = (4.0_PR*PI/3.0_PR)*rhosink*real(rads,PR)**3
  sink(s)%menc = sink(s)%mmax
#endif
#if defined(DEBUG_CREATE_SINK)
  write(6,*) "angmom : ",angmomsink(1:3)
#endif

! Only accrete particle itself with smooth-accretion.
! Otherwise, accrete all particles within new sink radius
! ----------------------------------------------------------------------------
#if defined(SINKS) && defined(SMOOTH_ACCRETION)
  ndead = 1
  deadlist(ndead) = psink
#else
  ndead = pp_tot + 1
  deadlist(1:pp_tot) = pp_templist(1:pp_tot)
  deadlist(ndead) = psink  
#endif

#if defined(DEBUG_CREATE_SINK)
  write(6,*) "stot  : ",stot
  write(6,*) "ndead : ",ndead
#endif

! Need to ensure lists of dead particles are sorted in ascending order
  call insertion_sort_int(ndead,deadlist(1:ndead))

! Record accretion history for new sink to file
  call write_accreted_particles(s,ndead,deadlist(1:ndead))

! Now shuffle all live particles down to fill in gaps
  call remove_from_list(ndead,deadlist(1:ndead))

! Free memory
  deallocate(deadlist)

! Need to re-build trees immediately on next timestep
  nbuild = nsteps 
  nstock = nsteps

! If a sink is created when using smooth accretion, set hmin = hp
#if defined(SMOOTH_ACCRETION) && defined(MINIMUM_H) && NDIM==3
  hmin = INVKERNRANGE*&
       &((3.0_PR*pp_gather*real(mgas,PR))/(4.0_PR*PI*pgas*rhosink))**(ONETHIRD)
!  hmin = INVKERNRANGE*sink(s)%radius
#endif

#if defined(DEBUG_CREATE_SINK)
  write(6,*) "Finished removing particles"
  write(6,*) "ptot : ",ptot
  write(6,*) "stot : ",stot
  call diagnostics
#endif

  return
END SUBROUTINE create_sink
