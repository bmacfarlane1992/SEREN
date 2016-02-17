! SINK_SEARCH.F90
! D. A. Hubber - 12/6/2007; 9/6/2010
! Searches through all gas particles for possible sink candidates.  
! Creates a sink particle once the following (optional) conditions are met:
! (i)   the density of the particle exceeds the sink density (rhosink).
! (ii)  the particle lies at the bottom of its local gravitational potential 
!       well (i.e. no SPH neighbour has a more negative potential).
! (iii) the candidate particle is not overlapping with any other sinks.
! (iv)  the Hill sphere of the particle does not overlap the Hill sphere of 
!       any existing sinks.
! (v)   the local acceleration divergence is negative, i.e. particles
!       are not accelerating away from local maximum.  This would otherwise
!       indicate a transient object that is 'bouncing' back out possibly
!       due to tidal forces of a nearby dense object/sink particle.
! (vi)  the local velocity divergence is negative, i.e. particles are on 
!       on average converging.
! (vii) the particle is unlikely to be accreted by an existing sink before 
!       it has time to condense out to a new object.
! Currently creates only one sink in any one timestep for simplicity.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_search
  use interface_module, only : BHhydrowalk_hgather,create_sink,&
       &distance2,distance3,insertion_sort_real
  use type_module
  use particle_module
  use hydro_module
  use sink_module
  use scaling_module
  use kernel_module
  implicit none

  integer :: i                           ! Auxilary loop counter
  integer :: iaux                        ! Aux. variable
  integer :: j                           ! ..
  integer :: kern_p                      ! Kernel table element
  integer :: p                           ! Particle counter
  integer :: p_counter                   ! Counts number of particles in list
  integer :: pp                          ! neighbouring particles (p')
  integer :: pp_max                      ! Max. no. of potential neighbours
  integer :: pp_pot                      ! No. of potential neighbours
  integer :: psink                       ! i.d. of candidate sink particle
  integer :: s                           ! Sink counter  
  integer, allocatable :: pp_potlist(:)  ! List of potential neighbours
  integer, allocatable :: rholist(:)     ! List of particles with rho > rhosink
  real(kind=PR) :: agravp(1:NDIM)        ! Grav. accel of particle p
  real(kind=PR) :: ap(1:VDIM)            ! Total accel of particle p
  real(kind=PR) :: arad                  ! Radial acceleration
  real(kind=PR) :: curl_v_p(1:3)         ! Curl of velocity of particle p
  real(kind=PR) :: da(1:NDIM)            ! Relative acceleration
  real(kind=PR) :: div_a_p               ! Divergence of acceleration for p
  real(kind=PR) :: div_v_p               ! Divergence of velocity for p
  real(kind=PR) :: drad                  ! Radial distance
  real(kind=PR) :: dr(1:NDIM)            ! Relative position vector 
  real(kind=PR) :: dr_unit(1:NDIM)       ! Relative unit vector
  real(kind=PR) :: drmag                 ! Distance
  real(kind=PR) :: drsqd                 ! Distance squared
  real(kind=PR) :: dv(1:NDIM)            ! Relative velocity
  real(kind=PR) :: dvdr                  ! Scalar product of dv and dr
  real(kind=PR) :: dvXdr(1:3)            ! Vector product of dv and dr
  real(kind=PR) :: hfactor_p             ! invhp^(NDIM)
  real(kind=PR) :: hrange                ! Tree search radius
  real(kind=PR) :: hrangesqd             ! hrange*hrange
  real(kind=PR) :: invdrmag              ! 1 / drmag
  real(kind=PR) :: invhp                 ! Inverse of smoothing length
  real(kind=PR) :: mpp                   ! Mass of neighbour pp
  real(kind=PR) :: raux                  ! Aux. rela variable
  real(kind=PR) :: rho_hill              ! Density of gas inside Hill sphere
  real(kind=PR) :: rhotemp               ! Aux. density summation variable
  real(kind=PR) :: rp(1:NDIM)            ! Position of particle p
  real(kind=PR) :: rs(1:NDIM)            ! Position of sink particle s
  real(kind=PR) :: skern_p               ! Aux. kernel variable
  real(kind=PR) :: vp(1:VDIM)            ! Velocity of particle p
  real(kind=PR) :: vrad                  ! Radial velocity
  real(kind=PR) :: vtansqd               ! Relative tangential velocity
  real(kind=PR), allocatable :: rhovalues(:)  ! rho for candidate particles

  type sinktest_node                     ! Sink candidate test data structure
     logical :: flag                     ! Flag if particle is viable sink
     integer :: p                        ! Id of sink candidate
     real(kind=PR) :: rho                ! SPH density of sink candidate p
     real(kind=PR) :: div_v              ! SPH velocity divergence of p
     real(kind=PR) :: div_a              ! SPH acceleration divergence of p
     real(kind=PR) :: curl_v(1:3)        ! Curl of velocity
     real(kind=PR) :: tacc1              ! Simple accretion timescale
     real(kind=PR) :: tacc2              ! Sophisticated accretion timescale
     real(kind=PR) :: tcond1             ! Simple condensation timescale
     real(kind=PR) :: tcond2             ! Sophisticated condensation timescale
  end type sinktest_node
  type(sinktest_node), allocatable :: stest(:)  ! Array containing info about 
                                                ! candidate sink particles

  debug2("Searching for candidate sink particles [sink_search.F90]")
  debug_timing("SINK_SEARCH")

! Allocate and initialize some values
  p_counter = 0
  psink = -1
  iaux = 0

! First, find no. of particles with rho > rhosink which lie at the bottom of 
! their local potential well.
  do p=pgasstart,pgasend
     if (rho_search .and. rho(p) < rhosink) cycle
     if (potmin_search .and. (.not. ispotmin(p))) cycle
     p_counter = p_counter + 1
  end do

! If there are no sink candidates, return from subrouinte immediatly
  if (p_counter == 0) return

! First, allocate required memory for sinks, and then create a list of all 
! particles with density greater than the sink density that also lie at the
! bottom of potential minima
  allocate(rholist(1:p_counter))
  allocate(rhovalues(1:p_counter))
  allocate(stest(1:p_counter))
  do p=pgasstart,pgasend
     if (rho_search .and. rho(p) < rhosink) cycle
     if (potmin_search .and. (.not. ispotmin(p))) cycle
     iaux = iaux + 1
     rholist(iaux) = p
     rhovalues(iaux) = rho(p)
  end do

! Order rholist and rhovalues in ascending order of rho
  call insertion_sort_real(p_counter,rholist(1:p_counter),rhovalues(1:p_counter))

#if defined(DEBUG_SINK_SEARCH)
  write(6,*) "Found ",p_counter," sink candidate particles"
  if (p_counter > 0) then
     do i=1,p_counter
        p = rholist(i)
        write(6,*) p,porig(p),rho(p)*rhoscale*rhocgs,ispotmin(p),gpot(p)
     end do
  end if
#endif

! Set values of sink candidate particles
  do i=1,p_counter
     p = rholist(i)
     stest(i)%flag  = .true.
     stest(i)%p     = p
     stest(i)%rho   = rho(p)
     stest(i)%div_v = div_v(p)
  end do


! Determine if any of the candidate sinks overlap existing sinks, and if not, 
! see if the Hills sphere of any particles overlap those of other nearby sinks.
! ============================================================================
  if (stot > 0) then

     ! Loop over list of candidate particles
     ! -----------------------------------------------------------------------
     do i=1,p_counter
        p = stest(i)%p

        ! Set variables for candidate sink
        agravp(1:NDIM) = a_grav(1:NDIM,p)
        ap(1:NDIM)     = a(1:NDIM,p)
        rp(1:NDIM)     = parray(1:NDIM,p)
        vp(1:NDIM)     = v(1:NDIM,p)

        ! Now loop over all existing sinks
        ! --------------------------------------------------------------------
        do s=1,stot
           rs(1:NDIM) = sink(s)%r(1:NDIM)
           call distance3(rp(1:NDIM),rs(1:NDIM),dr(1:NDIM),drsqd)

           ! Sink-overlap criterion
           if (drsqd < (NEW_SINK_RMAX*sink(s)%radius)**2) &
                & stest(i)%flag = .false.
#if defined(DEBUG_SINK_SEARCH)
           if (drsqd < (NEW_SINK_RMAX*sink(s)%radius)**2) &
                &write(6,*) "Overlap between sink ",s,&
                &" and particle ",p,porig(p)
#endif

           ! Hill-sphere search
           if (hill_sphere_search) then
              da(1:NDIM) = sink(s)%agrav(1:NDIM) - agravp(1:NDIM)
              rho_hill = -(3.0_PR*dot_product(da(1:NDIM),dr(1:NDIM))) &
                   & / (4.0_PR*PI*drsqd)
              if (stest(i)%rho < 3.0_PR*rho_hill) stest(i)%flag = .false.
#if defined(DEBUG_SINK_SEARCH)
              if (stest(i)%rho < 3.0_PR*rho_hill) write(6,*) &
                   &"Hill criteria violated : ",&
                   &p,porig(p),s,stest(i)%rho,3.0_PR*rho_hill
#endif
           end if

        end do
        ! --------------------------------------------------------------------

     end do
     ! -----------------------------------------------------------------------

  end if
! ============================================================================


! If there are no sink candidates, return from subrouinte immediatly
  if (p_counter == 0) then
     if (allocated(stest)) deallocate(stest)
     if (allocated(rhovalues)) deallocate(rhovalues)
     if (allocated(rholist)) deallocate(rholist)
     return
  end if

  pp_max = ptot
  allocate(pp_potlist(1:pp_max))


! Now loop over all candidate particles and calculate relevant SPH quantities 
! with enlarged smoothing kernel (if selected in params file).
! ============================================================================
  do i=1,p_counter
     if (.not. stest(i)%flag) cycle
     p = stest(i)%p
     rp(1:NDIM) = parray(1:NDIM,p)
     vp(1:NDIM) = v(1:NDIM,p)
     ap(1:NDIM) = a(1:NDIM,p)

     ! Find all potential gather neighbours within required radius
#if defined(HMULT_SINKRAD)
     hrange = max(sinkrad*parray(SMOO,p),KERNRANGE*parray(SMOO,p))
#else
     hrange = max(sinkrad,KERNRANGE*parray(SMOO,p))
#endif
     hrangesqd = hrange*hrange
     invhp = KERNRANGE / hrange
     hfactor_p = invhp**(NDIMPR)
#if defined(BH_TREE)
     call BHhydrowalk_hgather(rp(1:NDIM),hrange,pp_pot,pp_max,pp_potlist)
#else
     pp_pot = 0
     do pp=1,ptot
        if (pp == p) cycle
        call distance2(rp(1:NDIM),pp,dr(1:NDIM),drmag)
        if (drmag < hrangesqd) then
           pp_pot = pp_pot + 1
           pp_potlist(pp_pot) = pp
        end if
     end do
#endif

     ! First, add self-contribution to various SPH quantities
     rhotemp = parray(MASS,p)*w0(0)*hfactor_p
     div_v_p = 0.0_PR
     div_a_p = 0.0_PR
     curl_v_p(1:3) = 0.0_PR

#if defined(DEBUG_SINK_SEARCH)
     write(6,*) i,p,porig(p),hrange,hfactor_p
#endif

     ! Now loop over all gather neighbours
     ! -----------------------------------------------------------------------
     do j=1,pp_pot
        pp = pp_potlist(j)
        if (pp == p) cycle
        call distance2(rp(1:NDIM),pp,dr(1:NDIM),drmag)
        if (drmag > hrangesqd) cycle
        drmag = sqrt(drmag) + SMALL_NUMBER
        invdrmag = 1.0_PR / drmag
        skern_p = HALFKERNTOT * drmag * invhp
        kern_p  = int(skern_p)
        if (kern_p > KERNTOT) cycle
        mpp = parray(MASS,pp)
        dv(1:NDIM) = v(1:NDIM,pp) - vp(1:NDIM)
        dvdr = dot_product(dv(1:NDIM),dr(1:NDIM))
        rhotemp = rhotemp + mpp*hfactor_p*w0(kern_p)
        div_v_p = div_v_p - mpp*dvdr*w4(kern_p)*hfactor_p*invhp*invdrmag
        div_a_p = div_a_p - mpp*w4(kern_p)*hfactor_p*invhp/&
             &dot_product(a(1:NDIM,pp) - ap(1:NDIM),dr(1:NDIM))*invdrmag
#if NDIM==3
        dvXdr(1) = dv(2)*dr(3) - dv(3)*dr(2)
        dvXdr(2) = dv(3)*dr(1) - dv(1)*dr(3)
#endif
#if NDIM==2 || NDIM==3
        dvXdr(3) = dv(1)*dr(2) - dv(2)*dr(1)
#endif
        curl_v_p(1:3) = curl_v_p(1:3) - &
             &mpp*dvXdr(1:3)*invdrmag*w1(kern_p)*hfactor_p*invhp
     end do
     ! -----------------------------------------------------------------------

     stest(i)%rho         = rhotemp
     stest(i)%div_v       = div_v_p / rhotemp
     stest(i)%div_a       = div_a_p / rhotemp
     stest(i)%curl_v(1:3) = curl_v_p(1:3) / rhotemp

     ! Perform div and div_a checks here if selected
     if (div_a_search .and. stest(i)%div_a > 0.0_PR) stest(i)%flag = .false.
     if (div_v_search .and. stest(i)%div_v > 0.0_PR) stest(i)%flag = .false.

#if defined(DEBUG_SINK_SEARCH)
     write(6,*) "Smoothed quantities for : ",i,p,porig(p)
     write(6,*) "rho                     : ",stest(i)%rho*rhoscale*rhocgs,&
          &rho(p)*rhoscale*rhocgs,rhosink*rhoscale*rhocgs
     write(6,*) "div_v                   : ",stest(i)%div_v
     write(6,*) "div_a                   : ",stest(i)%div_a
     write(6,*) "curl_v                  : ",stest(i)%curl_v(1:3)
#endif

  end do
! ============================================================================


! If there are no sink candidates, return from subrouinte immediatly
  if (p_counter == 0) then
     if (allocated(pp_potlist)) deallocate(pp_potlist)
     if (allocated(stest)) deallocate(stest)
     if (allocated(rhovalues)) deallocate(rhovalues)
     if (allocated(rholist)) deallocate(rholist)
     return
  end if


! Determine if any of the candidate sinks are overlap existing sinks, 
! and if not, do the Hills sphere of any particles overlap those of other 
! nearby sinks.
! ============================================================================
  if (stot > 0) then

     ! -----------------------------------------------------------------------
     do i=1,p_counter
        if (.not. stest(i)%flag) cycle
        p = stest(i)%p
        rp(1:NDIM) = parray(1:NDIM,p)
        vp(1:NDIM) = v(1:NDIM,p)
        ap(1:NDIM) = a(1:NDIM,p)
        agravp(1:NDIM) = a_grav(1:NDIM,p)


        ! Now loop over all existing sinks
        ! --------------------------------------------------------------------
        do s=1,stot
           rs(1:NDIM) = sink(s)%r(1:NDIM)
           call distance3(rp(1:NDIM),rs(1:NDIM),dr(1:NDIM),drsqd)
           drmag = sqrt(drsqd) + SMALL_NUMBER
           dr_unit(1:NDIM) = dr(1:NDIM) / drmag
           dv(1:NDIM) = sink(s)%v(1:NDIM) - vp(1:NDIM)
           da(1:NDIM) = sink(s)%a(1:NDIM) - ap(1:NDIM)

           ! Simple accretion timescale
           drad = drmag - 2.0*sink(s)%radius
           vrad = dot_product(dv(1:NDIM),dr_unit(1:NDIM))
           stest(i)%tacc1 = - drad / vrad

           ! More sophisticated accretion timescale
           vtansqd = dot_product(dv(1:NDIM),dv(1:NDIM)) - vrad*vrad
           arad = dot_product(da(1:NDIM),dr_unit(1:NDIM)) + vtansqd/drmag
           raux = (vrad/arad)**2 - 2.0_PR*drad/arad
           if (raux > 0.0_PR) then
              stest(i)%tacc2 = -vrad/arad + sqrt(raux)
           else
              stest(i)%tacc2 = BIG_NUMBER
           end if

           ! Simple condensation timescale
           stest(i)%tcond1 = - 1.0_PR / stest(i)%div_v

           ! More sophisticated condensation timescale
           div_a_p = stest(i)%div_a + &
                & dot_product(stest(i)%curl_v(1:3),stest(i)%curl_v(1:3))
           raux = (stest(i)%div_v/div_a_p)**2 - 2.0_PR/div_a_p
           if (raux > 0.0_PR) then
              stest(i)%tcond2 = -stest(i)%div_v/div_a_p + sqrt(raux)
           else
              stest(i)%tcond2 = BIG_NUMBER
           end if

           ! Now look at tests
           if (timescale_search .and. stest(i)%tcond1 > stest(i)%tacc1) &
                &stest(i)%flag = .false.
           !if (timescale_search .and. stest(i)%tcond2 > stest(i)%tacc2) &
           !&stest(i)%flag = .false.

#if defined(DEBUG_SINK_SEARCH)
           write(6,*) "Timescales for : ",i,p,porig(p)
           write(6,*) "tacc1          : ",stest(i)%tacc1*tscale
           write(6,*) "tacc2          : ",stest(i)%tacc2*tscale
           write(6,*) "tcond1         : ",stest(i)%tcond1*tscale
           write(6,*) "tcond2         : ",stest(i)%tcond2*tscale
#endif

        end do
        ! --------------------------------------------------------------------

     end do
     ! -----------------------------------------------------------------------

  end if
! ============================================================================


! Find the densest particle that has passed all the above tests, if any 
! exist at all.
! ----------------------------------------------------------------------------
  do i=p_counter,1,-1
     if (.not. stest(i)%flag) cycle

#if defined(DEBUG_SINK_SEARCH)
     write(6,*) "All tests passed!! Create sink for particle ",i,p,porig(p)
#endif

     ! If all tests passed, record id of candidate and exit loop
     psink = p
     exit
  end do


! If all conditions are met, then create sink particle 
! ----------------------------------------------------------------------------
  if (psink > 0 .and. stot < SMAX) then
     call create_sink(psink)
  else if (psink > 0 .and. stot >= SMAX) then
     stop 'Maximum number of sinks exceeded.  Need to increase SMAX.'
  end if

  if (allocated(pp_potlist)) deallocate(pp_potlist)
  if (allocated(stest)) deallocate(stest)
  if (allocated(rhovalues)) deallocate(rhovalues)
  if (allocated(rholist)) deallocate(rholist)


  return
END SUBROUTINE sink_search
