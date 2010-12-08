! BHGRAV_ACCEL.F90
! D. A. Hubber - 23/01/2008
! Computes gravitational force exerted on particle p due to all other 
! particles by walking the BH tree.  Starting on the first level, each 
! cell is checked with one the opening criterion
! i) Opening angle
! ii) Octupole moment error term (GADGET MAC)
! iii) Maximum absolute multipole moment error (Salmon & Warren MAC)
! iv) Maximum quadrupole error using eigenvalues
! If the criterion is satisfied, the cell is not opened and the 
! gravitational force is given by the centre of mass of the cell plus 
! quadrupole moment correction terms.  If the above criterion is not 
! satisfied, then the cell is opened and the sub-cells are checked.  
! If a leaf cell is opened, then the contributions due to the particles 
! within the cell are added directly. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHgrav_accel(p,invhp,rp,agravp,potp)
  use interface_module, only : distance3,ewald_force,&
       &gravity_gradh,gravity_nbody,gravity_sph
  use definitions
  use tree_module
  use particle_module
#if defined(SINKS)
  use sink_module
#endif
#if defined(GRAD_H_SPH) && !defined(N_BODY)
  use hydro_module, only : omega
#endif
  implicit none

  integer, intent(in) :: p                      ! Id of current particle
  real(kind=PR), intent(in) :: invhp            ! Smoothing length of p
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(out) :: agravp(1:NDIM)  ! Gravitational accelertation
  real(kind=DP), intent(out) :: potp            ! Gravitational potential

  integer :: c                      ! Cell counter
  integer :: i                      ! Auxilary counter
  integer :: nlist                  ! Number of particles in treelist
  integer :: pp                     ! Second particle identifier
  integer :: treelist(1:GLISTSIZE)  ! List of nearby particle ids
  real(kind=PR) :: atemp(1:NDIM)    ! Grav acceleration contribution vector
  real(kind=PR) :: dr(1:NDIM)       ! Relative displacement vector
  real(kind=PR) :: dpot             ! Grav. potential contribution
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: invdrmag         ! ( 1 / drmag )
  real(kind=PR) :: invdrsqd         ! ( 1 / drsqd )
#if defined(QUADRUPOLE) || defined(OCTUPOLE)
  real(kind=PR) :: invdr5           ! ( 1 / drmag^5 )
#endif
#if defined(QUADRUPOLE)
  real(kind=PR) :: qscalar          ! Inner product of quad tensor
#endif
#if defined(OCTUPOLE)
  real(kind=PR) :: sscalar          ! Scalar component of octupole terms
#endif
#if defined(GRAD_H_SPH) && !defined(N_BODY)
  real(kind=PR) :: zo_p             ! local copy of zeta/omega for p
#endif
#if defined(SINKS)
  integer :: s                      ! Sink counter
  integer :: ss                     ! Sink counter
#endif
#if defined(EWALD)
  real(kind=PR) :: eaccel(1:NDIM)   ! Ewald grav. acceleration
#endif
#ifndef GEOMETRIC_MAC
  real(kind=PR) :: afactor          ! 'old' acceleration factor for MACs
#endif

  debug3("Calculating gravity forces [BHtreegravity.F90] for particle ", p)

! Zero arrays and variables
  agravp(1:NDIM) = 0.0_DP
  potp = 0.0_DP
  nlist = 0
#if defined(SINKS)
  if (p < 0) s = -p
  if (p >= 0) s = 0
#endif

! Make local copies of zeta/omega if required
#if defined(GRAD_H_SPH) && !defined(N_BODY)
  if (p < 0) then
     zo_p = 0.0_PR
  else if (p <= ptot) then
     zo_p = parray(ZETA,p)
  end if
#endif

! Prepare various MAC variables for tree walk
#ifndef GEOMETRIC_MAC
  afactor = 0.0_PR
#if defined(SINKS)
  if (p < 0) afactor = sink(s)%agravmag
#endif
  if (p > 0) afactor = agravmag(p)
  if (afactor > SMALL_NUMBER) then
#if defined(GADGET_MAC)
     afactor = afactor**(-ONETHIRD)
#elif defined(GADGET2_MAC)
     afactor = 1.0_PR/sqrt(afactor)
#elif defined(EIGEN_MAC)
     afactor = (1.0_PR/gpot(p))**(2.0_PR*ONETHIRD)
#endif
  else
     afactor = 0.0_PR
  end if
#endif

! Start on cell 0 and work our way down
  c = 0

! Walk gravity tree until we reach end cell pointer
! ============================================================================
  do

     ! If cell is leaf cell with just one particle, more efficient to 
     ! simply record particle id in list straight away
     ! -----------------------------------------------------------------------
     if (BHgrav(c)%leaf == 1) then
        pp = BHgrav(c)%plist(1)
        if (pp /= p) then
           nlist = nlist + 1
           treelist(nlist) = pp
        end if

        ! Point to next cell in list
        c = BHgrav(c)%nextcell
        
     ! Else check the opening angle of the cell
     ! -----------------------------------------------------------------------
     else

        ! Calculate distance depending on whether we use Ewald gravity
#if defined(EWALD)
        call distance3(BHgrav(c)%r(1:NDIM),rp(1:NDIM),dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = rp(1:NDIM) - BHgrav(c)%r(1:NDIM)
        drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif

        ! If distance between p and c is greater than min distance, and the 
        ! cell is not a leaf cell containing one particle, then calculate 
        ! gravitational acceleration due to cell
        ! --------------------------------------------------------------------
#if defined(GEOMETRIC_MAC)
        if (drsqd > BHgrav(c)%dminsqd) then
#else 
        if (drsqd > BHgrav(c)%mac*afactor .and. drsqd > BHgrav(c)%dminsqd) then
#endif        
           invdrsqd = 1.0_PR / (drsqd + SMALL_NUMBER)
           invdrmag = sqrt(invdrsqd)
           atemp(1:NDIM) = -BHgrav(c)%m*invdrsqd*invdrmag*dr(1:NDIM)
           dpot = BHgrav(c)%m*invdrmag

           ! Add quadrupole moment correction terms
#if defined(QUADRUPOLE)
           invdr5 = invdrsqd*invdrsqd*invdrmag
#if NDIM==2
           qscalar = BHgrav(c)%q(1)*dr(1)*dr(1) + BHgrav(c)%q(3)*dr(2)*dr(2) &
                    & + 2.0_PR*BHgrav(c)%q(2)*dr(1)*dr(2)
           atemp(1) = atemp(1) + (BHgrav(c)%q(1)*dr(1) + BHgrav(c)%q(2)*dr(2))&
                     &*invdr5 - 2.5_PR*qscalar*dr(1)*invdr5*invdrsqd
           atemp(2) = atemp(2) + (BHgrav(c)%q(2)*dr(1) + BHgrav(c)%q(3)*dr(2))&
                     &*invdr5 - 2.5_PR*qscalar*dr(2)*invdr5*invdrsqd
#elif NDIM==3
           qscalar = BHgrav(c)%q(1)*dr(1)*dr(1) + BHgrav(c)%q(3)*dr(2)*dr(2) &
                & - (BHgrav(c)%q(1) + BHgrav(c)%q(3))*dr(3)*dr(3) &
                & + 2.0_PR*(BHgrav(c)%q(2)*dr(1)*dr(2) &
                & + BHgrav(c)%q(4)*dr(1)*dr(3) + BHgrav(c)%q(5)*dr(2)*dr(3))
           atemp(1) = atemp(1) + (BHgrav(c)%q(1)*dr(1) + BHgrav(c)%q(2)*dr(2)&
                & + BHgrav(c)%q(4)*dr(3))*invdr5 &
                & - 2.5_PR*qscalar*dr(1)*invdr5*invdrsqd
           atemp(2) = atemp(2) + (BHgrav(c)%q(2)*dr(1)+BHgrav(c)%q(3)*dr(2)&
                & + BHgrav(c)%q(5)*dr(3))*invdr5 &
                & - 2.5_PR*qscalar*dr(2)*invdr5*invdrsqd
           atemp(3) = atemp(3) + (BHgrav(c)%q(4)*dr(1) + BHgrav(c)%q(5)*dr(2)&
                & - (BHgrav(c)%q(1)+BHgrav(c)%q(3))*dr(3))*invdr5 &
                & - 2.5_PR*qscalar*dr(3)*invdr5*invdrsqd
#endif
           dpot = dpot + 0.5_PR*qscalar*invdr5 
#endif
           ! Add octupole moment correction terms
#if defined(OCTUPOLE)
           sscalar = BHgrav(c)%s(1)*dr(1)*dr(1)*dr(1) + &
                   & BHgrav(c)%s(2)*dr(1)*dr(1)*dr(2) + &
                   & BHgrav(c)%s(3)*dr(2)*dr(2)*dr(1) + &
                   & BHgrav(c)%s(4)*dr(2)*dr(2)*dr(2) + &
                   & BHgrav(c)%s(5)*dr(3)*dr(3)*dr(1) + &
                   & BHgrav(c)%s(6)*dr(3)*dr(3)*dr(2) + &
                   & BHgrav(c)%s(7)*dr(3)*dr(3)*dr(3) + &
                   & BHgrav(c)%s(8)*dr(1)*dr(1)*dr(3) + &
                   & BHgrav(c)%s(9)*dr(2)*dr(2)*dr(3) + &
                   & BHgrav(c)%s(10)*dr(1)*dr(2)*dr(3)
           atemp(1) = atemp(1) + 0.5_PR*(3.0_PR*BHgrav(c)%s(1)*dr(1)*dr(1) + &
                & 2.0_PR*BHgrav(c)%s(2)*dr(1)*dr(2) + &
                & 2.0_PR*BHgrav(c)%s(8)*dr(1)*dr(3) + &
                & BHgrav(c)%s(3)*dr(2)*dr(2) + BHgrav(c)%s(5)*dr(3)*dr(3) + &
                & BHgrav(c)%s(10)*dr(2)*dr(3))*invdr5*invdrsqd - &
                & 3.5_PR*sscalar*dr(1)*invdr5*invdrsqd*invdrsqd
           atemp(2) = atemp(2) + 0.5_PR*(3.*BHgrav(c)%s(4)*dr(2)*dr(2) + &
                & 2.0_PR*BHgrav(c)%s(3)*dr(1)*dr(2) + &
                & 2.0_PR*BHgrav(c)%s(9)*dr(2)*dr(3) + &
                & BHgrav(c)%s(2)*dr(1)*dr(1) + BHgrav(c)%s(6)*dr(3)*dr(3) + &
                & BHgrav(c)%s(10)*dr(1)*dr(3))*invdr5*invdrsqd - &
                & 3.5_PR*sscalar*dr(2)*invdr5*invdrsqd*invdrsqd
           atemp(3) = atemp(3) + 0.5_PR*(3.*BHgrav(c)%s(7)*dr(3)*dr(3) + &
                & 2.0_PR*BHgrav(c)%s(5)*dr(1)*dr(3) + &
                & 2.0_PR*BHgrav(c)%s(6)*dr(2)*dr(3) + &
                & BHgrav(c)%s(8)*dr(1)*dr(1) + BHgrav(c)%s(9)*dr(2)*dr(2) + &
                & BHgrav(c)%s(10)*dr(1)*dr(2))*invdr5*invdrsqd - &
                & 3.5_PR*sscalar*dr(3)*invdr5*invdrsqd*invdrsqd
           dpot = dpot + 0.5_PR*sscalar*invdr5*invdrsqd
#endif

           ! Add contribution due to cell c to summation vector
           agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
           potp = potp + real(dpot,DP)

           ! Add Ewald correction force if required
#if defined(EWALD)
           call ewald_force(dr,BHgrav(c)%m,eaccel)
           agravp(1:NDIM) = agravp(1:NDIM) + real(eaccel(1:NDIM),DP)
#endif

           ! Move to next cell
           c = BHgrav(c)%nextcell

        ! If distance is smaller, open cell
        ! --------------------------------------------------------------------
        else
	
           ! Check if cell is actually a leaf cell
           if (BHgrav(c)%leaf > 0) then

              ! If so, loop over all particles in leaf cell
              do i=1,BHgrav(c)%leaf
                 pp = BHgrav(c)%plist(i)
                 if (pp == p) cycle
                 nlist = nlist + 1
                 treelist(nlist) = pp
              end do

              ! Point to next cell in list
              c = BHgrav(c)%nextcell

           ! If it's not a leaf cell, open cell.  If it's a dead cell 
           ! (e.g. due to particle accretion), point to next cell in list.
           else if (BHgrav(c)%leaf == 0) then
              c = BHgrav(c)%ifopen
           else
              c = BHgrav(c)%nextcell
           end if
        end if
        ! --------------------------------------------------------------------

     end if
     ! -----------------------------------------------------------------------

     ! Loop over particle list and add contributions to grav. accel.
     ! This is done when i) the particle id buffer has been filled up, 
     ! and/or ii) we have finished traversing the tree. 
     ! -----------------------------------------------------------------------
     if ((nlist > GLISTSIZE - LEAFMAX .or. (nlist > 0 .and. c > ctot_grav)) &
          & .and. p > 0) then
        do i=1,nlist
           pp = treelist(i)
#if defined(N_BODY)
           call gravity_nbody(parray(MASS,pp),rp(1:NDIM),&
                &parray(1:NDIM,pp),atemp(1:NDIM),dpot)
#elif defined(GRAD_H_SPH)
           call gravity_gradh(invhp,parray(SMOO,pp),parray(MASS,pp),rp(1:NDIM)&
                &,parray(1:NDIM,pp),zo_p,parray(ZETA,pp),atemp(1:NDIM),dpot)
#else
           call gravity_sph(invhp,parray(SMOO,pp),parray(MASS,pp),&
                &rp(1:NDIM),parray(1:NDIM,pp),atemp(1:NDIM),dpot)
#endif
           agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
           potp = potp + real(dpot,DP)
        end do
        nlist = 0
     end if
     ! -----------------------------------------------------------------------


     ! Loop over particle list and add contributions to grav. accel.
     ! This is done when i) the particle id buffer has been filled up, 
     ! and/or ii) we have finished traversing the tree. 
     ! -----------------------------------------------------------------------
#if defined(SINKS)
     if ((nlist > GLISTSIZE - LEAFMAX .or. (nlist > 0 .and. c > ctot_grav)) &
          & .and. p < 0) then
        do i=1,nlist
           pp = treelist(i)
#if defined(N_BODY)
           call gravity_nbody(parray(MASS,pp),rp(1:NDIM),&
                &parray(1:NDIM,pp),atemp(1:NDIM),dpot)
#else
           call gravity_sph(invhp,parray(SMOO,pp),parray(MASS,pp),&
                &rp(1:NDIM),parray(1:NDIM,pp),atemp(1:NDIM),dpot)
#endif
           agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
           potp = potp + real(dpot,DP)
        end do
        nlist = 0
     end if
#endif
     ! -----------------------------------------------------------------------


     ! Exit loop if we have finished traversing the tree
     if (c > ctot_grav) exit

  end do
! ============================================================================

! Record potential due to only SPH particles for polytropic cooling method
#if defined(RAD_WS)
  if (p > 0) sphgpot(p) = real(potp,PR)
#endif

! Include contribution due to sink particles        
! ----------------------------------------------------------------------------
#if defined(SINKS)
  do ss=1,stot
     if (s == ss) cycle
#if defined(N_BODY)
     call gravity_nbody(sink(ss)%m,rp,sink(ss)%r(1:NDIM),atemp,dpot)
#else
     call gravity_sph(invhp,sink(ss)%h,sink(ss)%m,rp(1:NDIM),&
          &sink(ss)%r(1:NDIM),atemp,dpot)
#endif
     agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
     potp = potp + real(dpot,DP)
  end do
#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE BHgrav_accel
