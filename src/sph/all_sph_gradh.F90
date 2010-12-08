! ALL_SPH_GRADH.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Calculates density and smoothing length of particle p using the 
! grad-h method of Price & Monaghan (2004).  Also calculates all other 
! gather-only SPH properties (e.g. velocity divergence, velocity curl, etc).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE all_sph_gradh(p)
  use interface_module, only : BHhydrowalk_hgather,distance2
  use particle_module
  use hydro_module
  use mhd_module
  use neighbour_module
  use kernel_module
  use time_module
  use sink_module
#if defined(IDEAL_MHD) && defined(INDUCTION_EQN)
  use mhd_module
#endif
  implicit none

  integer, intent(in) :: p              ! particle id
integer :: k

  integer :: i                          ! counter in neighbour search
  integer :: iteration                  ! Number of iterations
  integer :: kern                       ! kernel table value
  integer :: pp                         ! neighbouring particles (p')
  integer :: pp_max                     ! length of pot list array
  integer :: pp_pot                     ! number of potential neighbours
  integer, allocatable :: pp_potlist(:) ! list of potential neighbours
  real(kind=PR) :: div_v_p              ! Velocity divergence of particle p
  real(kind=PR) :: dr(1:NDIM)           ! vector displacements (p'-p)
!  real(kind=PR) :: drho_dt              ! Rate of change of density
  real(kind=PR) :: drmag                ! Distance
  real(kind=PR) :: drsqd                ! p'-p separation squared
!  real(kind=PR) :: dt                   ! Time since last force update 
!  real(kind=PR) :: dt_new               ! Latest timestep
  real(kind=PR) :: dt_old               ! Previous timestep
  real(kind=PR) :: dv(1:VDIM)           ! Relative velocity vector
  real(kind=PR) :: dvdr                 ! Scalar product of dv and dr
  real(kind=PR) :: fbisection           ! Aux variable for bisection method
  real(kind=PR) :: hrangesqd            ! particle radius (2*h_new) squared
  real(kind=PR) :: hfactor_p            ! invhp**(NDIM)
  real(kind=PR) :: h_high               ! high value estimate of h
  real(kind=PR) :: h_lower              ! Lower bound on h for bisection method
  real(kind=PR) :: h_new                ! New value of h
  real(kind=PR) :: h_old                ! old value of smoothing length
  real(kind=PR) :: hp                   ! Smoothing length of particle p
  real(kind=PR) :: h_upper              ! Upper bound on h for bisection method
  real(kind=PR) :: invhp                ! (1 / hp)
  real(kind=PR) :: mp                   ! Mass of particle p
  real(kind=PR) :: mpp                  ! Mass of particle pp
  real(kind=PR) :: omega_p              ! Omega correction factor for p
  real(kind=PR) :: rhotemp              ! Auxilary summation variable for rho
  real(kind=PR) :: rp(1:NDIM)           ! position of particle p
  real(kind=PR) :: skern                ! Kernel tabulation variable
  real(kind=PR) :: vp(1:NDIM)           ! Local copy of velocity of particle p
#if defined(GRAVITY)
  real(kind=PR) :: zeta_p               ! local value of zeta
#endif
#if defined(SINKS) && defined(GRAVITY)
  real(kind=PR) :: gpotmin              ! Minimum potential of neighbours 
#endif
#if !defined(BH_TREE)
  real(kind=PR) :: hrangesqd_high       ! particle radius (2*h_high) squared
#endif
#if defined(DIV_A)
  real(kind=PR) :: ap(1:NDIM)           ! accel. of particle p
  real(kind=PR) :: div_a_p              ! local value of div_a
#endif
#if defined(SMOOTHED_VELOCITY)
  real(kind=PR) :: v_smooth_p(1:VDIM)   ! smoothed velocity of particle p
#endif
#if defined(IONIZING_UV_RADIATION)
  real(kind=PR) :: gradrhotemp(1:NDIM)  ! Gradient of density field
#endif
#if defined(VISC_BALSARA)
  real(kind=PR) :: curl_v_p(1:3)        ! Curl of velocity vector (always 3-D)
  real(kind=PR) :: dvXdr(1:3)           ! Cross product of dv and dr
  real(kind=PR) :: mag_curl_v           ! Magnitude iof curl v
#endif
#if defined(VISC_PATTERN_REC)
  real(kind=PR) :: mv2, negcc, vneib, rneib, dotprod
  real(kind=PR) :: ratio, modr, modv, rcount
#endif
#if defined(IDEAL_MHD)
  real(kind=PR) :: Bp(1:BDIM)           ! B-field of particle B
  real(kind=PR) :: dBdt_temp(1:BDIM)    ! Rate of change of B
  real(kind=PR) :: dB(1:BDIM)           ! Difference in B-field
  real(kind=PR) :: div_B_p              ! Divergence of B-field for p
#endif
#if defined(DEBUG_GRAD_H_SPH)
  integer :: ngather                    ! Number of particles by gather
  real(kind=PR) :: mgather              ! Net mass by gather
#endif

  debug3("Calculating h and all gather SPH quantities [all_sph_gradh.F90] for particle ", p)

! Make local copy of old smoothing length and particle position
  rp(1:NDIM) = parray(1:NDIM,p)
  mp         = parray(MASS,p)
  hp         = parray(SMOO,p)
  vp(1:VDIM) = v(1:VDIM,p)
  rhotemp    = rho(p)
  omega_p    = omega(p)
  div_v_p    = div_v(p)
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
  Bp(1:BDIM) = B(1:BDIM,p)
#endif
#if defined(DIV_A)
  ap(1:NDIM) = a(1:NDIM,p)
#endif

! Estimate the new smoothing length h
! Determine time interval since last acceleration calculation
!  dt_old = real(laststep(p),PR)
!  dt_new = real(n - nlast(p),PR)*real(timestep,PR)
!#if defined(EULER) || defined(LEAPFROG_KDK)
!!  dt = dt_new
!  dt = dt_old
!#elif defined(RUNGE_KUTTA)
!  dt = 0.5_DP*dt_new
!#elif defined(LEAPFROG_DKD) || defined(PREDICTOR_CORRECTOR)
!  dt = 0.5_DP*(dt_old + dt_new)
!#endif
!  drho_dt = div_v_p*rhotemp/omega_p
!!  hp = hp * (1. - ((drho_dt / (NDIMPR * rhotemp))*dt))
  invhp = 1.0_PR / hp
  hfactor_p = invhp**(NDIMPR)
  hrangesqd  = KERNRANGESQD*hp*hp

! Gather neibs over a larger volume
  iteration = 0
  h_old     = hp
  h_new     = hp
  h_high    = 0.0_PR
  h_upper   = rextent
  h_lower   = 0.0_PR
  pp_pot    = 0
  pp_max    = min(LISTSIZE,ptot)
  allocate(pp_potlist(1:pp_max))


! Main iteration loop
! ============================================================================
  do

     ! Obtain new potential neighbour list either by direct sum or the tree
     ! -----------------------------------------------------------------------
     do 
        if (hp > h_high) then
           pp_pot = 0        
           h_high = hp*HMULT
#if defined(BH_TREE)
           call BHhydrowalk_hgather(rp(1:NDIM),KERNRANGE*h_high,&
                &pp_pot,pp_max,pp_potlist)
#else
           hrangesqd_high = KERNRANGESQD*h_high*h_high
           do pp=1,ptot
              if (pp == p) cycle
              call distance2(rp(1:NDIM),pp,dr(1:NDIM),drmag)
              if (drmag < hrangesqd_high) then
                 pp_pot = pp_pot + 1
                 pp_potlist(pp_pot) = pp
              end if
              if (pp_pot == pp_max) then
                 pp_pot = -1
                 exit
              end if
           end do
#endif
        end if

        ! Precaution against ridiculously high neighbour numbers
        if (pp_pot < 0) then
           h_high = 0.0_PR
           pp_pot = 0
           pp_max = ptot
           if (allocated(pp_potlist)) deallocate(pp_potlist)
           allocate(pp_potlist(1:pp_max))
        else if (pp_pot < pp_gather) then
           h_high = 0.0_PR
           pp_pot = 0
           hp     = hp*HMULT
        end if
        if (pp_pot > 0) exit
     end do 
     ! -----------------------------------------------------------------------

#if defined(DEBUG_GRAD_H_SPH)
     write(6,*) "Found ",pp_pot," potential neighbours for ",p,pp_max
     write(6,*) hp,h_new,KERNRANGE*h_high,hmin!,insinkflag
#endif

     h_new     = hp
     invhp     = 1.0 / hp
     hfactor_p = invhp**(NDIMPR)
     hrangesqd = KERNRANGESQD*hp*hp
#if defined(DEBUG_GRAD_H_SPH)
     ngather   = 0
     mgather   = 0.0_PR
#endif

     ! First, include self contributions
     rhotemp = mp*w0(0)*hfactor_p
     omega_p = mp*hfactor_p*invhp*wh(0)

     ! Now loop over all potential neighbours and calculate SPH quantities
     ! -----------------------------------------------------------------------
     do i=1,pp_pot
        pp = pp_potlist(i)
        if (pp == p) cycle
#if defined(PERIODIC)
        call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = parray(1:NDIM,pp) - rp(1:NDIM)
        drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
        if (drsqd >= hrangesqd) cycle
        drmag = sqrt(drsqd) + SMALL_NUMBER
        skern  = HALFKERNTOT * drmag * invhp
        kern  = int(skern)
        if (kern >= KERNTOT) cycle
        mpp = parray(MASS,pp)
        rhotemp = rhotemp + mpp*hfactor_p*w0(kern)
        omega_p = omega_p + mpp*hfactor_p*invhp*wh(kern)
#if defined(DEBUG_GRAD_H_SPH)
        ngather = ngather + 1
        mgather = mgather + mpp
#endif
     end do
     ! -----------------------------------------------------------------------

     omega_p = 1.0_PR + ((hp * omega_p) / (NDIMPR * rhotemp))   

#if defined(DEBUG_GRAD_H_SPH)
       write(6,*) "Before Newton-Raphson : ",hp,rhotemp,omega_p
       write(6,*) ngather,mgather
#endif

     ! First convergence check here
     if (abs(hp - h_fac*(mp/rhotemp)**(INVNDIM))/hp < H_CON &
          .and. iteration >= 1) exit


     ! Calculate new value of h_new depending on no. of iterations performed
     ! -----------------------------------------------------------------------
     if (iteration < GRADH_ITERATION_MAX) then
        h_new = h_fac*(mp / rhotemp)**(INVNDIM)
        
     ! If too many iterations have been attempted, prepare for bisection
     else if (iteration == GRADH_ITERATION_MAX) then
        h_upper = rextent
        h_lower = 0.0_PR
        h_new   = 0.5_PR*(h_lower + h_upper)
 
     ! Bisection method (after GRADH_ITERATION_MAX iterations)
     else if (iteration < 2*GRADH_ITERATION_MAX) then
        if (rhotemp > SMALL_NUMBER) then
           fbisection = hp**(NDIMPR) - ((h_fac)**(NDIMPR))*mp/rhotemp
        else
           fbisection = 1.0_PR
        end if
        if (fbisection > 0.0_PR) then
           h_upper = hp
        else
           h_lower = hp
        end if
        h_new = 0.5_PR*(h_lower + h_upper)

     ! If bisection is not converging, simply use the old value
     else
        h_new = h_old
        hp = h_new
        exit

     end if
     ! -----------------------------------------------------------------------

     hp = h_new
     iteration = iteration + 1

#if defined(MINIMUM_H)
     if (h_new < hp .and. h_new <= hmin) then
        h_new = hmin
        hp = hmin
        exit
     end if
#endif
     
#if defined(DEBUG_GRAD_H_SPH)
     write(6,*) "No of neighbours for ", p, "=",ngather,"(",pp_pot,")"
     write(6,*) "Total contained mass ",mgather,"  (",mgather/mp," mp)"
     write(6,*) "Smoothing length =", h_new," / Density =", rhotemp, &
          &" / Omega =", omega_p
     write(6,*) "Beginning iteration ", iteration
#endif

  end do
! ============================================================================


! Now we have converged on the correct value of h, calculate all other 
! SPH gather quantities.
! ============================================================================

! First initialize all quantities (and add self-contributions) depending 
! on compiler options
#if defined(MINIMUM_H)
  hp = h_fac*(mp/rhotemp)**(INVNDIM)
  if (hp < hmin) hp = hmin
#else
  hp        = h_fac*(mp/rhotemp)**(INVNDIM)
#endif
  h_new     = hp
  invhp     = 1.0_PR / hp
  hfactor_p = invhp**(NDIMPR)
  hrangesqd = KERNRANGESQD*hp*hp
  rhotemp   = mp*w0(0)*hfactor_p
  omega_p   = mp*hfactor_p*invhp*wh(0)
  div_v_p   = 0.0_PR
#if defined(DIV_A)
  div_a_p   = 0.0_PR
#endif
#if defined(SMOOTHED_VELOCITY)
  v_smooth_p(1:VDIM) = v(1:VDIM,p)*mp*hfactor_p*w0(0)
#endif
#if defined(IONIZING_UV_RADIATION)
  gradrhotemp(1:NDIM) = 0.0_PR
#endif
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = 0.0_PR
#endif
#if defined(VISC_PATTERN_REC)
  modr   = sqrt(parray(1,p)**2 + parray(2,p)**2)
  modv   = sqrt(v(1,p)**2 + v(2,p)**2)
  mv2    = modr*modv**2
  negcc  = 0.0_PR
  rcount = 0.0_PR
#endif
#if defined(GRAVITY)
  zeta_p = mp*invhp*invhp*wg(0)
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQN)
  dBdt_temp(1:BDIM) = 0.0_PR
#endif
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
  div_B_p = 0.0_PR
#endif
#if defined(SINKS)
  gpotmin = 0.0_PR
#endif

! Now loop over all gather neighbours
! ----------------------------------------------------------------------------
  do i=1,pp_pot
     pp = pp_potlist(i)
     if (pp == p) cycle
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drmag)
     if (drmag > hrangesqd) cycle
     drmag = sqrt(drmag) + SMALL_NUMBER
     skern = HALFKERNTOT * drmag * invhp
     kern  = int(skern)
     if (kern > KERNTOT) cycle
     mpp = parray(MASS,pp)
     dv(1:NDIM) = v(1:NDIM,pp) - vp(1:NDIM)
     dvdr = dot_product(dv(1:NDIM),dr(1:NDIM))
     rhotemp = rhotemp + mpp*hfactor_p*w0(kern)
     omega_p = omega_p + mpp*hfactor_p*invhp*wh(kern)
     div_v_p = div_v_p - mpp*dvdr*w4(kern)*hfactor_p*invhp/drmag
#if defined(DIV_A)
     div_a_p = div_a_p - mpp*w4(kern)*hfactor_p*invhp/&
          &dot_product(a(1:NDIM,pp)-ap(1:NDIM),dr(1:NDIM))/drmag
#endif
#if defined(SMOOTHED_VELOCITY)
     v_smooth_p(1:VDIM) = v_smooth_p(1:VDIM) + &
          &v(1:VDIM,pp)*mpp*hfactor_p*w0(kern)
#endif
#if defined(GRAVITY)
     zeta_p = zeta_p + mpp*invhp*invhp*wg(kern)
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQN)
     dBdt_temp(1:VDIM) = dBdt_temp(1:VDIM) + mpp*w4(kern)*hfactor_p*invhp* &
          & (dv(1:VDIM)*dot_product(Bp(1:NDIM),dr(1:NDIM)) &
          & - Bp(1:VDIM)*dvdr)/drmag
#endif
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
     dB(1:NDIM) = B(1:NDIM,pp) - Bp(1:NDIM)
     div_B_p = div_B_p - mpp*dot_product(dB(1:NDIM),dr(1:NDIM))&
          &*w4(kern)*hfactor_p*invhp/drmag
#endif
#if defined(VISC_BALSARA)
#if NDIM==3
     dvXdr(1) = dv(2)*dr(3) - dv(3)*dr(2)
     dvXdr(2) = dv(3)*dr(1) - dv(1)*dr(3)
#endif
#if NDIM==2 || NDIM==3
     dvXdr(3) = dv(1)*dr(2) - dv(2)*dr(1)
#endif
     curl_v_p(1:3) = curl_v_p(1:3) - &
          &(mpp * (dvXdr(1:3) / drmag) * (w1(kern)*hfactor_p*invhp))
#endif
#if defined(VISC_PATTERN_REC)
     rcount = rcount + 1.0_PR
     vneib  = sqrt(v(1,pp)**2 + v(2,pp)**2)
     rneib  = sqrt(parray(1,pp)**2 + parray(2,pp)**2)
     
     ! Two checks. If either fail, neighbour is anomalous
     dotprod = (v(1,pp)*parray(1,pp) + v(2,pp)*parray(2,pp))/vneib/rneib
     ratio   = mv2/vneib/vneib/rneib    
     if (dotprod > 0.03_PR .or. dotprod < -0.03_PR .or. &
          &ratio > 1.03_PR .or. ratio < 0.97_PR) negcc = negcc + 1
#endif
#if defined(IONIZING_UV_RADIATION)
     gradrhotemp(1:NDIM) = gradrhotemp(1:NDIM) - &
          &mpp*invhp*hfactor_p*w4(kern)*dr(1:NDIM)/drmag
#endif
#if defined(SINKS) && defined(GRAVITY)
     gpotmin = min(-gpot(pp),gpotmin)
#endif
  end do
! ----------------------------------------------------------------------------

  omega_p = 1.0_PR + ((hp * omega_p) / (NDIMPR * rhotemp))   
  div_v_p = div_v_p / rhotemp


! Store various SPH quantities in main arrays.  For the Balsara switch, 
! store magnitude of curl in balsara array for now, and calculate 
! balsara factor after thermal once we know sound speed
! ----------------------------------------------------------------------------
  rho(p)    = rhotemp
  omega(p)  = omega_p
  div_v(p)  = div_v_p
  drhodt(p) = - rhotemp * div_v_p / omega_p
  parray(SMOO,p) = hp
#if defined(GRAVITY)
  parray(ZETA,p) = -(hp*zeta_p) / (NDIMPR*rhotemp*omega_p)
#endif
#if defined(DIV_A)
  div_a(p) = div_a_p / rhotemp
#endif
#if defined(SMOOTHED_VELOCITY)
  v_smooth(1:VDIM,p) = v_smooth_p(1:VDIM) / rhotemp
#endif
#if defined(IONIZING_UV_RADIATION)
  gradrho(1:NDIM,p) = gradrhotemp(1:NDIM) / omega_p
#endif
#if defined(IDEAL_MHD) && defined(INDUCTION_EQN)
  dB_dt(1:BDIM,p) = dBdt_temp(1:BDIM) / omega_p / rhotemp
#endif
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
  div_B(p) = div_B_p / rhotemp
#endif
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = curl_v_p(1:3) / rhotemp
  mag_curl_v = sqrt(dot_product(curl_v_p(1:3),curl_v_p(1:3)))
  balsara(p) = mag_curl_v
#endif
#if defined(VISC_PATTERN_REC)
  if (rcount > 0) pattrec(p) = negcc/rcount
#endif
#if defined(SINKS) && defined(GRAVITY)
  if (gpotmin > -gpot(p)) then
     ispotmin(p) = .true.
  else
     ispotmin(p) = .false.
  endif
#endif

#if defined(DEBUG_GRAD_H_SPH)
  write(6,*) "No. of neighbours for ", p, "=",ngather,"(",pp_pot,")"
  write(6,*) "Total contained mass ",mgather,mp,mgather/mp
  write(6,*) "Smoothing length =", hp,h_old,hmin," / Density =", rhotemp, &
       &" / Omega =", omega_p
  if ((abs(hp - h_old)/h_old) > 0.1_PR) &
       & write(6,*) "Smoothing length changed by more than 10%"
#endif

  if (allocated(pp_potlist)) deallocate(pp_potlist)

  if (parray(SMOO,p) < 0.0_PR) then
     write(6,*) "Negative smoothing lengths : ",parray(SMOO,p)
     stop
  end if

  return
END SUBROUTINE all_sph_gradh
