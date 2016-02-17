! ALL_SPH.F90
! C. P. Batty & D. A. Hubber - 12/12/2006
! Calculates all SPH quantities (density, velocity divergence and velocity 
! curl if needed) for particle p using the gather method.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE all_sph(p)
  use interface_module, only : distance2,gather_neib_on_fly
  use neighbour_module
  use hydro_module
  use kernel_module
  use particle_module
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
  use mhd_module, only : B,div_B
#endif
  implicit none

  integer, intent(in) :: p              ! id of current particle

  integer :: i                          ! auxilary neighbour counter 
  integer :: kern                       ! kernel table element for p
  integer :: pp                         ! neighbouring particles (p')
  integer :: pp_numb                    ! number of neighbours 
  integer :: pp_templist(1:LISTSIZE)    ! temp. list of neighbours
  real(kind=PR) :: div_v_p              ! Local copy of velocity divergence
  real(kind=PR) :: dr(1:NDIM)           ! relative position vector
  real(kind=PR) :: dv(1:VDIM)           ! Relative velocity vector
  real(kind=PR) :: dvdr                 ! Scalar product of dv and dr
  real(kind=PR) :: drmag                ! magnitude of separation
  real(kind=PR) :: drsqd                ! separation squared
  real(kind=PR) :: hfactor              ! invhp ^ NDIM
  real(kind=PR) :: hp                   ! Smoothing length of particle p
  real(kind=PR) :: invdrmag             ! ( 1 / drmag )
  real(kind=PR) :: invhp                ! ( 1 / hp )
  real(kind=PR) :: mpp                  ! mass of neighbour pp
  real(kind=PR) :: rhotemp              ! local value of density
  real(kind=PR) :: rp(1:NDIM)           ! position of particle p
  real(kind=PR) :: skern                ! 0.5 * (r/h) * KERNTOT for p
  real(kind=PR) :: vp(1:VDIM)           ! Local copy of velocity of p
#if defined(IONIZING_UV_RADIATION)
  real(kind=PR) :: gradrhotemp(1:NDIM)  ! Gradient of density field
#endif
#if defined(SINKS) && defined(GRAVITY)
  real(kind=PR) :: gpotmin              ! Minimum potential of neighbours 
#endif
#if defined(VISC_BALSARA)
  real(kind=PR) :: curl_v_p(1:3)        ! Curl of velocity vector (always 3-D)
  real(kind=PR) :: dvXdr(1:3)           ! Cross product of dv and dr
#endif
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
  real(kind=PR) :: Bp(1:BDIM)           ! B-field of particle p
  real(kind=PR) :: dB(1:BDIM)           ! Difference in B-field
  real(kind=PR) :: div_B_p              ! Divergence of B-field for p
#endif

  debug3("Calculating SPH quantities [all_sph.F90] for particle ",p)

! Store local copies and initialize all other quantities
  rp(1:NDIM) = parray(1:NDIM,p)
  hp = parray(SMOO,p)
  invhp = 1.0_PR / hp
  hfactor = invhp**(NDIMPR)
  rhotemp = parray(MASS,p)*w0(0)*hfactor
  vp(1:VDIM) = v(1:VDIM,p)
  div_v_p = 0.0_PR
#if defined(IONIZING_UV_RADIATION)
  gradrhotemp(1:NDIM) = 0.0_PR
#endif
#if defined(SINKS)
  gpotmin = 0.0_PR
#endif
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = 0.0_PR
  dvXdr(1:3) = 0.0_PR
#endif
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
  div_B_p = 0.0_PR
  Bp(1:BDIM) = B(1:BDIM,p)
#endif

#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     pp_templist(1:pp_numb) = pplist(1:pp_numb,p)
  else
     call gather_neib_on_fly(p,parray(SMOO,p),pp_numb,pp_templist)
  end if
#else
  call gather_neib_on_fly(p,parray(SMOO,p),pp_numb,pp_templist)
#endif


! Now loop over all neighbours to find contributions
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     dv(1:VDIM) = v(1:VDIM,pp) - vp(1:VDIM)
     mpp = parray(MASS,pp)
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     dvdr    = dot_product(dv(1:NDIM),dr(1:NDIM))
     drmag   = sqrt(drsqd) + SMALL_NUMBER
     if (drmag >= KERNRANGE*hp) cycle
     invdrmag = 1.0_PR / drmag
     skern = HALFKERNTOT*drmag*invhp
     kern  = int(skern)
     kern  = min(kern,KERNTOT)
     rhotemp = rhotemp + mpp*hfactor*w0(kern)
     div_v_p = div_v_p - mpp*dvdr*w4(kern)*hfactor*invhp*invdrmag
#if defined(IONIZING_UV_RADIATION)
     gradrhotemp(1:NDIM) = gradrhotemp(1:NDIM) - &
          &mpp*invhp*hfactor*w4(kern)*invdrmag*dr(1:NDIM)
#endif
#if defined(SINKS) && defined(GRAVITY)
     gpotmin = min(-gpot(pp),gpotmin)
#endif
#if defined(VISC_BALSARA)
#if NDIM==3
     dvXdr(1) = dv(2)*dr(3) - dv(3)*dr(2)
     dvXdr(2) = dv(3)*dr(1) - dv(1)*dr(3)
#endif
#if NDIM==2 || NDIM==3
     dvXdr(3) = dv(1)*dr(2) - dv(2)*dr(1)
     curl_v_p(1:3) = curl_v_p(1:3) - mpp*dvXdr(1:3)*invdrmag* &
              & w1(kern)*hfactor*invhp
#endif
#endif
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
     dB(1:NDIM) = B(1:NDIM,pp) - Bp(1:NDIM)
     div_B_p = div_B_p - mpp*dot_product(dB(1:NDIM),dr(1:NDIM))&
          &*w1(kern)*hfactor*invhp*invdrmag
#endif

  end do
! ----------------------------------------------------------------------------


! Normalise and store SPH quantities in main arrays
  rho(p)    = rhotemp
  div_v_p   = div_v_p / rhotemp
  div_v(p)  = div_v_p
  drhodt(p) = -rhotemp * div_v_p
#if defined(IONIZING_UV_RADIATION)
  gradrho(1:NDIM,p) = gradrhotemp(1:NDIM)
#endif
#if defined(SINKS) && defined(GRAVITY)
  if (gpotmin > -gpot(p)) then
     ispotmin(p) = .true.
  else
     ispotmin(p) = .false.
  endif
#endif
 
! Store magnitude of curl in balsara array for now, and calculate 
! balsara factor after thermal once we know sound speed.
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = curl_v_p(1:3) / rhotemp
  balsara(p) = sqrt(dot_product(curl_v_p(1:3),curl_v_p(1:3)))
#endif  
#if defined(IDEAL_MHD) && defined(DEBUG_MHD)
  div_B(p) = div_B_p / rhotemp
#endif

  return
END SUBROUTINE all_sph
