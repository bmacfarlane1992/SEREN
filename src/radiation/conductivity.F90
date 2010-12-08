! CONDUCTIVITY.F90
! Dimitris Stamatellos 2/09/09
! Calculates thermal conductivity of p
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE conductivity(p)
  use interface_module, only : distance2
  use particle_module
  use neighbour_module
  use kernel_module
  use hydro_module
  use constant_module
  use eos_module
  use scaling_module
  implicit none

  integer, intent(in) :: p               ! particle identifier

  integer :: i                           ! neighbour counter
  integer :: kern_p                      ! kernel table element for p
  integer :: kern_pp                     ! kernel table element for p
  integer :: pp                          ! neighbouring particle number
  integer :: pp_numb                     ! number of neighbours
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: dr(1:NDIM)            ! Relative displacement vector
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: drsqd                 ! separation squared
  real(kind=PR) :: dr_unit(1:NDIM)       ! Unit displacement vector
  real(kind=PR) :: hfactor_p             ! invhp ^ NDIM
  real(kind=PR) :: hfactor_pp            ! invhpp ^ NDIM
  real(kind=PR) :: hp                    ! smoothing length of particle p
  real(kind=PR) :: hpp                   ! smoothing length of neighbour pp
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: invhpp                ! ( 1 / hpp )
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: rho_p                 ! density of particle p
  real(kind=PR) :: rho_pp                ! density of particle pp
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: skern                 ! 0.5 * (r/h) * KERNTOT
  real(kind=PR) :: wmean                 ! (W(p) + W(pp)) / 2
  real(kind=PR) :: radenergy             ! radiation energy 
  real(kind=PR) :: radenergygrad(1:NDIM) ! radiation energy gradient
  real(kind=PR) :: R_diff                ! R diffusion term 
  real(kind=PR) :: kappaT                ! Pseudo opacity


! Create local copies of important properties of particle p
  rp(1:NDIM) = parray(1:NDIM,p)
  hp         = parray(SMOO,p)
  invhp      = 1.0_PR / hp
  hfactor_p  = invhp**(NDIMPR)
  rho_p      = rho(p)
#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_numb))
     do i=1,pp_numb
        pp_templist(i) = pplist(i,p)
     end do
  else
     allocate(pp_templist(1:ptot))
     call get_neib_on_fly(p,hp,pp_numb,ptot,pp_templist)
  end if
#else
  allocate(pp_templist(1:ptot))
  call get_neib_on_fly(p,hp,pp_numb,ptot,pp_templist)
#endif


! Initialise summation variables
  radenergygrad(1:NDIM) = 0.0_PR
  radenergy = 0.0_PR


! Loop over all neighbours and sum terms
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     
     ! Create local copies for neighbour pp
     rho_pp  = rho(pp)
     hpp     = parray(SMOO,pp)
     invhpp  = 1.0_PR / hpp
     hfactor_pp = invhpp**(NDIMPR)
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     invdrmag = 1.0_PR / drmag
     dr_unit(1:NDIM) = dr(1:NDIM)*invdrmag
     skern   = HALFKERNTOT * drmag * invhp
     kern_p  = int(skern)
     kern_p  = min(kern_p,KERNTOT)
     skern   = HALFKERNTOT * drmag * invhpp 
     kern_pp = int(skern)
     kern_pp = min(kern_pp,KERNTOT)
     
     wmean   = 0.5_PR*(hfactor_p*invhp*w1(kern_p) + &
          & hfactor_pp*invhpp*w1(kern_pp))
     
     radenergygrad(1:NDIM) = radenergygrad(1:NDIM) + (parray(MASS,pp)/rho_pp)&
       & * (temp(pp)**4 - temp(p)**4)*wmean*dr_unit(1:NDIM)    
     radenergy = radenergy + &
          &(parray(MASS,pp)/rho_pp)*temp(pp)**4*hfactor_p*w0(kern_p)
  end do
! ----------------------------------------------------------------------------

  call getkappa(rho(p),temp(p),idens(p),kappaT)
  
  R_diff = sqrt(dot_product(radenergygrad,radenergygrad)) &
       & / radenergy/rho_p/kappaT
  
  k_cond(p) = 16.0_PR*rad_const*((2.0_PR + R_diff) / &
       & (6.0_PR + 3.0_PR*R_diff + R_diff**2))*temp(p)**3/rho_p/kappaT
  
  deallocate(pp_templist)

  return
END SUBROUTINE conductivity
