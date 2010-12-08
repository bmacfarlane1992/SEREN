! EWALD_INIT.F90
! A. McLeod & D. A. Hubber - 21/01/2008
! Calculates look-up tables for Ewald periodic forces
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE ewald_init
  use definitions
  use ewald_module
  use periodic_module
  implicit none

  INTERFACE
     real(kind=PR) FUNCTION erfcc(x)
       use definitions
       real(kind=PR), intent(in) :: x
     END FUNCTION erfcc
  END INTERFACE

  integer :: grid(1:NDIM)                ! ..
  integer :: h2(1:NDIM)                  ! ..
  integer :: hx                          ! ..
  integer :: hy                          ! ..
  integer :: hz                          ! ..
  integer :: i                           ! ..
  integer :: j                           ! ..
  integer :: k                           ! ..
  integer :: n(1:NDIM)                   ! ..
  integer :: nx                          ! ..
  integer :: ny                          ! ..
  integer :: nz                          ! ..
  integer, parameter :: nrange=5         ! ..
  integer, parameter :: hrange=5         ! ..

  real(kind=PR) :: alphae                ! ..
  real(kind=PR) :: dist                  ! ..
  real(kind=PR) :: ds(1:NDIM)            ! ..
  real(kind=PR) :: dr(1:NDIM)            ! ..
  real(kind=PR) :: dr2(1:NDIM)           ! ..
  real(kind=PR) :: f(1:NDIM)             ! ..
  real(kind=PR) :: k2(1:NDIM)            ! ..
  real(kind=PR) :: kr                    ! ..
  real(kind=PR) :: ksqrd                 ! ..
  real(kind=PR) :: x                     ! ..

  debug2("Initializing Ewald tables [ewald_init.F90]")

#if NDIM==1
  write (6,*) "One dimensional Ewald gravity not supported in this release"
  stop
#endif
#if NDIM==2
  allocate(fcorr(1:2,1:ewsize(1),1:ewsize(2),1))
  L = (/periodic_size(1),periodic_size(2)/)
#endif
#if NDIM==3
  allocate(fcorr(1:3,1:ewsize(1),1:ewsize(2),1:ewsize(3)))
  L = (/periodic_size(1),periodic_size(2),periodic_size(3)/)
#endif

! For a grid -L/2 < x,y,z < L/2 - we have quasi-periodic gravity
  ewsizeil = 2.0 * real(ewsize-1,PR) / L    

!    open(1, file="ewald.dat", status="old", form="formatted")
!    do i=1,ewsize(1)
!      do j=1,ewsize(2)
!        do k=1,ewsize(3)
!          read (1,*) dr(1:3), fcorr(1:3,i,j,k)
!        end do
!      end do
!    end do
!    close (1)
!
!    return

  debug1("Creating Ewald correction force table")

  alphae = 2.0 / minval(L) ! crazy scaling factor between real (close) and fourier (far) components

  ds = 1.0 / ewsizeil ! Grid cell size/spacing

! ============================================================================
  do i=1,ewsize(1)

     ! =======================================================================
#if NDIM==2 || NDIM==3
     do j=1,ewsize(2)
#endif
#if NDIM==3
        ! Iterate over all grid POINTS (which define the corners of grid CELLS
        ! =====================================================================
        !$OMP PARALLEL DO PRIVATE(grid,x,f,nx,ny,nz,hx,hy,hz,dist) &
        !$OMP PRIVATE (dr,dr2,k2,ksqrd,kr)
        do k=1,ewsize(3) 
           
           ! If quasi-periodic
           if (i == 1 .AND. j == 1 .AND. k == 1) then
#endif
#if NDIM==2
           ! quasi-periodic
           if (i == 1 .AND. j == 1) then 
#endif
           ! We are at dr=0, force is zero, skip
              cycle
           end if

#if NDIM==3
           grid = (/i,j,k/)
#endif
#if NDIM==2
           grid = (/i,j/)
#endif
           
           ! Zero correction term before recalculating for this grid point
           f = 0.0_PR 
           
           ! Grid POINT (top left corner of a grid cell, 
           ! except the last points)
           dr = real(grid - 1,PR)*ds 
           
           ! Compute first sum in ewald summation
           do nx=-nrange,nrange
#if NDIM==2 || NDIM==3
              do ny=-nrange,nrange
#endif
#if NDIM==3
                 do nz=-nrange,nrange
                    n(1:3)=(/nx,ny,nz/)
#endif
#if NDIM==2
                    n(1:2)=(/nx,ny/)
#endif
                    ! Implicit in this is that dr is actually -dr
                    dr2 = dr + L*n 
                    
                    ! Distance to periodic replica
                    dist = sqrt(sum(dr2**2))           
                    x = erfcc(alphae*dist) + (2*alphae*dist * SQRT(INVPI) &
                         & * EXP((-1.0_PR)*(alphae**2)*(dist**2)))
                    
                    ! First part of correction term
                    f = f - dr2*x/(dist**3)       
#if NDIM==3
                 end do
#endif
#if NDIM==3 || NDIM==2
              end do
#endif
           end do
           
           ! Now compute second sum in k-space
           do hx=-hrange,hrange
#if NDIM==2 || NDIM==3
              do hy=-hrange,hrange
#endif
#if NDIM==3
                 do hz=-hrange,hrange
                    if (hx == 0 .AND. hy == 0 .AND. hz == 0) cycle
                    h2(1:3) = (/hx,hy,hz/)
#endif
#if NDIM==2
                    if (hx == 0 .AND. hy == 0) cycle
                    h2(1:2) = (/hx,hy/)
#endif
                    ! Convert to k-space (k=2*PI*h/L)
                    k2 = REAL(h2,PR)*2.0*PI/L               
                    ksqrd = sum(k2**2)
                    
                    ! Because dr is actually -dr
                    kr = dot_product(k2,-dr)           
                    x = 4.0_PR*PI*exp((-1.0_PR)*ksqrd/(4.0_PR*(alphae**2)))&
                         &*sin(kr)/ksqrd
                    
                    ! Second part of correction term
                    f = f - k2*x/product(2*L)      
#if NDIM==3
                 end do
#endif
#if NDIM==3 || NDIM==2
              end do
#endif
           end do
           
           ! Subtract force from original particle
           dist = sqrt(sum(dr**2))
           f = f + dr/(dist**3)
           !f = f + dr/((dist**2 + 0.0025)**1.5) ! Temporary fixed smoothing length
           
           ! Record force in table
#if NDIM==3
           fcorr(1:3,i,j,k)=real(f(1:NDIM),PR)
#endif
#if NDIM==2
           fcorr(1:2,i,j,1)=real(f(1:NDIM),PR)
#endif
#if NDIM==3
        end do
        !$OMP END PARALLEL DO
#endif
#if NDIM==3 || NDIM==2
     end do
#endif
  end do
  
! At dr=0, fcorr = 0
  fcorr(1:NDIM,1,1,1) = 0.0_PR 
  
  open(1, file="ewald.dat", status="replace", form="formatted")
  
#if NDIM==2
  do i=1,ewsize(1)
     do j=1,ewsize(2)
        grid = (/i,j/)
        dr = real(grid - 1,PR)*ds
        write (1,*) dr(1:2), fcorr(1:2,i,j,1)
     end do
  end do
#endif
  
#if NDIM==3
  do i=1,ewsize(1)
     do j=1,ewsize(2)
        do k=1,ewsize(3)
           grid = (/i,j,k/)
           dr = real(grid-1,PR)*ds
           write (1,*) dr(1:3), fcorr(1:3,i,j,k)
        end do
     end do
  end do
#endif
  
  close (1)
     
  return
END SUBROUTINE ewald_init



! ----------------------------------------------------------------------------
function erfcc(x)
! Press et al. (1986) function to generate erfc using
! Chebyshev polynomial approximation
! (LIFTED SHAMELESSLY FROM DRAGON CODE)
  use definitions
  implicit none
  real (kind=PR), intent(in) :: x
  real (kind=PR)             :: erfcc,t,z
  z = abs(x)
  t = 1.0_PR/(1.0_PR + 0.5_PR*z)
  erfcc=t*exp(-z*z - 1.26551223_PR + t*(1.00002368_PR + t*(0.3740916_PR + &
           &t*(0.09678418_PR + t*(-0.18628806_PR + t*(0.27886807_PR + &
           &t*(-1.13520398_PR + t*(1.48851587_PR + t*(-0.82215223_PR + &
           &t*0.17087277_PR)))))))))
  if (x < 0.0_PR) erfcc = 2.0_PR - erfcc

  return
end function erfcc
