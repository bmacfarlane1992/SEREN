! REORDER_ARRAY.F90
! D. A. Hubber - 30/7/2008
! Reorders elements in all particle arrays to the order given in array aorder.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE reorder_particle_arrays(pstart,plast,dummylist)
  use interface_module, only : reorder_array_double,reorder_array_int,&
       &reorder_array_int_2D,reorder_array_logical,reorder_array_long_int,&
       &reorder_array_real,reorder_array_real_2D,reorder_inverse_array_int
  use definitions
  use particle_module
  use hydro_module
  use time_module
  use type_module
  use tree_module
  use Eos_module
  use HP_module
  implicit none

  integer, intent(in) :: pstart             ! id of first particle
  integer, intent(in) :: plast              ! id of last particle
  integer, intent(in) :: dummylist(1:ptot)  ! Array to be reordered

  debug2("Reorder all particle arrays to given order [reorder_particle_arrays]")

! Now re-order each array separately.  Uses OpenMP sections to 
! parallelize otherwise serial code.
! ============================================================================
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 

! Main particle data
! ----------------------------------------------------------------------------
!$OMP SECTION 
  call reorder_array_int(pstart,plast,porig,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(DATATOT,pstart,plast,parray,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,v,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,a,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(NDIM,pstart,plast,r_old,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,v_old,dummylist)

#if defined(RUNGE_KUTTA) || defined(LEAPFROG_KDK)
!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,v_half,dummylist)
#endif

#if defined(PREDICTOR_CORRECTOR)
!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,a_old,dummylist)
#endif

#if defined(GRAVITY)
!$OMP SECTION
  call reorder_array_real(pstart,plast,gpot,dummylist)

#if !defined(GEOMETRIC_MAC)
!$OMP SECTION
  call reorder_array_real(pstart,plast,agravmag,dummylist)
#endif

#if defined(RAD_WS)
!$OMP SECTION
  call reorder_array_real(pstart,plast,sphgpot,dummylist)
#endif
#endif

#if defined(DEBUG_FORCES) && defined(GRAVITY)
!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,a_grav,dummylist)
#endif

#if defined(DEBUG_FORCES) && defined(HYDRO)
!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,a_hydro,dummylist)
#endif

#if defined(DEBUG_FORCES) && defined(IDEAL_MHD)
!$OMP SECTION
  call reorder_array_real_2D(VDIM,pstart,plast,a_mag,dummylist)
#endif


! Internal energy arrays
! ----------------------------------------------------------------------------
#if defined(INTERNAL_ENERGY)
!$OMP SECTION
  call reorder_array_real(pstart,plast,u,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,du_dt,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,u_old,dummylist)

#if defined(DIFFUSION)
!$OMP SECTION
  call reorder_array_real(pstart,plast,du_dt_diff,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,k_cond,dummylist)
#endif
#endif


! Entropic function arrays
! ----------------------------------------------------------------------------
#if defined(ENTROPIC_FUNCTION)
!$OMP SECTION
  call reorder_array_real(pstart,plast,Aent,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,dA_dt,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,Aold,dummylist)
#endif


! Hydrodynamical arrays
! ----------------------------------------------------------------------------
!$OMP SECTION
  call reorder_array_real(pstart,plast,rho,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,rho_old,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,drhodt,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,temp,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,press,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,sound,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,div_v,dummylist)

#if defined(GRAD_H_SPH)
!$OMP SECTION
  call reorder_array_real(pstart,plast,omega,dummylist)
#endif

#if defined(VISC_TD)
!$OMP SECTION
  call reorder_array_real(pstart,plast,talpha,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,talpha_old,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,dalpha_dt,dummylist)
#endif

#if defined(VISC_BALSARA)
!$OMP SECTION
  call reorder_array_real(pstart,plast,balsara,dummylist)
#endif

#if defined(VISC_PATTERN_REC)
!$OMP SECTION
  call reorder_array_real(pstart,plast,pattrec,dummylist)
#endif

#if defined(IONIZING_UV_RADIATION)
!$OMP SECTION
  call reorder_array_logical(pstart,plast,ionizedo,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(NDIM,pstart,plast,gradrho,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,temp_min,dummylist)

!$OMP SECTION
  call reorder_array_int(pstart,plast,newtemp,dummylist)
#endif

#if defined(SINKS) && defined(GRAVITY)
!$OMP SECTION
  call reorder_array_logical(pstart,plast,ispotmin,dummylist)
#endif


! MHD arrays
! ----------------------------------------------------------------------------
#if defined(IDEAL_MHD)
!$OMP SECTION
  call reorder_array_real_2D(BDIM,pfirst,plast,B,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(BDIM,pfirst,plast,B_old,dummylist)

!$OMP SECTION
  call reorder_array_real_2D(BDIM,pfirst,plast,dB_dt,dummylist)

#if defined(DEBUG_MHD)
!$OMP SECTION
  call reorder_array_real_2D(BDIM,pfirst,plast,div_B,dummylist)
#endif
#endif


! Timestepping arrays
! ----------------------------------------------------------------------------
!$OMP SECTION 
  call reorder_array_logical(pstart,plast,accdo,dummylist)

!$OMP SECTION 
  call reorder_array_long_int(pstart,plast,nlevel,dummylist)

!$OMP SECTION
  call reorder_array_long_int(pstart,plast,nlast,dummylist)

!$OMP SECTION
  call reorder_array_double(pstart,plast,laststep,dummylist)

#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
!$OMP SECTION
  call reorder_array_long_int(pstart,plast,nminneib,dummylist)
#endif


! Radiative cooling arrays
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
!$OMP SECTION
  call reorder_array_int(pstart,plast,idens,dummylist)

!$OMP SECTION
  call reorder_array_int(pstart,plast,itemp,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,column2,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,dt_therm,dummylist)

!$OMP SECTION
  call reorder_array_real(pstart,plast,ueq,dummylist)

!$OMP SECTION
#if defined(DEBUG_DUDTRAD)
  call reorder_array_real(pstart,plast,dudt_rad,dummylist)
#endif

#endif

!$OMP END SECTIONS
!$OMP END PARALLEL
! ============================================================================

  return
END SUBROUTINE reorder_particle_arrays




! ============================================================================
SUBROUTINE reorder_array_int_2D(nsize,pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: nsize                   ! No. of elements per part
  integer, intent(in)    :: pstart                  ! id of first particle
  integer, intent(in)    :: ptot                    ! No. of particles
  integer, intent(inout) :: iarray(1:nsize,1:ptot)  ! Array to be reordered
  integer, intent(in)    :: aorder(1:ptot)          ! New order of array

  integer :: k                                      ! Element counter
  integer :: p                                      ! Particle counter
  integer :: pold                                   ! Old p
  integer, allocatable :: itemp(:,:)                ! Aux storage array

  debug2("[reorder_array_int_2D.F90]")

  allocate(itemp(1:nsize,1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
    do k=1,nsize
      itemp(k,p) = iarray(k,p)
    end do
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
    pold = aorder(p)
    do k=1,nsize
       iarray(k,p) = itemp(k,pold)
    end do
  end do

  deallocate(itemp)

  return
END SUBROUTINE reorder_array_int_2D



! ============================================================================
SUBROUTINE reorder_array_int(pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart            ! id of first particle
  integer, intent(in)    :: ptot              ! No. of particles
  integer, intent(inout) :: iarray(1:ptot)    ! Array to be reordered
  integer, intent(in)    :: aorder(1:ptot)    ! New array order

  integer :: p                                ! Particle counter
  integer :: pold                             ! old id
  integer, allocatable :: itemp(:)            ! Aux. storage array

  debug2("[reorder_array_int.F90]")

  allocate(itemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     itemp(p) = iarray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     iarray(p) = itemp(pold)
  end do

  deallocate(itemp)

  return
END SUBROUTINE reorder_array_int



! ============================================================================
SUBROUTINE reorder_array_long_int(pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart                    ! id of first particle
  integer, intent(in)    :: ptot                      ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)            ! New array order
  integer(kind=ILP), intent(inout) :: iarray(1:ptot)  ! Array to be reordered

  integer :: p                                        ! Particle counter
  integer :: pold                                     ! old id
  integer(kind=ILP), allocatable :: itemp(:)          ! Aux. storage array

  debug2("[reorder_array_long_int.F90]")

  allocate(itemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     itemp(p) = iarray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     iarray(p) = itemp(pold)
  end do

  deallocate(itemp)

  return
END SUBROUTINE reorder_array_long_int



! ============================================================================
SUBROUTINE reorder_array_real_2D(nsize,pstart,ptot,rmain,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: nsize                       ! # elements in array
  integer, intent(in)    :: pstart                      ! id of first particle
  integer, intent(in)    :: ptot                        ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)              ! New array order
  real(kind=PR), intent(inout) :: rmain(1:nsize,1:ptot) ! Array to reorder

  integer :: k                                          ! Dimension counter
  integer :: p                                          ! Particle counter
  integer :: pold                                       ! Old id
  real(kind=PR), allocatable :: rtemp(:,:)              ! Aux. storage array

  debug2("[reorder_array_real_2D.F90]")

  allocate(rtemp(1:nsize,1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
    do k=1,nsize
      rtemp(k,p) = rmain(k,p)
    end do
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
    pold = aorder(p)
    do k=1,nsize
       rmain(k,p) = rtemp(k,pold)
    end do
  end do

  deallocate(rtemp)

  return
END SUBROUTINE reorder_array_real_2D



! ============================================================================
SUBROUTINE reorder_array_real(pstart,ptot,rmain,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart               ! id of first particle
  integer, intent(in)    :: ptot                 ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)       ! New array order
  real(kind=PR), intent(inout) :: rmain(1:ptot)  ! Array to reorder

  integer :: p                                   ! Particle counter
  integer :: pold                                ! Old id
  real(kind=PR), allocatable :: rtemp(:)         ! Aux. storage array

  debug2("[reorder_array_real.F90]")

  allocate(rtemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     rtemp(p) = rmain(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     rmain(p) = rtemp(pold)
  end do

  deallocate(rtemp)

  return
END SUBROUTINE reorder_array_real



! ============================================================================
SUBROUTINE reorder_array_double(pstart,ptot,rmain,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart               ! id of first particle
  integer, intent(in)    :: ptot                 ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)       ! New array order
  real(kind=DP), intent(inout) :: rmain(1:ptot)  ! Array to reorder

  integer :: p                                   ! Particle counter
  integer :: pold                                ! Old id
  real(kind=DP), allocatable :: rtemp(:)         ! Aux. storage array

  debug2("[reorder_array_double.F90]")

  allocate(rtemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     rtemp(p) = rmain(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     rmain(p) = rtemp(pold)
  end do

  deallocate(rtemp)

  return
END SUBROUTINE reorder_array_double



! ============================================================================
SUBROUTINE reorder_array_logical(pstart,ptot,larray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart          ! id of first particle
  integer, intent(in)    :: ptot            ! No. of particles
  logical, intent(inout) :: larray(1:ptot)  ! Array to reorder
  integer, intent(in)    :: aorder(1:ptot)  ! New array order

  integer :: p                              ! Particle counter
  integer :: pold                           ! Old id
  logical, allocatable :: ltemp(:)          ! Aux. storage array

  debug2("[reorder_array_logical.F90]")

  allocate(ltemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     ltemp(p) = larray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     larray(p) = ltemp(pold)
  end do

  deallocate(ltemp)

  return
END SUBROUTINE reorder_array_logical



! ============================================================================
SUBROUTINE reorder_inverse_array_int(pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart            ! id of first particle
  integer, intent(in)    :: ptot              ! No. of particles
  integer, intent(inout) :: iarray(1:ptot)    ! Array to be reordered
  integer, intent(in)    :: aorder(1:ptot)    ! New array order

  integer :: p                                ! Particle counter
  integer :: pold                             ! old id
  integer, allocatable :: itemp(:)            ! Aux. storage array

  debug2("[reorder_inverse_array_int.F90]")

  allocate(itemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     itemp(p) = iarray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     iarray(p) = itemp(pold)
  end do

  deallocate(itemp)

  return
END SUBROUTINE reorder_inverse_array_int
