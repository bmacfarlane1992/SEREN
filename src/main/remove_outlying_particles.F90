! REMOVE_OUTLYING_PARTICLES
! D. A. Hubber - 8/11/2009
! Remove any outlying particles that have escaped the system depending on 
! various criteria selected in the Makefile.  Current options are : 
! 1. - Energy of particle relative to COM is positive (i.e. unbound)
! 2. - Density of particle is below some minimum threshold
! 3. - Particle has passed through the 'removal' sphere (i.e. r > rremove)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE remove_outlying_particles
  use particle_module
  use hydro_module
  use type_module
  implicit none

  logical :: remove                     ! Remove particle or not?
  integer :: ndead                      ! No. of 'dead' (removed) particles
  integer :: p                          ! Particle counter
  integer, allocatable :: deadlist(:)   ! List of dead particles
  real(kind=DP) :: dr(1:NDIM)           ! Relative position
  real(kind=DP) :: drsqd                ! Distance squared
  real(kind=DP) :: dv(1:VDIM)           ! Relative velocities
  real(kind=DP) :: gpe                  ! Grav. potential energy
  real(kind=DP) :: ke                   ! Kinetic energy
  real(kind=DP) :: mdead                ! Mass of dead particles
  real(kind=DP) :: mp                   ! Mass of particle p
  real(kind=DP) :: momdead(1:NDIM)      ! Momentum of dead particles
  real(kind=DP) :: angmomdead(1:NDIM)   ! Ang. mom of dead particles
  real(kind=DP) :: rcomdead(1:NDIM)     ! Centre of mass of dead particles
  real(kind=DP) :: r0(1:NDIM)           ! Position of origin

  debug2("Remove inconsequential particles [remove_outlying_particles.F90]")

  allocate(deadlist(1:ptot))

  ndead = 0
  mdead = 0.0_DP
  r0(1:NDIM)         = 0.0_DP
  rcomdead(1:NDIM)   = 0.0_DP
  momdead(1:NDIM)    = 0.0_DP
  angmomdead(1:NDIM) = 0.0_DP


! Loop over all particles and find outlying particles
! ----------------------------------------------------------------------------
  do p=1,ptot

     ! Calculate any relevant variables here needed for removal conditions
     remove = .true.
     mp = real(parray(MASS,p))
     call distance2_dp(r0(1:NDIM),p,dr(1:NDIM),drsqd)
     dv(1:NDIM) = real(v(1:NDIM,p),DP) - vcom(1:NDIM)
     ke = 0.5_DP*mp*dot_product(dv(1:NDIM),dv(1:NDIM))
     gpe = mp*real(gpot(p),DP)

     ! Process removal criteria here
     if (rho_remove .and. rho(p) >= rholost) remove = .false.
     if (energy_remove .and. ke - gpe < 0.0_DP) remove = .false.
     if (rad_remove .and. drsqd <= real(rad_lost*rad_lost,DP)) remove = .false.
     if (.not. remove) cycle

     ndead = ndead + 1
     deadlist(ndead) = p

     mdead = mdead + mp
     momdead(1:NDIM) = momdead(1:NDIM) + mp*real(v(1:NDIM,p),DP)
     rcomdead(1:NDIM) = rcomdead(1:NDIM) + mp*real(parray(1:NDIM,p),DP)

#if defined(DEBUG_REMOVE_OUTLIERS)
     write(6,*) "REMOVING OUTLYING PARTICLE : ",p,rho(p)/rholost,(ke-gpe)/gpe
#endif

  end do
! ----------------------------------------------------------------------------


! If particles have been removed, need to reorder arrays
  if (ndead > 0) then
     call remove_from_list(ndead,deadlist(1:ndead))
     rlost(1:NDIM)   = mlost*rlost(1:NDIM) + rcomdead(1:NDIM)
     mlost           = mlost + mdead
     momlost(1:NDIM) = momlost(1:NDIM) + momdead(1:NDIM)
     rlost(1:NDIM)   = rlost(1:NDIM)/mlost
     angmomlost(1:3) = 0.0_DP
  end if

  deallocate(deadlist)


  return
END SUBROUTINE remove_outlying_particles
