! DIAGNOSTICS.F90
! C. P. Batty & D. A. Hubber - 19/1/2007
! Computes various diagnostic quantities during simulation, i.e.
! position and velocity of the centre of mass, total linear and angular
! momentum,  kinetic, gravitational, thermal and total energies, net
! hydrodynamical and gravitational forces on all particles
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE diagnostics
  use seren_sim_module
  use particle_module
  use scaling_module
  use type_module
  use hydro_module
#if defined(SINKS)
  use sink_module
  use Nbody_module
#endif
#if defined(DEBUG_DIAGNOSTICS)
  use time_module, only : nsteps,time
  use filename_module, only : run_dir,run_id
#endif
  implicit none

  integer       :: p                   ! particle counter
  real(kind=DP) :: agravp(1:VDIM)      ! grav. acceleration of particle p
  real(kind=DP) :: ap(1:VDIM)          ! acceleration of particle p
  real(kind=DP) :: angmom(1:3)         ! angular momentum
  real(kind=DP) :: Atot                ! entropic function summation
  real(kind=DP) :: etot                ! total energy
  real(kind=DP) :: force(1:VDIM)       ! net force on all particles
  real(kind=DP) :: gpetot              ! total grav. potential energy
  real(kind=DP) :: ketot               ! total kinetic energy
  real(kind=DP) :: mom(1:VDIM)         ! linear momentum
  real(kind=DP) :: mp                  ! mass of particle p
  real(kind=DP) :: rp(1:NDIM)          ! position of particle p
  real(kind=DP) :: utot                ! total internal energy
  real(kind=DP) :: vp(1:VDIM)          ! velocity of particle p
#if defined(IDEAL_MHD)
  real(kind=DP) :: betot               ! total magnetic energy
#endif
#if defined(SINKS)
  integer       :: s                   ! sink counter
  real(kind=DP) :: angmomp(1:3)        ! angular momentum of p
  real(kind=DP) :: gpep                ! gpe of particle p
#endif
#if defined(DEBUG_DIAGNOSTICS)
  logical :: ex                        ! Does file exist already?
  character(len=40) :: out_file        ! Name of outputted debug file
  integer :: itemp                     ! ..
  real(kind=DP) :: ttemp               ! ..
#endif
#if defined(DEBUG_FORCES)
  real(kind=DP) :: force_grav(1:VDIM)  ! net grav. force
  real(kind=DP) :: force_hydro(1:VDIM) ! net hydro force
  real(kind=DP) :: force_mag(1:VDIM)   ! net magnetic force
#endif

  debug2("Calculating total energy and momentum [diagnostics.F90]")

! Zero all accumulation variables
! ----------------------------------------------------------------------------
  mtot   = 0.0_DP
  etot   = 0.0_DP
  ketot  = 0.0_DP
  gpetot = 0.0_DP
  utot   = 0.0_DP
  Atot   = 0.0_DP
  force(1:VDIM) = 0.0_DP
  mom(1:VDIM)   = 0.0_DP
  rcom(1:NDIM)  = 0.0_DP
  vcom(1:VDIM)  = 0.0_DP
  angmom(1:3)   = 0.0_DP
#if defined(IDEAL_MHD)
  betot  = 0.0_DP
#endif
#if defined(DEBUG_FORCES)
  force_grav(1:VDIM)  = 0.0_DP
  force_hydro(1:VDIM) = 0.0_DP
  force_mag(1:VDIM)   = 0.0_DP
#endif


! Sum up the contributions from all SPH particles
! ============================================================================
  if (ptot > 0) then

     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ap,mp,rp,vp) &
     !$OMP REDUCTION(+:angmom,force,ketot,mom,mtot,rcom,vcom)
     do p=1,ptot

        ! Make local copies of r, m, v and a
        rp(1:NDIM) = real(parray(1:NDIM,p),DP)
        mp         = real(parray(MASS,p),DP)
        vp(1:VDIM) = real(v(1:VDIM,p),DP)
        ap(1:VDIM) = real(a(1:VDIM,p),DP)

        ! Sum up the contributions to net force, momentum, energy and mass
        mtot = mtot + mp
        mom(1:VDIM) = mom(1:VDIM) + mp*vp(1:VDIM)
        force(1:VDIM) = force(1:VDIM) + mp*ap(1:VDIM)
        rcom(1:NDIM) = rcom(1:NDIM) + mp*rp(1:NDIM)
        vcom(1:VDIM) = vcom(1:VDIM) + mp*vp(1:VDIM)
#if NDIM==3
        angmom(1) = angmom(1) + mp*(rp(2)*vp(3) - rp(3)*vp(2))
        angmom(2) = angmom(2) + mp*(rp(3)*vp(1) - rp(1)*vp(3))
        angmom(3) = angmom(3) + mp*(rp(1)*vp(2) - rp(2)*vp(1))
#elif NDIM==2
        angmom(3) = angmom(3) + mp*(rp(1)*vp(2) - rp(2)*vp(1))
#endif
        ketot = ketot + mp*dot_product(vp(1:VDIM),vp(1:VDIM))
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

#if defined(GRAVITY)
     do p=pgravitystart,pgravityend
        gpetot = gpetot + real(parray(MASS,p)*gpot(p),DP)
     end do
#endif
#if defined(HYDRO)
     do p=phydrostart,phydroend
#if defined(BAROTROPIC) || defined(ISOTHERMAL) || defined(POLYTROPIC)
        utot = utot + real(parray(MASS,p)*press(p)/rho(p),DP)
#elif defined(INTERNAL_ENERGY)
        utot = utot + real(parray(MASS,p)*u(p),DP)
#endif
        Atot = Atot + real(parray(MASS,p)*press(p)/(rho(p)**gamma),DP)
     end do
#endif
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     do p=pgravitystart,pgravityend
        force_grav(1:VDIM) = force_grav(1:VDIM) + &
             &real(parray(MASS,p),DP)*real(a_grav(1:VDIM,p),DP)
     end do
#endif
#if defined(DEBUG_FORCES) && defined(HYDRO)
     do p=phydrostart,phydroend
        force_hydro(1:VDIM) = force_hydro(1:VDIM) + &
             &real(parray(MASS,p),DP)*real(a_hydro(1:VDIM,p),DP)
     end do
#endif
#if defined(DEBUG_FORCES) && defined(IDEAL_MHD)
     do p=phydrostart,phydroend
        force_mag(1:VDIM) = force_mag(1:VDIM) + &
             &real(parray(MASS,p),DP)*real(a_mag(1:VDIM,p),DP)
     end do
#endif
  end if
! ============================================================================


! Sum up the contributions from all sinks/stars
! ============================================================================
#if defined(SINKS)
  do s=1,stot

     if (nbody_sph_sim .or. nbody_sim) then
        mp = star(s)%m
        rp(1:NDIM) = star(s)%r(1:NDIM)
        vp(1:VDIM) = star(s)%v(1:VDIM)
        ap(1:VDIM) = star(s)%a(1:VDIM)
        agravp(1:NDIM) = star(s)%a(1:NDIM)
        angmomp(1:3) = star(s)%angmom(1:3)
        gpep = star(s)%gpe
     else
        mp = real(sink(s)%m,DP)
        rp(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
        vp(1:VDIM) = real(sink(s)%v(1:VDIM),DP)
        ap(1:VDIM) = real(sink(s)%a(1:VDIM),DP)
        agravp(1:NDIM) = real(sink(s)%agrav(1:NDIM),DP)
        angmomp(1:3) = real(sink(s)%angmom(1:3),DP)
        gpep = real(sink(s)%gpe,DP)
#if defined(DEBUG_FORCES)
        force_hydro(1:VDIM) = &
             &force_hydro(1:VDIM) + mp*real(sink(s)%ahydro(1:VDIM),DP)
#endif
     end if
     
     ! Sum up the contributions to net force, momentum, energy and mass
     mtot = mtot + mp
     mom(1:VDIM)   = mom(1:VDIM) + mp*vp(1:VDIM)
     force(1:VDIM) = force(1:VDIM) + mp*ap(1:VDIM)
     angmom(1:3)   = angmom(1:3)   + angmomp(1:3)
     rcom(1:NDIM) = rcom(1:NDIM) + mp*rp(1:NDIM)
     vcom(1:VDIM) = vcom(1:VDIM) + mp*vp(1:VDIM)
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     force_grav(1:NDIM) = force_grav(1:NDIM) + mp*agravp(1:NDIM)
#endif
     
     ! Now add orbital angular momentum of sink COM
#if NDIM==3
     angmom(1) = angmom(1) + mp*(rp(2)*vp(3) - rp(3)*vp(2))
     angmom(2) = angmom(2) + mp*(rp(3)*vp(1) - rp(1)*vp(3))
     angmom(3) = angmom(3) + mp*(rp(1)*vp(2) - rp(2)*vp(1))
#elif NDIM==2
     angmom(3) = angmom(3) + mp*(rp(1)*vp(2) - rp(2)*vp(1))
#endif
     ketot  = ketot + mp*dot_product(vp(1:VDIM),vp(1:VDIM))
#if defined(GRAVITY)
     gpetot = gpetot + gpep
#endif
  end do
#endif
! ============================================================================


! Normalise centre of mass variables
  rcom(1:NDIM) = rcom(1:NDIM) / mtot
  vcom(1:VDIM) = vcom(1:VDIM) / mtot


! Account for constant multipliers, double summation (gravity),
! degrees of freedom (thermal) etc.
! ----------------------------------------------------------------------------
  ketot  = 0.5_DP*ketot
  etot   = ketot
#if defined(GRAVITY)
  gpetot = 0.5_DP*gpetot
  etot   = etot - gpetot
#endif
#if defined(HYDRO)
#if defined(BAROTROPIC) || defined(POLYTROPIC)
  utot = real(1.0_PR/(gamma - 1.0_PR),DP)*utot
  etot = etot + utot
#elif defined(ISOTHERMAL)
  utot = 1.5_DP*utot
  etot = etot + utot
#elif defined(INTERNAL_ENERGY)
  etot = etot + utot
#endif
#endif
#if defined(IDEAL_MHD)
  etot = etot + betot
#endif


! Output to screen
! ----------------------------------------------------------------------------
5 format(1X,A27,1I10)
10 format(1X,A8,2X,1G18.10)
15 format(1X,A8,3X,1G18.10,3X,1G18.10,3X,1G18.10)
#if NDIM==1
20 format(1X,A8,3X,1G18.10)
#elif NDIM==2
20 format(1X,A8,3X,1G18.10,3X,1G18.10)
#elif NDIM==3
20 format(1X,A8,3X,1G18.10,3X,1G18.10,3X,1G18.10)
#endif

  write(6,5) "Total number of particles : ", ptot
#if defined(SINKS)
  write(6,5) "Total number of sinks     : ", stot
#endif
  write(6,10) "mtot   :", mtot*mscale
  write(6,20) "rcom   :", rcom(1:NDIM)*rscale
  write(6,20) "vcom   :", vcom(1:VDIM)*vscale
  write(6,20) "mom    :", mom(1:VDIM)*momscale
  write(6,15) "ang    :", angmom(1:3)*angmomscale
  write(6,20) "force  :", force(1:VDIM)
#if defined(DEBUG_FORCES)
#if defined(GRAVITY)
  write(6,20) "forceG :", force_grav(1:VDIM)
#endif
#if defined(HYDRO)
  write(6,20) "forceH :", force_hydro(1:VDIM)
#endif
#if defined(IDEAL_MHD)
  write(6,20) "forceB :", force_mag(1:VDIM)
#endif
#endif
  write(6,10) "etot   :", (ketot + utot - gpetot)*Escale
  write(6,10) "ketot  :", ketot*Escale
#if defined(HYDRO)
  write(6,10) "utot   :", utot*Escale
#endif
#if defined(GRAVITY)
  write(6,10) "gpetot :", -gpetot*Escale
#endif
#if defined(IDEAL_MHD)
  write(6,10) "betot :", betot*Escale
#endif


! Write to file
! ----------------------------------------------------------------------------
#if defined(DEBUG_DIAGNOSTICS)
  out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".diag"

  ! Check if file exists
  inquire(file=out_file,exist=ex)

  ! Read through file to synchronise files with time of simulation.
  if (ex) then
     open(1,file=out_file,status="unknown",&
          &form="formatted",position="rewind")
     do
        read(1,'(I8,E18.10)',end=50,err=50) itemp,ttemp
        if (ttemp > time*tscale) exit
     end do
50   backspace (1,err=100)
  else
     open(1,file=out_file,status="unknown",form="formatted")
  end if
#if defined(DEBUG_DIAGNOSTICS) && defined(DEBUG_FORCES)
100 write(1,'(I8,60E18.10)') nsteps,time*tscale,mtot*mscale,&
         & rcom(1:NDIM)*rscale,vcom(1:VDIM)*vscale,mom(1:VDIM)*momscale,&
         & angmom(1:3)*angmomscale,force(1:VDIM),force_grav(1:VDIM),&
         & force_hydro(1:VDIM),force_mag(1:VDIM),utot*Escale,gpetot*Escale,&
         & ketot*Escale,etot*Escale,Atot
  close(1)
#else
100 write(1,'(I8,60E18.10)') nsteps,time*tscale,mtot*mscale,&
         & rcom(1:NDIM)*rscale,vcom(1:VDIM)*vscale,mom(1:VDIM)*momscale,&
         & angmom(1:3)*angmomscale,force(1:VDIM),utot*Escale,gpetot*Escale,&
         & ketot*Escale,etot*Escale,Atot
  close(1)
#endif
#endif

  return
END SUBROUTINE diagnostics
