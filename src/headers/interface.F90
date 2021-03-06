! INTERFACE.F90
! D. A. Hubber - 6/9/2010
! Interfaces for all subroutines with an argument list.  All subroutines 
! without arguments are excluded (but perhaps should be included in the 
! future for safety).
! ============================================================================

#include "macros.h"

! ============================================================================
MODULE interface_module
  use definitions

  INTERFACE

     SUBROUTINE add_external_gravitational_force(rp,atemp,pottemp)
       use definitions
       real(kind=PR), intent(in) :: rp(1:NDIM)
       real(kind=DP), intent(out) :: atemp(1:NDIM)
       real(kind=DP), intent(out) :: pottemp
     END SUBROUTINE add_external_gravitational_force

     SUBROUTINE advance_boundary_particle(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE advance_boundary_particle

     SUBROUTINE advance_euler(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE advance_euler

     SUBROUTINE advance_leapfrog_dkd(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE advance_leapfrog_dkd

     SUBROUTINE advance_leapfrog_kdk(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE advance_leapfrog_kdk

     SUBROUTINE advance_predictor_corrector(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE advance_predictor_corrector

     SUBROUTINE advance_runge_kutta(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE advance_runge_kutta

     SUBROUTINE all_sph(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE all_sph

     SUBROUTINE all_sph_gradh(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE all_sph_gradh

     SUBROUTINE ambient_temp(p,atemp)
       !DEC$ ATTRIBUTES FORCEINLINE :: ambient_temp
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(out) :: atemp
     END SUBROUTINE ambient_temp

     SUBROUTINE ang2pix_nest(nside, theta, phi, ipix, x2pix, y2pix)
       use healpix_types
       INTEGER(KIND=I4B), INTENT(IN) :: nside
       INTEGER(KIND=I4B), INTENT(OUT) :: ipix
       INTEGER(KIND=I4B), INTENT(IN), DIMENSION(128) :: x2pix, y2pix
       REAL(KIND=DP), INTENT(IN) ::  theta, phi
     END SUBROUTINE ang2pix_nest

     SUBROUTINE BHgrav_accel(p,invhp,rp,agravp,potp)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(in) :: invhp
       real(kind=PR), intent(in) :: rp(1:NDIM)
       real(kind=DP), intent(out) :: agravp(1:NDIM)
       real(kind=DP), intent(out) :: potp
     END SUBROUTINE BHgrav_accel

     SUBROUTINE BHgrav_accel_jerk(s,invhs,rs,vs,agravp,adotp,potp)
       use definitions
       integer, intent(in) :: s
       real(kind=DP), intent(in) :: invhs
       real(kind=DP), intent(in) :: rs(1:NDIM)
       real(kind=DP), intent(in) :: vs(1:NDIM)
       real(kind=DP), intent(out) :: agravp(1:NDIM)
       real(kind=DP), intent(out) :: adotp(1:NDIM)
       real(kind=DP), intent(out) :: potp
     END SUBROUTINE BHgrav_accel_jerk

     SUBROUTINE BHhydro_walk(rp,hrange,pp_pot,pp_max,pp_list)
       use definitions
       integer, intent(in)  :: pp_max
       integer, intent(out) :: pp_pot
       integer, intent(out) :: pp_list(1:pp_max)
       real(kind=PR),intent(in) :: hrange
       real(kind=PR), intent(in) :: rp(1:NDIM)
     END SUBROUTINE BHhydro_walk

     SUBROUTINE BHhydrowalk_hgather(rp,hrange,pp_pot,pp_max,pp_list)
       use definitions
       integer, intent(in)  :: pp_max
       integer, intent(out) :: pp_pot
       integer, intent(out) :: pp_list(1:pp_max)
       real(kind=PR),intent(in) :: hrange
       real(kind=PR),intent(in) :: rp(1:NDIM)
     END SUBROUTINE BHhydrowalk_hgather

     SUBROUTINE BH_remove_particles(newid)
       use definitions
       use particle_module, only : ptot
       integer, intent(in) :: newid(1:ptot)
     END SUBROUTINE BH_remove_particles

     SUBROUTINE binary_energy(s1,s2,h1,h2,m1,m2,r1,r2,v1,v2,binen)
       use definitions
       integer,intent(in) :: s1
       integer,intent(in) :: s2
       real(kind=DP),intent(in) :: h1
       real(kind=DP),intent(in) :: h2
       real(kind=DP),intent(in) :: m1
       real(kind=DP),intent(in) :: m2
       real(kind=DP),intent(in) :: r1(1:NDIM)
       real(kind=DP),intent(in) :: r2(1:NDIM)
       real(kind=DP),intent(in) :: v1(1:NDIM)
       real(kind=DP),intent(in) :: v2(1:NDIM)
       real(kind=DP),intent(out) :: binen
     END SUBROUTINE binary_energy

     SUBROUTINE binary_gather_walk(rp,hp,pp_pot,pp_max,pp_list)
       use definitions
       integer, intent(in)  :: pp_max
       integer, intent(out) :: pp_pot
       integer, intent(out) :: pp_list(1:pp_max)
       real(kind=PR), intent(in) :: hp
       real(kind=PR), intent(in) :: rp(1:NDIM)
     END SUBROUTINE binary_gather_walk

     SUBROUTINE binary_neibfind(rp,hp,pp_pot,pp_max,pp_list)
       use definitions
       integer, intent(in)  :: pp_max
       integer, intent(out) :: pp_pot
       integer, intent(out) :: pp_list(1:pp_max)
       real(kind=PR), intent(in) :: hp
       real(kind=PR), intent(in) :: rp(1:NDIM)
     END SUBROUTINE binary_neibfind

     SUBROUTINE binary_properties(s1,s2,binen,m1,m2,r1,r2,v1,v2)
       use definitions
       integer, intent(in) :: s1
       integer, intent(in) :: s2
       real(kind=DP),intent(in) :: binen
       real(kind=DP),intent(in) :: m1
       real(kind=DP),intent(in) :: m2
       real(kind=DP),intent(in) :: r1(1:NDIM)
       real(kind=DP),intent(in) :: r2(1:NDIM)
       real(kind=DP),intent(in) :: v1(1:NDIM)
       real(kind=DP),intent(in) :: v2(1:NDIM)
     END SUBROUTINE binary_properties

     SUBROUTINE bounding_box(pstart,pend,rmax,rmin)
       !DEC$ ATTRIBUTES FORCEINLINE :: bounding_box
       use definitions
       integer, intent(in) :: pstart
       integer, intent(in) :: pend
       real(kind=PR), intent(out) :: rmax(1:NDIM)
       real(kind=PR), intent(out) :: rmin(1:NDIM)
     END SUBROUTINE bounding_box

     SUBROUTINE check_boundary_conditions(rp,vp)
       use definitions
       real(kind=PR), intent(inout) :: rp(1:NDIM)
       real(kind=PR), intent(inout) :: vp(1:VDIM)
     END SUBROUTINE check_boundary_conditions

     SUBROUTINE comperror(errmsg)
       use definitions
       character(len=*), intent(in) :: errmsg
     END SUBROUTINE comperror

     SUBROUTINE conductivity(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE conductivity

     SUBROUTINE create_HP_source(ssink,rsource)
       use definitions
       integer, intent(in) :: ssink
       real(kind=PR), intent(in) :: rsource(1:NDIM)
     END SUBROUTINE create_HP_source

     SUBROUTINE create_sink(psink)
       use definitions
       integer, intent(in) :: psink
     END SUBROUTINE create_sink

     SUBROUTINE diffusion(p,dt)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(in) :: dt
     END SUBROUTINE diffusion

     SUBROUTINE direct_sink_gravity(p,invhp,rp,agravp,potp)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(in) :: invhp
       real(kind=PR), intent(in) :: rp(1:NDIM)
       real(kind=DP), intent(out) :: agravp(1:NDIM)
       real(kind=DP), intent(out) :: potp
     END SUBROUTINE direct_sink_gravity

     SUBROUTINE direct_sph_gravity(p,invhp,rp,agravp,potp)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(in) :: invhp
       real(kind=PR), intent(in) :: rp(1:NDIM)
       real(kind=DP), intent(out) :: agravp(1:NDIM)
       real(kind=DP), intent(out) :: potp
     END SUBROUTINE direct_sph_gravity

     SUBROUTINE distance(p,pp,dr,drsqd)
       !DEC$ ATTRIBUTES FORCEINLINE :: distance
       use definitions
       integer, intent(in) :: p
       integer, intent(in) :: pp
       real(kind=PR), intent(out) :: dr(1:NDIM)
       real(kind=PR), intent(out) :: drsqd
     END SUBROUTINE distance

     SUBROUTINE distance2(rp,pp,dr,drsqd)
       !DEC$ ATTRIBUTES FORCEINLINE :: distance2
       use definitions
       integer, intent(in) :: pp
       real(kind=PR), intent(in)  :: rp(1:NDIM)
       real(kind=PR), intent(out) :: dr(1:NDIM)
       real(kind=PR), intent(out) :: drsqd
     END SUBROUTINE distance2

     SUBROUTINE distance2_dp(rp,pp,dr,drsqd)
       !DEC$ ATTRIBUTES FORCEINLINE :: distance2_dp
       use definitions
       integer, intent(in) :: pp
       real(kind=DP), intent(in)  :: rp(1:NDIM)
       real(kind=DP), intent(out) :: dr(1:NDIM)
       real(kind=DP), intent(out) :: drsqd
     END SUBROUTINE distance2_dp

     SUBROUTINE distance3(rp,rpp,dr,drsqd)
       !DEC$ ATTRIBUTES FORCEINLINE :: distance3
       use definitions
       real(kind=PR), intent(in)  :: rp(1:NDIM)
       real(kind=PR), intent(in)  :: rpp(1:NDIM)
       real(kind=PR), intent(out) :: dr(1:NDIM)
       real(kind=PR), intent(out) :: drsqd
     END SUBROUTINE distance3

     SUBROUTINE distance3_dp(rp,rpp,dr,drsqd)
       !DEC$ ATTRIBUTES FORCEINLINE :: distance3_dp
       use definitions
       real(kind=DP), intent(in)  :: rp(1:NDIM)
       real(kind=DP), intent(in)  :: rpp(1:NDIM)
       real(kind=DP), intent(out) :: dr(1:NDIM)
       real(kind=DP), intent(out) :: drsqd
     END SUBROUTINE distance3_dp

     SUBROUTINE ebalance(energybalance,du_dt,To,T,kappapT,kappaT,column2)
       use definitions
       real(kind=PR), intent(in)  ::du_dt,To,T,kappapT,kappaT,column2
       real(kind=PR), intent(out) ::energybalance
     END SUBROUTINE ebalance

     SUBROUTINE effective_gamma(p1,p2,gamma_eff)
       !DEC$ ATTRIBUTES FORCEINLINE :: effective_gamma
       use definitions
       integer, intent(in) :: p1
       integer, intent(in) :: p2
       real(kind=PR), intent(out) :: gamma_eff
     END SUBROUTINE effective_gamma

     SUBROUTINE eigenvalue_mac(qc,mac)
       !DEC$ ATTRIBUTES FORCEINLINE :: eigenvalue_mac
       use definitions
       real(kind=DP), intent(in) :: qc(1:NQUAD)
       real(kind=PR), intent(out) :: mac
     END SUBROUTINE eigenvalue_mac

     FUNCTION eosmu(dens,temp,idens,itemp)
       !DEC$ ATTRIBUTES FORCEINLINE :: eosmu
       use definitions
       real(kind=PR) :: eosmu
       real(kind=PR), intent(in)  :: dens,temp
       integer, intent(in)  :: idens,itemp
     END FUNCTION eosmu
     
     FUNCTION eosenergy(dens,temp,idens,itemp)
       !DEC$ ATTRIBUTES FORCEINLINE :: eosenergy
       use definitions
       real(kind=PR) :: eosenergy
       real(kind=PR), intent(in)  :: dens,temp
       integer, intent(in)  :: idens,itemp
     END FUNCTION eosenergy

     SUBROUTINE ewald_force(dr,m,eaccel)
       use definitions
       real(kind=PR), intent(in) :: dr(1:NDIM)
       real(kind=PR), intent(in) :: m
       real(kind=PR), intent(out):: eaccel(1:NDIM)
     END SUBROUTINE ewald_force

     SUBROUTINE fatal_error (msg)
       character(len=*), intent(in) :: msg
     END SUBROUTINE fatal_error

     SUBROUTINE find_equilibrium_temp_ws(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE find_equilibrium_temp_ws

     SUBROUTINE find_idens(dens,idens_index)
       !DEC$ ATTRIBUTES FORCEINLINE :: find_idens
       use definitions
       integer, intent(out)  :: idens_index
       real(kind=PR), intent(in)   :: dens
     END SUBROUTINE find_idens

     SUBROUTINE find_itemp(temp,itemp_index)
       !DEC$ ATTRIBUTES FORCEINLINE :: find_itemp
       use definitions
       integer, intent(out)  :: itemp_index
       real(kind=PR), intent(in)   :: temp
     END SUBROUTINE find_itemp

     SUBROUTINE find_temp_from_energy(idens,energy,itemp,temp)
       use definitions
       integer, intent(in)  :: idens
       real(kind=PR), intent(in)     :: energy
       integer, intent(out) :: itemp
       real(kind=PR), intent(out)    ::temp
     END SUBROUTINE find_temp_from_energy

     SUBROUTINE gather_neib(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE gather_neib

     SUBROUTINE gather_neib_on_fly(p,hp,pp_tot,pp_list)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(in) :: hp
       integer, intent(out) :: pp_tot
       integer, intent(out) :: pp_list(1:LISTSIZE)
     END SUBROUTINE gather_neib_on_fly

     SUBROUTINE getkappa(dens,temp,idens,kappa_mean,kappa_rosseland,kappa_planck)
       !DEC$ ATTRIBUTES FORCEINLINE :: getkappa
       use definitions
       integer, intent(in) :: idens
       real(kind=PR), intent(in)  :: temp
       real(kind=PR), intent(in) :: dens
       real(kind=PR), intent(out) :: kappa_mean
       real(kind=PR), intent(out) :: kappa_rosseland
       real(kind=PR), intent(out) :: kappa_planck
     END SUBROUTINE getkappa

     SUBROUTINE get_neib(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE get_neib

     SUBROUTINE get_neib_on_fly(p,hp,pp_tot,pp_totmax,pp_list)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(in) :: hp
       integer, intent(out) :: pp_tot
       integer, intent(in) :: pp_totmax
       integer, intent(out) :: pp_list(1:pp_totmax)
     END SUBROUTINE get_neib_on_fly

     SUBROUTINE gravity_gradh(invhp,hpp,mpp,rp,rpp,zo_p,zo_pp,atemp,dpotp)
       !DEC$ ATTRIBUTES FORCEINLINE :: gravity_gradh
       use definitions
       real(kind=PR), intent(in)  :: hpp
       real(kind=PR), intent(in)  :: invhp
       real(kind=PR), intent(in)  :: mpp
       real(kind=PR), intent(in)  :: rp(1:NDIM)
       real(kind=PR), intent(in)  :: rpp(1:NDIM)
       real(kind=PR), intent(out) :: atemp(1:NDIM)
       real(kind=PR), intent(out) :: dpotp
       real(kind=PR), intent(in)  :: zo_p
       real(kind=PR), intent(in)  :: zo_pp
     END SUBROUTINE gravity_gradh

     SUBROUTINE gravity_hermite4(invhp,hpp,mpp,rp,rpp,vp,vpp,&
          &atemp,adottemp,dpotp)
       !DEC$ ATTRIBUTES FORCEINLINE :: gravity_hermite4
       use definitions
       real(kind=DP), intent(in)  :: hpp
       real(kind=DP), intent(in)  :: invhp
       real(kind=DP), intent(in)  :: mpp
       real(kind=DP), intent(in)  :: rp(1:NDIM)
       real(kind=DP), intent(in)  :: rpp(1:NDIM)
       real(kind=DP), intent(in)  :: vp(1:NDIM)
       real(kind=DP), intent(in)  :: vpp(1:NDIM)
       real(kind=DP), intent(out) :: atemp(1:NDIM)
       real(kind=DP), intent(out) :: adottemp(1:NDIM)
       real(kind=DP), intent(out) :: dpotp
     END SUBROUTINE gravity_hermite4

     SUBROUTINE gravity_hermite4_meanh(hmean,mpp,rp,rpp,vp,vpp,&
          &atemp,adottemp,dpotp)
       !DEC$ ATTRIBUTES FORCEINLINE :: gravity_hermite4
       use definitions
       real(kind=DP), intent(in)  :: hmean
       real(kind=DP), intent(in)  :: mpp
       real(kind=DP), intent(in)  :: rp(1:NDIM)
       real(kind=DP), intent(in)  :: rpp(1:NDIM)
       real(kind=DP), intent(in)  :: vp(1:NDIM)
       real(kind=DP), intent(in)  :: vpp(1:NDIM)
       real(kind=DP), intent(out) :: atemp(1:NDIM)
       real(kind=DP), intent(out) :: adottemp(1:NDIM)
       real(kind=DP), intent(out) :: dpotp
     END SUBROUTINE gravity_hermite4_meanh

     SUBROUTINE gravity_meanh(hmean,mpp,rp,rpp,atemp,dpotp)
       !DEC$ ATTRIBUTES FORCEINLINE :: gravity_gradh
       use definitions
       real(kind=PR), intent(in)  :: hmean
       real(kind=PR), intent(in)  :: mpp
       real(kind=PR), intent(in)  :: rp(1:NDIM)
       real(kind=PR), intent(in)  :: rpp(1:NDIM)
       real(kind=PR), intent(out) :: atemp(1:NDIM)
       real(kind=PR), intent(out) :: dpotp
     END SUBROUTINE gravity_meanh

     SUBROUTINE gravity_nbody(mpp,rp,rpp,atemp,dpotp)
       !DEC$ ATTRIBUTES FORCEINLINE :: gravity_nbody
       use definitions
       real(kind=PR), intent(in)  :: mpp
       real(kind=PR), intent(in)  :: rp(1:NDIM)
       real(kind=PR), intent(in)  :: rpp(1:NDIM)
       real(kind=PR), intent(out) :: atemp(1:NDIM)
       real(kind=PR), intent(out) :: dpotp
     END SUBROUTINE gravity_nbody

     SUBROUTINE gravity_sph(invhp,hpp,mpp,rp,rpp,atemp,dpotp)
       !DEC$ ATTRIBUTES FORCEINLINE :: gravity_sph
       use definitions
       real(kind=PR), intent(in)  :: hpp
       real(kind=PR), intent(in)  :: invhp
       real(kind=PR), intent(in)  :: mpp
       real(kind=PR), intent(in)  :: rp(1:NDIM)
       real(kind=PR), intent(in)  :: rpp(1:NDIM)
       real(kind=PR), intent(out) :: atemp(1:NDIM)
       real(kind=PR), intent(out) :: dpotp
     END SUBROUTINE gravity_sph

     SUBROUTINE h_density_evaluation_point(rep,rmag_ep,rsqd_ep,rsource,pprev,&
          &p,plast,nextptcl,rhoep,hep,gradrhoep,rtree,htree,rhotree,&
          &gradrhotree,treeep)
       use definitions
       use particle_module, only : ptot
       logical, intent(inout) :: treeep
       integer, intent(in) :: nextptcl(1:ptot)
       integer, intent(in) :: p
       integer, intent(in) :: plast
       integer, intent(in) :: pprev
       real(kind=PR), intent(in) :: rep(1:NDIM)
       real(kind=PR), intent(in) :: rsource(1:NDIM)
       real(kind=PR), intent(in) :: rmag_ep
       real(kind=PR), intent(in) :: rsqd_ep
       real(kind=PR), intent(inout) :: hep
       real(kind=PR), intent(inout) :: rhoep
       real(kind=PR), intent(out) :: gradrhoep(1:NDIM)
       real(kind=PR), intent(inout) :: rtree(1:NDIM)
       real(kind=PR), intent(inout) :: htree
       real(kind=PR), intent(inout) :: rhotree
       real(kind=PR), intent(inout) :: gradrhotree(1:NDIM)
     END SUBROUTINE h_density_evaluation_point

     SUBROUTINE vec2ang(vector,theta,phi)
       !DEC$ ATTRIBUTES FORCEINLINE :: vec2ang
       use healpix_types
       real(kind=DP), intent(in), dimension(1:3) :: vector
       real(kind=DP), intent(out) :: theta
       real(kind=DP), intent(out) :: phi
     END SUBROUTINE vec2ang
 
     SUBROUTINE heapsort_real(nsort,rarray,iarray)
       !DEC$ ATTRIBUTES FORCEINLINE :: heapsort_real
       use definitions
       integer, intent(in)  :: nsort
       integer, intent(inout) :: iarray(1:nsort)
       real(kind=PR), intent(inout) :: rarray(1:nsort)
     END SUBROUTINE heapsort_real

     SUBROUTINE h_gather(p,hguess,rp)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(inout) :: hguess
       real(kind=PR), intent(in) :: rp(1:NDIM)
     END SUBROUTINE h_gather

     SUBROUTINE HP_calculate_basis_vector(level,ipix,rbasis,isource)
       !DEC$ ATTRIBUTES FORCEINLINE :: HP_calculate_basis_vector
       use definitions
       use healpix_types
       integer, intent(in) :: isource
       integer, intent(in) :: level
       integer(kind=I4B), intent(in) :: ipix
       real(kind=PR), intent(out) :: rbasis(1:NDIM)
     END SUBROUTINE HP_calculate_basis_vector

     SUBROUTINE HP_evaluation_point(rep,rmag_ep,rsqd_ep,rsource,pprev,&
          &p,plast,nextptcl,rhoep,hep,gradrhoep,rtree,htree,rhotree,&
          &gradrhotree,treeep)
       use definitions
       use particle_module, only : ptot
       logical, intent(inout) :: treeep
       integer, intent(in) :: nextptcl(1:ptot)
       integer, intent(in) :: p
       integer, intent(in) :: plast
       integer, intent(in) :: pprev
       real(kind=PR), intent(in) :: rep(1:NDIM)
       real(kind=PR), intent(in) :: rsource(1:NDIM)
       real(kind=PR), intent(in) :: rmag_ep
       real(kind=PR), intent(in) :: rsqd_ep
       real(kind=PR), intent(inout) :: hep
       real(kind=PR), intent(inout) :: rhoep
       real(kind=PR), intent(out) :: gradrhoep(1:NDIM)
       real(kind=PR), intent(inout) :: rtree(1:NDIM)
       real(kind=PR), intent(inout) :: htree
       real(kind=PR), intent(inout) :: rhotree
       real(kind=PR), intent(inout) :: gradrhotree(1:NDIM)
     END SUBROUTINE HP_evaluation_point

     SUBROUTINE HP_heapsort(qtot,Q,qorig,qq)
       use definitions
       integer, intent(in)  :: qtot
       real(kind=PR),intent(in)  :: Q(1:qtot)
       integer, intent(in)  :: qorig(1:qtot)
       integer, intent(out) :: qq(1:qtot)
     END SUBROUTINE HP_heapsort

     SUBROUTINE HP_initialize_source(isource,distids,nextptcl,&
          &nextptclaux,prevptcl,prevptclaux,raylist,distsqd)
       use definitions
       use healpix_types
       use particle_module, only : ptot
       integer, intent(in) :: isource
       integer, intent(inout) :: distids(1:ptot)
       integer, intent(inout) :: nextptcl(1:ptot)
       integer, intent(inout) :: nextptclaux(1:ptot)
       integer, intent(inout) :: prevptcl(1:ptot)
       integer, intent(inout) :: prevptclaux(1:ptot)
       integer(kind=I4B), intent(inout) :: raylist(1:ptot)
       real(kind=PR), intent(inout) :: distsqd(1:ptot)
     END SUBROUTINE HP_initialize_source

     SUBROUTINE HP_inverse_positions(rvec,isource)
       !DEC$ ATTRIBUTES FORCEINLINE :: HP_inverse_positions
       use definitions
       integer, intent(in) :: isource
       real(kind=DP), intent(inout) :: rvec(1:NDIM)
     END SUBROUTINE HP_inverse_positions

     SUBROUTINE HP_propagate_UV_radiation(i,isource,level,p,pprev,fh,hep,&
          &oldfh,oldstep,ray_integral,rep,rep_prev,rhoep,rsource,rsqd_ep,&
          &rtree,htree,rhotree,gradrhotree,treeep,nextptcl,prevptcl)
       use definitions
       use particle_module, only : ptot
       logical, intent(inout) :: treeep
       integer, intent(in) :: i
       integer, intent(in) :: isource
       integer, intent(in) :: level
       integer, intent(inout) :: p
       integer, intent(inout) :: pprev
       real(kind=PR), intent(inout) :: fh
       real(kind=PR), intent(inout) :: hep
       real(kind=PR), intent(inout) :: oldfh
       real(kind=PR), intent(in) :: oldstep
       real(kind=PR), intent(inout) :: ray_integral
       real(kind=PR), intent(inout) :: rep(1:NDIM)
       real(kind=PR), intent(in) :: rep_prev(1:NDIM)
       real(kind=PR), intent(in) :: rhoep
       real(kind=PR), intent(in) :: rsource(1:NDIM)
       real(kind=PR), intent(inout) :: rsqd_ep
       real(kind=PR), intent(inout) :: rtree(1:NDIM)
       real(kind=PR), intent(inout) :: htree
       real(kind=PR), intent(inout) :: rhotree
       real(kind=PR), intent(inout) :: gradrhotree(1:NDIM)
       integer, intent(inout) :: nextptcl(1:ptot)
       integer, intent(inout) :: prevptcl(1:ptot)
     END SUBROUTINE HP_propagate_UV_radiation

     SUBROUTINE HP_reorder_lists(newid)
       use definitions
       use particle_module, only : ptot
       integer, intent(in) :: newid(1:ptot)
     END SUBROUTINE HP_reorder_lists

     SUBROUTINE HP_rhoh_ep(p,hguess,rp,rhop,gradrhop)
       use definitions
       integer, intent(in) :: p
       real(kind=PR), intent(inout) :: hguess
       real(kind=PR), intent(in) :: rp(1:NDIM)
       real(kind=PR), intent(inout) :: rhop
       real(kind=PR), intent(out) :: gradrhop(1:NDIM)
     END SUBROUTINE HP_rhoh_ep

     SUBROUTINE HP_split_active_rays(isource,itot,level,nliverays,rsource,&
          &livelist,nextptcl,nextptclaux,prevptcl,prevptclaux,raylist)
       use definitions
       use healpix_types
       use particle_module, only : ptot
       integer, intent(in) :: isource
       integer, intent(inout) :: itot
       integer, intent(inout) :: level
       integer, intent(inout) :: nliverays
       integer, intent(inout) :: livelist(1:ptot)
       integer, intent(inout) :: nextptcl(1:ptot)
       integer, intent(inout) :: nextptclaux(1:ptot)
       integer, intent(inout) :: prevptcl(1:ptot)
       integer, intent(inout) :: prevptclaux(1:ptot)
       integer(kind=I4B), intent(inout) :: raylist(1:ptot)
       real(kind=PR), intent(in) :: rsource(1:NDIM)
     END SUBROUTINE HP_split_active_rays
  
     SUBROUTINE HP_walk_all_rays(isource)
       use definitions
       integer, intent(in) :: isource
     END SUBROUTINE HP_walk_all_rays

     SUBROUTINE HP_walk_ray(i,isource,level,resolution,&
          &nep,rsource,nextptcl,prevptcl)
       use definitions
       use particle_module, only : ptot
       integer, intent(in) :: i
       integer, intent(in) :: isource
       integer, intent(in) :: level
       integer, intent(inout) :: nep
       integer, intent(inout) :: nextptcl(1:ptot)
       integer, intent(inout) :: prevptcl(1:ptot)
       real(kind=PR), intent(out) :: resolution
       real(kind=PR), intent(in) :: rsource(1:NDIM)
     END SUBROUTINE HP_walk_ray

     SUBROUTINE hydro(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE hydro

     SUBROUTINE hydro_gradh(p)
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE hydro_gradh

     SUBROUTINE IFbinarychop(rayint_max,rep,rep_prev,hep,hprev,rhoep,&
          &rc,ray_integral,oldstep,pprev,p,nextptcl,plast,&
          &rtree,htree,rhotree,gradrhotree,treeep)
       !DEC$ ATTRIBUTES FORCEINLINE :: IFbinarychop
       use definitions
       use particle_module, only : ptot
       integer, intent(in) :: p
       integer, intent(in) :: plast
       integer, intent(in) :: pprev
       integer, intent(in) :: nextptcl(1:ptot)
       real(kind=PR), intent(in)    :: rayint_max
       real(kind=PR), intent(inout) :: rep(1:NDIM)
       real(kind=PR), intent(in)    :: rep_prev(1:NDIM)
       real(kind=PR), intent(inout) :: hep
       real(kind=PR), intent(in)    :: hprev
       real(kind=PR), intent(inout) :: rhoep
       real(kind=PR), intent(in)    :: rc(1:NDIM)
       real(kind=PR), intent(inout) :: ray_integral
       real(kind=PR), intent(inout)    :: oldstep
       real(kind=PR), intent(inout) :: rtree(1:NDIM)
       real(kind=PR), intent(inout) :: htree
       real(kind=PR), intent(inout) :: rhotree
       real(kind=PR), intent(inout) :: gradrhotree(1:NDIM)
       logical, intent(inout) :: treeep
     END SUBROUTINE IFbinarychop

     SUBROUTINE insertion_sort_real(nsort,iarray,rarray)
       !DEC$ ATTRIBUTES FORCEINLINE :: insertion_sort_real
       use definitions
       integer, intent(in) :: nsort
       integer, intent(inout) :: iarray(1:nsort)
       real(kind=PR), intent(inout) :: rarray(1:nsort)
     END SUBROUTINE insertion_sort_real

     SUBROUTINE insertion_sort_dp(nsort,iarray,rarray)
       !DEC$ ATTRIBUTES FORCEINLINE :: insertion_sort_dp
       use definitions
       integer, intent(in) :: nsort
       integer, intent(inout) :: iarray(1:nsort)
       real(kind=DP), intent(inout) :: rarray(1:nsort)
     END SUBROUTINE insertion_sort_dp

     SUBROUTINE insertion_sort_int(nsort,ilist)
       !DEC$ ATTRIBUTES FORCEINLINE :: insertion_sort_int
       use definitions
       integer, intent(in)    :: nsort
       integer, intent(inout) :: ilist(1:nsort)
     END SUBROUTINE insertion_sort_int

     SUBROUTINE isothermal_riemann_solver(vp,vpp,dr_unit,&
          &rho_p,rho_pp,sound_p,Pstar,vstar)
       !DEC$ ATTRIBUTES FORCEINLINE :: isothermal_riemann_solver
       use definitions
       real(kind=PR), intent(in) :: dr_unit(1:NDIM)
       real(kind=PR), intent(in) :: rho_p
       real(kind=PR), intent(in) :: rho_pp
       real(kind=PR), intent(in) :: sound_p
       real(kind=PR), intent(in) :: vp(1:NDIM)
       real(kind=PR), intent(in) :: vpp(1:NDIM)
       real(kind=PR), intent(out) :: Pstar
       real(kind=PR), intent(out) :: vstar
     END SUBROUTINE isothermal_riemann_solver

     SUBROUTINE mk_pix2xy(pix2x,pix2y)
       use healpix_types
       integer(kind=I4B), intent(inout), dimension(1:1023) :: pix2x
       integer(kind=I4B), intent(inout), dimension(1:1023) :: pix2y
     END SUBROUTINE mk_pix2xy

     SUBROUTINE mk_xy2pix(x2pix,y2pix)
       use healpix_types
       integer(kind=I4B), intent(inout), dimension(1:1023) :: x2pix
       integer(kind=I4B), intent(inout), dimension(1:1023) :: y2pix
     END SUBROUTINE mk_xy2pix

     SUBROUTINE nbody_hermite4_direct_gravity(s,invhs,rs,vs,agravs,adots,potp)
       use definitions
       integer, intent(in) :: s
       real(kind=DP), intent(in) :: invhs
       real(kind=DP), intent(in) :: rs(1:NDIM)
       real(kind=DP), intent(in) :: vs(1:NDIM)
       real(kind=DP), intent(out) :: agravs(1:NDIM)
       real(kind=DP), intent(out) :: adots(1:NDIM)
       real(kind=DP), intent(out) :: potp
     END SUBROUTINE nbody_hermite4_direct_gravity

     SUBROUTINE paramerror(errmsg)
       use definitions
       character(len=*), intent(in) :: errmsg
     END SUBROUTINE paramerror

     SUBROUTINE paramstore(store_file)
       use definitions
       character(len=256), intent(in) :: store_file
     END SUBROUTINE paramstore
     
     SUBROUTINE pix2vec_nest(nside,ipix,pix2x,pix2y,vector)
       use healpix_types
       integer(kind=I4B), intent(in) :: nside,ipix
       real(kind=DP), intent(out), dimension(1:3) :: vector
       integer(kind=I4B), intent(inout), dimension(0:1023) :: pix2x,pix2y
     END SUBROUTINE pix2vec_nest

     SUBROUTINE read_data(filename,file_form)
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: file_form
     END SUBROUTINE read_data

     SUBROUTINE read_parameters(filename)
       use definitions
       character(len=*), intent(in) :: filename
     END SUBROUTINE read_parameters

     SUBROUTINE reduce_particle_timestep(p,lnew)
       use definitions
       integer, intent(in) :: p
       integer(kind=ILP), intent(in) :: lnew
     END SUBROUTINE reduce_particle_timestep

     SUBROUTINE remove_from_list(ndead,deadlist)
       use definitions
       use particle_module, only : ptot
       integer, intent(in) :: ndead
       integer, intent(in) :: deadlist(1:ptot)
     END SUBROUTINE remove_from_list

     SUBROUTINE reorder_array_double(pstart,ptot,rmain,aorder)
       use definitions
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       integer, intent(in)    :: aorder(1:ptot)
       real(kind=DP), intent(inout) :: rmain(1:ptot)
     END SUBROUTINE reorder_array_double

     SUBROUTINE reorder_array_int(pstart,ptot,iarray,aorder)
       use definitions
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       integer, intent(inout) :: iarray(1:ptot)
       integer, intent(in)    :: aorder(1:ptot)
     END SUBROUTINE reorder_array_int

     SUBROUTINE reorder_array_int_2D(nsize,pstart,ptot,iarray,aorder)
       use definitions
       integer, intent(in)    :: nsize
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       integer, intent(inout) :: iarray(1:nsize,1:ptot)
       integer, intent(in)    :: aorder(1:ptot)
     END SUBROUTINE reorder_array_int_2D

     SUBROUTINE reorder_array_logical(pstart,ptot,larray,aorder)
       use definitions
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       logical, intent(inout) :: larray(1:ptot)
       integer, intent(in)    :: aorder(1:ptot)
     END SUBROUTINE reorder_array_logical

     SUBROUTINE reorder_array_long_int(pstart,ptot,iarray,aorder)
       use definitions
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       integer, intent(in)    :: aorder(1:ptot)
       integer(kind=ILP), intent(inout) :: iarray(1:ptot)
     END SUBROUTINE reorder_array_long_int

     SUBROUTINE reorder_array_real(pstart,ptot,rmain,aorder)
       use definitions
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       integer, intent(in)    :: aorder(1:ptot)
       real(kind=PR), intent(inout) :: rmain(1:ptot)
     END SUBROUTINE reorder_array_real

     SUBROUTINE reorder_array_real_2D(nsize,pstart,ptot,rmain,aorder)
       use definitions
       integer, intent(in)    :: nsize
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       integer, intent(in)    :: aorder(1:ptot)
       real(kind=PR), intent(inout) :: rmain(1:nsize,1:ptot)
     END SUBROUTINE reorder_array_real_2D

     SUBROUTINE reorder_inverse_array_int(pstart,ptot,iarray,aorder)
       use definitions
       integer, intent(in)    :: pstart
       integer, intent(in)    :: ptot
       integer, intent(inout) :: iarray(1:ptot)
       integer, intent(in)    :: aorder(1:ptot)
     END SUBROUTINE reorder_inverse_array_int

     SUBROUTINE reorder_particle_arrays(pstart,plast,dummylist)
       use definitions
       use particle_module, only : ptot
       integer, intent(in) :: pstart
       integer, intent(in) :: plast
       integer, intent(in) :: dummylist(1:ptot)
     END SUBROUTINE reorder_particle_arrays

     SUBROUTINE riemann_solver(gamma_eff,vp,vpp,dr_unit,&
          &rho_p,rho_pp,press_p,press_pp,Pstar,vstar)
       !DEC$ ATTRIBUTES FORCEINLINE :: riemann_solver
       use definitions
       real(kind=PR), intent(in) :: dr_unit(1:NDIM)
       real(kind=PR), intent(in) :: gamma_eff
       real(kind=PR), intent(in) :: press_p
       real(kind=PR), intent(in) :: press_pp
       real(kind=PR), intent(in) :: rho_p
       real(kind=PR), intent(in) :: rho_pp
       real(kind=PR), intent(in) :: vp(1:NDIM)
       real(kind=PR), intent(in) :: vpp(1:NDIM)
       real(kind=PR), intent(out) :: Pstar
       real(kind=PR), intent(out) :: vstar
     END SUBROUTINE riemann_solver

     SUBROUTINE sink_accretion_properties(s,maccreted)
       use definitions
       integer, intent(in) :: s
       real(kind=DP), intent(in) :: maccreted
     END SUBROUTINE sink_accretion_properties

     SUBROUTINE sink_timestep(s,dt_ideal)
       use definitions
       integer, intent(in) :: s
       real(kind=DP), intent(out) :: dt_ideal
     END SUBROUTINE sink_timestep

     SUBROUTINE thermal(p)
       !DEC$ ATTRIBUTES FORCEINLINE :: thermal
       use definitions
       integer, intent(in) :: p
     END SUBROUTINE thermal

     SUBROUTINE timestep_size(p,dt)
       use definitions
       integer, intent(in)        :: p
       real(kind=DP), intent(out) :: dt
     END SUBROUTINE timestep_size

     SUBROUTINE timing(current_mark)
       use definitions
       character(len=*), intent(in) :: current_mark
     END SUBROUTINE timing

     SUBROUTINE track_particle(p,rorigin)
       use definitions
       integer, intent(in)      :: p
       real(kind=PR), intent(in) :: rorigin(1:NDIM)
     END SUBROUTINE track_particle

     SUBROUTINE write_accreted_particles(s,pp_tot,pp_templist)
       use definitions
       use particle_module, only : ptot
       integer, intent(in) :: s
       integer, intent(in) :: pp_tot
       integer, intent(in) :: pp_templist(1:pp_tot)
     END SUBROUTINE write_accreted_particles

     SUBROUTINE write_data(filename,file_form)
       use definitions
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: file_form
     END SUBROUTINE write_data

     SUBROUTINE write_data_debug(out_file,rorigin)
       use definitions
       character(len=*), intent(in) :: out_file
       real(kind=PR), intent(inout) :: rorigin(1:NDIM)
     END SUBROUTINE write_data_debug

     SUBROUTINE write_debug_column_info(unitno)
       use definitions
       integer, intent(in) :: unitno
     END SUBROUTINE write_debug_column_info

     SUBROUTINE write_makefile_options(unitno)
       use definitions
       integer, intent(in) :: unitno
     END SUBROUTINE write_makefile_options


     SUBROUTINE add_integer_parameter(param_name,ipointer,idefault)
       use definitions
       character(len=*), intent(in) :: param_name
       integer, target, intent(in) :: ipointer
       integer, intent(in) :: idefault
     END SUBROUTINE add_integer_parameter


  SUBROUTINE add_long_integer_parameter(param_name,ipointer,idefault)
       use definitions
       character(len=*), intent(in) :: param_name
       integer(kind=ILP), target, intent(in) :: ipointer
       integer(kind=ILP), intent(in) :: idefault
     END SUBROUTINE add_long_integer_parameter

     SUBROUTINE add_real_parameter(param_name,rpointer,rdefault)
       use definitions
       character(len=*), intent(in) :: param_name
       real(kind=PR), target, intent(in) :: rpointer
       real(kind=PR), intent(in) :: rdefault
     END SUBROUTINE add_real_parameter

     SUBROUTINE add_double_parameter(param_name,dpointer,ddefault)
       use definitions
       character(len=*), intent(in) :: param_name
       real(kind=DP), target, intent(in) :: dpointer
       real(kind=DP), intent(in) :: ddefault
     END SUBROUTINE add_double_parameter

     SUBROUTINE add_string_parameter(param_name,cpointer,cdefault)
       use definitions
       character(len=*), intent(in) :: param_name
       character(len=*), target, intent(in) :: cpointer
       character(len=*), intent(in) :: cdefault
     END SUBROUTINE add_string_parameter

     SUBROUTINE add_unit_parameter(param_name,upointer,udefault)
       use definitions
       character(len=*), intent(in) :: param_name
       character(len=*), target, intent(in) :: upointer
       character(len=*), intent(in) :: udefault
     END SUBROUTINE add_unit_parameter

     SUBROUTINE add_logical_parameter(param_name,lpointer,ldefault)
       use definitions
       character(len=*), intent(in) :: param_name
       logical, target, intent(in) :: lpointer
       logical, intent(in) :: ldefault
     END SUBROUTINE add_logical_parameter

  END INTERFACE

END MODULE interface_module

