! CONVERT_TO_CODE_UNITS.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Converts all particle variables from physical units to dimensionless code 
! units using scaling factors calculated in units.F90.  Variables are 
! converted to dimensionless form by dividing by Xscale 
! (where X is r,t,m etc..). 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE convert_to_code_units
  use constant_module
  use scaling_module
  use particle_module
  use periodic_module
  use sink_module
  use hydro_module
  use time_module
  use type_module
  use neighbour_module
  use Nbody_module
  implicit none

  integer :: p        ! counter to loop over particles
#if defined(SINKS)
  integer :: s        ! sink counter
#endif

  debug1("Converting to dimensionless code units [convert_to_code_units.F90]")

! Scale main particle data to code units 
! ----------------------------------------------------------------------------
  do p=1,ptot
     parray(1:NDIM,p) = parray(1:NDIM,p) / real(rscale,PR)
     parray(MASS,p)   = parray(MASS,p) / real(mscale,PR)
     parray(SMOO,p)   = parray(SMOO,p) / real(rscale,PR)
     v(1:VDIM,p)      = v(1:VDIM,p) / real(vscale,PR)
     rho(p)           = rho(p) / real(rhoscale,PR)
#if defined(INTERNAL_ENERGY)
     u(p)             = u(p) / real(uscale,PR)
#endif
#if defined(IDEAL_MHD)
     B(1:BDIM,p)      = B(1:BDIM,p) / real(Bscale,PR)
#endif
  end do

! Scale sink variables
! ----------------------------------------------------------------------------
#if defined(SINKS)
  do s=1,stot
     sink(s)%tcreate     = sink(s)%tcreate / tscale
     sink(s)%r(1:NDIM)   = sink(s)%r(1:NDIM) / real(rscale,PR)
     sink(s)%v(1:VDIM)   = sink(s)%v(1:VDIM) / real(vscale,PR)
     sink(s)%m           = sink(s)%m / real(mscale,PR)
     sink(s)%h           = sink(s)%h / real(rscale,PR)
     sink(s)%radius      = sink(s)%radius / real(rscale,PR)
     sink(s)%angmom(1:3) = sink(s)%angmom(1:3) / angmomscale
     sink(s)%dmdt        = sink(s)%dmdt / dmdtscale
     sink(s)%star_radius = sink(s)%star_radius / rscale
     sink(s)%luminosity  = sink(s)%luminosity / Lscale
     sink(s)%macc(1:DMDT_RANGE) = sink(s)%macc(1:DMDT_RANGE) / mscale
     sink(s)%tacc(1:DMDT_RANGE) = sink(s)%tacc(1:DMDT_RANGE) / tscale
  end do
#endif

! Scale periodic variables
  periodic_min(1:3)  = periodic_min(1:3) / real(rscale,PR)
  periodic_max(1:3)  = periodic_max(1:3) / real(rscale,PR)
  periodic_half(1:3) = periodic_half(1:3) / real(rscale,PR)
  periodic_size(1:3) = periodic_size(1:3) / real(rscale,PR)

! Scale SPH simulation time variables 
  dt_fixed    = dt_fixed  / tscale
  firstsnap   = firstsnap / tscale
  lastsnap    = lastsnap / tscale
  snaptime    = snaptime  / tscale
  sph_endtime = sph_endtime / tscale
  time        = time / tscale

! Scale N-body/SPH simulation time variables
  nbody_sph_endtime = nbody_sph_endtime / tscale

! Scale N-body simulation time variables
  nbody_endtime  = nbody_endtime / tscale

! Scale density variables 
  rhobary = rhobary / real(rhoscale*rhocgs,PR)
  rhosink = rhosink / real(rhoscale*rhocgs,PR)

! Other misc. variables
  Pext       = Pext / real(Pscale,PR)
  mgas_orig  = mgas_orig / mscale
  hmin       = hmin / real(rscale,PR)
  rspheremax = rspheremax / real(rscale,PR)
  rholost    = rholost / real(rhoscale*rhocgs,PR)
  rad_lost   = rad_lost / real(rscale,PR)

! If using fixed absolute sink radius, scale variable from au to code units
#if defined(FIXED_ABSOLUTE_SINKRAD)
  sinkrad = sinkrad * real(r_au/(rscale*r_SI),PR)
#endif


  return
END SUBROUTINE convert_to_code_units
