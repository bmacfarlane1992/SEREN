! CALCULATE_SINK_PROPERTIES.F90
! D. A. Hubber - 30/11/2009
! Calculate the stellar, and stellar evolution properties of all sinks.
! For now, compares the sink mass to a look up table (c.f. Crowther et al.) 
! and interpolates the number of UV ionizing photons per second.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE calculate_sink_properties
  use sink_module
  use scaling_module
  use constant_module
  implicit none

  integer :: s                  ! Sink counter
  integer :: UVid               ! UV source id
  real(kind=DP) :: intmax_s     ! Maximum ionization integral
  real(kind=PR) :: ms           ! Mass of sink s
  real(kind=DP) :: N_LyC_s      ! No. of ionizing photons per second

  debug2("[calculate_sink_properties.F90]")

! Loop over all sinks
! ============================================================================
  do s=1,stot
     ms = sink(s)%m
        
     ! -----------------------------------------------------------------------
#if defined(MULTIPLE_SINK_SOURCES)

     ! Check where sink lies in look-up table
     itable = 0
     do
        itable = itable + 1
        if (ms > stellar_table(itable)%mass .or. itable = Ntable - 1) exit
     end do
     if (itable >= Ntable) itable = Ntable - 1
     intfactor = (ms - stellar_table(itable)%mass) / &
          (stellar_table(itable + 1)%mass - stellar_table(itable)%mass)
     if (intfactor > 1.0) intfactor = 1.0
     
     ! Calculated (interpolated) value of no. of UV photons per second
     N_LyC_s = (1 - intfactor)*stellar_table(itable)%log_N_LyC &
          & + intfactor*stellar_table(itable + 1)%log_N_LyC
     N_LyC_s = 10.D0**(real(N_LyC_s,DP))
     
     
     ! If there is some ionizing flux from the star, create a UV source 
     ! from the sink (if required) and record intmax for sources.
     if (ms > 0.0 .and. N_LyC_s > SMALL_NUMBER) then
        UVid = sink(s)%UVid
        
        ! Create a new source if it doesn't exist
        if (UVid == -1) then
           call create_UV_source(s,sink(s)%r(1:NDIM))
           UVid = UVtot
        end if
        
        N_LyC_s = N_LyC_s*tscale*t_SI
        intmax_s = N_LyC_s*(m_hydrogen/(mscale*m_SI*Xfrac))**2 / &
             &(4.*PI*a_star) 
        UVsource(UVid)%N_LyC = N_LyC_s
        UVsource(UVid)%intmax = intmax_s
     end if
     
#endif
     ! -----------------------------------------------------------------------
          
  end do
! ============================================================================


  return
END SUBROUTINE calculate_sink_properties
