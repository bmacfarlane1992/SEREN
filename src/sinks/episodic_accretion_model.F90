! EPISODIC_ACCRETION_MODEL
! D. Stamatellos - 08/06/
! Calculates ....
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE episodic_accretion_model(s)
  use sink_module
  use time_module
  use scaling_module
  use constant_module
  use Eos_module, only : rad_const
  implicit none

  integer, intent(in) :: s              ! sink id

  real(kind=PR) :: dt                 	! sink timestep
  real(kind=PR) :: dmdt_max	      	! maximum accretion rate
  real(kind=PR) :: DMDT_episodic      	! episodic accretion rate
  real(kind=PR) :: DM_EA		! mass delivered in one accretion event (from Zhu et al. 2010, ApJ 713)
  real(kind=PR) :: DT_EA		! duraton of accretion event
  real(kind=PR) :: DMDT_EA		! accretion rate during the event
  real(kind=PR) :: dmdt_star            ! aux accretion rate

!  dt=time-episodic_time
  dmdt_star=0.0_DP
  DMDT_EA=0.0_DP
  DT_EA=0.0_DP
  DM_EA=0.0_DP
  dt= sink_dt 
  dmdt_max=1e-2/dmdtscale  
  DMDT_episodic=0.0_DP
 
  episodic_time=time
! when there is no accretion onto the star

 if (sink(s)%accretion_flag==0)  then 

! calculate mass threshold for initiating an episodic accretion event

  DT_EA=(960*1e-6/tscale)*(0.1/alpha_EA)*(sink(s)%dmdt*dmdtscale/1e-4)**(1./9.)*(sink(s)%Mstar*mscale)**(2./3.)
  DMDT_EA=(5e-4/dmdtscale)*(alpha_EA/0.1)
  DM_EA=DMDT_EA*DT_EA

       	if (sink(s)%Mdisc>DM_EA .and. real(sink(s)%m,DP)>feedback_minmass) then
        	sink(s)%accretion_flag=1
         	sink(s)%t_episode_start=time
         	sink(s)%t_episode_duration=DT_EA
         	sink(s)%dmdt_0=min(1.58*sink(s)%Mdisc/sink(s)%t_episode_duration,dmdt_max)    
        else

       		if (dmdt_regular>sink(s)%dmdt) then
         		dmdt_star=sink(s)%dmdt
       		else 
         		dmdt_star=dmdt_regular
       		endif
! ensure that the drop in accretion rate cannot be more than 10% in each step
         	if (dmdt_star<0.9*sink(s)%dmdt_star) then
          		sink(s)%dmdt_star=0.9*sink(s)%dmdt_star
         	else
          		sink(s)%dmdt_star=dmdt_star
         	endif

         	sink(s)%Mdisc=sink(s)%Mdisc-sink(s)%dmdt_star*dt
         
		if (sink(s)%Mdisc<0) sink(s)%Mdisc=0.0

         	sink(s)%Mstar=sink(s)%m-sink(s)%Mdisc 
       endif
 
 endif

! when there is episodic accretion onto the star

  if (sink(s)%accretion_flag==1)  then	

! check if 'restart' of EA event is needed
      DT_EA=(960*1e-6/tscale)*(0.1/alpha_EA)*(sink(s)%dmdt*dmdtscale/1e-4)**(1./9.)*(sink(s)%Mstar*mscale)**(2./3.)
      DMDT_EA=(5e-4/dmdtscale)*(alpha_EA/0.1)
      DM_EA=DMDT_EA*DT_EA

       if (sink(s)%Mdisc>DM_EA) then
         sink(s)%t_episode_start=time
         sink(s)%t_episode_duration=DT_EA
         sink(s)%dmdt_0=min(1.58*sink(s)%Mdisc/sink(s)%t_episode_duration,dmdt_max)
       endif
 
       DMDT_episodic=sink(s)%dmdt_0*exp(-(time-sink(s)%t_episode_start)/sink(s)%t_episode_duration)

        dmdt_star=dmdt_regular+DMDT_episodic
! ensure the increase in accretion rate cannot be more than 100%
        if (dmdt_star>2*sink(s)%dmdt_star) then 
           sink(s)%dmdt_star=2*sink(s)%dmdt_star
        else
           sink(s)%dmdt_star=dmdt_star
        endif
        !if (sink(s)%dmdt>=sink(s)%dmdt_star) sink(s)%dmdt_star=sink(s)%dmdt
        sink(s)%Mdisc=sink(s)%Mdisc-sink(s)%dmdt_star*dt
        if (sink(s)%Mdisc<0) sink(s)%Mdisc=0
        sink(s)%Mstar=sink(s)%m-sink(s)%Mdisc 
     
        if (time-sink(s)%t_episode_start>sink(s)%t_episode_duration.OR. sink(s)%Mdisc<=0.0_DP) sink(s)%accretion_flag=0

 endif





  return
END SUBROUTINE episodic_accretion_model
