! REDISTRIBUTE_SINK_ANGMOM.F90
! D. A. Hubber & S. K. Walch - 1/6/2009
! Redistribute the angular momentum accured by a sink due to accretion onto 
! all neighbouring (r < ANGMOMRAD*rads) SPH particles.  
! Currently, the angular momentum is transfered to a particle at the end of 
! its timestep in the form of an instantaneous impulse at the end of the 
! particles timestep.  This allows neighbour particles to retain their 
! ideal timestep rather than being forced onto the smallest timestep.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE redistribute_sink_angmom
  use interface_module, only : distance3
  use particle_module
  use sink_module
  use type_module
  use hydro_module
  use kernel_module
  use time_module
  use scaling_module
  use filename_module
  use constant_module
  implicit none

  character(len=*),parameter :: string_jdistr = 'uniform'
  character(len=7),parameter :: deposit = 'shakura' ! 'instant' or 'shakura'

  real(kind=PR), parameter ::thicknesstoradius=0.1 ! h/r ratio
  integer :: i                              ! Aux. particle counter
  integer :: p                              ! Particle counter
  integer :: pp_near                        ! Particles near to sinks
  integer :: s                              ! Sink counter
  integer, allocatable :: pp_templist(:)    ! Temp. list for particle ids
  real(kind=PR) :: angle                    ! For angle dependent weighting
  real(kind=PR) :: checkLp(1:3)             ! Sum of added ang. mom. for check
  real(kind=PR) :: dr(1:3)                  ! Distance particle sink	
  real(kind=PR) :: drsqd                    ! Distance squared	
  real(kind=PR) :: Lfunc                    ! Ang. mom. redistribution funciton
  real(kind=PR) :: Lpart(1:3)               ! Ang. mom added to all particles
  real(kind=PR) :: Lpartmag                 ! Part of sink ang. mom to be given
  real(kind=PR) :: Lsinkmag                 ! Ang. mom of sink
  real(kind=PR) :: Lnetmag                  ! Net ang. mom of sink
  real(kind=PR) :: mp                       ! Particle mass
  real(kind=PR) :: pextratot(1:3)           ! Sum of added momentum
  real(kind=PR) :: rotfreq(1:3)             ! Rotational frequ Omega
  real(kind=PR) :: rp(1:NDIM)               ! Position of particle p
  real(kind=PR) :: rperp(1:3)               ! Shortest dist vector connecting 
                                            ! part and rrot 	
  real(kind=PR) :: rrot(1:3)                ! Rotation axis
  real(kind=PR) :: tvisc                    ! viscous time scale
  real(kind=PR) :: weightsum                ! Sum of weighting functions
  real(kind=PR) :: vextra(1:3)              ! Added velocity components
  real(kind=PR), allocatable :: rperpp(:)   ! Length of shortest dist vector 
  real(kind=PR), allocatable :: weightpp(:) ! Weighting for all neib. particles
#if defined(DEBUG_SINK_REMOVE_ANGMOM)
  character(len=8)  :: file_ext             ! filename ext for data output
  character(len=40) :: out_file             ! filename ext for data output
  logical :: ex                             ! Does file exist already?
  integer :: idummy                         ! ..
  integer :: ierr                           ! Aux. error variable
  integer :: ndead                          ! Total number of dead particles
  real(kind=PR) :: alpha1                   ! ..
  real(kind=PR) :: alpha2                   ! ..
  real(kind=PR) :: alpha3                   ! ..
  real(kind=PR) :: checkLptot               ! Tot ang mom removed from sink
  real(kind=PR) :: rdummy(1:13)             ! ..
  real(kind=PR) :: rrot2(1:3)               ! Rotation axis
  real(kind=PR) :: zaxis(1:3)		    ! z-axis unit vector
#endif

  debug2("[redistribute_sink_angmom.F90]")
  
  allocate(pp_templist(1:ptot))
  allocate(rperpp(1:ptot))
  allocate(weightpp(1:ptot))


! Loop over all sink particles
! ============================================================================
  do s=1,stot

     ! Zero variables
     weightsum        = 0.0_PR
     Lpart(1:3)       = 0.0_PR
     pextratot(1:3)   = 0.0_PR
     rperpp(1:ptot)   = 0.0_PR
     weightpp(1:ptot) = 0.0_PR

     ! Magnitude of total angular momentum ever accreted by sink
     Lnetmag = sqrt(dot_product(sink(s)%angmomnet(1:3),&
          &sink(s)%angmomnet(1:3)))
     Lsinkmag = sqrt(dot_product(sink(s)%angmom(1:3),sink(s)%angmom(1:3)))

#if defined(DEBUG_SINK_REMOVE_ANGMOM)     
     print*,'Lsinkmag : ',Lsinkmag,Lsinkmag*angmomscale
     print*,'Lnetmag  : ',Lnetmag,Lnetmag*angmomscale
     print*,'sinkmass : ',sink(s)%m*mscale
#endif

     ! If the sink contains any angular momentum, calculate the weightings  
     ! of how it will be distributed amongst the neighbouring particles.
     ! -----------------------------------------------------------------------
     if (Lsinkmag > SMALL_NUMBER .and. Lnetmag > SMALL_NUMBER) then    

        pp_near = 0
        rrot(1:3) = sink(s)%angmom(1:3) / Lsinkmag
        
        ! Assume Shakura-Sunyaev type viscosity to find viscous timescale
    sink(s)%menc=sink(s)%m

    tvisc =thicknesstoradius**(-2)*sqrt((sink(s)%radius)**3/(sink(s)%menc))/alpha_ss
    sink(s)%tvisc=tvisc

        if (deposit .eq. 'shakura') then
           Lpartmag = Lsinkmag*(1.0_PR - &
                 &exp(-timestep*real(2**(level_step - nlevel_sinks),DP)/tvisc))
      else if (deposit .eq. 'instant') then
           Lpartmag = Lsinkmag
        else
           print *,'No valid amount of angmom specified'
        end if

print*,Lpartmag, Lsinkmag, Lpartmag/Lsinkmag,tvisc,alpha_ss, sink(s)%radius, sink(s)%menc
print*, 'here:', timestep, level_step, nlevel_sinks

#if defined(DEBUG_SINK_REMOVE_ANGMOM)
        print*,'Lpartmag : ',Lpartmag*Lscale
        print*,'tvisc    : ',tvisc*tscale,timestep*real(2**(level_step - nlevel_sinks),DP)/tvisc
        print*,'rrot     : ',rrot(1:3)
#endif     

        ! Loop over all near-sink particles and calculate the sum of the 
        ! weighting functions to normalize later.
        ! --------------------------------------------------------------------
        do p=pgasstart,pgasend
           
           rp(1:NDIM) = parray(1:NDIM,p)
           mp = parray(MASS,p)
           call distance3(sink(s)%r(1:NDIM),rp(1:NDIM),dr(1:NDIM),drsqd)
           
           ! Only include gas particles within certain distance of sink
           ! -----------------------------------------------------------------
           if (drsqd <= (ANGMOMRAD*sink(s)%radius)**2 .and. &
                &p >= pgasstart .and. p <= pgasend) then

              rperp(1:3) = dr(1:3) - &
                   &dot_product(dr(1:3),rrot(1:3))*rrot(1:3)

              ! Store for further usage:
              pp_near = pp_near + 1
              pp_templist(pp_near) = p
              rperpp(pp_near) = sqrt(rperp(1)**2 + rperp(2)**2 + rperp(3)**2)

              ! Calculate the weighting:
              ! --------------------------------------------------------------
              select case (string_jdistr)
              case ('uniform')
                 weightpp(pp_near) = mp*laststep(p)/tvisc
                 weightsum = weightsum + mp
                 
              ! Redistribute according to Keplerian profile
              case ('kepler')
                 weightpp(pp_near) = mp/sqrt(rperpp(pp_near))
                 
              ! Redistribute according to v_phi ~ R
              case ('solidb')
                 weightpp(pp_near) = mp*(rperpp(pp_near)**2)*laststep(p)/tvisc
                 weightsum = weightsum + mp*(rperpp(pp_near)**2)
                 
              ! Angle dependent redistribution:
              case ('kepler_angle')
                 call get_angle(angle,dr(1:NDIM),rrot)
                 
              case ('soldib_angle')
                 call get_angle(angle,dr(1:NDIM),rrot)
   
              case default
                 print*, 'No valid redistribution-method defined!'
                 stop
              end select

              
           end if
           ! -----------------------------------------------------------------
           
        end do
        ! --------------------------------------------------------------------

#if defined(DEBUG_SINK_REMOVE_ANGMOM)
        print*,'pp_near   :',pp_near	
        print*,'weightsum :',weightsum 
#endif
    	
        ! Skip to the next sink if there are no neighbouring particles
        if (pp_near == 0) cycle
        checkLp(1:3) = 0.0_PR

        ! Get extra velocities by redistributing Lpartmag 
        ! (only nearest neighbours)
        ! --------------------------------------------------------------------
        do i=1,pp_near
           p = pp_templist(i)
           
           ! Skip particle unless it is at the end of its timestep
           if (n /= nlast(p)) cycle

           mp = parray(MASS,p)
           rp(1:NDIM) = parray(1:NDIM,p)
           
           ! Calculate relative displacement in plane of ang. mom. vector
           call distance3(sink(s)%r(1:NDIM),rp(1:NDIM),dr(1:NDIM),drsqd)
           
           ! Angular momentum prescription for particle
           ! (for now, uniform specific angular momentum)
           Lfunc = Lpartmag*weightpp(i)/weightsum
           rotfreq(1:NDIM) = rrot(1:NDIM)*Lfunc/(mp*rperpp(i)**2)

           ! v = rotfreq x dr
           vextra(1) = rotfreq(2)*dr(3) - rotfreq(3)*dr(2)
           vextra(2) = rotfreq(3)*dr(1) - rotfreq(1)*dr(3)
           vextra(3) = rotfreq(1)*dr(2) - rotfreq(2)*dr(1)

           ! checkLp calculated with dr
           checkLp(1) = checkLp(1) + mp*(dr(2)*vextra(3) - dr(3)*vextra(2))
           checkLp(2) = checkLp(2) + mp*(dr(3)*vextra(1) - dr(1)*vextra(3))
           checkLp(3) = checkLp(3) + mp*(dr(1)*vextra(2) - dr(2)*vextra(1))

           ! Linear momentum and mass
           pextratot(1:3) = pextratot(1:3) + vextra(1:3)*mp

           ! Update particle velocities : (For leapfrog-kdk scheme, we
           ! need to add extra velocity to v_half due to the way v_old is 
           ! recalculated at the beginning of advance_leapfrog_kdk.F90).
           v(1:VDIM,p) = v(1:VDIM,p) + vextra(1:3)
           v_old(1:VDIM,p) = v(1:3,p)
#if defined(LEAPFROG_KDK)
           v_half(1:VDIM,p) = v_half(1:VDIM,p)  + vextra(1:3)
#endif

!#if defined(DEBUG_SINK_REMOVE_ANGMOM)
!           write(6,*) "i :",i,"  p :",p,rp(1:NDIM)
!           write(6,*) "Lfunc     :",Lfunc,rperpp(i)
!           write(6,*) "rotfreq   :",rotfreq(1:NDIM)
!           write(6,*) "rrot      :",rrot(1:NDIM)
!           write(6,*) "vextra    :",vextra(1:3)
!           write(6,*) "checkLp   :",checkLp(1:3)
!           write(6,*) "pextratot :",pextratot(1:3)
!#endif

        end do        
        ! --------------------------------------------------------------------

        ! Adjust sink properties to conserve linear and angular momentum
        sink(s)%angmom(1:3)  = sink(s)%angmom(1:3) - checkLp(1:3)
        sink(s)%v(1:NDIM)    = sink(s)%v(1:NDIM) - pextratot(1:NDIM)/sink(s)%m
        sink(s)%vold(1:NDIM) = sink(s)%v(1:NDIM)
#if defined(LEAPFROG_KDK)
        sink(s)%vhalf(1:VDIM) = sink(s)%vhalf(1:VDIM) &
             & - pextratot(1:NDIM)/sink(s)%m
#endif


        ! --------------------------------------------------------------------
#if defined(DEBUG_SINK_REMOVE_ANGMOM)
        ! Calculate total angular momentum taken away from sink
        checkLptot = sqrt(checkLp(1)**2 + checkLp(2)**2 + checkLp(3)**2)
        zaxis(1)   = 0.0_PR
        zaxis(2)   = 0.0_PR
        zaxis(3)   = 1.0_PR
        alpha1     = acos(dot_product(rrot(1:3),checkLp(1:3)/checkLptot))
        alpha2     = acos(dot_product(rrot(1:3),zaxis(1:3)))
        rrot2(1:3) = sink(s)%angmomnet(1:3) / Lnetmag
        alpha3     = acos(dot_product(rrot2(1:3),zaxis(1:3)))
        Lsinkmag   = sqrt(dot_product(sink(s)%angmom(1:3),sink(s)%angmom(1:3)))

        print*,'========== TOTAL J TAKEN AWAY FROM SINK ===='
        print*,'checkLptot   :',checkLptot*angmomscale
        print*,'Lpartmag     :',Lpartmag*angmomscale
!        print*,'alpha1       :',alpha1
!        print*,'alpha2       :',alpha2
!        print*,'alpha3       :',alpha3
        print*,'checkLp      :',checkLp(1:3)*angmomscale
        print*,'angmomnet    :',sink(s)%angmomnet(1:3)*angmomscale
        print*,'angmom       :',sink(s)%angmom(1:3)*angmomscale

        ! Write information to file
        
        ! Make filename from runid and sink number
        if (s >= 100) then
           write(file_ext,"(I3)") s
        else if (s>=10) then
           write(file_ext,"(I2)") s
        else
           write(file_ext,"(I1)") s
        end if
        out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))&
             &//".angmom"//trim(adjustl(file_ext))
        
        ! Check if file exists
        inquire(file=out_file,exist=ex)
        
        ! Read through file to synchronise files with time of simulation.
        if (ex) then
           open(1,file=out_file,status="unknown",&
                &form="formatted",position="append")
           do
              read(1,'(I8,8E18.10)',end=10,err=10) idummy,rdummy(1:13)
              if (rdummy(1) > time*tscale) exit
           end do
10         backspace (1,err=20)
        else
           open(1,file=out_file,status="unknown",form="formatted")
        end if
20      write(1,'(1I8,13E18.10)') nsteps,time*tscale,Lsinkmag*angmomscale,&
             &Lnetmag*angmomscale,alpha1,alpha2,alpha3,&
             &sink(s)%m*mscale,sink(s)%dmdt*dmdtscale,&
             &sink(s)%menc*mscale,tvisc*tscale,sink(s)%trot*tscale,&
             &weightsum,sum(weightpp(1:ptot))/weightsum
        close(1)
#endif
        ! --------------------------------------------------------------------

     end if
     ! -----------------------------------------------------------------------

  end do
! ============================================================================
  
  deallocate(weightpp)
  deallocate(rperpp)
  deallocate(pp_templist)

  return
END SUBROUTINE redistribute_sink_angmom



! ============================================================================
! GET_ANGLE
! ..
! ============================================================================
SUBROUTINE get_angle(angle,dr, rrot)
  use definitions
  
  real(kind=PR) :: angle
  real(kind=PR) :: dr(1:NDIM)
  real(kind=PR) :: rrot(1:NDIM)
  
  ! if 0< angle< 180 then angle=acos(dr rrot/(|dr| |rrot|))
  ! otherwise make sure to change direction of rrot!
  
  
  return
END SUBROUTINE get_angle
! ============================================================================
