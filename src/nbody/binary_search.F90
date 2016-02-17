! BINARY_SEARCH.F90
! D. A. Hubber - 28/01/2008
! Identifies bound binary systems from the ensemble of sink particles 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE binary_search
  use interface_module, only : binary_energy,binary_properties,&
       &distance3_dp,insertion_sort_dp
  use definitions
  use sink_module
  use Nbody_module
  use scaling_module
  use constant_module
  use time_module, only : n,nsteps,time
  use filename_module, only : run_dir,run_id
  implicit none

  character(len=256) :: filename1             ! Output filename
  character(len=256) :: filename2             ! Output filename
  character(len=256) :: filename3             ! Output filename
  logical :: ex                               ! Does file exist?
  integer :: binflag(1:2*SMAX)                ! Status flag of stars
  integer :: i                                ! Auxilary loop counter
  integer :: most_bound(1:2*SMAX,1:2*SMAX)    ! List of most bound members
  integer :: nsingles                         ! Number of single stars
  integer :: s                                ! Sink counter
  integer :: ss                               ! Secondary sink counter
  integer :: sinkid(1:2*SMAX)                 ! Sink id list
  real(kind=DP) :: binen                      ! Two-body energy
  real(kind=DP) :: energy(1:2*SMAX,1:2*SMAX)  ! Energy of all system pairs
  real(kind=DP) :: energy_s(1:2*SMAX)         ! Energy of pairs for system s
  real(kind=DP) :: htemp(1:2*SMAX)            ! Smoothing lengths of systems
  real(kind=DP) :: mtemp(1:2*SMAX)            ! Masses of systems
  real(kind=DP) :: rtemp(1:NDIM,1:2*SMAX)     ! Positions of systems
  real(kind=DP) :: vtemp(1:NDIM,1:2*SMAX)     ! Velocities of systems

! Immediately return if there are not enough sinks
  if (stot <= 1) return

  debug2("Search for bound multiple systems [binary_search.F90]")

! Initialize variables and arrays
  nbin = 0
  nsingles = stot
  do s=1,2*SMAX
     binflag(s) = 0
     do ss=1,2*SMAX
        energy(s,ss) = BIG_NUMBER_DP
        most_bound(s,ss) = -1
     end do
  end do

! Store mechanical data of sinks
  do s=1,stot
     rtemp(1:NDIM,s) = star(s)%r(1:NDIM)
     vtemp(1:NDIM,s) = star(s)%v(1:NDIM)
     htemp(s) = star(s)%h
     mtemp(s) = star(s)%m
  end do

#if defined(DEBUG_BINARY_SEARCH)
  write(6,*) "Finished initialising variables [binary_search.F90]"
#endif


! Calculate the energy of all sink pairs
! ----------------------------------------------------------------------------
  do s=1,stot
     do ss=1,stot
        if (s == ss) cycle
        call binary_energy(s,ss,htemp(s),htemp(ss),mtemp(s),mtemp(ss),&
             &rtemp(1:NDIM,s),rtemp(1:NDIM,ss),vtemp(1:NDIM,s),&
             &vtemp(1:NDIM,ss),binen)
        energy(s,ss) = binen

     end do
  end do

#if defined(DEBUG_BINARY_SEARCH)
  write(6,*) "Sort energies into ascending order [binary_search.F90]"
#endif


! Sort energies into ascending order and store for later analysis
! ----------------------------------------------------------------------------
  do s=1,stot
     i = 0
     do ss=1,stot
        if (s == ss) cycle        
        i = i + 1
        energy_s(i) = energy(s,ss)
        sinkid(i) = ss
     end do

     ! Arrange in ascending order and store in main array
     call insertion_sort_dp(i,sinkid(1:i),energy_s(1:i))
     most_bound(s,1:i) = sinkid(1:i)

  end do


! Identify all mutually bound pairs of stars as binaries and calculate 
! and store binary properties
! ----------------------------------------------------------------------------
  do s=1,stot-1
     do ss=s+1,stot
        if (most_bound(s,1) == ss .and.  most_bound(ss,1) == s .and. &
             &energy(s,ss) < 0.0_DP .and. energy(ss,s) < 0.0_DP) then

           ! Adjust binary/single counters, and flag stars as part of binary
           nbin        = nbin + 1
           nsingles    = nsingles - 2
           binflag(s)  = stot + nbin
           binflag(ss) = stot + nbin
           
           call binary_properties(s,ss,energy(s,ss),mtemp(s),&
                &mtemp(ss),rtemp(1:NDIM,s),rtemp(1:NDIM,ss),&
                &vtemp(1:NDIM,s),vtemp(1:NDIM,ss))

           rtemp(1:NDIM,stot+nbin) = binary(nbin)%r(1:NDIM)
           vtemp(1:NDIM,stot+nbin) = binary(nbin)%v(1:NDIM)
           mtemp(stot+nbin) = binary(nbin)%m
           htemp(stot+nbin) = 0.0_DP

        end if
     end do
  end do

#if defined(DEBUG_BINARY_SEARCH)
  write(6,*) "Search for multiple systerms [binary_search.F90]"
#endif


! If we've identified any binaries and have at least two systems, repeat 
! search using binary COM particles for hierarchical multiple systems
! ----------------------------------------------------------------------------
  if (nbin > 0 .and. nbin + nsingles > 1) then

     ! Initialize variables and arrays for second search
     do s=1,2*SMAX
        do ss=1,2*SMAX
           energy(s,ss) = BIG_NUMBER_DP
           most_bound(s,ss) = -1
        end do
     end do

     ! Calculate 2-body energies of all pairs
     ! -----------------------------------------------------------------------
     do s=1,stot+nbin
        if (binflag(s) > 0) cycle

        ! First loop over all remaining single particles and new binaries
        do ss=1,stot+nbin
           if (binflag(ss) > 0 .or. s == ss) cycle
           call binary_energy(s,ss,htemp(s),htemp(ss),mtemp(s),mtemp(ss),&
                &rtemp(1:NDIM,s),rtemp(1:NDIM,ss),vtemp(1:NDIM,s),&
                &vtemp(1:NDIM,ss),binen)
           energy(s,ss) = binen
        end do
     end do


     ! Sort energies into ascending order and store for later analysis
     ! -----------------------------------------------------------------------
     do s=1,stot+nbin
        if (binflag(s) > 0) cycle
        i = 0
        do ss=1,stot+nbin
           if (binflag(ss) > 0 .or. s == ss) cycle
           i = i + 1
           energy_s(i) = energy(s,ss)
           sinkid(i) = ss
        end do
        
        ! Arrange in ascending order and store in main array
        call insertion_sort_dp(i,sinkid(1:i),energy_s(1:i))
        most_bound(s,1:i) = sinkid(1:i)
        
     end do


     ! Calculate binary properties of identified hierarchical systems
     ! -----------------------------------------------------------------------
     do s=1,stot+nbin-1
        if (binflag(s) > 0) cycle
        do ss=s+1,stot+nbin
           if (binflag(ss) > 0) cycle

           if (most_bound(s,1) == ss .and.  most_bound(ss,1) == s .and. &
                &energy(s,ss) < 0.0_DP .and. energy(ss,s) < 0.0_DP) then
              
              nbin = nbin + 1
              if (s <= stot)  nsingles = nsingles - 1
              if (ss <= stot) nsingles = nsingles - 1

              binflag(s)  = stot + nbin
              binflag(ss) = stot + nbin

              call binary_properties(s,ss,energy(s,ss),mtemp(s),&
                   &mtemp(ss),rtemp(1:NDIM,s),rtemp(1:NDIM,ss),&
                   &vtemp(1:NDIM,s),vtemp(1:NDIM,ss))
              
           end if
        end do
     end do

  end if
! ----------------------------------------------------------------------------



! Write properties to file
! ----------------------------------------------------------------------------
  if (nbin > 0) then
     filename1 = trim(adjustl(run_dir))//trim(adjustl(run_id))//".binstats"
     filename2 = trim(adjustl(run_dir))//trim(adjustl(run_id))//".finbinstats"
     inquire(file=filename1,exist=ex)
     if (.NOT. ex) then
        filename3 = trim(adjustl(run_dir))//trim(adjustl(run_id))&
             &//".inibinstats"
        open(unit=3,file=filename3,status="unknown")
     end if
     open(unit=1,file=filename1,status="unknown",position="append")
     open(unit=2,file=filename2,status="unknown")
     
#if NDIM==2
     ! -----------------------------------------------------------------------
10   format(2i10,2X,E12.4,2X,3i10,11E12.4)
     
     do s=1,nbin
        write(1,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(3)*mscale*vscale*rscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale*t_SI/yr, &
             &binary(s)%q, &
             &binary(s)%sma*rscale*r_SI/r_au, &
             &binary(s)%drmag*rscale*r_SI/r_au
        write(2,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(3)*mscale*vscale*rscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale*t_SI/yr, &
             &binary(s)%q, &
             &binary(s)%sma*rscale*r_SI/r_au, &
             &binary(s)%drmag*rscale*r_SI/r_au
        if (.not. ex) then
           write(2,10) n,nsteps,time*tscale,&
                &binary(s)%id, &
                &binary(s)%s1, &
                &binary(s)%s2, &
                &binary(s)%r(1)*rscale, &
                &binary(s)%r(2)*rscale, &
                &binary(s)%v(1)*vscale, &
                &binary(s)%v(2)*vscale, &
                &binary(s)%m*mscale, &
                &binary(s)%angmom(3)*mscale*vscale*rscale, &
                &binary(s)%ecc, &
                &binary(s)%period*tscale*t_SI/yr, &
                &binary(s)%q, &
                &binary(s)%sma*rscale*r_SI/r_au, &
                &binary(s)%drmag*rscale*r_SI/r_au
        end if
     end do

#elif NDIM==3
     ! -----------------------------------------------------------------------
10   format(2i10,2X,E12.4,2X,3i10,15E12.4)
     
     do s=1,nbin
        write(1,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%r(3)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%v(3)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(1)*mscale*vscale*rscale, &
             &binary(s)%angmom(2)*mscale*vscale*rscale, &
             &binary(s)%angmom(3)*mscale*vscale*rscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale*t_SI/yr, &
             &binary(s)%q, &
             &binary(s)%sma*rscale*r_SI/r_au, &
             &binary(s)%drmag*rscale*r_SI/r_au
        write(2,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%r(3)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%v(3)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(1)*mscale*vscale*rscale, &
             &binary(s)%angmom(2)*mscale*vscale*rscale, &
             &binary(s)%angmom(3)*mscale*vscale*rscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale*t_SI/yr, &
             &binary(s)%q, &
             &binary(s)%sma*rscale*r_SI/r_au, &
             &binary(s)%drmag*rscale*r_SI/r_au
        if (.not. ex) then
           write(2,10) n,nsteps,time*tscale,&
                &binary(s)%id, &
                &binary(s)%s1, &
                &binary(s)%s2, &
                &binary(s)%r(1)*rscale, &
                &binary(s)%r(2)*rscale, &
                &binary(s)%r(3)*rscale, &
                &binary(s)%v(1)*vscale, &
                &binary(s)%v(2)*vscale, &
                &binary(s)%v(3)*vscale, &
                &binary(s)%m*mscale, &
                &binary(s)%angmom(1)*mscale*vscale*rscale, &
                &binary(s)%angmom(2)*mscale*vscale*rscale, &
                &binary(s)%angmom(3)*mscale*vscale*rscale, &
                &binary(s)%ecc, &
                &binary(s)%period*tscale*t_SI/yr, &
                &binary(s)%q, &
                &binary(s)%sma*rscale*r_SI/r_au, &
                &binary(s)%drmag*rscale*r_SI/r_au
        end if
     end do
#endif
     ! -----------------------------------------------------------------------

     if (.not. ex) close(3)
     close(2)
     close(1)
  end if

  return
END SUBROUTINE binary_search
