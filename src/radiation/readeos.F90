! READEOS.F90
! D. Stamatellos - 3/1/2007
! Reads the EOS table (dens, temp, energy, mu)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_eos
  use hydro_module
  use Eos_module
  use scaling_module
  use constant_module, only : m_sun,kappa_const
#if defined(SPIEGEL_TEST)
  use Tprof_module
#endif
  implicit none
  
  logical :: ex                  ! Does 'eos.dat' file exist?
  integer :: m                   ! Density grid counter 
  integer :: n                   ! Temperature grid counter
  real(kind=DP) :: auxscale      ! Aux. scaling variable
  real(kind=DP) :: rdummy6(1:6)  ! Dummy array for reading arrays
  
  debug1("Reading in eos tables [read_eos.F90]")

! Auxilary scaling variable for specific energy
  auxscale = (Escale*Ecgs) / (mscale*mcgs)

! Open and read opacity tables (if the file exists)
! ----------------------------------------------------------------------------
  inquire(file='eos.dat',exist=ex)
  if (ex) then
     open(2,file='eos.dat',status='old',form='formatted')
     
     ! Read-in dimensions of EOS table
     read (2,*) dim_dens, dim_temp
     
     ! Now allocate array sizes according to table dimensions
     allocate(eos_dens(1:dim_dens))
     allocate(eos_temp(1:dim_temp))
     allocate(eos_energy(dim_dens,dim_temp))
     allocate(eos_mu(dim_dens,dim_temp))  
     allocate(kappa(dim_dens,dim_temp))
     allocate(kappap(dim_dens,dim_temp))
     
     ! Read-in data from eos.dat
     do m=1,dim_dens
        do n=1,dim_temp
           read (2,*) rdummy6(1:6)
           eos_dens(m)     = real(rdummy6(1),PR)
           eos_temp(n)     = real(rdummy6(2),PR)
           eos_energy(m,n) = real(rdummy6(3),PR)
           eos_mu(m,n)     = real(rdummy6(4),PR)
           kappa(m,n)      = real(rdummy6(5),PR)
           kappap(m,n)     = real(rdummy6(6),PR)
        end do
     end do

     close(2)
  else
     STOP "eos.dat file not found"
  end if
! ----------------------------------------------------------------------------

#ifdef SPIEGEL_TEST 

#ifndef SPIEGEL_DISPERSION
  write(*,*) "SPIEGEL_TEST: multiplying opacities by tau=",ptemp_q,&
       &"(defined as ptemp_q in params.dat)"
!assumes an input sphere of density 1.41ee-19 g/cm3 
! (from Masunaga & Inutsuka test)
    do m=1,dim_dens
        do n=1,dim_temp
        kappa(m,n)=kappa(m,n)*ptemp_q/4.259e-3
        kappap(m,n)=kappap(m,n)*ptemp_q/4.259e-3
        end do
     end do

#else

write(*,*) "SPIEGEL_TEST: multiplying opacities by 2.952e6*10**(-0.2*ptemp_q)(defined as ptemp_q in params.dat)" 

! Assumes an input sphere of density 1.41ee-19 g/cm3 
! (from Masunaga & Inutsuka test)
   do m=1,dim_dens
        do n=1,dim_temp
         kappa(m,n)=kappa(m,n)*2.952e6*10**(-0.2*ptemp_q)
         kappap(m,n)=kappap(m,n)*2.952e6*10**(-0.2*ptemp_q)
        end do
   end do

#endif 

#endif

! Convert arrays to code units 
  do m=1,dim_dens
     eos_dens(m) = eos_dens(m) / (rhoscale*rhocgs)
     do n=1,dim_temp
        eos_energy(m,n) = eos_energy(m,n) / auxscale
        kappa(m,n) = kappa(m,n) / (kappascale*kappacgs)
        kappap(m,n) = kappap(m,n) / (kappascale*kappacgs)
     end do
  end do

! Calculate log factors for fast table referencing
  densmin = eos_dens(1)
  densmax = eos_dens(dim_dens)
  bdens = (real(dim_dens,PR) - 1.0_PR) / (log10(densmax/densmin))

  tempmin = eos_temp(1)
  tempmax = eos_temp(dim_temp)
  btemp = (real(dim_temp,PR) - 1.0_PR) / (log10(tempmax/tempmin))
  
  return
END SUBROUTINE read_eos
