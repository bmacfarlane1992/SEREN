! COPY_BHGRAV_TO_BHHYDRO.F90
! D. A. Hubber - 19/7/2008
! Copies tree structure of BHgrav tree to BHhydro tree.  Used when all 
! of the SPH particles are self-gravitating gas particles therefore 
! avoiding duplicating an identical tree build for both hydro and gravity.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE copy_BHgrav_to_BHhydro
  use tree_module
  implicit none

  integer :: c                     ! Cell counter

  debug2("Copying gravity to hydro tree [copy_BHgrav_to_BHhydro.F90]")

  ctot_hydro = ctot_grav
  ltot_hydro = ltot_grav

! Copy level information for stocking
  first_cell_hydro(0:ltot_hydro) = first_cell_grav(0:ltot_hydro)
  last_cell_hydro(0:ltot_hydro)  = last_cell_grav(0:ltot_hydro)

! Now copy all necessary lists and variables
!$OMP PARALLEL DO DEFAULT(SHARED)
  do c=0,ctot_hydro
     BHhydro(c)%leaf = BHgrav(c)%leaf
     BHhydro(c)%nextcell = BHgrav(c)%nextcell
     BHhydro(c)%ifopen = BHgrav(c)%ifopen
     BHhydro(c)%plist(1:LEAFMAX) = BHgrav(c)%plist(1:LEAFMAX)
  end do
!$OMP END PARALLEL DO

  return
END SUBROUTINE copy_BHgrav_to_BHhydro
