! DEFINITIONS.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Definition of precision kind variables
! ============================================================================

#include "macros.h"

! ============================================================================
MODULE definitions

  integer, parameter :: DP = selected_real_kind(p=15) ! double precision
  integer, parameter :: QP = selected_real_kind(p=33) ! quadruple precision
  integer, parameter :: SP = selected_real_kind(p=6)  ! single precision
  integer, parameter :: ILP = selected_int_kind(r=15) ! integer long precision

#if defined(QUADRUPLE_PRECISION)
  integer, parameter :: PR = QP
#elif defined(DOUBLE_PRECISION)
  integer, parameter :: PR = DP
#else
  integer, parameter :: PR = SP                       ! default = single
#endif

END MODULE definitions


! ============================================================================
