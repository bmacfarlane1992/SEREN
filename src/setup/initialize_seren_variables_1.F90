! INITIALIZE_SEREN_VARIABLES_1.F90
! D. A. Hubber - 1/10/2007
! Initializes various variables in Seren before the main simulation begins.  
! Should be called immediatly after the parameters have been set/read-in.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_seren_variables_1
  use neighbour_module
  use filename_module
  use type_module
  use particle_module
  use sink_module
  use hydro_module
  use timing_module
  use periodic_module
  implicit none

  logical :: flag                       ! flag if restart file exists
  character(len=256) :: restart_file    ! snapshot from restart log

  debug1("First initialisation of seren variables [initialize_seren_variables_1.F90]")


! Calculate average neighbour number for 'grad-h' SPH scheme
! ----------------------------------------------------------------------------
#if defined(GRAD_H_SPH)
#if NDIM==1
  pp_gather = int(2.0_PR*KERNRANGE*h_fac)
#elif NDIM==2
  pp_gather = int(PI*(KERNRANGE*h_fac)**2)
#elif NDIM==3
  pp_gather = int((4.0_PR*PI/3.0_PR)*(KERNRANGE*h_fac)**3)
#endif
#endif
  pp_limit = pp_gather*2


! Set filename variables
! ----------------------------------------------------------------------------
  if (out_file_form=="dragon_form" .or. out_file_form=="df") then
     fileform_ext = ".df"
  else if (out_file_form=="dragon_unform" .or. out_file_form=="du") then
     fileform_ext = ".du"
  else if (out_file_form=="seren_form" .or. out_file_form=="sf") then
     fileform_ext = ".sf"
  else if (out_file_form=="seren_unform" .or. out_file_form=="su") then
     fileform_ext = ".su"
  end if

  run_dir = trim(adjustl(run_dir))//"/"
  out_init  = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".ini"
  out_final = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".fin"
  out_temp1 = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".tmp1"
  out_temp2 = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".tmp2"
  restart_log = trim(adjustl(run_dir))//trim(adjustl(run_id))//".restart"


! Check if a temp restart file has been generated.  If yes, use this 
! instead of the file named in the params file.  
! ----------------------------------------------------------------------------
  inquire(file=restart_log,exist=flag)
  if (flag) then
     open(1,file=restart_log,status='unknown',form='formatted')
     read(1,'(a)') restart_file
     read(1,'(a)') in_file_form
     in_file = restart_file
     write(6,*) "Restart file   : ",trim(restart_file)
     write(6,*) "Restart format : ",trim(in_file_form)
     if (restart_file == out_temp1) then
        ntemp = 2
     else if (restart_file == out_temp2) then
        ntemp = 1
     else
        ntemp = 1
     end if
     close(1)
     restart = .true.
     inifile = .true.
  else
     ntemp = 1
     inifile = .false.
  end if


! Initialise particle counters
! ----------------------------------------------------------------------------
  pboundary = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 0


! Initialize all timing variables
! ----------------------------------------------------------------------------
#if defined(TIMING)
  itime    = 0
  rtime    = 0.0_DP
  last_id  = 0
  mark_tot = 0
  marker_id(1:NBLOCKS) = " "
  ngravcomp = 0_ILP
  nhydrocomp = 0_ILP
  nsphcomp = 0_ILP
  iblock(1:NBLOCKS) = 0
  rblock(1:NBLOCKS) = 0.0_DP
#endif


! Initialize periodic variables
! ----------------------------------------------------------------------------
  periodic_size(1:3) = periodic_max(1:3) - periodic_min(1:3)
  periodic_half(1:3) = 0.5_PR*periodic_size(1:3)


! Initialize misc. variables
! ----------------------------------------------------------------------------
  call random_seed(rseed)
#if defined(DEBUG_PLOT_DATA)
  rzero(1:NDIM) = 0.0_PR
#endif
#if defined(ARTIFICIAL_CONDUCTIVITY)
  alpha_cond = 1.0_PR
#endif


  return
END SUBROUTINE initialize_seren_variables_1
