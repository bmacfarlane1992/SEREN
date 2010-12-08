! READ_ARGUMENTS.F90
! D. A. Hubber - 21/02/2010
! Process command line arguments for the parameters file.  
! If no argument is given, then use the default 'params.dat' file.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_arguments
  use interface_module, only : write_debug_column_info,write_makefile_options
  use filename_module
  implicit none

  integer :: narguments            ! No. of command line arguments
  character(len=256) :: auxstring  ! Aux. string variable

  narguments = command_argument_count()
  write(6,*) "No. of command line arguments : ",narguments

! Terminate program if there are too many command line arguments
! ----------------------------------------------------------------------------
  if (narguments > 1) then 
     stop 'Too many command line arguments!!'

! Use default "params.dat" file if there is no argument
! ----------------------------------------------------------------------------
  else if (narguments == 0) then
     param_file = "params.dat"

! Process individual argument
! ----------------------------------------------------------------------------
  else if (narguments == 1) then
     call get_command_argument(1,auxstring)
     select case(auxstring)
     case ("-H","-h","-help")
        write(6,*)
        write(6,*) "Available command line options in Seren"
        write(6,*) "---------------------------------------"
        write(6,*) "-d, -D, -debug         : Data columns in debug files"
        write(6,*) "-diag                  : Data columns in diagnostic file"
        write(6,*) "-h, -H, -help          : help"
        write(6,*) "-m, -M, -makefile      : &
             &Makefile options used to compile Seren"
        write(6,*) "-s, -S, -sinks, -stars : Data columns in sink files"
        write(6,*) "-v, -V, -version       : version number"
        write(6,*) "'paramsfile'           : Reads named parameter file instead of default"
        stop
     case ("-V","-v","-version")
        write(6,*) "Seren version 0.9.9.9"
        stop
     case ("-M","-m","-makefile")
        call write_makefile_options(6)
        stop
     case ("-D","-d","-debug")
        call write_debug_column_info(6)
        stop
     case ("-S","-s","-sinks","-stars")
        call write_sink_column_info(6)
        stop
     case ("-diag")
        call write_diagnostic_column_info(6)
        stop
     case default
        param_file = trim(adjustl(auxstring))
     end select

  end if
! ----------------------------------------------------------------------------

  write(6,*) "param_file : ",trim(adjustl(param_file))

  return
END SUBROUTINE read_arguments
