!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz 
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!=================================================================================================================================
!> Module to handle commandline arguments 
!=================================================================================================================================
MODULE MOD_Commandline_Arguments
IMPLICIT NONE

! Global variables for command line argument parsing
INTEGER                              :: nArgs              !< number of command line argumens
CHARACTER(LEN=255),ALLOCATABLE       :: Args(:)            !< command line arguments

INTERFACE ParseCommandlineArguments
  MODULE PROCEDURE ParseCommandlineArguments
END INTERFACE ParseCommandlineArguments

!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Reads all commandline arguments
!==================================================================================================================================
SUBROUTINE ParseCommandlineArguments()
! MODULES
USE MOD_Globals
USE MOD_StringTools     ,ONLY: STRICMP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iArg,nArgs_tmp
CHARACTER(LEN=255)      :: tmp
LOGICAL,ALLOCATABLE     :: alreadyRead(:)
!==================================================================================================================================
! Get number of command line arguments
nArgs_tmp=COMMAND_ARGUMENT_COUNT()
ALLOCATE(alreadyRead(nArgs_tmp))
alreadyRead = .FALSE.

! Get keyword arguments (arbitrary order)
doGenerateUnittestReferenceData = .FALSE.
doPrintHelp = 0
nArgs = nArgs_tmp
DO iArg = 1, nArgs_tmp
  CALL GET_COMMAND_ARGUMENT(iArg,tmp)
  IF (STRICMP(tmp, "--generateUnittestReferenceData")) THEN
    doGenerateUnittestReferenceData = .TRUE.
    alreadyRead(iArg) = .TRUE.
    nArgs = nArgs - 1
  END IF
  IF (STRICMP(tmp, "--help").OR.STRICMP(tmp,"-h")) THEN
    doPrintHelp = 1
    alreadyRead(iArg) = .TRUE.
    nArgs = nArgs - 1
  END IF
  IF (STRICMP(tmp, "--markdown")) THEN
    doPrintHelp = 2
    alreadyRead(iArg) = .TRUE.
    nArgs = nArgs - 1
  END IF
END DO ! iArg = 1, nArgs

! Get all remaining parameters
nArgs = MAX(1,nArgs) ! at least one argument is generated (empty)
ALLOCATE(Args(nArgs))
Args(1) = ""
nArgs = 0
DO iArg = 1,nArgs_tmp
  IF (.NOT.alreadyRead(iArg)) THEN
    nArgs = nArgs + 1
    CALL GET_COMMAND_ARGUMENT(iArg,Args(nArgs))
    alreadyRead(iArg) = .TRUE.
  END IF
END DO

DEALLOCATE(alreadyRead)
END SUBROUTINE ParseCommandlineArguments

!==================================================================================================================================
!> Finalizes commandline arguments
!==================================================================================================================================
SUBROUTINE FinalizeCommandlineArguments() 
IMPLICIT NONE
!===================================================================================================================================
SDEALLOCATE(Args)
END SUBROUTINE FinalizeCommandlineArguments

END MODULE MOD_Commandline_Arguments
