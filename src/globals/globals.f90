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

!==================================================================================================================================
!> Provides parameters, used globally (please use EXTREMLY carefully!)
!==================================================================================================================================
MODULE MOD_Globals
USE ISO_C_BINDING
! MODULES
#if USE_MPI
USE mpi
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)::ParameterFile                                             !< filename of the parameter file
INTEGER,PARAMETER ::UNIT_stdOut=6                                             !< unit for writing to standard output (e.g. terminal)
INTEGER,PARAMETER ::UNIT_logOut=133                                           !< unit for writing log files
INTEGER           ::UNIT_errOut=999                                           !< unit for writing error files
LOGICAL           ::Logging                                                   !< switch to turn log file writing on or of
LOGICAL           ::ErrorFiles                                                !< switch to turn error file writing on or of
CHARACTER(LEN=255)::ErrorFileName='NOT_SET'                                   !< file to write error data into
INTEGER           ::iError                                                    !< default error handle
REAL              ::StartTime                                                 !< start time of the simulation
INTEGER           ::myRank,myLocalRank,myLeaderRank,myWorkerRank
INTEGER           ::nProcessors,nLocalProcs,nLeaderProcs,nWorkerProcs
INTEGER           ::MPI_COMM_NODE                                             !< local node subgroup
INTEGER           ::MPI_COMM_LEADERS                                          !< all node masters
INTEGER           ::MPI_COMM_WORKERS                                          !< all non-master nodes
LOGICAL           ::MPIRoot                                                   !< flag whether process is MPI root process
LOGICAL           ::MPILocalRoot                                              !< flag whether process is root of MPI subgroup
#if USE_MPI
INTEGER           ::MPIStatus(MPI_STATUS_SIZE)
#endif

LOGICAL           :: doGenerateUnittestReferenceData                         
INTEGER           :: doPrintHelp ! 0: no help, 1: help, 2: markdown-help

INTERFACE Abort
  MODULE PROCEDURE Abort
END INTERFACE Abort

INTERFACE CollectiveStop
  MODULE PROCEDURE CollectiveStop
END INTERFACE CollectiveStop

INTERFACE PrintWarning
  MODULE PROCEDURE PrintWarning
END INTERFACE PrintWarning

INTERFACE INTSTAMP
  MODULE PROCEDURE INTSTAMP
END INTERFACE INTSTAMP

INTERFACE TIMESTAMP
  MODULE PROCEDURE TIMESTAMP
END INTERFACE

INTERFACE FLEXITIME
  MODULE PROCEDURE FLEXITIME
END INTERFACE

INTERFACE CreateErrFile
  MODULE PROCEDURE CreateErrFile
END INTERFACE CreateErrFile

!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> \brief Safely terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!> 
!> Safely terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!> Terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!> This routine can only be used if ALL processes are guaranteed to generate the same error at the same time!
!> Prime use is to exit FLEXI without MPI errors and with a single error message if some parameters are not set in the init
!> routines or a file is not found.
!>
!> Criteria where CollectiveStop may be used:
!> 0. In case of doubt stick with Abort, which is always safe!
!> 1. A routine is BY DESIGN (!) called by all processes, i.e. does not permit to be called by single processes or subgroups.
!> 2. The criteria for the CollectiveStop must be identical among all processors.
!> 3. The routine is only used during the init phase.
!> 4. The error must not originate from MPI errors (e.g. during MPI init)
!> 5. The error must not originate from checking roundof errors (e.g. accuracy of interpolation matrices)
!>
!==================================================================================================================================
SUBROUTINE CollectiveStop(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfo,RealInfo)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                  :: SourceFile      !< Source file where error has occurred
INTEGER                           :: SourceLine      !< Line in source file
CHARACTER(LEN=*)                  :: CompDate        !< Compilation date
CHARACTER(LEN=*)                  :: CompTime        !< Compilation time
CHARACTER(LEN=*)                  :: ErrorMessage    !< Error message
INTEGER,OPTIONAL                  :: IntInfo         !< Error info (integer)
REAL,OPTIONAL                     :: RealInfo        !< Error info (real)
!   There is no way back!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=50)                 :: IntString,RealString
!==================================================================================================================================
IntString = ""
RealString = ""

IF (PRESENT(IntInfo))  WRITE(IntString,"(A,I0)")  "\nIntInfo:  ", IntInfo
IF (PRESENT(RealInfo)) WRITE(RealString,"(A,F24.19)") "\nRealInfo: ", RealInfo

SWRITE(UNIT_stdOut,*) '_____________________________________________________________________________\n', &
                     'Program abort caused on Proc ',myRank, '\n', &
                     '  in File : ',TRIM(SourceFile),' Line ',SourceLine, '\n', &
                     '  This file was compiled at ',TRIM(CompDate),'  ',TRIM(CompTime), '\n', &
                     'Message: ',TRIM(ErrorMessage), &
                     TRIM(IntString), TRIM(RealString)

CALL FLUSH(UNIT_stdOut)
#if USE_MPI
CALL MPI_FINALIZE(iError)
#endif
ERROR STOP 1
END SUBROUTINE CollectiveStop


!==================================================================================================================================
!> Terminate program correctly if an error has occurred (important in MPI mode!).
!> Uses a MPI_ABORT which terminates FLEXI if a single proc calls this routine.
!==================================================================================================================================
SUBROUTINE Abort(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfo,RealInfo,ErrorCode)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                  :: SourceFile      !< Source file where error has occurred
INTEGER                           :: SourceLine      !< Line in source file
CHARACTER(LEN=*)                  :: CompDate        !< Compilation date
CHARACTER(LEN=*)                  :: CompTime        !< Compilation time
CHARACTER(LEN=*)                  :: ErrorMessage    !< Error message
INTEGER,OPTIONAL                  :: IntInfo         !< Error info (integer)
REAL,OPTIONAL                     :: RealInfo        !< Error info (real)
INTEGER,OPTIONAL                  :: ErrorCode       !< Error info (integer)
!   There is no way back!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=50)                 :: IntString,RealString
INTEGER                           :: errOut          ! Output of MPI_ABORT
INTEGER                           :: signalout       ! Output errorcode
!==================================================================================================================================
IntString = ""
RealString = ""

IF (PRESENT(IntInfo))  WRITE(IntString,"(A,I0)")  "\nIntInfo:  ", IntInfo
IF (PRESENT(RealInfo)) WRITE(RealString,"(A,F24.19)") "\nRealInfo: ", RealInfo

WRITE(UNIT_stdOut,*) '_____________________________________________________________________________\n', &
                     'Program abort caused on Proc ',myRank, '\n', &
                     '  in File : ',TRIM(SourceFile),' Line ',SourceLine, '\n', &
                     '  This file was compiled at ',TRIM(CompDate),'  ',TRIM(CompTime), '\n', &
                     'Message: ',TRIM(ErrorMessage), &
                     TRIM(IntString), TRIM(RealString)

CALL FLUSH(UNIT_stdOut)
#if USE_MPI
signalout=2 ! MPI_ABORT requires an output error-code /=0
IF(PRESENT(ErrorCode)) signalout=ErrorCode
CALL MPI_ABORT(MPI_COMM_WORLD,signalout,errOut)
#endif  
ERROR STOP 2
END SUBROUTINE Abort


!==================================================================================================================================
!> print a warning to the command line (only MPI root)
!==================================================================================================================================
SUBROUTINE PrintWarning(msg) 
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*) :: msg
!===================================================================================================================================
IF (myRank.EQ.0) THEN 
  WRITE(UNIT_stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(UNIT_stdOut,*) 'WARNING:'
  WRITE(UNIT_stdOut,*) TRIM(msg)
  WRITE(UNIT_stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
END IF 
END SUBROUTINE PrintWarning


!==================================================================================================================================
!> Open file for error output
!==================================================================================================================================
SUBROUTINE CreateErrFile()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
LOGICAL                        :: isOpen
!==================================================================================================================================
IF (ErrorFiles) THEN
  INQUIRE(UNIT=UNIT_errOut,OPENED=isOpen)
  IF(.NOT.isOpen)THEN
    OPEN(UNIT=UNIT_errOut,  &
        FILE=ErrorFileName,&
        STATUS='REPLACE',  &
        ACTION='WRITE',    &
        IOSTAT=OpenStat)
  END IF
END IF
END SUBROUTINE CreateErrFile


!==================================================================================================================================
!> Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!==================================================================================================================================
FUNCTION INTSTAMP(Nam,Num)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)   :: Nam      !< Name
INTEGER            :: Num      !< Number
CHARACTER(LEN=200) :: IntStamp !< The stamp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
WRITE(IntStamp,'(A,A5,I6.6)')TRIM(Nam),'_Proc',Num
END FUNCTION INTSTAMP


!==================================================================================================================================
!> Creates a timestamp, consistent of a filename (project name + processor) and current time niveau
!==================================================================================================================================
FUNCTION TIMESTAMP(Filename,Time)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)   :: Filename  !< (file)name
REAL               :: Time      !< physical time
CHARACTER(LEN=255) :: TimeStamp !< the complete timestamp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i         ! loop variable
!==================================================================================================================================
WRITE(TimeStamp,'(F17.9)')Time
! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
END FUNCTION TIMESTAMP


!==================================================================================================================================
!> Calculates current time (own function because of a laterMPI implementation)
!==================================================================================================================================
FUNCTION FLEXITIME(Comm)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN),OPTIONAL    :: Comm                                       !< global mpi communicator
REAL                            :: FlexiTime                                  !< output time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if USE_MPI
IF(PRESENT(Comm))THEN
  CALL MPI_BARRIER(Comm,iError)
ELSE
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END IF
#endif
GETTIME(FlexiTime)
END FUNCTION FLEXITIME


END MODULE MOD_Globals
