!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
! MODULES
USE ISO_C_BINDING
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: HP => INT16,  & ! half precision (only defined for integer)
                                         SP => REAL32, & ! single precision
                                         DP => REAL64, & ! double precision
                                         QP => REAL128   ! quadruple precision
! MODULES
#if USE_MPI
USE mpi_f08
#endif
USE ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255):: ParameterFile                                            !< filename of the parameter file
INTEGER,PARAMETER :: UNIT_stdOut=6                                            !< unit for writing to standard output (e.g. terminal)
INTEGER,PARAMETER :: UNIT_logOut=133                                          !< unit for writing log files
INTEGER           :: UNIT_errOut=999                                          !< unit for writing error files
LOGICAL           :: Logging                                                  !< switch to turn log file writing on or of
LOGICAL           :: ErrorFiles                                               !< switch to turn error file writing on or of
CHARACTER(LEN=255):: ErrorFileName='NOT_SET'                                  !< file to write error data into
INTEGER           :: iError                                                   !< default error handle
INTEGER           :: myRank,myLocalRank,myLeaderRank,myWorkerRank
INTEGER           :: nProcessors,nLocalProcs,nLeaderProcs,nWorkerProcs
#if USE_MPI
TYPE(MPI_Comm)    :: MPI_COMM_FLEXI                                           !< Flexi MPI communicator
#endif
LOGICAL           :: MPIRoot                                                  !< flag whether process is MPI root process
LOGICAL           :: MPILocalRoot                                             !< flag whether process is root of MPI subgroup
#if USE_MPI
! TYPE(MPI_Status)  :: MPIStatus(MPI_STATUS_SIZE)
TYPE(MPI_Comm)    :: MPI_COMM_NODE   =MPI_COMM_NULL                           !< local node subgroup
TYPE(MPI_Comm)    :: MPI_COMM_LEADERS=MPI_COMM_NULL                           !< all node masters
TYPE(MPI_Comm)    :: MPI_COMM_WORKERS=MPI_COMM_NULL                           !< all non-master nodes
#endif

LOGICAL           :: doGenerateUnittestReferenceData
INTEGER           :: doPrintHelp ! 0: no help, 1: help, 2: markdown-help

LOGICAL           :: postiMode=.FALSE.                                        !< set TRUE if called from posti

! Overload the MPI interface because MPICH fails to provide it
! > https://github.com/pmodels/mpich/issues/2659
! > https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node263.htm
#if LIBS_MPICH_FIX_SHM_INTERFACE
INTERFACE MPI_WIN_ALLOCATE_SHARED
  SUBROUTINE PMPI_WIN_ALLOCATE_SHARED(SIZE, DISP_UNIT, INFO, COMM, BASEPTR, WIN, IERROR)
    ! MODULES
    USE, INTRINSIC ::  ISO_C_BINDING, ONLY : C_PTR
    IMPORT         ::  MPI_ADDRESS_KIND
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    INTEGER,INTENT(IN)     ::  DISP_UNIT, INFO, COMM, WIN, IERROR
    INTEGER(KIND=MPI_ADDRESS_KIND),INTENT(IN) :: SIZE
    TYPE(C_PTR),INTENT(IN) ::  BASEPTR
  END SUBROUTINE PMPI_WIN_ALLOCATE_SHARED
END INTERFACE

INTERFACE MPI_WIN_SHARED_QUERY
  SUBROUTINE PMPI_WIN_SHARED_QUERY(WIN, RANK, SIZE, DISP_UNIT, BASEPTR, IERROR)
    ! MODULES
    USE, INTRINSIC ::  ISO_C_BINDING, ONLY : C_PTR
    IMPORT         ::  MPI_ADDRESS_KIND
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    INTEGER,INTENT(IN)     :: WIN, RANK, DISP_UNIT, IERROR
    INTEGER(KIND=MPI_ADDRESS_KIND),INTENT(IN) :: SIZE
    TYPE(C_PTR),INTENT(IN) ::  BASEPTR
  END SUBROUTINE PMPI_WIN_SHARED_QUERY
END INTERFACE
#endif /*LIBS_MPICH_FIX_SHM_INTERFACE*/

INTERFACE
  SUBROUTINE setstacksizeunlimited() BIND(C)
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
  END SUBROUTINE setstacksizeunlimited
END INTERFACE

PUBLIC:: setstacksizeunlimited
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)       :: SourceFile      !< Source file where error has occurred
INTEGER,INTENT(IN)                :: SourceLine      !< Line in source file
CHARACTER(LEN=*),INTENT(IN)       :: CompDate        !< Compilation date
CHARACTER(LEN=*),INTENT(IN)       :: CompTime        !< Compilation time
CHARACTER(LEN=*),INTENT(IN)       :: ErrorMessage    !< Error message
INTEGER,OPTIONAL,INTENT(IN)       :: IntInfo         !< Error info (integer)
REAL,OPTIONAL,INTENT(IN)          :: RealInfo        !< Error info (real)
!   There is no way back!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=50)                 :: IntString,RealString
!==================================================================================================================================
IntString  = ''
RealString = ''

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
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
CALL MPI_FINALIZE(iError)
#endif
IF (MPIRoot) THEN
  ERROR STOP 1
ELSE
  STOP
END IF
END SUBROUTINE CollectiveStop


!==================================================================================================================================
!> Terminate program correctly if an error has occurred (important in MPI mode!).
!> Uses a MPI_ABORT which terminates FLEXI if a single proc calls this routine.
!==================================================================================================================================
SUBROUTINE Abort(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfo,RealInfo,ErrorCode)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)       :: SourceFile      !< Source file where error has occurred
INTEGER,INTENT(IN)                :: SourceLine      !< Line in source file
CHARACTER(LEN=*),INTENT(IN)       :: CompDate        !< Compilation date
CHARACTER(LEN=*),INTENT(IN)       :: CompTime        !< Compilation time
CHARACTER(LEN=*),INTENT(IN)       :: ErrorMessage    !< Error message
INTEGER,OPTIONAL,INTENT(IN)       :: IntInfo         !< Error info (integer)
REAL,OPTIONAL,INTENT(IN)          :: RealInfo        !< Error info (real)
INTEGER,OPTIONAL,INTENT(IN)       :: ErrorCode       !< MPI Error info (integer)
! There is no way back!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=50)                 :: IntString,RealString
#if USE_MPI
INTEGER                           :: errOut          ! Output of MPI_ABORT
INTEGER                           :: signalout       ! Output errorcode
#endif
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
CALL MPI_ABORT(MPI_COMM_FLEXI,signalout,errOut)
#endif
ERROR STOP 2
END SUBROUTINE Abort


!==================================================================================================================================
!> print a warning to the command line (only MPI root)
!==================================================================================================================================
SUBROUTINE PrintWarning(msg)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: msg  !< output message
!===================================================================================================================================
IF (myRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(UNIT_stdOut,'(A)') 'WARNING:'
  WRITE(UNIT_stdOut,'(A)') TRIM(msg)
  WRITE(UNIT_stdOut,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
END IF
END SUBROUTINE PrintWarning


!==================================================================================================================================
!> Convert a String to an Integer
!==================================================================================================================================
SUBROUTINE int2str(str,int_number,stat)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=255),INTENT(OUT) :: str
INTEGER,INTENT(IN)             :: int_number
INTEGER,INTENT(OUT)            :: stat
!===================================================================================================================================
WRITE(str,'(I0)',IOSTAT=stat)  int_number
END SUBROUTINE int2str


!==================================================================================================================================
!> Convert an Integer to a String
!==================================================================================================================================
!SUBROUTINE int2strf(str,int_number,stat)
FUNCTION int2strf(int_number)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=3) :: int2strf
INTEGER,INTENT(IN) :: int_number
!===================================================================================================================================
WRITE(int2strf,'(I0)')  int_number
int2strf = TRIM(ADJUSTL(int2strf))
END FUNCTION int2strf


!==================================================================================================================================
!> Convert a String to an Integer
!==================================================================================================================================
SUBROUTINE str2int(str,int_number,stat)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: str        !< input string
INTEGER,INTENT(OUT)         :: int_number !< output integer
INTEGER,INTENT(OUT)         :: stat       !< status
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  int_number
END SUBROUTINE str2int


!==================================================================================================================================
!> Convert a String to a REAL
!==================================================================================================================================
SUBROUTINE str2real(str,real_number,stat)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(len=*),INTENT(IN) :: str         !< input string
REAL,INTENT(OUT)            :: real_number !< output real
INTEGER,INTENT(OUT)         :: stat        !< status
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  real_number
END SUBROUTINE str2real


!==================================================================================================================================
!> Convert a String to a LOGICAL
!==================================================================================================================================
SUBROUTINE str2logical(str,logical_number,stat)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(len=*),INTENT(IN) :: str            !< input string
LOGICAL,INTENT(OUT)         :: logical_number !< output logical
INTEGER,INTENT(OUT)         :: stat           !< status
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  logical_number
END SUBROUTINE str2logical


! !==================================================================================================================================
! !> read compile flags from a specified file
! !> example line in "CMakeLists.txt": SET(FLEXI_EQNSYSNAME "navierstokes" CACHE STRING "Used equation system")
! !> ParameterName: timestep
! !> output: 0.1
! !> Type of Msg: [G]et[P]arameter[F]rom[File] -> GPFF: not ordinary read-in tool
! !==================================================================================================================================
! SUBROUTINE GetParameterFromFile(FileName,ParameterName,output,DelimiterSymbolIN,CommentSymbolIN,DoDisplayInfo)
! ! MODULES
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT/OUTPUT VARIABLES
! CHARACTER(LEN=*),INTENT(IN)          :: FileName          !< e.g. './../myfile'
! CHARACTER(LEN=*),INTENT(IN)          :: ParameterName     !< e.g. 'timestep'
! CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: DelimiterSymbolIN !< e.g. '=' (default is '=')
! CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: CommentSymbolIN   !< e.g. '#' (default is '!')
! CHARACTER(LEN=*),INTENT(INOUT)       :: output            !< e.g. '0.1'
! LOGICAL,OPTIONAL,INTENT(IN)          :: DoDisplayInfo     !< default is: TRUE
!                                                           !< display DefMsg or errors if the parameter or the file is not found
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! LOGICAL                              :: ExistFile         !> file exists=.true., file does not exist=.false.
! INTEGER                              :: iSTATUS           !> status
! CHARACTER(LEN=255)                   :: temp,temp2,temp3  !> temp variables for read in of file lines
! CHARACTER(LEN=255)                   :: DelimiterSymbol   !> symbol for commenting out code, e.g., "#" or "!"
! CHARACTER(LEN=255)                   :: CommentSymbol     !> symbol for commenting out code, e.g., "#" or "!"
! INTEGER                              :: ioUnit            !> field handler unit and ??
! INTEGER                              :: IndNum            !> Index Number
! CHARACTER(LEN=8)                     :: DefMsg            !> additional flag like "DEFAULT" or "*CUSTOM"
! !===================================================================================================================================
! IF(PRESENT(DelimiterSymbolIN))THEN
!   DelimiterSymbol=TRIM(ADJUSTL(DelimiterSymbolIN))
! ELSE
!   DelimiterSymbol='='
! END IF
! IF(PRESENT(CommentSymbolIN))THEN
!   CommentSymbol=TRIM(ADJUSTL(CommentSymbolIN))
! ELSE
!   CommentSymbol='!'
! END IF
! output=''
! ! read from file
! INQUIRE(File=TRIM(FileName),EXIST=ExistFile)
! IF(ExistFile) THEN
!   OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ')
!   DO
!     READ(ioUnit,'(A)',iostat=iSTATUS)temp
!     IF(ADJUSTL(temp(1:LEN(TRIM(CommentSymbol)))).EQ.TRIM(CommentSymbol)) CYCLE  ! complete line is commented out
!     IF(iSTATUS.EQ.-1)EXIT                           ! end of file is reached
!     IF(LEN(trim(temp)).GT.1)THEN                    ! exclude empty lines
!       IndNum=INDEX(temp,TRIM(ParameterName))        ! e.g. 'timestep'
!       IF(IndNum.GT.0)THEN
!         IF(IndNum-1.GT.0)THEN                       ! check if the parameter name is contained within a substring of another
!           IF(temp(IndNum-1:IndNum-1).NE.' ')CYCLE   ! parameter, e.g., "timestep" within "fd_timestep" -> skip
!         END IF
!         temp2=TRIM(ADJUSTL(temp(IndNum+LEN(TRIM(ParameterName)):LEN(temp))))
!         IF(DelimiterSymbol.NE.'')THEN               ! demiliting symbol must not be empty
!           IndNum=INDEX(temp2,TRIM(DelimiterSymbol)) ! only use string FROM delimiting symbol +1
!           IF(IndNum.GT.0)THEN
!             temp3=TRIM(ADJUSTL(temp2(IndNum+1:LEN(temp2))))
!             temp2=temp3
!           END IF
!         ELSE
!           ! no nothing?
!         END IF
!         IndNum=INDEX(temp2,TRIM(CommentSymbol)) ! only use string UP TO commenting symbol
!         IF(IndNum.EQ.0)IndNum=LEN(TRIM(temp2))+1
!         output=TRIM(ADJUSTL(temp2(1:IndNum-1)))
!         DefMsg='GPFF'
!         SWRITE(UNIT_stdOut,'(a3,a30,a3,a33,a3,a7,a3)')' | ',TRIM(ParameterName),' | ', output,' | ',TRIM(DefMsg),' | '
!         EXIT ! found the parameter -> exit loop
!       END IF
!     END IF
!   END DO
!   CLOSE(ioUnit)
!   IF(output.EQ.'')THEN
!     IF(PRESENT(DoDisplayInfo))THEN
!       IF(DoDisplayInfo)THEN
!         SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: Parameter ['//TRIM(ParameterName)//'] not found.'
!       END IF
!     ELSE
!       SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: Parameter ['//TRIM(ParameterName)//'] not found.'
!     END IF
!     output='ParameterName does not exist'
!   END IF
! ELSE
!   IF(PRESENT(DoDisplayInfo))THEN
!     IF(DoDisplayInfo)THEN
!       SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: File ['//TRIM(FileName)//'] not found.'
!     END IF
!   ELSE
!     SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: File ['//TRIM(FileName)//'] not found.'
!   END IF
!   output='file does not exist'
! END IF
! END SUBROUTINE GetParameterFromFile


!==================================================================================================================================
!> Open file for error output
!==================================================================================================================================
SUBROUTINE CreateErrFile()
! MODULES
! IMPLICIT VARIABLE HANDLING
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
FUNCTION FILEEXISTS(filename)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename   !< current filename
LOGICAL                     :: FILEEXISTS !< logical indicating if file with current filename exists
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
INQUIRE(FILE=TRIM(filename), EXIST=FILEEXISTS)
END FUNCTION FILEEXISTS


!==================================================================================================================================
!> Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!==================================================================================================================================
FUNCTION INTSTAMP(Nam,Num)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: Nam      !< Name
INTEGER,INTENT(IN)          :: Num      !< Number
CHARACTER(LEN=200)          :: IntStamp !< The stamp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
WRITE(IntStamp,'(A,A5,I6.6)')TRIM(Nam),'_Proc',Num
END FUNCTION INTSTAMP


!==================================================================================================================================
!> Creates a timestamp, consistent of a filename (project name + processor) and current time niveau
!==================================================================================================================================
FUNCTION TIMESTAMP(Filename,Time,Time2)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: Filename  !< (file)name
REAL,INTENT(IN)             :: Time      !< physical time
REAL,OPTIONAL,INTENT(IN)    :: Time2     !< physical time (in case of range)
CHARACTER(LEN=255)          :: TimeStamp !< the complete timestamp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i         ! loop variable
CHARACTER(LEN=255) :: tmp
!==================================================================================================================================
WRITE(TimeStamp,'(F17.9)')Time
! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
IF(PRESENT(Time2))THEN
  WRITE(tmp,'(F17.9)')Time2
  DO i=1,LEN(TRIM(tmp))
    IF(tmp(i:i).EQ.' ') tmp(i:i)='0'
  END DO
  TimeStamp=TRIM(tmp)//'-'//TRIM(TimeStamp)
END IF
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
END FUNCTION TIMESTAMP


SUBROUTINE DisplaySimulationTime(Time, StartTime, Message)
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: Message         !< Output message
REAL,INTENT(IN)             :: Time, StartTime !< Current simulation time and beginning of simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: SimulationTime,mins,secs,hours,days
CHARACTER(LEN=60) :: hilf
!===================================================================================================================================
! Return with all procs except root if not called during abort
IF(.NOT.MPIRoot.AND.(Message.NE.'ABORTED')) RETURN

! Output particle info
WRITE(UNIT_stdOut,'(132("="))')

! Calculate simulation time
SimulationTime = Time-StartTime

! Get secs, mins, hours and days
secs = MOD(SimulationTime,60.)
SimulationTime = SimulationTime / 60.
mins = MOD(SimulationTime,60.)
SimulationTime = SimulationTime / 60.
hours = MOD(SimulationTime,24.)
SimulationTime = SimulationTime / 24.
!days = MOD(SimulationTime,365.) ! Use this if years are also to be displayed
days = SimulationTime

! Output message with all procs, as root might not be the calling process during abort
WRITE(hilf,'(F16.2)') Time-StartTime
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')  ' FLEXI '//TRIM(Message)//'! [ '//TRIM(ADJUSTL(hilf))//' sec ]'
WRITE(UNIT_stdOut,'(A3,I0,A1,I0.2,A1,I0.2,A1,I0.2,A2)') ' [ ',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),' ]'
IF(MPIRoot.AND.(Message.NE.'ABORTED')) WRITE(UNIT_stdOut,'(132("="))')
END SUBROUTINE DisplaySimulationTime


!===================================================================================================================================
! Output message to UNIT_stdOut and an elapsed time is seconds as well as min/hour/day format
!===================================================================================================================================
SUBROUTINE DisplayMessageAndTime(ElapsedTimeIn, Message, DisplayLine, rank)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: Message          !< Output message
REAL,INTENT(IN)             :: ElapsedTimeIn    !< Time difference
LOGICAL,INTENT(IN),OPTIONAL :: DisplayLine      !< Display 132*"-" (default is TRUE)
INTEGER,INTENT(IN),OPTIONAL :: rank             !< if 0, some kind of root is assumed, every other processor return this routine
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: ElapsedTime,mins,secs,hours,days
LOGICAL           :: DisplayLineLoc,LocalRoot
CHARACTER(LEN=60) :: hilf
!===================================================================================================================================

! Define who returns and who does the output (default is MPIRoot)
LocalRoot = .FALSE. ! default
IF(PRESENT(rank))THEN
  IF(rank.EQ.0) LocalRoot = .TRUE.
ELSE
  IF(MPIRoot) LocalRoot = .TRUE.
END IF ! PRESENT(rank)

! Return with all procs except LocalRoot
IF(.NOT.LocalRoot) RETURN

! Check if 132*"-" is required
IF(PRESENT(DisplayLine))THEN
  DisplayLineLoc = DisplayLine
ELSE
  DisplayLineLoc = .TRUE.
END IF ! PRESENT(DisplayLine)

! Aux variable
ElapsedTime=ElapsedTimeIn

! Get secs, mins, hours and days
secs = MOD(ElapsedTime,60.)
ElapsedTime = ElapsedTime / 60.
mins = MOD(ElapsedTime,60.)
ElapsedTime = ElapsedTime / 60.
hours = MOD(ElapsedTime,24.)
ElapsedTime = ElapsedTime / 24.
!days = MOD(ElapsedTime,365.) ! Use this if years are also to be displayed
days = ElapsedTime

! Output message
IF(LocalRoot)THEN
  WRITE(hilf,'(F16.2)')  ElapsedTimeIn
  ! Only output the second if it is actually useful
  IF (ElapsedTimeIn.GT.60) THEN
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')  ' '//TRIM(Message)//' [ '//TRIM(ADJUSTL(hilf))//' sec ]'
    WRITE(UNIT_stdOut,'(A3,I0,A1,I0.2,A1,I0.2,A1,I0.2,A2)') ' [ ',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),' ]'
  ELSE
    WRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' '//TRIM(Message)//' [ '//TRIM(ADJUSTL(hilf))//' sec ]'
  END IF
  IF(DisplayLineLoc) WRITE(UNIT_StdOut,'(132("-"))')
END IF ! LocalRoot)

END SUBROUTINE DisplayMessageAndTime


!==================================================================================================================================
!> Calculates current time (own function because of a laterMPI implementation)
!==================================================================================================================================
FUNCTION FLEXITIME(     &
#if USE_MPI
                   Comm &
#endif /*USE_MPI*/
                  )
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
#if USE_MPI
TYPE(mpi_Comm),INTENT(IN),OPTIONAL :: Comm                                       !< global mpi communicator
#endif /*USE_MPI*/
REAL                               :: FlexiTime                                  !< output time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if USE_MPI
IF (PRESENT(Comm)) THEN
  CALL MPI_BARRIER(Comm          ,iError)
ELSE
  CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
END IF
#endif /*USE_MPI*/
GETTIME(FlexiTime)
END FUNCTION FLEXITIME

END MODULE MOD_Globals
