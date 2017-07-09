#include "flexi.h"

!==================================================================================================================================
!> Contains the utilize routines of the regressioncheck
!> -GetExampleList extracts the examples which are subfolders in examples
!> -InitExample reads in the parameter_reggie.ini file
!==================================================================================================================================
MODULE MOD_RegressionCheck_Tools
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE GetCommandLineOption
  MODULE PROCEDURE GetCommandLineOption
END INTERFACE

INTERFACE GetExampleList
  MODULE PROCEDURE GetExampleList
END INTERFACE

INTERFACE InitExample
  MODULE PROCEDURE InitExample
END INTERFACE

INTERFACE CheckForExecutable
  MODULE PROCEDURE CheckForExecutable
END INTERFACE

INTERFACE SummaryOfErrors
  MODULE PROCEDURE SummaryOfErrors
END INTERFACE

INTERFACE AddError
  MODULE PROCEDURE AddError
END INTERFACE


INTERFACE CheckFileForString
  MODULE PROCEDURE CheckFileForString
END INTERFACE

INTERFACE REGGIETIME
  MODULE PROCEDURE REGGIETIME
END INTERFACE

INTERFACE CalcOrder
  MODULE PROCEDURE CalcOrder
END INTERFACE

INTERFACE GetConfigurationFile
  MODULE PROCEDURE GetConfigurationFile
END INTERFACE

INTERFACE ConfigurationCounter
  MODULE PROCEDURE ConfigurationCounter
END INTERFACE


PUBLIC::GetExampleList,InitExample,CheckForExecutable,GetCommandLineOption
PUBLIC::SummaryOfErrors
PUBLIC::AddError
PUBLIC::CheckFileForString
PUBLIC::REGGIETIME
PUBLIC::CalcOrder
PUBLIC::GetConfigurationFile
PUBLIC::ConfigurationCounter
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> reads the command line options for the regressioncheck
!> options are:
!> run [default]:  - runs only the regressioncheck
!> build           - builds all valid compiler flag combinations (default uses the configuration.reggie from run_basic) and
!>                   performs the tests
!>
!> ./regressioncheck [RuntimeOption] [RuntimeOptionType]
!>
!> ./regressioncheck                -> uses default "run" and runs the current compiler build and all "run_" examples
!> ./regressioncheck
!> ./regressioncheck build          -> runs "run_basic" for numerous builds
!> ./regressioncheck build convtest -> runs "feature_convtest" for numerous builds def. "feature_convtest/configuration.reggie"
!> ./regressioncheck build all      -> runs all examples for numerous builds defined in "run_basic/configuration.reggie"
!==================================================================================================================================
SUBROUTINE GetCommandLineOption()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars, ONLY: RuntimeOption,BuildNoDebug,BuildDebug,DoFullReggie
USE MOD_RegressionCheck_Vars, ONLY: BuildContinue,BuildSolver,NumberOfProcs,RunContinue
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: I
INTEGER                        :: nArgs                   !> Number of supplied command line arguments
INTEGER                        :: iSTATUS                 !> I/O status
!===================================================================================================================================
RuntimeOption(1)='run'        ! only run pre-built executable (no building of new cmake compiler flag combinations)
RuntimeOption(2)='run_basic'  ! set to standard case folder 'run_basic'
RuntimeOption(3)=''           ! default
RuntimeOption(4)=''           ! default
BuildDebug=.FALSE.            ! default: display the compelte compilation output
BuildNoDebug=.FALSE.          ! default: don't display the compilation output
BuildSolver=.FALSE.           ! default: no compiling, just pre-built binary execution
DoFullReggie=.FALSE.          ! default: don't run reggie recursively using gitlab-ci.yml file
! Get number of command line arguments and read in runtime option of regressioncheck
nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs.GE.7)THEN
  SWRITE(UNIT_stdOut,'(A)') ' ERROR: too many arguments for regressioncheck!'
  ERROR STOP 2
END IF


DO I=1,nArgs
  CALL GET_COMMAND_ARGUMENT(I,RuntimeOption(I),STATUS=iSTATUS)
  ! check failure 
  IF(iSTATUS.NE.0)THEN
    SWRITE(UNIT_stdOut,'(A)') ' ERROR: wrong argument for regressioncheck!'
    ERROR STOP 2
  END IF
  ! check help output request
  IF((TRIM(RuntimeOption(I)).EQ.'--help').OR.(TRIM(RuntimeOption(I)).EQ.'help').OR.(TRIM(RuntimeOption(I)).EQ.'HELP')) THEN
    CALL Print_Help_Information()
    STOP 0
  END IF
  ! set the options
  SELECT CASE(I)
  CASE(1) ! RuntimeOption 1
    IF(TRIM(RuntimeOption(I)).EQ.'run') THEN
      BuildSolver=.FALSE.
    ELSEIF(TRIM(RuntimeOption(I)).EQ.'build') THEN
      BuildSolver=.TRUE.
    ELSEIF(TRIM(RuntimeOption(I)).EQ.'full') THEN
      BuildSolver=.FALSE.
      DoFullReggie=.TRUE.
    ELSEIF(TRIM(RuntimeOption(I)).EQ.'no-full') THEN ! needed for "regressioncheck full" because 'no-full' is added to the flags
      DoFullReggie=.FALSE.
      RuntimeOption(1)='run'
    ELSEIF(TRIM(RuntimeOption(I)).EQ.'tutorials') THEN 
      BuildSolver=.FALSE.
      RuntimeOption(1)='tutorials'
      RuntimeOption(2)=''
    ELSE
      SWRITE(UNIT_stdOut,'(A)') ' ERROR: wrong argument for regressioncheck!' 
      ERROR STOP 2
    END IF
  CASE DEFAULT ! RuntimeOption 2 - 4
    IF(BuildSolver.EQV..TRUE.)THEN
      IF(TRIM(RuntimeOption(I)).EQ.'build-continue') BuildContinue=.TRUE.
      IF(TRIM(RuntimeOption(I)).EQ.'run-continue') RunContinue=.TRUE.
      IF(TRIM(RuntimeOption(I)).EQ.'debug')THEN
        BuildDebug=.TRUE.             ! e.g. "./regressioncheck debug" or "./regressioncheck build feature_convtest debug"
      ELSEIF(TRIM(RuntimeOption(I)).EQ.'no-debug')THEN
        BuildNoDebug=.TRUE.           ! redirect std- and err-output channels to "/build_reggie/build_reggie.out"
        IF(I.EQ.2)RuntimeOption(2)='run_basic'  ! set to standard case folder 'run_basic'
      END IF
      IF((I.EQ.2).AND.((BuildDebug.EQV..TRUE.).OR.(BuildNoDebug.EQV..TRUE.)))RuntimeOption(I)='run_basic' ! debug uses
                                                                                            !"configuration.reggie" from "run_basic"
    ELSE
      IF(TRIM(RuntimeOption(I)).EQ.'no-debug')THEN ! needed for "regressioncheck full" because 'no-debug' is added to the flags
        IF(I.EQ.2)RuntimeOption(2)='run_basic'  ! set to standard case folder 'run_basic'
      END IF
    END IF
    ! [RuntimeOption] = all: run all example folders
    IF((I.EQ.2).AND.((TRIM(RuntimeOption(I)).EQ.'all').OR.(TRIM(RuntimeOption(I)).EQ.'ALL')))RuntimeOption(I)=''
    ! prevent infinite recursive loops
    IF(TRIM(RuntimeOption(I)).EQ.'no-full') THEN
      DoFullReggie=.FALSE.
      IF(I.EQ.2)RuntimeOption(2)='run_basic'  ! set to standard case folder 'run_basic'
    END IF
  END SELECT
END DO

! get number of threads/procs for parallel building
CALL GetNumberOfProcs(nArgs)

! display the resulting options
SWRITE(UNIT_stdOut,'(A,4(A1,A,A1),A21,I3,A20)')' Running with arguments: ',&
'[',TRIM(RuntimeOption(1)),']',&
'[',TRIM(RuntimeOption(2)),']',&
'[',TRIM(RuntimeOption(3)),']',&
'[',TRIM(RuntimeOption(4)),']',   ' and compiling with [',NumberOfProcs,'] threads (-1 is -j)'
SWRITE(UNIT_stdOut,'(A,L1,A1)')      ' BuildSolver  :  [',BuildSolver,']'
SWRITE(UNIT_stdOut,'(A,L1,A1)')      ' BuildContinue:  [',BuildContinue,']'
SWRITE(UNIT_stdOut,'(A,L1,A1)')      ' RunContinue  :  [',RunContinue,']'
SWRITE(UNIT_stdOut,'(A,L1,A1)')      ' BuildDebug   :  [',BuildDebug,']'
SWRITE(UNIT_stdOut,'(A,L1,A1)')      ' BuildNoDebug :  [',BuildNoDebug,']'
SWRITE(UNIT_stdOut,'(A,L1,A1)')      ' DoFullReggie :  [',DoFullReggie,']'

END SUBROUTINE GetCommandLineOption


!==================================================================================================================================
!> Check if examples exist. Next, scan the folder for all examples. The routine returns the number of examples, their name
!> and nullifies the parameter entries for each example
!> If the regression check is run within an example folder, only said example is executed
!==================================================================================================================================
SUBROUTINE GetExampleList()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: nExamples,ExampleNames,Examples,ExamplesDir,BuildDir,RuntimeOption
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLESzz
CHARACTER(LEN=500)            :: SYSCOMMAND           ! string to fit the system command
CHARACTER(LEN=500)            :: SYSCOMMANDTWO        ! string to fit the system command
CHARACTER(LEN=255)            :: FilePathName,FileName! file and path strings
INTEGER                       :: ioUnit               ! io-unit
INTEGER                       :: iSTATUS              ! status
INTEGER                       :: iExample             ! loop index for example
CHARACTER(len=255)            :: cwd                  ! current cworking directory CALL getcwd(cwd)
LOGICAL                       :: ExistFile            ! T=file exists, F=file does not exist
!==================================================================================================================================
iSTATUS=0 ! nullify
! check if regressioncheck is executed von /bin/regressioncheck directory, otherwise use the current working directory
CALL getcwd(cwd)                                      ! get path of current directory
FileName=cwd(INDEX(cwd,'/',BACK = .TRUE.)+1:LEN(cwd)) ! cut the last folder name from current diectory path
FilePathName=TRIM(cwd)//'/parameter_reggie.ini'       ! check if parameter_reggie.ini is located within cwd
INQUIRE(File=FilePathName,EXIST=ExistFile)            ! inquire

IF(ExistFile.EQV..FALSE.)THEN ! use existing example folder
  IF(RuntimeOption(1).EQ."tutorials")THEN
    ExamplesDir=TRIM(BASEDIR(2:LEN(BASEDIR)-1))//'../tutorials/'
  ELSE
    ExamplesDir=TRIM(BASEDIR(2:LEN(BASEDIR)-1))//'../regressioncheck/examples/'
  END IF
  SYSCOMMAND='cd '//TRIM(ExamplesDir)//' && ls -d */ > tmp.txt'
ELSE ! run regressioncheck for a single folder located anywhere from which the reggie is executed
  ExamplesDir='./../'
  SYSCOMMAND='cd '//TRIM(ExamplesDir)//' && ls -d '//TRIM(FileName)//'/ > tmp.txt'
  RuntimeOption(2)=TRIM(FileName) ! override RuntimeOption in order to select only this directory
END IF
BuildDir=TRIM(BASEDIR(2:LEN(BASEDIR)-1))! use basedir because one cannot use: TRIM(cwd)//'/'
                                        ! because the checked out code source
                                        ! code is needed for building new binaries

! sanity check: change directory into the example folder, if it does not exist this check fails
SYSCOMMANDTWO='cd '//TRIM(ExamplesDir)
CALL EXECUTE_COMMAND_LINE(SYSCOMMANDTWO, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Error: Example folder does not work!'
  ERROR STOP 66
END IF

! get number of examples by complicated fortran hack:
! ls is called and the output is piped into the file tmp.txt. The number of lines is the number of available examples.
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Could not create tmp.txt to get number of examples'
  ERROR STOP 99
END IF

! read tmp.txt | list of directories if regressioncheck/examples
FileName=TRIM(ExamplesDir)//'tmp.txt'
OPEN(NEWUNIT = ioUnit, FILE = FileName, STATUS ="OLD", IOSTAT = iSTATUS )

nExamples=0
DO
  READ(ioUnit,FMT='(A)',IOSTAT=iSTATUS) FileName
  IF (iSTATUS.NE.0) EXIT
  nExamples=nExamples+1
END DO
SWRITE(UNIT_stdOut,'(A,I3)')  ' Number of Examples: ', nExamples

! read in the directory name for each example and initialization of default values a.k.a. nullify
ALLOCATE(ExampleNames(nExamples))
ALLOCATE(Examples(nExamples))
REWIND(ioUnit)
DO iExample=1,nExamples
  READ(ioUnit,FMT='(A)') ExampleNames(iExample)
  SWRITE(UNIT_stdOut,'(A,I3.3,3x,A)')  ' Example-',iExample, ExampleNames(iExample)
  ! fill PATH of examples
  Examples(iExample)%PATH = TRIM(ExamplesDir)//TRIM(ExampleNames(iExample))
  Examples(iExample)%ReferenceFile=''
  Examples(iExample)%ReferenceNormFile=''
  Examples(iExample)%H5DIFFCheckedStateFile=''
  Examples(iExample)%H5DIFFReferenceStateFile=''
  Examples(iExample)%H5DIFFReferenceDataSetName=''
  Examples(iExample)%H5diffToleranceType='absolute'
  Examples(iExample)%H5diffTolerance=-1.
  Examples(iExample)%RestartFileName=''
  Examples(iExample)%ErrorStatus=0
END DO
CLOSE(ioUnit)

! and remove tmp.txt || cleanup of ls
! clean tmp.txt
SYSCOMMAND='cd '//TRIM(ExamplesDir)//' && rm tmp.txt'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Could not remove tmp.txt'
  ERROR STOP 99
END IF

END SUBROUTINE GetExampleList


!==================================================================================================================================
!> Read the parameter_reggie.ini file of an example given by its relative path. It
!> contains information for the computation of the example:
!>  MPI - a mpi or serial example
!>  optional reference files for error-norms, reference state file and tested dataset and name of the checked state file
!>  optional a restart filename
!==================================================================================================================================
SUBROUTINE InitExample(FilePath,Example,SkipExample)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: tExample,readRHS
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)               :: FilePath
TYPE(tExample),INTENT(INOUT)              :: Example
LOGICAL,INTENT(OUT)                       :: SkipExample
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ioUnit
INTEGER                                   :: iSTATUS,IndNum,IndNum2,IndNum3,MaxNum
CHARACTER(LEN=255)                        :: FileName,temp1,temp2
LOGICAL                                   :: ExistFile
!==================================================================================================================================
SkipExample=.FALSE.
! test if file exists and open
FileName=TRIM(FilePath)//'parameter_reggie.ini'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SkipExample=.TRUE.
  SWRITE(UNIT_stdOut,'(A16,A,A1)') '   FileName  : [', TRIM(FileName),']'
  SWRITE(UNIT_stdOut,'(A16,L,A1)') '   ExistFile : [', ExistFile,']'
  SWRITE(UNIT_stdOut,'(A16,A)')    '   ERROR     : ','no parameter_reggie.ini found.                        ...skipping'
  RETURN
ELSE
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF

! init ALWAYS
Example%MPIrun                  = .FALSE. ! don't use "mpirun" n default
Example%MPIcommand              = 'mpirun'! use mpirun for running parallel simulations as default
Example%MPIthreadsStr           = '1'     ! run with 1 MPI thread on default
Example%MPIthreadsN             = 1       ! minimum
Example%nRuns                   = 1       ! minimum
Example%nVar                    = 0
Example%ReferenceTolerance      = -1.
Example%SubExample              = '-'     ! init
Example%SubExampleNumber        = 0       ! init total number of subexamples
Example%SubExampleOption(1:100) = '-'     ! default option is nothing
Example%IntegrateLine           = .FALSE. ! nullify
Example%IntegrateLineMultiplier = 1.0     ! default option for IntegrateLine
Example%IntegrateLineOption     = 'default' ! default value is "default"
Example%ConvergenceTest         = .FALSE. ! nullify
Example%CompareDatafileRow      = .FALSE. ! nullify
Example%CompareHDF5ArrayBounds  = .FALSE. ! nullify
DO ! extract reggie information
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1 ! get first line assuming it is something like "nVar= 5"
  IF(iSTATUS.EQ.-1) EXIT ! end of file (EOF) reached
  IF(INDEX(temp1,'!').GT.0)temp1=temp1(1:INDEX(temp1,'!')-1) ! if temp1 contains a "!", remove it and the following characters
  temp1=TRIM(ADJUSTL(temp1)) ! remove trailing and leading white spaces
  temp2=''
  IndNum=INDEX(temp1,'=')
  MaxNum=LEN(TRIM(temp1))
  IF(IndNum.GT.0)THEN ! definition found
    readRHS(1)=temp1(1:IndNum-1)      ! parameter name
    readRHS(2)=temp1(IndNum+1:MaxNum) ! parameter setting
    ! number of DG cells in one direction
    IF(TRIM(readRHS(1)).EQ.'NumberOfCells')CALL GetParameterList(ParameterList   = Example%NumberOfCellsStr,&
                                                                 nParameters     = 100,               &
                                                                 ParameterNumber = Example%NumberOfCellsN)
    ! get size of EQNSYS (deprecated)
    IF(TRIM(readRHS(1)).EQ.'nVar')CALL str2int(readRHS(2),Example%nVar,iSTATUS)
    ! single or parallel
    IF(TRIM(readRHS(1)).EQ.'MPIrun')    CALL str2logical(readRHS(2),Example%MPIrun,iSTATUS) ! True/False
    IF(TRIM(readRHS(1)).EQ.'MPIcommand')            Example%MPIcommand        =TRIM(ADJUSTL(readRHS(2)))
    !IF(TRIM(readRHS(1)).EQ.'MPIthreads')CALL str2int(    readRHS(2),Example%MPIthreadsStr,iSTATUS)!number of threads
    IF(TRIM(readRHS(1)).EQ.'MPIthreads')CALL GetParameterList(ParameterList   = Example%MPIthreadsStr,&
                                                              nParameters     = 100,               &
                                                              ParameterNumber = Example%MPIthreadsN)
    ! number of scaling runs
    IF(TRIM(readRHS(1)).EQ.'nRuns')CALL str2int(readRHS(2),Example%nRuns,iSTATUS)
    ! Reference Norm/State
    IF(TRIM(readRHS(1)).EQ.'ReferenceTolerance')CALL str2real(readRHS(2),Example%ReferenceTolerance,iSTATUS)
    IF(TRIM(readRHS(1)).EQ.'ReferenceFile')                Example%ReferenceFile         =TRIM(ADJUSTL(readRHS(2)))
    IF(TRIM(readRHS(1)).EQ.'ReferenceNormFile')            Example%ReferenceNormFile     =TRIM(ADJUSTL(readRHS(2)))
    IF(TRIM(readRHS(1)).EQ.'H5DIFFReferenceStateFile')           Example%H5DIFFReferenceStateFile    =TRIM(ADJUSTL(readRHS(2)))
    IF(TRIM(readRHS(1)).EQ.'H5DIFFCheckedStateFile')       Example%H5DIFFCheckedStateFile      =TRIM(ADJUSTL(readRHS(2)))
    IF(TRIM(readRHS(1)).EQ.'H5DIFFReferenceDataSetName')   Example%H5DIFFReferenceDataSetName  =TRIM(ADJUSTL(readRHS(2)))
    IF(TRIM(readRHS(1)).EQ.'H5diffToleranceType')          Example%H5diffToleranceType   =TRIM(ADJUSTL(readRHS(2)))
    IF(TRIM(readRHS(1)).EQ.'H5diffTolerance')   CALL str2real(readRHS(2),Example%H5diffTolerance,iSTATUS)
    IF(TRIM(readRHS(1)).EQ.'RestartFileName')              Example%RestartFileName       =TRIM(ADJUSTL(readRHS(2)))
    ! test if restart file exists
    IF(TRIM(Example%RestartFileName).NE.'')THEN
      FileName=TRIM(FilePath)//TRIM(Example%RestartFileName)
      INQUIRE(File=FileName,EXIST=ExistFile)
      IF(.NOT.ExistFile) THEN
        SkipExample=.TRUE.
        SWRITE(UNIT_stdOut,'(A16,A,A1)') '   FileName  : [', TRIM(FileName),']'
        SWRITE(UNIT_stdOut,'(A16,L,A1)') '   ExistFile : [', ExistFile,']'
        ERROR STOP 'Restart file supplied, but not found.'
      END IF
    END IF
    ! SubExamples - currently one subexample class is allowed with multiple options
    IF(TRIM(readRHS(1)).EQ.'SubExample') CALL GetParameterList(ParameterName   = Example%SubExample,       &
                                                               ParameterList   = Example%SubExampleOption, &
                                                               nParameters     = 100,                       &
                                                               ParameterNumber = Example%SubExampleNumber)
    ! Line integration (e.g. integrate a property over time, the data is read from a .csv or .dat file)
    ! e.g. in parameter_reggie.ini:
    ! 
    ! IntegrateLine= Database.csv   , 1           , ','      , 1:5       , 44.00         , 2.5e-2   , DivideByTimeStep , 5.340e-03 
    !                data file name , header lines, delimiter, colums x:y, integral value, tolerance, option           , multiplier
    IF(TRIM(readRHS(1)).EQ.'IntegrateLine')THEN
       Example%IntegrateLine            = .TRUE.
       Example%IntegrateLineRange(1:2)  = 0         ! init
       Example%IntegrateLineHeaderLines = 0         ! init
       Example%IntegrateLineDelimiter   = '999'     ! init
       Example%IntegrateLineValue       = 0         ! init
       Example%IntegrateLineTolerance   = 5E-2      ! init 5% tolerance
       Example%IntegrateLineOption      = 'default' ! init 
       Example%IntegrateLineMultiplier  = 1.0       ! init: 1.0 changed nothing
       IndNum2=INDEX(readRHS(2),',')
       IF(IndNum2.GT.0)THEN ! get the name of the data file
         temp2                     = readRHS(2)
         Example%IntegrateLineFile = TRIM(ADJUSTL(temp2(1:IndNum2-1))) ! data file name
         temp2                     = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
         IndNum2                   = INDEX(temp2,',')

         IF(IndNum2.GT.0)THEN ! get number of header lines in data file (they are ignored on reading the file)
           CALL str2int(temp2(1:IndNum2-1),Example%IntegrateLineHeaderLines,iSTATUS)
           temp2                   = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
           IndNum2                 = INDEX(temp2,"'")
           IF(IndNum2.GT.0)THEN ! get delimiter for separating the columns in the data file
             IndNum3=INDEX(temp2(IndNum2+1:LEN(TRIM(temp2))),"'")+IndNum2
             Example%IntegrateLineDelimiter=temp2(IndNum2+1:IndNum3-1)
             temp2                 = temp2(IndNum3+1:LEN(TRIM(temp2))) ! next
             IndNum2               = INDEX(temp2,',')
             IF(IndNum2.GT.0)THEN ! get column ranges
               temp2               = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
               IndNum2             = INDEX(temp2,',')
               IndNum3=INDEX(temp2(1:IndNum2),':')
               IF(IndNum3.GT.0)THEN ! check range
                 CALL str2int(temp2(1        :IndNum3-1),Example%IntegrateLineRange(1),iSTATUS) ! column number 1
                 CALL str2int(temp2(IndNum3+1:IndNum2-1),Example%IntegrateLineRange(2),iSTATUS) ! column number 2
                 temp2             = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
                 IndNum2           = INDEX(temp2,',')
                 IF(IndNum2.GT.0)THEN ! get integral value
                   CALL str2real(temp2(1:IndNum2-1),Example%IntegrateLineValue,iSTATUS)
                   temp2           = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
                   IndNum2         = INDEX(temp2,',')
                   IF(IndNum2.GT.0)THEN ! get tolerance - only if supplied, if not then use the default value (see top)
                     CALL str2real(temp2(1:IndNum2-1),Example%IntegrateLineTolerance,iSTATUS)
                     temp2           = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
                     IndNum2         = INDEX(temp2,',')
                     IF(IndNum2.GT.0)THEN ! get option
                       Example%IntegrateLineOption=TRIM(temp2(1:IndNum2-1))
                       IF(temp2(IndNum2+1:LEN(TRIM(temp2))).NE.'')THEN ! get multiplier
                         CALL str2real(temp2(IndNum2+1:LEN(TRIM(temp2))),Example%IntegrateLineMultiplier,iSTATUS)
                       END IF ! get multiplier
                     END IF ! get option
                   ELSE ! try and get tolerance
                     CALL str2real(temp2(1:LEN(TRIM(temp2))),Example%IntegrateLineTolerance,iSTATUS)
                   END IF ! get tolerance
                 END IF ! get integral value
               END IF ! check range
             END IF ! get column ranges
           END IF ! get delimiter
         END IF !  get number of header lines
       END IF ! get file name
      IF(ANY( (/iSTATUS.NE.0                                       ,&
                Example%IntegrateLineRange(1:2).EQ.0               ,&
                Example%IntegrateLineFile.EQ.''                    ,&
                Example%IntegrateLineDelimiter(1:3).EQ.'999' /)     &
       )) THEN
         SWRITE(UNIT_stdOut,'(A,I5)')      ' iSTATUS                         : ',iSTATUS
         SWRITE(UNIT_stdOut,'(A,I5)')      ' Example%IntegrateLineRange(1)   : ',Example%IntegrateLineRange(1)
         SWRITE(UNIT_stdOut,'(A,I5)')      ' Example%IntegrateLineRange(2)   : ',Example%IntegrateLineRange(2)
         SWRITE(UNIT_stdOut,'(A,A)')       ' Example%IntegrateLineFile       : ',TRIM(Example%IntegrateLineFile)
         SWRITE(UNIT_stdOut,'(A,A)')       ' Example%IntegrateLineDelimiter  : ',TRIM(Example%IntegrateLineDelimiter)
         SWRITE(UNIT_stdOut,'(A,I5)')      ' Example%IntegrateLineHeaderLines: ',Example%IntegrateLineHeaderLines
         SWRITE(UNIT_stdOut,'(A,E25.14E3)')' Example%IntegrateLineTolerance  : ',Example%IntegrateLineTolerance
         SWRITE(UNIT_stdOut,'(A,E25.14E3)')' Example%IntegrateLineValue      : ',Example%IntegrateLineValue
         SWRITE(UNIT_stdOut,'(A,A)')       ' Example%IntegrateLineOption     : ',TRIM(Example%IntegrateLineOption)
         SWRITE(UNIT_stdOut,'(A,E25.14E3)')' Example%IntegrateLineMultiplier : ',Example%IntegrateLineMultiplier
         SWRITE(UNIT_stdOut,'(A)')         ' Setting Example%IntegrateLine=.FALSE.'
         Example%IntegrateLine=.FALSE.
       END IF
    END IF ! 'IntegrateLine'
    ! Line comparison (compare one complete line in, e.g., a *.csv or *.dat file)
    ! in parameter_reggie.ini define:
    !            CompareDatafileRow  =  Database.csv   ,Database.csv_ref, 1e-3      , 1            , ','       , last
    !                                   data file name ,ref data file   , Tolerance , header lines , delimiter , line number
    IF(TRIM(readRHS(1)).EQ.'CompareDatafileRow')THEN
       Example%CompareDatafileRow            = .TRUE.
       Example%CompareDatafileRowHeaderLines = 0       ! init
       Example%CompareDatafileRowTolerance   = -1      ! init
       Example%CompareDatafileRowReadHeader  = .FALSE. ! init
       Example%CompareDatafileRowDelimiter   = '999'   ! init
       IndNum2=INDEX(readRHS(2),',')
       IF(IndNum2.GT.0)THEN ! get the name of the data file
         temp2                     = readRHS(2)
         Example%CompareDatafileRowFile      = TRIM(ADJUSTL(temp2(1:IndNum2-1))) ! data file name
         temp2                     = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
         IndNum2                   = INDEX(temp2,',')
         IF(IndNum2.GT.0)THEN ! get the name of the reference data file
           Example%CompareDatafileRowRefFile = TRIM(ADJUSTL(temp2(1:IndNum2-1))) ! ref data file name
           temp2                   = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
           IndNum2                 = INDEX(temp2,',')
           IF(IndNum2.GT.0)THEN ! get tolerance value for comparison
             CALL str2real(temp2(1:IndNum2-1),Example%CompareDatafileRowTolerance,iSTATUS)
             temp2                 = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
             IndNum2               = INDEX(temp2,',')
             IF(IndNum2.GT.0)THEN ! get number of header lines in data file (they are ignored on reading the file)
               CALL str2int(temp2(1:IndNum2-1),Example%CompareDatafileRowHeaderLines,iSTATUS)
               temp2               = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
               IndNum2             = INDEX(temp2,"'")
               IF(IndNum2.GT.0)THEN ! get delimiter for separating the columns in the data file
                 IndNum3=INDEX(temp2(IndNum2+1:LEN(TRIM(temp2))),"'")+IndNum2
                 Example%CompareDatafileRowDelimiter=temp2(IndNum2+1:IndNum3-1)
                 temp2             = temp2(IndNum3+1:LEN(TRIM(temp2))) ! next
                 IndNum2           = INDEX(temp2,',')
                 IF(IndNum2.GT.0)THEN ! get row number
                   temp2           = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
                   IF(ADJUSTL(TRIM(temp2)).EQ.'last')THEN ! use the 'last' line number in file for comparison
                     Example%CompareDatafileRowNumber=HUGE(1)
                     iSTATUS=0
                   ELSE
                     CALL str2int(temp2,Example%CompareDatafileRowNumber,iSTATUS)
                   END IF
                 END IF ! get row number
               END IF ! get delimiter
             END IF !  get number of header lines
           END IF ! get tolerance value for comparison
         END IF ! get the name of the reference data file
       END IF ! get file name
      ! set CompareDatafileRow to false if any of the following cases is true
      IF(ANY( (/iSTATUS.NE.0                                       ,&
                Example%CompareDatafileRowFile.EQ.''               ,&
                Example%CompareDatafileRowRefFile.EQ.''            ,&
                Example%CompareDatafileRowTolerance.LT.0.          ,&
                Example%CompareDatafileRowHeaderLines.EQ.0         ,&
                Example%CompareDatafileRowDelimiter(1:3).EQ.'999'/) )) THEN
         SWRITE(UNIT_stdOut,'(A,I5)')      ' iSTATUS                                : ',iSTATUS
         SWRITE(UNIT_stdOut,'(A,A)')       ' Example%CompareDatafileRowFile         : ',TRIM(Example%CompareDatafileRowFile)
         SWRITE(UNIT_stdOut,'(A,A)')       ' Example%CompareDatafileRowRefFile      : ',TRIM(Example%CompareDatafileRowRefFile)
         SWRITE(UNIT_stdOut,'(A,E25.14E3)')' Example%CompareDatafileRowTolerance    : ',Example%CompareDatafileRowTolerance
         SWRITE(UNIT_stdOut,'(A,I5)')      ' Example%CompareDatafileRowHeaderLines  : ',Example%CompareDatafileRowHeaderLines
         SWRITE(UNIT_stdOut,'(A,A)')       ' Example%CompareDatafileRowDelimiter    : ',TRIM(Example%CompareDatafileRowDelimiter)
         SWRITE(UNIT_stdOut,'(A)')         ' Setting Example%CompareDatafileRow=.FALSE.'
         Example%CompareDatafileRow=.FALSE.
       ELSE
         IF(Example%CompareDatafileRowHeaderLines.GT.0)Example%CompareDatafileRowReadHeader=.TRUE.
       END IF
    END IF ! 'CompareDatafileRow'

    ! Convergence Test (h- or p-convergence for increasing the number of cells or increasing the polynomial degree N)
    ! in parameter_reggie.ini define:
    !  for p: ConvergenceTest =       p     ,                   IntegrateLine                      , 0.12434232          , 1e-2
    !  for h: ConvergenceTest =       h     ,                   Constant                           , 3.99                , 1e-2
    !                          type (h or p), comparison type (IntegrateLine or power law exponent), value for comparison, Tolerance
    IF(TRIM(readRHS(1)).EQ.'ConvergenceTest')THEN
      Example%ConvergenceTest           = .TRUE.
      Example%ConvergenceTestType       = ''     ! init
      Example%ConvergenceTestDomainSize = -999.0 ! init
      Example%ConvergenceTestValue      = -999.0 ! init
      Example%ConvergenceTestTolerance  = -1.     ! init
      Example%ConvergenceTestSuccessRate= 0.5    ! default: 50 % of tests of all nVar must hold
      IndNum2=INDEX(readRHS(2),',')
      IF(IndNum2.GT.0)THEN ! get the type of the convergence test (h- or p-convergence)
        temp2 = readRHS(2)
        Example%ConvergenceTestType= TRIM(ADJUSTL(temp2(1:IndNum2-1))) ! type
        temp2 = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
        IndNum2 = INDEX(temp2,',')
        IF(IndNum2.GT.0)THEN ! get the size of the domain
          CALL str2real(temp2(1:IndNum2-1),Example%ConvergenceTestDomainSize,iSTATUS)
          temp2 = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
          IndNum2 = INDEX(temp2,',')
          IF((IndNum2.GT.0).AND.(iSTATUS.EQ.0))THEN ! get value for comparison
            CALL str2real(temp2(1:IndNum2-1),Example%ConvergenceTestValue,iSTATUS)
            temp2 = TRIM(ADJUSTL(temp2(IndNum2+1:LEN(TRIM(temp2))))) ! next
            IndNum2 = INDEX(temp2,',')
            IF(IndNum2.LE.0) IndNum2 = LEN(temp2)+1 ! no SuccessRate given -> use default
            IF((IndNum2.GT.0).AND.(iSTATUS.EQ.0))THEN ! get tolerance value for comparison
              CALL str2real(temp2(1:IndNum2-1),Example%ConvergenceTestTolerance,iSTATUS)
              IF((IndNum2.LT.LEN(TRIM(temp2))).AND.(iSTATUS.EQ.0))THEN
                temp2 = TRIM(ADJUSTL(temp2(IndNum2+1:LEN(TRIM(temp2))))) ! next
                CALL str2real(TRIM(ADJUSTL(temp2)),Example%ConvergenceTestSuccessRate,iSTATUS)
              END IF
            END IF ! get tolerance value for comparison
          END IF ! get value for comparison
        END IF ! get the comparison type
      END IF ! get the type of the convergence test (h- or p-convergence)
      ! set ConvergenceTest to false if any of the following cases is true
      IF(ANY( (/iSTATUS.NE.0                             ,&
                Example%ConvergenceTestType.EQ.''        ,&
                Example%ConvergenceTestDomainSize.LT.-1. ,&
                Example%ConvergenceTestValue.LT.-1.      ,&
                Example%ConvergenceTestTolerance.EQ.-1.  ,&
                Example%ConvergenceTestSuccessRate.LT.0. ,&
                Example%ConvergenceTestSuccessRate.GT.1.  &
                /) )) Example%ConvergenceTest=.FALSE.
      SWRITE(UNIT_stdOut,'(A)') ' ERROR in Convergence-Test Read-in. Deactivating Convergence-Test!'
      IF(Example%ConvergenceTest.EQV..FALSE.)THEN
        SWRITE(UNIT_stdOut,'(A,A)')        ' Example%ConvergenceTestType       : ',Example%ConvergenceTestType
        SWRITE(UNIT_stdOut,'(A,E25.14E3)') ' Example%ConvergenceTestDomainSize : ',Example%ConvergenceTestDomainSize
        SWRITE(UNIT_stdOut,'(A,E25.14E3)') ' Example%ConvergenceTestValue      : ',Example%ConvergenceTestValue
        SWRITE(UNIT_stdOut,'(A,E25.14E3)') ' Example%ConvergenceTestTolerance  : ',Example%ConvergenceTestTolerance
        SWRITE(UNIT_stdOut,'(A,E25.14E3)') ' Example%ConvergenceTestSuccessRate: ',Example%ConvergenceTestSuccessRate
      END IF
    END IF ! 'ConvergenceTest'
    ! Check the bounds of an array in a HDF5 file, if they are outside the supplied ranges -> fail
    IF(TRIM(readRHS(1)).EQ.'CompareHDF5ArrayBounds')THEN
      Example%CompareHDF5ArrayBounds           = .TRUE. ! read an array from a HDF5 file and compare certain entry
      Example%CompareHDF5ArrayBoundsValue(1:2) = 0.     ! value ranges for comparison
      Example%CompareHDF5ArrayBoundsRange(1:2) = -1     ! HDF5 array dim ranges
      Example%CompareHDF5ArrayBoundsName       = '-1'   ! array name in HDF5 file
      Example%CompareHDF5ArrayBoundsFile       = '-1'   ! name of HDF5 file
      IndNum2=INDEX(readRHS(2),',')
      IF(IndNum2.GT.0)THEN ! get name of array in HDF5 file
        temp2                              = readRHS(2)
        Example%CompareHDF5ArrayBoundsFile = TRIM(ADJUSTL(temp2(1:IndNum2-1))) ! type
        temp2                              = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
        IndNum2                            = INDEX(temp2,',')
        IF(IndNum2.GT.0)THEN ! get name of array in HDF5 file
          Example%CompareHDF5ArrayBoundsName = TRIM(ADJUSTL(temp2(1:IndNum2-1))) ! type
          temp2                              = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
          IndNum2                            = INDEX(temp2,',')
          IF(IndNum2.GT.0)THEN ! HDF5 array dim ranges
            IndNum2             = INDEX(temp2,',')
            IndNum3=INDEX(temp2(1:IndNum2),':')
            IF(IndNum3.GT.0)THEN ! check range
              CALL str2int(temp2(1        :IndNum3-1),Example%CompareHDF5ArrayBoundsRange(1),iSTATUS) ! column number 1
              CALL str2int(temp2(IndNum3+1:IndNum2-1),Example%CompareHDF5ArrayBoundsRange(2),iSTATUS) ! column number 2
              temp2             = temp2(IndNum2+1:LEN(TRIM(temp2))) ! next
              IndNum2           = LEN(TRIM(temp2))
              IF(IndNum2.GT.0)THEN ! value ranges for comparison
                IndNum2           = LEN(temp2)
                IndNum3=INDEX(temp2(1:IndNum2),':')
                IF(IndNum3.GT.0)THEN ! check range
                  CALL str2real(temp2(1        :IndNum3-1),Example%CompareHDF5ArrayBoundsValue(1),iSTATUS) ! column number 1
                  CALL str2real(temp2(IndNum3+1:IndNum2-1),Example%CompareHDF5ArrayBoundsValue(2),iSTATUS) ! column number 2
                END IF ! check range
              END IF ! value ranges for comparison
            END IF ! check range
          END IF ! HDF5 array dim ranges
        END IF ! get name of array in HDF5 file
      END IF ! get name of HDF5 file
      ! set "CompareHDF5ArrayBounds" to false if any of the following cases is true
      IF(ANY( (/iSTATUS.NE.0                                                                           ,&
                ANY(Example%CompareHDF5ArrayBoundsRange(1:2).EQ.-1)                                    ,&
                    Example%CompareHDF5ArrayBoundsName.EQ.'-1'                                         ,&
                    Example%CompareHDF5ArrayBoundsFile.EQ.'-1'                                         ,&
                    Example%CompareHDF5ArrayBoundsValue(1).GT.Example%CompareHDF5ArrayBoundsValue(2)   ,&
                    Example%CompareHDF5ArrayBoundsRange(1).GT.Example%CompareHDF5ArrayBoundsRange(2)  /)&
                  ))Example%CompareHDF5ArrayBounds=.FALSE.
      IF(Example%CompareHDF5ArrayBounds.EQV..FALSE.)THEN
        SWRITE(UNIT_stdOut,'(A,E25.14,A)') 'Example%CompareHDF5ArrayBoundsValue(1) : '&
                                           ,Example%CompareHDF5ArrayBoundsValue(1),' (lower)'
        SWRITE(UNIT_stdOut,'(A,E25.14,A)') 'Example%CompareHDF5ArrayBoundsValue(2) : '&
                                           ,Example%CompareHDF5ArrayBoundsValue(2),' (upper)'
        SWRITE(UNIT_stdOut,'(A,I6,A)') 'Example%CompareHDF5ArrayBoundsRange(1) : '&
                                           ,Example%CompareHDF5ArrayBoundsRange(1),' (lower)'
        SWRITE(UNIT_stdOut,'(A,I6,A)') 'Example%CompareHDF5ArrayBoundsRange(2) : '&
                                           ,Example%CompareHDF5ArrayBoundsRange(2),' (upper)'
        SWRITE(UNIT_stdOut,'(A,A)')        'Example%CompareHDF5ArrayBoundsName     : ',Example%CompareHDF5ArrayBoundsName
        SWRITE(UNIT_stdOut,'(A,A)')        'Example%CompareHDF5ArrayBoundsFile     : ',Example%CompareHDF5ArrayBoundsFile
      END IF
    END IF ! 'CompareHDF5ArrayBounds'
    ! Next feature
    !IF(TRIM(readRHS(1)).EQ.'NextFeature')
  END IF ! IndNum.GT.0 -> definition found
END DO
CLOSE(ioUnit)

END SUBROUTINE InitExample


!==================================================================================================================================
!> reads a list of parameters from parameter_reggie: Option=Name{parameter1,parameter2,parameter3,.....,parameterN}
!> e.g. SubExample=TimeDiscMethod{standardrk3-3,niegemannrk4-14,toulorgerk4-8c}
!>                 ParameterName   = TimeDiscMethod
!>                 ParameterList   = [standardrk3-3,niegemannrk4-14,toulorgerk4-8c]
!>                 ParameterNumber = 3
!==================================================================================================================================
SUBROUTINE GetParameterList(ParameterName,ParameterList,nParameters,ParameterNumber)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: tExample,readRHS
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                       :: nParameters
CHARACTER(LEN=*),INTENT(INOUT)           :: ParameterList(1:nParameters)
CHARACTER(LEN=*),INTENT(INOUT),OPTIONAL  :: ParameterName
INTEGER,INTENT(INOUT),OPTIONAL           :: ParameterNumber
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: IndNum2,MaxNum
INTEGER                                   :: iParameter
CHARACTER(LEN=255)                        :: temp1
!==================================================================================================================================
temp1=readRHS(2)
MaxNum=LEN(TRIM(temp1))
IndNum2=INDEX(temp1(1:MaxNum),'{') ! check for bracket "{" -> indicates a list, otherwise only single argument is expected
IF(IndNum2.GT.0)THEN
  IF(PRESENT(ParameterName))THEN
    ParameterName  = TRIM(ADJUSTL(temp1(1:IndNum2-1))) ! e.g. SubExample=[TimeDiscMethod]{sta...
  END IF
  temp1=TRIM(ADJUSTL(temp1(IndNum2+1:MaxNum))) ! e.g. SubExample=TimeDiscMethod{[sta...  ...enterrk4-5}]
  IndNum2=INDEX(temp1(1:MaxNum),'}')
  IF(IndNum2.GT.0)THEN
  ParameterList=''  ! set default
  ParameterNumber=0 ! set default
    temp1=temp1(1:IndNum2-1) ! e.g. SubExample=TimeDiscMethod{[sta...  ...enterrk4-5]}
    DO iParameter=1,19
      IndNum2=INDEX(temp1(1:LEN(TRIM(temp1))),',')
      ParameterNumber = ParameterNumber + 1 ! get total number of parameters in list {1,2,3,4,5,....,N}
      IF(IndNum2.GT.0)THEN
        ParameterList(iParameter)=TRIM(ADJUSTL(temp1(1:IndNum2-1)))
        temp1=temp1(IndNum2+1:LEN(TRIM(temp1)))
      ELSE
        ParameterList(iParameter)=TRIM(ADJUSTL(temp1(1:LEN(TRIM(temp1)))))
        EXIT
      END IF
    END DO
  ELSE
    CALL abort(&
      __STAMP__&
      ,'GetParameterList: expecting closing bracket "}" in string')
  END IF
ELSE ! no opening bracket "{" found, single argument is expected
  IndNum2=INDEX(temp1(1:MaxNum),',') ! look for "," which would indicate an error
  IF(IndNum2.GT.0)THEN
    CALL abort(&
      __STAMP__&
      ,'GetParameterList(): found more than one parameter argument -> expecting brackets "{ }" in string')
  ELSE ! single parameter argument expected
    ParameterList(1)=TRIM(ADJUSTL(temp1(IndNum2+1:MaxNum)))
    ParameterList(2:nParameters)=''
    IF(PRESENT(ParameterNumber))ParameterNumber=1
  END IF
END IF
END SUBROUTINE GetParameterList



!==================================================================================================================================
!>  Check if executable exists
!==================================================================================================================================
SUBROUTINE CheckForExecutable(Mode)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: EXECPATH,CodeNameLowCase
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: Mode             ! which mode (1 or 2)
                                                  ! 1: pre-compiled binary
                                                  ! 2: code compiled via reggie with defined "EXECPATH", e.g.,
                                                  ! "~/code/code/build_reggie.dev/build_reggie/bin/code"
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: ExistSolver      ! logical to flag solver
CHARACTER(len=255)            :: cwd             ! current cworking directory CALL getcwd(cwd)
!===================================================================================================================================
IF(Mode.EQ.1) THEN
  ! 1.) check if binary is located at the position of the regressioncheck binary (e.g. in a specific case folder)
  CALL getcwd(cwd) ! get path of current directory
  EXECPATH=TRIM(cwd)//'/'//CodeNameLowCase
  INQUIRE(File=EXECPATH,EXIST=ExistSolver)
  IF(.NOT.ExistSolver) THEN
    ! 2.) binary not found in location from where the reggie was executed -> check standard directory instead
    EXECPATH=BASEDIR(2:LEN(BASEDIR)-1)//'bin/'//CodeNameLowCase
  END IF
END IF
! EXECPATH has been set, inquire the existence of the binary
INQUIRE(File=EXECPATH,EXIST=ExistSolver)
IF(.NOT.ExistSolver) THEN
  SWRITE(UNIT_stdOut,'(A38,I1,A)') ' CALL CheckForExecutable() with Mode=',&
                                       Mode,': no executable found. Error during compilation or linking?'
  SWRITE(UNIT_stdOut,'(A41,A)')                         ' EXECPATH: ', TRIM(EXECPATH)
  SWRITE(UNIT_stdOut,'(A41,L)')                      ' ExistSolver: ', ExistSolver
  ERROR STOP 77
END IF
END SUBROUTINE CheckForExecutable


!==================================================================================================================================
!>  Check if configuration.cmake exists and set path
!==================================================================================================================================
SUBROUTINE GetConfigurationFile(Proposal)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: EXECPATH,CodeNameLowCase,RunContinue
USE MOD_RegressionCheck_Vars,  ONLY: BuildSolver,configuration_cmake
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN),OPTIONAL  :: Proposal
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: ExistFile      ! logical to flag file
!===================================================================================================================================
IF(PRESENT(Proposal))THEN ! use path supplied by user
  !configuration_cmake=TRIM(BuildDir)//'build_reggie/bin/configuration.cmake'
  configuration_cmake=TRIM(Proposal)
ELSE
  IF(BuildSolver)THEN ! use pre-existing executables and cmake configrations files which are labelled 0001, 0002, etc.
    configuration_cmake=TRIM(EXECPATH)//'_configuration.cmake'
  ELSE
    configuration_cmake=EXECPATH(1:INDEX(EXECPATH,'/',BACK = .TRUE.))//'configuration.cmake'
  END IF
END IF

! 'configuration_cmake' has been set, inquire the existence of the configuration.cmake file
INQUIRE(File=configuration_cmake,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A)')   ' CALL GetConfigurationFile() : configuration.cmake not found.'
  SWRITE(UNIT_stdOut,'(A,A)') ' configuration_cmake         : ', TRIM(configuration_cmake)
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile                   : ', ExistFile
  ERROR STOP 77
ELSE
  SWRITE(UNIT_stdOut,'(A)')'configuration_cmake=['//TRIM(configuration_cmake)//']'
END IF
SWRITE(UNIT_stdOut,'(A)') ' '
!SWRITE(UNIT_stdOut,'(A)') '...    Continue running?'
!read*
END SUBROUTINE GetConfigurationFile


!==================================================================================================================================
!> Iterate the Configuration permutation counter 
!==================================================================================================================================
SUBROUTINE ConfigurationCounter(N_compile_flags)
! MODULES
USE MOD_RegressionCheck_Vars,  ONLY: BuildCounter,BuildIndex
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                 :: N_compile_flags
! LOCAL VARIABLES
LOGICAL                            :: ExistSolver      ! logical to flag solver
INTEGER                            :: J
!===================================================================================================================================
DO J=1,N_compile_flags
  BuildCounter(J)=BuildCounter(J)+1
  IF(BuildCounter(J).GT.BuildIndex(J))THEN
    BuildCounter(J)=1
  ELSE
    EXIT
  END IF
END DO
END SUBROUTINE ConfigurationCounter


!==================================================================================================================================
!> Get the number of threads/procs for a parallel compilation
!> Check each input argument for being an integer and use it for the number of mpi ranks when compiling the code
!==================================================================================================================================
SUBROUTINE GetNumberOfProcs(nArgs)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars, ONLY: RuntimeOption
USE MOD_RegressionCheck_Vars, ONLY: NumberOfProcs,NumberOfProcsStr
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nArgs
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSTATUS                           !> Error status
INTEGER                        :: I                                 !> loop variable
!===================================================================================================================================
NumberOfProcs=0
IF(nArgs.GE.1)THEN ! first input argument must be "build"
  DO I=1,nArgs
    ! use -j notation like in make (i.e. -j x uses x procs for compiling)
    IF(TRIM(RuntimeOption(I)).EQ.'-j')THEN ! use all available threads for compiling
      NumberOfProcs=-1
      IF(I.EQ.2)THEN
        RuntimeOption(I)='run_basic' ! set back to default
      !ELSE
        !RuntimeOption(I)=''        ! return to default
      END IF
      IF(nArgs.GT.I)THEN ! check if next argument is a number as in "-j 5" -> use 5 threads for compiling
        CALL str2int(RuntimeOption(I+1),NumberOfProcs,iSTATUS) ! if the argument is not a number then "-1" is returned
        !IF(iSTATUS.EQ.0)THEN
          !RuntimeOption(I+1)=''           ! return to default
        !END IF
      END IF
    END IF
  END DO
  IF(NumberOfProcs.EQ.0)THEN
    DO I=1,nArgs
      CALL str2int(RuntimeOption(I),NumberOfProcs,iSTATUS)
      IF(iSTATUS.EQ.0)THEN
        IF(I.EQ.1)THEN
          RuntimeOption(1)='run'        ! return to default -> needed for setting it to 'run_basic'
          SWRITE(UNIT_stdOut,'(A)') ' First argument cannot be a number! Specify "run" or "build"'
          STOP 1
        ELSEIF(I.EQ.2)THEN
          RuntimeOption(2)='run_basic'  ! set to standard case folder 'run_basic'
          IF(TRIM(RuntimeOption(1)).EQ.'run')THEN
            SWRITE(UNIT_stdOut,'(A)') ' An argument with a specific number for building can only be used in "build" mode'
            STOP 1
          END IF
        ELSE
          RuntimeOption(I)=''           ! return to default
        END IF
        EXIT ! first integer input argument found -> use as number of procs for compilation
      ELSE
        NumberOfProcs=1
      END IF
    END DO
  ELSE IF(NumberOfProcs.EQ.-1) THEN
    SWRITE(UNIT_stdOut,'(A)')        ' Building regression checks with [ -j] threads/processors'
  END IF
    ! sanity check
    IF((iSTATUS.EQ.0).AND.(NumberOfProcs.GT.0))THEN
      SWRITE(UNIT_stdOut,'(A,I3,A)') ' Building regression checks with [',NumberOfProcs,'] threads/processors'
    END IF
ELSE
  NumberOfProcs=1
END IF

! set the number of procs INTEGER/CHARACTER
IF(NumberOfProcs.GT.0 .OR. NumberOfProcs.EQ.-1)THEN
  WRITE(UNIT=NumberOfProcsStr,FMT='(I5)') NumberOfProcs
ELSE
  NumberOfProcsStr='fail'
  SWRITE(UNIT_stdOut,'(A)') ' The number of MPI ranks for building must be [>= 1] or [-j]'
  STOP 1
END IF
END SUBROUTINE GetNumberOfProcs


!==================================================================================================================================
!> Print a table containing the information of the error code pointer list
!==================================================================================================================================
SUBROUTINE SummaryOfErrors(EndTime)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars, ONLY: firstError,aError,nErrors,BuildDir,BuildSolver
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT),OPTIONAL    :: EndTime            ! Used to track computation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time
CHARACTER(LEN=255)             :: SpacingBetweenExamples(2) ! set row spacing between examples in table
!===================================================================================================================================
IF(.NOT.PRESENT(EndTime))THEN
  Time=REGGIETIME() ! Measure processing duration -> when e.g. failed to compile occurs
ELSE
  Time=EndTime ! when reggie terminates normally
END IF
SWRITE(UNIT_stdOut,'(132("="))')
nErrors=0
IF(.NOT.ASSOCIATED(aError))THEN
  SWRITE(UNIT_stdOut,'(A)') ' No Examples were executed'
ELSE
  NULLIFY(aError%nextError) ! nullify unused next error pointer
  SWRITE(UNIT_stdOut,'(A)') ' Summary of Errors (0=no Error): '
  SWRITE(UNIT_stdOut,'(A)') ' '
  aError=>firstError ! set aError to first error in list
  SpacingBetweenExamples(:)=''
  SWRITE(UNIT_stdOut,'(A5,2x,A45,2x,A30,2x,A10,2x,A15,2x,A12,2x,A35,2x)')&
                       'run#','Example','SubExample','ErrorCode','build','MPI threads','Information'
  DO WHILE (ASSOCIATED(aError))
    IF( (TRIM(SpacingBetweenExamples(1)).NE.TRIM(aError%Example)) .OR. & ! IF the exmaple changes, include blank line
        (TRIM(SpacingBetweenExamples(2)).NE.TRIM(aError%Build  )) )THEN  ! IF the binary  changes, include blank line
      SWRITE(UNIT_stdOut,'(A)') '' ! include empty line
    END IF
    SpacingBetweenExamples(1)=TRIM(aError%Example)
    SpacingBetweenExamples(2)=TRIM(aError%Build)
    SWRITE(UNIT_stdOut,'(I5,2x)',ADVANCE='no') aError%RunNumber                                               ! run#
    SWRITE(UNIT_stdOut,'(A45,2x)',ADVANCE='no') TRIM(aError%Example)                                          ! Example
    IF(TRIM(aError%SubExampleOption).EQ.'-')THEN
      SWRITE(UNIT_stdOut,'(A30,2x)',ADVANCE='no') '-'
    ELSE
      SWRITE(UNIT_stdOut,'(A30,2x)',ADVANCE='no') TRIM(aError%SubExample)//'='//TRIM(aError%SubExampleOption) ! SubExample
    END IF
    SWRITE(UNIT_stdOut,'(I10,2x)',ADVANCE='no') aError%ErrorCode                                              ! ErrorCode
    SWRITE(UNIT_stdOut,'(A15,2x)',ADVANCE='no') TRIM(aError%Build)                                            ! build
    SWRITE(UNIT_stdOut,'(A12,2x)',ADVANCE='no') TRIM(aError%MPIthreadsStr)                                    ! MPI threads
    SWRITE(UNIT_stdOut,'(A35,2x)',ADVANCE='no') TRIM(aError%Info)                                             ! Information
    SWRITE(UNIT_stdOut,'(A)') ' '
    IF(aError%ErrorCode.NE.0) nErrors=nErrors+1
    aError=>aError%nextError
  END DO
  IF(BuildSolver.AND.(nErrors.GT.0))THEN
    SWRITE(UNIT_stdOut,'(A)') ' '
    SWRITE(UNIT_stdOut,'(A)') ' Run [./regressioncheck build] with the flag [build-continue] in order to skip the successful builds'
    SWRITE(UNIT_stdOut,'(A)') ' and continue with the failed one. The building process can be started at an arbitrary position by'
    SWRITE(UNIT_stdOut,'(A)') ' modifying the file ['//TRIM(BuildDir)//'build_reggie/BuildContinue.reggie] '
    SWRITE(UNIT_stdOut,'(A)') ' where the number corresponds to the build executable'
    SWRITE(UNIT_stdOut,'(A)') ' '
  END IF
  SWRITE(UNIT_stdOut,'(A,I4)') ' Number of errors:  ', nErrors
END IF

SWRITE(UNIT_stdOut,'(132("-"))')
IF(nErrors.GT.0)THEN
  SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' RegressionCheck FAILED! [',Time-StartTime,' sec ]'
ELSE
  SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' RegressionCheck SUCCESSFUL! [',Time-StartTime,' sec ]'
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(132("="))')
END SUBROUTINE SummaryOfErrors


!==================================================================================================================================
!> Add an Error entry to the list of pointers containing information regarding the compilation process, execution process,
!> Error codes and example info
!==================================================================================================================================
SUBROUTINE AddError(MPIthreadsStr,Info,iExample,iSubExample,ErrorStatus,ErrorCode)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: ExampleNames,Examples,EXECPATH,firstError,aError,GlobalRunNumber
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(len=*),INTENT(IN) :: Info,MPIthreadsStr
INTEGER,INTENT(IN)          :: iExample,iSubExample,ErrorStatus,ErrorCode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!    !Debugging:
!    SWRITE(UNIT_stdOut,*)"ErrorCode                                        = [",ErrorCode,"]"
!    SWRITE(UNIT_stdOut,*)"ExampleNames(iExample)                           = [",Trim(ExampleNames(iExample)),"]"
!    SWRITE(UNIT_stdOut,*)"Examples(iExample)%SubExample                    = [",Trim(Examples(iExample)%SubExample),"]"
!    SWRITE(UNIT_stdOut,*)"Examples(iExample)%SubExampleOption(iSubExample) = [",&
!    Trim(Examples(iExample)%SubExampleOption(iSubExample)),"]"
!    SWRITE(UNIT_stdOut,*)"Info                                             = [",Trim(Info),"]"
!    SWRITE(UNIT_stdOut,*)"MPIthreadsStr                                    = [",Trim(MPIthreadsStr),"]"
!    SWRITE(UNIT_stdOut,*)"Build                                            = [",&
!    TRIM(EXECPATH(INDEX(EXECPATH,'/',BACK=.TRUE.)+1:LEN(EXECPATH))),"]"
Examples(iExample)%ErrorStatus=ErrorStatus
IF(firstError%ErrorCode.EQ.-1)THEN ! first error pointer
  firstError%ErrorCode              =ErrorCode ! no error
  firstError%RunNumber              =GlobalRunNumber !
  firstError%Example                =TRIM(ExampleNames(iExample))
  firstError%SubExample             =TRIM(Examples(iExample)%SubExample)
  firstError%SubExampleOption       =TRIM(Examples(iExample)%SubExampleOption(iSubExample))
  firstError%Info                   =TRIM(Info)
  firstError%MPIthreadsStr          =TRIM(MPIthreadsStr)
  firstError%Build                  =TRIM(EXECPATH(INDEX(EXECPATH,'/',BACK=.TRUE.)+1:LEN(EXECPATH)))
  !ALLOCATE(aError) WRONG
  aError=>firstError
ELSE ! next error pointer
  ALLOCATE(aError%nextError)
  aError%nextError%ErrorCode        =ErrorCode ! no error
  aError%nextError%RunNumber        =GlobalRunNumber !
  aError%nextError%Example          =TRIM(ExampleNames(iExample))
  aError%nextError%SubExample       =TRIM(Examples(iExample)%SubExample)
  aError%nextError%SubExampleOption =TRIM(Examples(iExample)%SubExampleOption(iSubExample))
  aError%nextError%Info             =TRIM(Info)
  aError%nextError%MPIthreadsStr    =TRIM(MPIthreadsStr)
  aError%nextError%Build            =TRIM(EXECPATH(INDEX(EXECPATH,'/',BACK=.TRUE.)+1:LEN(EXECPATH)))
  aError=>aError%nextError
END IF

END SUBROUTINE AddError



!==================================================================================================================================
!> search a file for a specific string, return .TRUE. if it is found
!==================================================================================================================================
SUBROUTINE CheckFileForString(FileName,ParameterName,ExistString)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileName       !> e.g. './../build_reggie/bin/configuration.cmake'
CHARACTER(LEN=*),INTENT(IN)  :: ParameterName  !> e.g. 'XX_EQNSYSNAME'
LOGICAL,INTENT(OUT) :: ExistString ! e.g. 'navierstokes'
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: ExistFile    !> file exists=.true., file does not exist=.false.
INTEGER                        :: iSTATUS      !> status
CHARACTER(LEN=255)             :: temp,temp2   !> temp variables for read in of file lines
INTEGER                        :: ioUnit       !> field handler unit and ??
INTEGER                        :: IndNum       !> Index Number
!===================================================================================================================================
!print*,"FileName      = ",TRIM(FileName)
!print*,"ParameterName = ",TRIM(ParameterName)
!print*,"Continue? ... "
!read*
ExistString=.FALSE.
INQUIRE(File=TRIM(FileName),EXIST=ExistFile)
IF(ExistFile) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ')
  DO
    READ(ioUnit,'(A)',iostat=iSTATUS)temp
    temp2=ADJUSTL(temp)
    IF(ADJUSTL(temp2(1:1)).EQ.'!') CYCLE  ! complete line is commented out
    IF(iSTATUS.EQ.-1)EXIT           ! end of file is reached
    IF(LEN(trim(temp)).GT.1)THEN    ! exclude empty lines
      IndNum=INDEX(temp,TRIM(ParameterName)) ! e.g. 'XX_EQNSYSNAME'
      IF(IndNum.GT.0)THEN
        !temp2=TRIM(ADJUSTL(temp(IndNum+LEN(TRIM(ParameterName)):LEN(temp))))
        !IndNum=INDEX(temp2, '=')
        !IF(IndNum.GT.0)THEN
          !temp2=temp2(IndNum+1:LEN(TRIM(temp2)))
          !IndNum=INDEX(temp2, '!')
          !IF(IndNum.GT.0)THEN
            !temp2=temp2(1:IndNum-1)
          !END IF
        !END IF
        ExistString=.TRUE.!TRIM(ADJUSTL(temp2))
        EXIT
      END IF
    END IF
  END DO
  CLOSE(ioUnit)
  !IF(ExistString.EQ.'')ExistString=.FALSE.
ELSE
  ExistString=.FALSE.
END IF
END SUBROUTINE CheckFileForString


!==================================================================================================================================
!> Calculates current time (own function because of a laterMPI implementation)
!==================================================================================================================================
#if USE_MPI
FUNCTION REGGIETIME(Comm)
USE MOD_Globals, ONLY:iError,MPI_COMM_WORLD
USE mpi
#else
FUNCTION REGGIETIME()
#endif
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if USE_MPI
INTEGER, INTENT(IN),OPTIONAL    :: Comm
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: REGGIETIME
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
IF(PRESENT(Comm))THEN
  CALL MPI_BARRIER(Comm,iError)
ELSE
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END IF
#endif
GETTIME(REGGIETIME)
END FUNCTION REGGIETIME


!==================================================================================================================================
!> Calculate the average convergence order for vectors h and E
!==================================================================================================================================
SUBROUTINE CalcOrder(DimOfVectors,h,E,order)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)       :: DimOfVectors
REAL,INTENT(IN)          :: h(DimOfVectors)
REAL,INTENT(IN)          :: E(DimOfVectors)
REAL,INTENT(OUT)         :: order
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: I
REAL                     :: m
!===================================================================================================================================
IF(DimOfVectors.LT.2)THEN
  SWRITE(UNIT_stdOut,'(A)') 'SUBROUTINE CalcOrder(DimOfVectors,h,E,order): input dimension must be larger than 1!'
  order=-999.
END IF

!print*,h
!print*,E

order=0.
DO I=1,DimOfVectors-1
  m=LOG(E(i+1)/E(i))/LOG(h(i+1)/h(i))
!print*,"m=",m
!read*
  order=order+m
END DO
order=order/(DimOfVectors-1)
!print*,"order=",order
!read*
END SUBROUTINE CalcOrder


!==================================================================================================================================
!> print regression check information
!> 1.) give regressioncheck option parameters:
!>     ./regressioncheck [RuntimeOption] [RuntimeOptionType]
!>
!>     ./regressioncheck                -> uses default "run" and runs the current compiler build and all "run_" examples
!>     ./regressioncheck
!>     ./regressioncheck build          -> runs "run_basic" for numerous builds
!>     ./regressioncheck build convtest -> runs "feature_convtest" for numerous builds defined in
!>                                              "feature_convtest/configuration.reggie"
!>     ./regressioncheck build all     -> runs all examples for numerous builds defined in "run_basic/configuration.reggie"
!> 2.) information on input files, e.g., comfiguration.reggie and parameter_reggie.ini
!> 3.) give information on error codes for builing/compiling the source code and running the code
!==================================================================================================================================
SUBROUTINE Print_Help_Information()
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' 1.) How do I run the regression check?'
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' Regression Check Execution: ./regressioncheck [RuntimeOption] [RuntimeOptionType]                     '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' First input argument [RuntimeOption]                                                                  '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' [RuntimeOption]            | mode                                                                     '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' help                       | prints this information output                                           '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' run (default)              | runs all examples with prefix "run", e.g., "run_1" or                    '
SWRITE(UNIT_stdOut,'(A)') '                            | "run_basic". These tests all inlcude very                                '
SWRITE(UNIT_stdOut,'(A)') '                            | short simulations with <1sec execution time                              '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for on-check-in tests                                          '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' build                      | runs the example "run_basic" on default and                              '
SWRITE(UNIT_stdOut,'(A)') '                            | (requires locally built HDF5 or loaded HDF5 paths)                       '
SWRITE(UNIT_stdOut,'(A)') '                            | compiles all possible compiler flag combinations                         '
SWRITE(UNIT_stdOut,'(A)') '                            | specified in "configuration.reggie" and considers                        '
SWRITE(UNIT_stdOut,'(A)') '                            | the specified exclude list for invalid combinations                      '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for nightly tests                                              '
SWRITE(UNIT_stdOut,'(A)') '                            | The cmake output may be shown (disabled) on-screen by adding "debug"     '
SWRITE(UNIT_stdOut,'(A)') '                            | ("no-debug") after possible [RuntimeOptionType] commands                 '
SWRITE(UNIT_stdOut,'(A)') '                            | Multi-processor compilation is supported by adding "XX"  after possible  '
SWRITE(UNIT_stdOut,'(A)') '                            | [RuntimeOptionType] commands, where "XX" is the number of processors     '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' build-continue             | run ./regressioncheck build build-continue in order to continue the      '
SWRITE(UNIT_stdOut,'(A)') '                            | building process at the last possibly failed compilation step            '
SWRITE(UNIT_stdOut,'(A)') '                            | does the same as "build", but starts the process at the last failed      '
SWRITE(UNIT_stdOut,'(A)') '                            | building position of a previous run, hence, all previously successful    '
SWRITE(UNIT_stdOut,'(A)') '                            | builds are skipped (-> debugging purposes)                               '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' run-continue               | run ./regressioncheck build run-continue in order to skip already built  '
SWRITE(UNIT_stdOut,'(A)') '                            | executables and simply run the regression check tests. Continue building '
SWRITE(UNIT_stdOut,'(A)') '                            | as soon as the corresponding executables are missing in the              '
SWRITE(UNIT_stdOut,'(A)') '                            | /build_reggie_bin directory                                              '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' conv_test                  | specific feature test: runs the "conv_test" example                      '
SWRITE(UNIT_stdOut,'(A)') '                            | runs two modes: p-convergence and h-convergence                          '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for weekly tests                                               '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' performance (ToDo!)        | specific feature test: runs the "performance" example                    '
SWRITE(UNIT_stdOut,'(A)') '                            | automatically checks out specified code version tag                      '
SWRITE(UNIT_stdOut,'(A)') '                            | and run the example to acquire the reference                             '
SWRITE(UNIT_stdOut,'(A)') '                            | performance wall time                                                    '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for weekly tests                                               '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' Second input argument [RuntimeOptionType] depends on the first input argument                         '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' [RuntimeOptionType]        | [RuntimeOption] | mode                                                   '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' run (default)              |                 |                                                        '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' cavity (example)           | run             | specific feature test: runs all the the examples that  '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | begin with "cavity", e.g., "cavity_1" or "cavity_xyz"  '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | the example "cavity", e.g., tests long time stability, '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | some BC and ExactFunc ( e.g. used for weekly tests)    '
SWRITE(UNIT_stdOut,'(A)') '                            |                 |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            | build           | use the "configurations.reggie" within the "cavity"    '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | directory and builds all compiler flag combinations    '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | that are specified there and runs them all on "cavity" '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'


! needed/optional files used by the regression check
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' 2.) Which input files are needed, e.g., comfiguration.reggie and parameter_reggie.ini'
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' The following parameter files are supported (within each example folder):           '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' parameter file             | Description                                            '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' configuration.reggie       | needed for "./regressioncheck build"                   '
SWRITE(UNIT_stdOut,'(A)') '                            | contains all required compilation flags needed to      '
SWRITE(UNIT_stdOut,'(A)') '                            | create a minimum of one combinations. If multiple      '
SWRITE(UNIT_stdOut,'(A)') '                            | compiler flag cmake combinations are specified, all    '
SWRITE(UNIT_stdOut,'(A)') '                            | possible combinations are created if they to not       '
SWRITE(UNIT_stdOut,'(A)') '                            | violate the listed exlclude combinations. An example is'
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' parameter_reggie.ini       | contains parameters for the specific regression check  '
SWRITE(UNIT_stdOut,'(A)') '                            |                                                        '
SWRITE(UNIT_stdOut,'(A)') '        number of variables | nVar= 5 (depricated)                                   '
SWRITE(UNIT_stdOut,'(A)') '                 MPI on/off | MPI= T                                                 '
SWRITE(UNIT_stdOut,'(A)') '     L2/Linf reference file | ReferenceNormFile= referencenorm.txt                   '
SWRITE(UNIT_stdOut,'(A)') '  ref state file for h5diff | H5DIFFReferenceStateFile= cavity_refe_State_0.20.h5    '
SWRITE(UNIT_stdOut,'(A)') '      state file for h5diff | H5DIFFCheckedStateFile= cavity_State_0000000.200.h5    '
SWRITE(UNIT_stdOut,'(A)') '      array name for h5diff | H5DIFFReferenceDataSetName= DG_Solution                '
SWRITE(UNIT_stdOut,'(A)') '       if restart is wanted | RestartFileName=                                       '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'



! error code description
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' 3.) What do the Error Codes tell me?                                                '
SWRITE(UNIT_stdOut,'(A)') '                                                                                     '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' Error Code                 | Description                                            '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') '                          0 | no error                                               '
SWRITE(UNIT_stdOut,'(A)') '                          1 | failed during build                                    '
SWRITE(UNIT_stdOut,'(A)') '                          2 | computation of example failed                          '
SWRITE(UNIT_stdOut,'(A)') '                          3 | mismatch in norms                                      '
SWRITE(UNIT_stdOut,'(A)') '                          4 | mismatch in dataset                                    '
SWRITE(UNIT_stdOut,'(A)') '                          5 | fail during comparison                                 '
SWRITE(UNIT_stdOut,'(A)') '                          6 | mismatch in IntegrateLine                              '
SWRITE(UNIT_stdOut,'(A)') '                         77 | no executable found for option run                     '
SWRITE(UNIT_stdOut,'(A)') '                         99 | fail of execute_system_command                         '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
END SUBROUTINE Print_Help_Information


END MODULE MOD_RegressionCheck_Tools
