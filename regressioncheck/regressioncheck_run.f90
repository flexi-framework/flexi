#include "flexi.h"

!==================================================================================================================================
!> Contains the routines to 
!> - perform the actual regressioncheck
!==================================================================================================================================
MODULE MOD_RegressionCheck_Run
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE PerformRegressionCheck
  MODULE PROCEDURE PerformRegressionCheck
END INTERFACE

INTERFACE PerformFullRegressionCheck
  MODULE PROCEDURE PerformFullRegressionCheck
END INTERFACE

PUBLIC::PerformRegressionCheck,PerformFullRegressionCheck
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Routine which performs the actual regressioncheck. It triggers the builds and execute commands. Additionally, it performs
!> the checks for L2-error norms, h5-diff and runtime
!==================================================================================================================================
SUBROUTINE PerformRegressionCheck()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Compare, ONLY: CompareResults,CompareConvergence
USE MOD_RegressionCheck_Tools,   ONLY: InitExample
USE MOD_RegressionCheck_Vars,    ONLY: nExamples,ExampleNames,Examples,EXECPATH,RuntimeOption,GlobalRunNumber 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: TESTCASE                          !> compilation flag in FLEXI
CHARACTER(LEN=255)             :: TIMEDISCMETHOD                    !> compilation flag in PICLas
CHARACTER(LEN=255)             :: ReggieBuildExe                    !> cache executables when doing "BuildSolver":
                                                                    !> this means don't build the same cmake configuration twice
CHARACTER(LEN=255)             :: parameter_ini                     !> input parameter file depending on EQNSYS
CHARACTER(LEN=255)             :: parameter_ini2                    !> input parameter file depending on TIMEDISC
CHARACTER(LEN=15)              :: MPIthreadsStr                     !> string for the number of MPI threads for execution
                                                                    !> it may differ from "Examples(iExample)%MPIthreadsStr"
INTEGER                        :: iExample                          !> loop index for example
INTEGER                        :: N_compile_flags                   !> number of compile-flags
INTEGER                        :: iReggieBuild,nReggieBuilds ! field handler unit and ??
INTEGER                        :: iSubExample,iScaling,iRun
LOGICAL                        :: SkipExample,SkipBuild,ExitBuild,SkipFolder,SkipComparison
LOGICAL                        :: UseFV,Use2D,UsePARABOLIC      !> compiler flags currently used for ConvergenceTest
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Performing tests ...'
ReggieBuildExe='' ! init
GlobalRunNumber=0      ! init
!==================================================================================================================================
DO iExample = 1, nExamples ! loop level 1 of 5
!==================================================================================================================================
  CALL CheckExampleName(ExampleNames(iExample),RuntimeOption(2),SkipExample)
  IF(SkipExample)CYCLE ! ignore the example folder and continue with the next

  ! set the build configuration environment when BuildSolver=.TRUE.
  CALL GetnReggieBuilds(iExample,ReggieBuildExe,N_compile_flags,nReggieBuilds)

  ! delete pre-existing data files at the beginning of the reggie
  CALL CleanFolder(iExample,MODE=0) ! MODE=0: INITIAL -> delete pre-existing files and folders
  
!==================================================================================================================================
  DO iReggieBuild = 1, nReggieBuilds ! loop level 2 of 5: cycle the number of build configurations (no configuration = only 1 run)
!==================================================================================================================================
    ! read the parameters for the current example (parameter_reggie.ini), must be read separately for every iReggieBuild
    ! because changes to specific parameters are made depending on the cmake compilation flags, e.g., in "CALL SetParameters(...)"
    CALL InitExample(Examples(iExample)%PATH,Examples(iExample),SkipExample)
    IF(SkipExample)CYCLE ! skip if "parameter_reggie.ini" file is missing

    ! Get code binary (build or find it)
    CALL GetCodeBinary(iExample,iReggieBuild,nReggieBuilds,ReggieBuildExe,SkipBuild,ExitBuild)
    IF(SkipBuild)CYCLE ! invalid reggie build but not last reggie build
    IF(ExitBuild)EXIT  ! last reggie build -> exit ("cycle" would start an infinite loop)

    ! depending on the equation system -> get different Nvar 
    CALL GetNvar(iExample,iReggieBuild)

    ! check if executable is compiled with correct TESTCASE (e.g. for tylorgreenvortex)
    CALL CheckCompilerFlags(iExample,iReggieBuild,TESTCASE,TIMEDISCMETHOD,UseFV,Use2D,UsePARABOLIC)

    ! remove subexample (before printing the case overview) for certain configurations: e.g. Preconditioner when running explicitly
    CALL CheckSubExample(iExample,iReggieBuild,TIMEDISCMETHOD)

    ! check folder name and decide whether it can be executed with the current binary (e.g. testcases ...)
    CALL CheckFolderName(iExample,TESTCASE,SkipFolder)
    IF(SkipFolder)CYCLE ! e.g. TESTCASE folder and non-TESTCASE binary or vice versa

    ! get list of parameter files for running the simulation
    CALL GetParameterFiles(iExample,TIMEDISCMETHOD,parameter_ini,parameter_ini2)

    ! Output settings (before going into subexamples): display table with settings
    CALL PrintExampleInfo(iExample,EXECPATH,parameter_ini,parameter_ini2)
 
    ! Set options in parameter.ini file
    CALL SetParameters(iExample,parameter_ini,UseFV,Use2D,UsePARABOLIC,SkipFolder)
    IF(SkipFolder)CYCLE ! e.g. p-convergence folder and FV subcells (p-convergence not meaningful)

!==================================================================================================================================
    DO iSubExample = 1, MAX(1,Examples(iExample)%SubExampleNumber) ! loop level 3 of 5: SubExamples (e.g. different TimeDiscMethods)
!==================================================================================================================================
      ! Set the SubExample in the parameter.ini file
      CALL SetSubExample(iExample,iSubExample,parameter_ini)

      ! delete pre-existing data files before running the code (e.g. "TGVAnalysis.dat" or "Database.csv")
      CALL CleanFolder(iExample,MODE=1) ! MODE=1: delete pre-existing files and folders

!==================================================================================================================================
      DO iScaling = 1, Examples(iExample)%MPIthreadsN ! loop level 4 of 5: multiple MPI runs with different MPI threads
!==================================================================================================================================
!==================================================================================================================================
        DO iRun = 1, Examples(iExample)%nRuns ! loop level 5 of 5: repeat the same run multiple times
!==================================================================================================================================
          CALL RunTheCode(iExample,iSubExample,iScaling,iRun,MPIthreadsStr,EXECPATH,parameter_ini,parameter_ini2,SkipComparison)
          IF(SkipComparison)CYCLE ! the execution has failed, no comparisons are needed

          ! compare the results and write error messages for the current case
          CALL CompareResults(iExample,iSubExample,MPIthreadsStr)

          ! IF all comparisons are successful the error status is 0 -> delete created files in CleanFolder(iExample)
          IF(Examples(iExample)%ErrorStatus.EQ.0) CALL CleanFolder(iExample,MODE=2) ! MODE=2: delete files after simulation
        END DO ! iScalingRuns = 1, Examples(iExample)%nRuns
      END DO ! iScaling = 1, Examples(iExample)%MPIthreadsN
    END DO ! iSubExample = 1, MAX(1,SubExampleNumber) (for cases without specified SubExamples: SubExampleNumber=0)
  END DO ! iReggieBuild = 1, nReggieBuilds
END DO ! iExample=1,nExamples

END SUBROUTINE PerformRegressionCheck


!==================================================================================================================================
!> Routine which performs the actual regressioncheck. It triggers the builds and execute commands. Additionally, it performs
!> the checks for L2-error norms, h5-diff and runtime
!==================================================================================================================================
SUBROUTINE PerformFullRegressionCheck()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Compare, ONLY: CompareResults,CompareConvergence
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=500)            :: SYSCOMMAND                         !> string to fit the system command
CHARACTER(LEN=500)            :: FileName,OutputFile                !> file and path strings
INTEGER                       :: ioUnit                             !> io-unit
INTEGER                       :: iSTATUS                            !> status
LOGICAL                       :: ExistFile                          !> T=file exists, F=file does not exist
CHARACTER(len=500)            :: temp,temp2,TmpStr                  !> auxiliary variables
CHARACTER(len=500)            :: reggie(50)                         !> consider 50 cases maximum
CHARACTER(len=500)            :: reggieUnique(50)                   !> consider 50 cases maximum
LOGICAL                       :: isNotUnique                        !> don't repeat the same case twice
INTEGER                       :: IndNum                             !> index
INTEGER                       :: iReggie,nReggie,iReggieDONE        !> number of cases
INTEGER                       :: iReggieUnique,nReggieUnique        !> number of cases
CHARACTER(len=50)             :: PreFix                             !> auxiliary variable
CHARACTER(LEN=50)             :: FileSuffix                         !> auxiliary vars for file name endings
INTEGER                       :: IndNum1,IndNum2,Ind1,Ind2          !> index numbers
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Performing FULL regressioncheck ...'
nReggie=0
iSTATUS=0 ! nullify

! check for the '.gitlab-ci.yml' file
FileName=TRIM(BASEDIR(2:LEN(BASEDIR)-1))//'../.gitlab-ci.yml'
INQUIRE(File=FileName,EXIST=ExistFile)              ! inquire existence
IF(ExistFile.EQV..FALSE.)THEN ! check a second possibility
  FileName=TRIM(BASEDIR(2:LEN(BASEDIR)-1))//'.gitlab-ci.yml'
  INQUIRE(File=FileName,EXIST=ExistFile)            ! inquire existence
END IF

! Read the file if is exists
IF(ExistFile) THEN
  print*,"FileName=",TRIM(FileName)
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ')
  IF(iSTATUS.NE.0)THEN
    SWRITE(UNIT_stdOut,'(A)') ' Could not open .gitlab-ci.yml'
    ERROR STOP 1
  END IF
  DO
    READ(ioUnit,'(A)',iostat=iSTATUS)temp
    IF(iSTATUS.GT.0)THEN
      SWRITE(UNIT_stdOut,'(A,I5)') " Read-in failed: iSTATUS=",iSTATUS
      SWRITE(UNIT_stdOut,'(A)') temp
      ERROR STOP 1
    ELSE IF(iSTATUS.LT.0)THEN
      EXIT
    ELSE ! iSTATUS = 0
      IndNum=INDEX(temp,'./regressioncheck')
      IF(IndNum.GT.0)THEN
        nReggie=nReggie+1
        temp2=TRIM(ADJUSTL(temp(IndNum:LEN(temp)))) ! +17
        IndNum=INDEX(temp2,';')
        IF(IndNum.GT.0)THEN
          temp=TRIM(ADJUSTL(temp2(1:IndNum-1)))
          temp2=temp
        END IF
        !print*,'[',trim(temp2),']   length=',LEN(TRIM(ADJUSTL(temp2)))
        reggie(nReggie)=TRIM(ADJUSTL(temp2))//' no-full'//' no-debug'
      END IF
    END IF
  END DO
  CLOSE(ioUnit)
ELSE ! could not find '.gitlab-ci.yml'
  SWRITE(UNIT_stdOut,'(A13,A)') ' ERROR     : ','no ".gitlab-ci.yml" found.'
  SWRITE(UNIT_stdOut,'(A13,A)') ' FileName  : ', TRIM(FileName)
  SWRITE(UNIT_stdOut,'(A13,L)') ' ExistFile : ', ExistFile
  ERROR STOP 1
END IF

FileName=TRIM(BASEDIR(2:LEN(BASEDIR)-1))//'bin/regressioncheck'
INQUIRE(File=FileName,EXIST=ExistFile)              ! inquire existence
IF(ExistFile.EQV..FALSE.)THEN ! check a second possibility
  FileName=TRIM(BASEDIR(2:LEN(BASEDIR)-1))//'regressioncheck'
  INQUIRE(File=FileName,EXIST=ExistFile)            ! inquire existence
  PreFix=''
ELSE
  PreFix='bin/'
END IF

IF(ExistFile.EQV..FALSE.)THEN
  SWRITE(UNIT_stdOut,'(A13,A)') ' ERROR     : ','no "regressioncheck" found.'
  SWRITE(UNIT_stdOut,'(A13,A)') ' FileName  : ', TRIM(FileName)
  SWRITE(UNIT_stdOut,'(A13,L)') ' ExistFile : ', ExistFile
  ERROR STOP 1
END IF
SWRITE(UNIT_stdOut,'(A)')' '

! Display an overview of all cases that were found
SWRITE(UNIT_stdOut,'(A)') " FULL list: The following cases were found in .gitlab-ci.yml"
SWRITE(UNIT_stdOut,'(132("-"))')
DO iReggie=1,nReggie
  SWRITE(UNIT_stdOut,'(A)') ' cd '//TRIM(BASEDIR)//' && '//TRIM(PreFix)//TRIM(reggie(iReggie))
END DO
SWRITE(UNIT_stdOut,'(132("-"))')

! sort out duplicate cases
iReggieUnique=0
reggieUnique='-1'
DO iReggie=1,nReggie
  isNotUnique=.FALSE.
  DO iReggieDONE=1,iReggie
    IF(TRIM(reggieUnique(iReggieDONE)).EQ.TRIM(reggie(iReggie)))THEN ! if any old case (reggieUnique) matches 
                                                                     ! the current case (reggie)
      isNotUnique=.TRUE.
    END IF
  END DO
  IF(isNotUnique)CYCLE ! don't repeat the same case twice
  WRITE(FileSuffix,'(I4.4)')iReggieUnique
  reggieUnique(iReggie)=TRIM(reggie(iReggie))
  iReggieUnique=iReggieUnique+1
END DO
nReggieUnique=iReggieUnique
iReggieUnique=1
! map unique cases to reggie()
DO iReggie=1,nReggie
  IF(reggieUnique(iReggie).EQ.'-1')CYCLE
  reggie(iReggieUnique)=reggieUnique(iReggie)
  iReggieUnique=iReggieUnique+1
END DO

! Display an overview of all cases that will be run
SWRITE(UNIT_stdOut,'(A)') " UNIQUE list: The following cases will be executed"
SWRITE(UNIT_stdOut,'(132("-"))')
DO iReggie=1,nReggieUnique
  SWRITE(UNIT_stdOut,'(A)') ' cd '//TRIM(BASEDIR)//' && '//TRIM(PreFix)//TRIM(reggie(iReggie))
END DO
SWRITE(UNIT_stdOut,'(132("-"))')

! loop the cases that will be executed
DO iReggie=1,nReggieUnique
  WRITE(FileSuffix,'(I4.4)')iReggie
  SYSCOMMAND='cd '//TRIM(BASEDIR)//' && '//TRIM(PreFix)//TRIM(reggie(iReggie))//' > '& !' | tee '&
             //TRIM(PreFix)//'full-reggie-'//TRIM(FileSuffix)//'.out'
  SWRITE(UNIT_stdOut,'(A)')' '
  ! display the file where the output will be written
  IndNum=INDEX(BASEDIR,'/')
  IndNum1=INDEX(BASEDIR,'"') ! remove " from file path
  Ind1=1
  IF((IndNum1.GT.0).AND.(IndNum1.LT.IndNum))  Ind1=MAX(1,IndNum1+1) ! only valid if ["] comes before [/]
  IndNum2=INDEX(BASEDIR,'"',BACK = .TRUE.) ! remove " from file path
  Ind2=LEN(BASEDIR)
  IF((IndNum2.GT.1).AND.(IndNum1.NE.IndNum2)) Ind2=IndNum2-1        ! only valid if a second ["] exists 
  
  OutputFile=TRIM(BASEDIR(Ind1:Ind2))//TRIM(PreFix)//'full-reggie-'//TRIM(FileSuffix)//'.out'
  SWRITE(UNIT_stdOut,'(A,A,A,A)') ' Running case [',TRIM(FileSuffix),'] and writing the output to ',TRIM(OutputFile)
                                   !TRIM(BASEDIR)//TRIM(PreFix)//'full-reggie-'//TRIM(FileSuffix)//'.out'
  ! display the command that will be executed
  SWRITE(UNIT_stdOut,'(A)')'   SYSCOMMAND: '//TRIM(SYSCOMMAND)
  ! run the recursive regressioncheck
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  IF(iSTATUS.NE.0)THEN
    SWRITE(UNIT_stdOut,'(A)')    '   Failed running recursive regressioncheck'
    SWRITE(UNIT_stdOut,'(A,I5)') '     iSTATUS: ',iSTATUS
    ERROR STOP 1
  ELSE
    CALL GetParameterFromFile(OutputFile,'RegressionCheck SUCCESSFUL!',TmpStr)
    SWRITE(UNIT_stdOut,'(A,A)') '   RegressionCheck SUCCESSFUL! ',ADJUSTL(TRIM(TmpStr))
  END IF
END DO

END SUBROUTINE PerformFullRegressionCheck


!==================================================================================================================================
!> depending on the equation system -> get different Nvar (number of variables in the equation system) for the current example
!> currently supports: - navierstokes             ->    Examples(iExample)%Nvar=5
!>                     - linearscalaradvection    ->    Examples(iExample)%Nvar=1
!>                     - maxwell                  ->    Examples(iExample)%Nvar=8
!> Add equation systems at will!!!!
!==================================================================================================================================
SUBROUTINE GetNvar(iExample,iReggieBuild)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples,EXECPATH,BuildEQNSYS,BuildSolver,CodeNameUppCase,CodeNameLowCase,configuration_cmake
USE MOD_RegressionCheck_Vars,    ONLY: RunContinue
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iReggieBuild
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: IndNum
CHARACTER(LEN=255)             :: temp,temp2                        ! temp variables for read in of file lines
LOGICAL                        :: ExistFile                         ! file exists=.true., file does not exist=.false.
INTEGER                        :: iSTATUS                           ! status
INTEGER                        :: ioUnit                            ! IO channel
!===================================================================================================================================
iSTATUS=0 ! nullify
IndNum=INDEX(EXECPATH, '/')
IF(IndNum.GT.0)THEN
  IndNum=INDEX(EXECPATH,'/',BACK = .TRUE.) ! get path without binary
  !configuration_cmake=EXECPATH(1:IndNum)//'configuration.cmake'
  INQUIRE(File=configuration_cmake,EXIST=ExistFile)
  IF(ExistFile) THEN
    OPEN(NEWUNIT=ioUnit,FILE=TRIM(configuration_cmake),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ')
    DO
      READ(ioUnit,'(A)',iostat=iSTATUS)temp
      IF(ADJUSTL(temp).EQ.'!') CYCLE
      IF(iSTATUS.EQ.-1)EXIT
      IF(LEN(trim(temp)).GT.1)THEN
        IndNum=INDEX(temp,CodeNameUppCase//'_EQNSYSNAME')
        IF(IndNum.GT.0)THEN
          temp2=TRIM(ADJUSTL(temp(IndNum+LEN(CodeNameUppCase//'_EQNSYSNAME'):LEN(temp))))
          IndNum=INDEX(temp2, '"')
          IF(IndNum.GT.0)THEN
            temp2=temp2(IndNum+1:LEN(TRIM(temp2)))
            IF(INDEX(temp2(IndNum+1:LEN(TRIM(temp2))), '"')+IndNum.GT.IndNum)THEN ! get binary path up to 2nd '"' in name
              IndNum=INDEX(temp2(IndNum+1:LEN(TRIM(temp2))), '"')+IndNum
            END IF
          END IF
          Examples(iExample)%EQNSYSNAME=temp2(1:IndNum-1)
          EXIT
        END IF
      END IF
    END DO
    CLOSE(ioUnit)
  ELSE ! could not find 'configuration.cmake' at location of execution binary
    IF((BuildSolver).AND.(.NOT.RunContinue))THEN ! get EQNSYSNAME from cmake build configuration settings
      Examples(iExample)%EQNSYSNAME=BuildEQNSYS(iReggieBuild)
    ELSE ! stop 
      SWRITE(UNIT_stdOut,'(A24,A)') ' ERROR: ','no "configuration.cmake" found at the location of the '&
                                                                                                   //CodeNameLowCase//' binary.'
      SWRITE(UNIT_stdOut,'(A24,A)') ' configuration_cmake: ', TRIM(configuration_cmake)
      SWRITE(UNIT_stdOut,'(A24,L)') ' ExistFile: ', ExistFile
      ERROR STOP 1
    END IF
  END IF
END IF
! select the correct equation system: ADD NEW SYSTEMS WHEN NEEDED!!!
SELECT CASE (TRIM(Examples(iExample)%EQNSYSNAME))
  CASE ('navierstokes')  
    Examples(iExample)%Nvar=5
  CASE ('linearscalaradvection')  
    Examples(iExample)%Nvar=1
  CASE ('maxwell')  
    Examples(iExample)%Nvar=8
  CASE ('poisson')  
    Examples(iExample)%Nvar=1
  CASE DEFAULT
    Examples(iExample)%Nvar=-1
    SWRITE(UNIT_stdOut,'(A)') ' ERROR: missing case select for this '&
                                                                   //CodeNameUppCase//'_EQNSYSNAME with appropriate Nvar. Fix it by'
    SWRITE(UNIT_stdOut,'(A)') '        adding the correct line of code to ../regressioncheck/regressioncheck_run.f90'
    SWRITE(UNIT_stdOut,'(A)') '        Examples(iExample)%EQNSYSNAME=['//TRIM(Examples(iExample)%EQNSYSNAME)//']'
    ERROR STOP 77
END SELECT
END SUBROUTINE GetNvar


!===================================================================================================================================
!> Check the example folder name if it matches the supplied input variable (either the complete string or only the first characters)
!===================================================================================================================================
SUBROUTINE CheckExampleName(ExampleName,RuntimeOptionType,SkipExample)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: ExampleName,RuntimeOptionType
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: SkipExample
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: dummystr     ! temp string variable
!===================================================================================================================================
! check example: currently each parameter build configuration is only built and tested for the "run_basic" example
dummystr=TRIM(ADJUSTL(ExampleName)) ! e.g. "run_basic"
IF(dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType)))).EQ.RuntimeOptionType)THEN ! e.g. "run[_basic]" = "run"
  SWRITE(UNIT_stdOut,'(A)') ''
END IF
SWRITE(UNIT_stdOut,'(A65)',ADVANCE='no') ' Example-Name: '//dummystr
IF(dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType)))).NE.RuntimeOptionType)THEN
  SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...skipping'
  SkipExample=.TRUE.
ELSE
  SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...running'
  SkipExample=.FALSE.
END IF
END SUBROUTINE CheckExampleName


!===================================================================================================================================
!> IF BuildSolver=T: read the configurations.reggie and determine the number of builds
!> ELSE            : nReggieBuilds=1
!===================================================================================================================================
SUBROUTINE GetnReggieBuilds(iExample,ReggieBuildExe,N_compile_flags,nReggieBuilds)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Build,   ONLY: ReadConfiguration
USE MOD_RegressionCheck_Vars,    ONLY: BuildDir,BuildCounter,BuildSolver,BuildContinue,RunContinue
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
CHARACTER(LEN=*),INTENT(IN)    :: ReggieBuildExe
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)          :: nReggieBuilds,N_compile_flags
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=500)             :: SYSCOMMAND                        !> string to fit the system command
INTEGER                        :: iSTATUS                           !> status
!===================================================================================================================================
iSTATUS=0 ! nullify
! if "BuildSolver" is true, the complete (valid) compiler-flag parameter combination
! is tested (specified in "configuration.reggie", default example is "run_basic")
IF(BuildSolver)THEN
  IF(ReggieBuildExe.EQ.'')THEN
    CALL ReadConfiguration(iExample,nReggieBuilds,N_compile_flags)
    BuildCounter=1 ! reset the counter between read and build (used for selecting the build configuration for compilation)
    IF((.NOT.BuildContinue).AND.(.NOT.RunContinue))THEN ! only delete the executables folder if reggie starts at the beginning
      SYSCOMMAND='rm -rf '//TRIM(BuildDir)//'build_reggie_bin > /dev/null 2>&1'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      SYSCOMMAND= 'mkdir '//TRIM(BuildDir)//'build_reggie_bin'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
    END IF
  ELSE
    SWRITE(UNIT_stdOut,'(A)')'regressioncheck_run.f90: CALL ReadConfiguration() has already been performed, skipping...'
  END IF
ELSE ! pure run, no compilation
  nReggieBuilds=1
END IF
END SUBROUTINE GetnReggieBuilds


!===================================================================================================================================
!> set the correct path to the binary for running a simulation
!> 1.) build a new binary with cmake, e.g., XXXXX which will be moved to /build_reggie_bin/XXXXX0001 for re-use
!> 2.) use a previously built binary that was created with reggie, e.g., /build_reggie_bin/XXXXX0001
!> 3.) use an existing binary
!===================================================================================================================================
SUBROUTINE GetCodeBinary(iExample,iReggieBuild,nReggieBuilds,ReggieBuildExe,SkipBuild,ExitBuild)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: BuildDir,CodeNameLowCase,EXECPATH
USE MOD_RegressionCheck_Vars,    ONLY: BuildValid,BuildSolver
USE MOD_RegressionCheck_Tools,   ONLY: CheckForExecutable,GetConfigurationFile!,ConfigurationCounter
USE MOD_RegressionCheck_Build,   ONLY: BuildConfiguration
USE MOD_RegressionCheck_Vars,    ONLY: BuildContinue,BuildContinueNumber,RunContinue
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iReggieBuild,nReggieBuilds
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: ReggieBuildExe
LOGICAL,INTENT(OUT)            :: SkipBuild,ExitBuild
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: ExistFile                         !> file exists=.true., file does not exist=.false.
CHARACTER(LEN=255)             :: FileName                          !> path to a file or its name
INTEGER                        :: iSTATUS                           !> status
CHARACTER(LEN=500)             :: SYSCOMMAND                        !> string to fit the system command
!===================================================================================================================================
SkipBuild=.FALSE.
ExitBuild=.FALSE.
iSTATUS=0 ! nullify
! Get code binary (build or find it)
IF(BuildSolver)THEN
!print*,""
!print*,"iReggieBuild=[",iReggieBuild,"]"
!print*,"BuildValid(iReggieBuild)=[",BuildValid(iReggieBuild),"]"

  ! if build is not valid no binary has been built and the lopp can cycle here
  IF(.NOT.BuildValid(iReggieBuild))THEN ! invalid reggie build 
    WRITE (ReggieBuildExe, '(a, i4.4)') "[BUILD_is_invalid]"
  ELSE
    WRITE (ReggieBuildExe, '(a, i4.4)') CodeNameLowCase, COUNT(BuildValid(1:iReggieBuild)) ! e.g. XXXXX0001
  END IF
  ! check if build exists -> if it does, don't build a new executable with cmake
  FileName=TRIM(BuildDir)//'build_reggie_bin/'//ReggieBuildExe
  INQUIRE(File=FileName,EXIST=ExistFile)

  IF(BuildContinue)THEN ! even if the executable exists, re-build it when it equels BuildContinueNumber
    IF(COUNT(BuildValid(1:iReggieBuild)).EQ.BuildContinueNumber)THEN
      ExistFile=.FALSE. ! set in order to skip building
    END IF
  END IF

  IF(RunContinue)THEN ! info when re-using the built binaries
    SWRITE(UNIT_stdOut,'(A)') ' '
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A)') TRIM(FileName)
    SWRITE(UNIT_stdOut,'(A,L)')'Binary exists: ',ExistFile
    ! NEW -> NOT USED ANY MORE !  CALL ConfigurationCounter(N_compile_flags) ! call here because below BuildConfiguration(...) will be skipped
  END IF

  IF(ExistFile) THEN ! 1. build already exists (e.g. XXXX0001 located in ../build_reggie_bin/)
    EXECPATH=TRIM(FileName)
  ELSE ! 2. build does not exists -> create it
    CALL BuildConfiguration(iExample,iReggieBuild,nReggieBuilds)
    IF(BuildValid(iReggieBuild))THEN ! only move binary if it has been created (only for valid builds)
      SYSCOMMAND='cd '//TRIM(BuildDir)//' && mv build_reggie/bin/'//CodeNameLowCase//' build_reggie_bin/'//TRIM(ReggieBuildExe)
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      SYSCOMMAND='cd '//TRIM(BuildDir)//' && cp build_reggie/bin/configuration.cmake build_reggie_bin/'//TRIM(ReggieBuildExe)//&
                 '_configuration.cmake'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      EXECPATH=TRIM(BuildDir)//'build_reggie_bin/'//ReggieBuildExe
    END IF
  END IF

  ! if build is not valid no exe has been built and the lopp can cycle here
  IF((.NOT.BuildValid(iReggieBuild)).AND.(iReggieBuild.NE.nReggieBuilds))THEN ! invalid reggie build but not last reggie build
    SkipBuild=.TRUE. ! -> does a CYLCE
  ELSEIF((.NOT.BuildValid(iReggieBuild)).AND.(iReggieBuild.EQ.nReggieBuilds))THEN ! last reggie build -> exit 
                                                                                  ! ("cycle" would start an infinite loop)
    ExitBuild=.TRUE. ! -> does an EXIT
  END IF


  IF(BuildContinue)THEN ! skip build if the build is smaller than BuildContinueNumber
    IF(COUNT(BuildValid(1:iReggieBuild)).LT.BuildContinueNumber)THEN ! skip the first N builds when BuildContinue is used
    !IF(iReggieBuild.LT.BuildContinueNumber)THEN ! skip the first N builds when BuildContinue is used
!print*,"COUNT(BuildValid(1:iReggieBuild))",COUNT(BuildValid(1:iReggieBuild))
!print*,"BuildContinueNumber              ",BuildContinueNumber,"----> RETURN"
!read*
      SkipBuild=.TRUE. ! set in order to skip running
      !RETURN
    END IF
  END IF

  CALL CheckForExecutable(Mode=2) ! check if executable was created correctly
END IF

! set the configuration file path info depending on reggie build/run setting
IF(.NOT.SkipBuild)THEN ! only set the path when the build is not to be skipped
  CALL GetConfigurationFile()
END IF

END SUBROUTINE GetCodeBinary



!===================================================================================================================================
!> check if executable is compiled with correct TESTCASE (e.g. for tylorgreenvortex)
!===================================================================================================================================
SUBROUTINE CheckCompilerFlags(iExample,iReggieBuild,TESTCASE,TIMEDISCMETHOD,UseFV,Use2D,UsePARABOLIC)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: CodeNameLowCase,Examples,BuildFV,Build2D,BuildPARABOLIC
USE MOD_RegressionCheck_Vars,    ONLY: BuildSolver
USE MOD_RegressionCheck_Build,   ONLY: BuildConfiguration
USE MOD_RegressionCheck_Vars,    ONLY: BuildTESTCASE,BuildTIMEDISCMETHOD,BuildMPI,CodeNameUppCase,configuration_cmake
USE MOD_RegressionCheck_Build,   ONLY: GetFlagFromFile
USE MOD_RegressionCheck_Vars,    ONLY: RunContinue
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iReggieBuild
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: TESTCASE,TIMEDISCMETHOD
LOGICAL,INTENT(INOUT)          :: UseFV,Use2D,UsePARABOLIC      !> compiler flags currently used for ConvergenceTest
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: ExistFile                         !> file exists=.true., file does not exist=.false.
CHARACTER(LEN=255)             :: tempStr
LOGICAL                        :: UseMPI
!===================================================================================================================================
UseMPI       = .FALSE. ! default
UseFV        = .FALSE. ! default
Use2D    = .FALSE. ! default
UsePARABOLIC = .FALSE. ! default
IF((BuildSolver).AND.(.NOT.RunContinue))THEN ! only when freshly built files and folders are present
  TESTCASE       = BuildTESTCASE(iReggieBuild)
  TIMEDISCMETHOD = BuildTIMEDISCMETHOD(iReggieBuild)
  IF(ADJUSTL(TRIM(BuildMPI(       iReggieBuild))).EQ.'ON')UseMPI      =.TRUE.
  IF(ADJUSTL(TRIM(BuildFV(        iReggieBuild))).EQ.'ON')UseFV       =.TRUE.
  IF(ADJUSTL(TRIM(Build2D(        iReggieBuild))).EQ.'ON')Use2D       =.TRUE.
  IF(ADJUSTL(TRIM(BuildPARABOLIC( iReggieBuild))).EQ.'ON')UsePARABOLIC=.TRUE.
ELSE ! pre-compiled binary
  ! check if "configuration.cmake" exists and read specific flags from it
  !configuration_cmake=EXECPATH(1:INDEX(EXECPATH,'/',BACK = .TRUE.))//'configuration.cmake'
  !print*,"configuration_cmake=["//TRIM(configuration_cmake)//"]"
  INQUIRE(File=configuration_cmake,EXIST=ExistFile)
  IF(ExistFile) THEN
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! check if binary was compiled with MPI
    CALL  GetFlagFromFile(configuration_cmake,CodeNameUppCase//'_MPI',tempStr,BACK=.TRUE.)
    IF(ADJUSTL(TRIM(tempStr)).EQ.'ON')UseMPI=.TRUE.
    IF(TRIM(tempStr).EQ.'flag does not exist')CALL abort(&
      __STAMP__&
      ,CodeNameUppCase//'_MPI flag not found in configuration.cmake!',999,999.)
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! check if binary was compiled for a certain testcase
    CALL  GetFlagFromFile(configuration_cmake,CodeNameUppCase//'_TESTCASE',TESTCASE)
    IF(CodeNameLowCase.EQ.'boltzplatz')TESTCASE='default'! set default for boltzplatz (currently no testcases are implemented)
    IF(TRIM(TESTCASE).EQ.'flag does not exist')CALL abort(&
      __STAMP__&
      ,CodeNameUppCase//'_TESTCASE flag not found in configuration.cmake!',999,999.)
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! check if binary was compiled with a certain time integration method
    CALL GetFlagFromFile(configuration_cmake,CodeNameUppCase//'_TIMEDISCMETHOD',TIMEDISCMETHOD)
    IF(CodeNameLowCase.EQ.'flexi')TIMEDISCMETHOD='default'! set default for flexi (TIMEDISCMETHOD is not a compile flag)
    IF(TRIM(TIMEDISCMETHOD).EQ.'flag does not exist')CALL abort(&
      __STAMP__&
      ,CodeNameUppCase//'_TIMEDISCMETHOD flag not found in configuration.cmake!',999,999.)
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! check if binary was compiled for specific convergence tests
    CALL  GetFlagFromFile(configuration_cmake,CodeNameUppCase//'_FV',tempStr,BACK=.TRUE.)
    IF(ADJUSTL(TRIM(tempStr)).EQ.'ON')UseFV=.TRUE.
    CALL  GetFlagFromFile(configuration_cmake,CodeNameUppCase//'_2D',tempStr,BACK=.TRUE.)
    IF(ADJUSTL(TRIM(tempStr)).EQ.'ON')Use2D=.TRUE.
    CALL  GetFlagFromFile(configuration_cmake,CodeNameUppCase//'_PARABOLIC',tempStr,BACK=.TRUE.)
    IF(ADJUSTL(TRIM(tempStr)).EQ.'ON')UsePARABOLIC=.TRUE.
    ! -----------------------------------------------------------------------------------------------------------------------------
  ELSE
    SWRITE(UNIT_stdOut,'(A24,A)') ' ERROR: ','no "configuration.cmake" found at the location of the'//CodeNameLowCase//' binary.'
    SWRITE(UNIT_stdOut,'(A24,A)') ' configuration_cmake: ', TRIM(configuration_cmake)
    SWRITE(UNIT_stdOut,'(A24,L)') ' ExistFile            : ', ExistFile
    ERROR STOP 1
  END IF
END IF

IF(UseMPI.EQV..FALSE.)THEN ! parameter_reggie.ini supplied with MPI and possibly multiple runs with different numbers of mpi ranks
  IF(Examples(iExample)%MPIrun)THEN
    SWRITE(UNIT_stdOut,'(A)') ' MPIrun is set .TRUE. but the supplied binary was compiled without MPI. Setting MPIrun=.FALSE.'
  END IF
  Examples(iExample)%MPIrun=.FALSE. ! deactivate MPI for running reggie
END IF
IF(Examples(iExample)%MPIrun.EQV..FALSE.)Examples(iExample)%MPIthreadsN=1  ! set number of mpi ranks to 1

END SUBROUTINE CheckCompilerFlags


!===================================================================================================================================
!> 1.) Exclude subexamples set in parameter_reggie.ini when the binary is not suited for its execution
!> 2.) Check if SubExample is set correctly for possible convergence test (only once when iReggieBuild==1)
!===================================================================================================================================
SUBROUTINE CheckSubExample(iExample,iReggieBuild,TIMEDISCMETHOD)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iReggieBuild
CHARACTER(LEN=*),INTENT(IN)    :: TIMEDISCMETHOD
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! remove subexample (before displaying the case overview) for certain configurations: e.g. Preconditioner when running explicitly
IF((TRIM(TIMEDISCMETHOD).NE.'ImplicitO3').AND.(TRIM(Examples(iExample)%SubExample).EQ.'PrecondType'))THEN
  Examples(iExample)%SubExample       = ''
  Examples(iExample)%SubExampleNumber = 0
  Examples(iExample)%SubExampleOption(1:20) = '-' ! default option is nothing
END IF

! ConvergenceTest: allocate the array in which the L2 error norms are saved for convergence calculation
IF(iReggieBuild.EQ.1)THEN ! DON'T allocate the field "ConvergenceTestError" multiple times
  ! Check if SubExample is set correctly for possible convergence test
  IF(Examples(iExample)%ConvergenceTest)THEN ! Do ConvergenceTest
    ! Check number of sub examples used
    IF(Examples(iExample)%SubExampleNumber.LT.1)THEN ! less than 2 SubExample are executed, hence, a minimum of 2 points is not used
      SWRITE(UNIT_stdOut,'(A)') 'SubExampleNumber<2: Cannot allocate array in "CheckSubExample" Deactivating convergence test'
      Examples(iExample)%ConvergenceTest=.FALSE.
    ELSE
      ! check which type of convergence test is used
      IF(   ((TRIM(Examples(iExample)%SubExample).EQ.'N'       ).AND.(TRIM(Examples(iExample)%ConvergenceTestType).EQ.'p')))THEN
        ! p-convergence: constant number of DG cells, varying polynomial degree
        ! 1.) NumberOfCellsN == 1
        ! 2.) NumberOfCellsN < SubExampleNumber
        IF( (Examples(iExample)%NumberOfCellsN.EQ.1).AND.&
            (Examples(iExample)%NumberOfCellsN.LT.Examples(iExample)%SubExampleNumber) )THEN
          SWRITE(UNIT_stdOut,'(A)') ' Selected p-convergence: constant number of DG cells, varying polynomial degree'
          !ALLOCATE(Examples(iExample)%ConvergenceTestGridSize(1))
        ELSE
          Examples(iExample)%ConvergenceTest=.FALSE.
        END IF
      ELSEIF((TRIM(Examples(iExample)%SubExample).EQ.'MeshFile').AND.(TRIM(Examples(iExample)%ConvergenceTestType).EQ.'h'))THEN
        ! h-convergence: constant polynomial degree, varying number of DG cells
        ! NumberOfCellsN == SubExampleNumber
        IF(Examples(iExample)%NumberOfCellsN.EQ.Examples(iExample)%SubExampleNumber)THEN
          SWRITE(UNIT_stdOut,'(A)') ' Selected h-convergence: constant polynomial degree, varying number of DG cells'
          !ALLOCATE(Examples(iExample)%ConvergenceTestGridSize(Examples(iExample)%SubExampleNumber))
        ELSE
          Examples(iExample)%ConvergenceTest=.FALSE.
        END IF
      ELSE
        SWRITE(UNIT_stdOut,'(A)') 'SubExample not suitable for convergence test. Deactivating convergence test'
        SWRITE(UNIT_stdOut,'(A)') '  Choose from ...'
        SWRITE(UNIT_stdOut,'(A)') '  1.) SubExample=N        + ConvergenceTestType=p    OR'
        SWRITE(UNIT_stdOut,'(A)') '  2.) SubExample=MeshFile + ConvergenceTestType=h'
        Examples(iExample)%ConvergenceTest=.FALSE.
      END IF
      ! Allocate the arrays for L2 Error and GridSize
      IF(Examples(iExample)%ConvergenceTest)THEN ! Do ConvergenceTest was not changed to false above
        ! finally: allocate the fields
        ALLOCATE(Examples(iExample)%ConvergenceTestError(Examples(iExample)%SubExampleNumber,Examples(iExample)%nVar))
        Examples(iExample)%ConvergenceTestError   =0 ! default value
        ALLOCATE(Examples(iExample)%ConvergenceTestGridSize(Examples(iExample)%SubExampleNumber))
        Examples(iExample)%ConvergenceTestGridSize=0 ! default value
      END IF ! IF(ConvergenceTest==TRUE) - ConvergenceTest was not changed to false above
    END IF ! Examples(iExample)%SubExampleNumber.LT.1 - less than 2 SubExamples cannot yield a convergence rate
  END IF ! IF(ConvergenceTest==TRUE)
END IF ! iReggieBuild.EQ.1
END SUBROUTINE CheckSubExample



!===================================================================================================================================
!> check folder name and decide whether it can be executed with the current binary (e.g. testcases ...)
!===================================================================================================================================
SUBROUTINE CheckFolderName(iExample,TESTCASE,SkipFolder)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: CodeNameLowCase,ExampleNames
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
CHARACTER(LEN=*),INTENT(IN)    :: TESTCASE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)            :: SkipFolder
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FolderName                          !> path to a file or its name
INTEGER                        :: IndNum                           !> index number
!===================================================================================================================================
SkipFolder=.False.
IndNum=INDEX(ExampleNames(iExample), 'TESTCASE') ! look for TESTCASE- flag in folder name of example
IF(IndNum.GT.0)THEN ! folder name contains 'TESTCASE'
  FolderName=ExampleNames(iExample)               ! e.g. run_TESTCASE-taylorgreenvortex/
  FolderName=FolderName(IndNum+9:LEN(FolderName)) ! e.g. taylorgreenvortex/
  IndNum=INDEX(FolderName, '_')                   ! e.g. taylorgreenvortex_full/
  IF(IndNum.GT.0)FolderName=FolderName(1:IndNum-1)! e.g. taylorgreenvortex
  IndNum=INDEX(FolderName, '/')                   ! e.g. taylorgreenvortex/
  IF(IndNum.GT.0)FolderName=FolderName(1:IndNum-1)! e.g. taylorgreenvortex
  IF(FolderName.NE.TESTCASE)THEN ! e.g. taylorgreenvortex .NE. phill
    ! TESTCASE folder and non-TESTCASE binary
    SWRITE(UNIT_stdOut,'(A,2x,A)') ' TESTCASE not found in '//CodeNameLowCase//' binary ...skipping'
    SkipFolder=.True.!CYCLE
  ELSE
    SWRITE(UNIT_stdOut,'(A,2x,A)') ' TESTCASE is correct ...running' ! TESTCASE folder and TESTCASE binary
  END IF
ELSE ! folder name does not contain 'TESTCASE'
  FolderName='default' ! for non-TESTCASE setups, set the default settings
  IF(FolderName.NE.TESTCASE)THEN ! e.g. default .NE. phill
    ! non-TESTCASE folder, but TESTCASE binary
    SWRITE(UNIT_stdOut,'(A)') ' TESTCASE "default" not found in '//CodeNameLowCase//&
                                                     ' binary: TESTCASE=['//TRIM(TESTCASE)//'] ...skipping'
    SkipFolder=.True.!CYCLE
  ELSE
    ! non-TESTCASE folder and non-TESTCASE binary
    SWRITE(UNIT_stdOut,'(A)') ' TESTCASE "default" is correct: TESTCASE=['//TRIM(TESTCASE)//'] ...running' 
  END IF
END IF
END SUBROUTINE CheckFolderName


!===================================================================================================================================
!> get list of parameter files for running the simulation [parameter_ini] and [parameter_ini2]
!===================================================================================================================================
SUBROUTINE GetParameterFiles(iExample,TIMEDISCMETHOD,parameter_ini,parameter_ini2)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: CodeNameLowCase,Examples
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
CHARACTER(LEN=*),INTENT(IN)    :: TIMEDISCMETHOD
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: parameter_ini,parameter_ini2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: ExistFile                         !> file exists=.true., file does not exist=.false.
LOGICAL                        :: UseDSMC
CHARACTER(LEN=255)             :: TempStr
INTEGER                        :: iSTATUS,IndNum
!===================================================================================================================================
iSTATUS=0 ! nullify
! get list of parameter files for running the simulation
parameter_ini2=''
IF(TRIM(TIMEDISCMETHOD).EQ.'DSMC')THEN
  ! set main parameter.ini file: e.g. parameter_XXX_EQNSYSNAME.ini
  parameter_ini='parameter_'//CodeNameLowCase//'_'//TRIM(ADJUSTL(Examples(iExample)%EQNSYSNAME))//'_DSMC.ini'
  parameter_ini2='parameter_DSMC.ini'
ELSE ! standard flexi or PIC related simulation
  ! set main parameter.ini file: e.g. parameter_XXX_EQNSYSNAME.ini
  parameter_ini='parameter_'//CodeNameLowCase//'_'//TRIM(ADJUSTL(Examples(iExample)%EQNSYSNAME))//'.ini'

  ! PIC-DSMC run: Check for DSMC in "parameter_ini" -> "UseDSMC=T"
  CALL GetParameterFromFile(TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini),'UseDSMC',TempStr)
  IF(TempStr.EQ.'ParameterName does not exist'.OR.Tempstr.EQ.'file does not exist')THEN
    UseDSMC=.FALSE.
  ELSE
    CALL str2logical(TempStr,UseDSMC,iSTATUS)
  END IF
  IF(UseDSMC)THEN ! if UseDSMC is true, a parameter file for DSMC info needs to bespecified (must be located where the mesh is)
                  ! and must be named 'parameter_DSMC.ini'
    CALL GetParameterFromFile(TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini),'MeshFile',TempStr) ! find mesh file lcoation
    IndNum=INDEX(TempStr,'/',BACK = .TRUE.) ! get path without mesh file name (*.h5)
    IF(IndNum.GT.0)THEN
      TempStr=TempStr(1:IndNum)                    ! e.g. "./poisson/turner2013_mesh.h5" -> "./poisson/"
      IF(TRIM(ADJUSTL(TempStr)).EQ.'./')TempStr='' ! e.g. "./turner2013_mesh.h5"         -> "./"
    ELSE ! parameter file not located within a different directory
      TempStr=''
    END IF
    !parameter_folder ! get folder where the mesh is
    parameter_ini2=TRIM(ADJUSTL(TempStr))//'parameter_DSMC.ini'
  END IF
END IF
INQUIRE(File=TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini),EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A)')   ' ERROR: no parameter file found under : ['//TRIM(Examples(iExample)%PATH)//']'
  SWRITE(UNIT_stdOut,'(A)')   ' parameter_ini                        : ['//TRIM(parameter_ini)//']'
  ERROR STOP 1
END IF
IF(parameter_ini2.NE.'')THEN
  INQUIRE(File=TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini2),EXIST=ExistFile)
  IF(.NOT.ExistFile) THEN
    SWRITE(UNIT_stdOut,'(A)') ' ERROR: no parameter file found under : ['//TRIM(Examples(iExample)%PATH)//']'
    SWRITE(UNIT_stdOut,'(A)') ' parameter_ini2                       : ['//TRIM(parameter_ini2)//']'
    ERROR STOP 1
  END IF
END IF
END SUBROUTINE


!===================================================================================================================================
!> Output settings (before going into subexamples)
!===================================================================================================================================
SUBROUTINE PrintExampleInfo(iExample,EXECPATH,parameter_ini,parameter_ini2)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
CHARACTER(LEN=*),INTENT(IN)    :: EXECPATH,parameter_ini,parameter_ini2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')      ' EXECPATH:                       ['//TRIM(EXECPATH)//']'
SWRITE(UNIT_stdOut,'(A)')      ' EQNSYSNAME:                     ['//TRIM(Examples(iExample)%EQNSYSNAME)//']'
SWRITE(UNIT_stdOut,'(A,I6,A1)')' nVar:                           [',      Examples(iExample)%Nvar,']'
SWRITE(UNIT_stdOut,'(A)')      ' PATH (to example):              ['//TRIM(Examples(iExample)%PATH)//']'
SWRITE(UNIT_stdOut,'(A,L1,A1)')' MPIrun (run with MPI):          [',      Examples(iExample)%MPIrun,']'
IF(Examples(iExample)%MPIthreadsN.GT.1)THEN
  SWRITE(UNIT_stdOut,'(A)')    ' MPIthreads (number of threads): ['//TRIM(Examples(iExample)%MPIthreadsStr(1          ))//&
                                                 ' to '//TRIM(Examples(iExample)%MPIthreadsStr(Examples(iExample)%MPIthreadsN))//']'
ELSE
  SWRITE(UNIT_stdOut,'(A)')    ' MPIthreads (number of threads): ['//TRIM(Examples(iExample)%MPIthreadsStr(1))//']'
END IF
SWRITE(UNIT_stdOut,'(A)')      ' Reference:                      ['//TRIM(Examples(iExample)%ReferenceFile)//']'
SWRITE(UNIT_stdOut,'(A)')      ' Reference Norm:                 ['//TRIM(Examples(iExample)%ReferenceNormFile)//']'
SWRITE(UNIT_stdOut,'(A)')      ' State:                          ['//TRIM(Examples(iExample)%H5DIFFReferenceStateFile)//']'
SWRITE(UNIT_stdOut,'(A)')      ' HDF5 dataset:                   ['//TRIM(Examples(iExample)%H5DIFFReferenceDataSetName)//']'
SWRITE(UNIT_stdOut,'(A)')      ' Restart:                        ['//TRIM(Examples(iExample)%RestartFileName)//']'
SWRITE(UNIT_stdOut,'(A)')      ' Example%SubExample:             ['//TRIM(Examples(iExample)%SubExample)//']'
SWRITE(UNIT_stdOut,'(A,I6,A1)')' Example%SubExampleNumber:       [',      Examples(iExample)%SubExampleNumber,']'
SWRITE(UNIT_stdOut,'(A)')      ' parameter files:                ['//TRIM(parameter_ini)//' '//TRIM(parameter_ini2)//']'
END SUBROUTINE PrintExampleInfo



!===================================================================================================================================
!> check if executable is compiled with correct TESTCASE (e.g. for tylorgreenvortex)
!===================================================================================================================================
SUBROUTINE SetParameters(iExample,parameter_ini,UseFV,Use2D,UsePARABOLIC,SkipFolder)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples,CodeNameUppCase
USE MOD_RegressionCheck_Vars,    ONLY: ExampleNames
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(OUT)            :: SkipFolder
INTEGER,INTENT(IN)             :: iExample
LOGICAL,INTENT(INOUT)          :: UseFV,Use2D,UsePARABOLIC      !> compiler flags currently used for ConvergenceTest
CHARACTER(LEN=*),INTENT(INOUT) :: parameter_ini
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!LOGICAL                        :: UseFV,Use2D,UsePARABOLIC      !> compiler flags currently used for ConvergenceTest
!LOGICAL                        :: ExistFile                         !> file exists=.true., file does not exist=.false.
!CHARACTER(LEN=255)             :: FileName                          !> path to a file or its name
!LOGICAL                        :: UseMPI
INTEGER                        :: IndNum
INTEGER                        :: iSubExample
CHARACTER(LEN=255)             :: MeshFile
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') " SUBROUTINE SetParameters ..."
SkipFolder=.FALSE.
! Check specific compile flags for ConvergenceTest
IF(Examples(iExample)%ConvergenceTest)THEN ! Do ConvergenceTest
  ! Check for UseFV (finite volume operator cells)
  !print*,"UseFV=",UseFV
  SWRITE(UNIT_stdOut,'(A15,L1,A2)',ADVANCE='NO')" UseFV       =[",UseFV,"] "
  IF(UseFV)THEN
    IF(    CodeNameUppCase.EQ.'FLEXI')THEN
      CALL SetSubExample(iExample,-1,parameter_ini,'IndicatorType','33')
      IF(    TRIM(Examples(iExample)%ConvergenceTestType).EQ.'h')THEN ! h-convergence    
        Examples(iExample)%ConvergenceTestValue=1.4 ! set redueced order of convergence due to FV cells
        SWRITE(UNIT_stdOut,'(A18,A)')"","Setting option in parameter_reggie.ini: h-ConvergenceTestValue=1.4 (reduced mean order)"
      ELSEIF(TRIM(Examples(iExample)%ConvergenceTestType).EQ.'p')THEN ! p-convergence
        SkipFolder=.TRUE.
        SWRITE(UNIT_stdOut,'(A18,A)')"","Skipping example because [p-convergence] and [FV=ON]"
        RETURN
      END IF
    ELSEIF(CodeNameUppCase.EQ.'BOLTZPLATZ')THEN
      SWRITE(UNIT_stdOut,'(A)')"Setting option in parameter.ini: NOTHING (not implemented yet)"
      ! DO NOTHING
    END IF
  ELSE
    IF(    CodeNameUppCase.EQ.'FLEXI')THEN
      CALL SetSubExample(iExample,-1,parameter_ini,'IndicatorType','0')
    ELSEIF(CodeNameUppCase.EQ.'BOLTZPLATZ')THEN
      SWRITE(UNIT_stdOut,'(A)')"Setting option in parameter.ini: NOTHING (not implemented yet)"
      ! DO NOTHING
    END IF
  END IF

  ! Check for UsePARABOLIC==OFF -> Euler simulation
  SWRITE(UNIT_stdOut,'(A15,L1,A2)',ADVANCE='NO')" UsePARABOLIC=[",UsePARABOLIC,"] "
  IF(UsePARABOLIC)THEN
    CALL SetSubExample(iExample,-1,parameter_ini,'IniExactFunc','4') !Sine wave in vel
  ELSE
    IF(    CodeNameUppCase.EQ.'FLEXI')THEN
      CALL SetSubExample(iExample,-1,parameter_ini,'IniExactFunc','2') !Sine wave in density
    ELSEIF(CodeNameUppCase.EQ.'BOLTZPLATZ')THEN
      CALL SetSubExample(iExample,-1,parameter_ini,'IniExactFunc','12')
    END IF
  END IF

 
  ! Check for 2D version of the code
  SWRITE(UNIT_stdOut,'(A15,L1,A2)',ADVANCE='NO')" Use2D   =[",Use2D,"] "
  SWRITE(UNIT_stdOut,'(A)')"Setting option in parameter.ini: NOTHING (not implemented yet)"
  IF(Use2D)THEN
      IndNum=INDEX(ExampleNames(iExample), 'mortar') ! look for mortar in folder name of example
      IF(IndNum.GT.0)THEN ! folder name contains 'mortar'
        SkipFolder=.TRUE.
        SWRITE(UNIT_stdOut,'(A18,A)')"","Skipping example because [mortar] and [2D=ON]"
        RETURN
      END IF
      CALL SetSubExample(iExample,-1,parameter_ini,'IniExactFunc','2') !Sine wave in density
      IF  (TRIM(Examples(iExample)%SubExample).EQ.'MeshFile') THEN
         DO iSubExample = 1, MAX(1,Examples(iExample)%SubExampleNumber) 
            IndNum=INDEX(TRIM(Examples(iExample)%SubExampleOption(iSubExample)),'3D')
            IF(IndNum.GT.0)THEN
              Examples(iExample)%SubExampleOption(iSubExample)(IndNum:IndNum)='2'
              SWRITE(UNIT_stdOut,'(A)')"Setting SubExampleOption for MeshFile from 3D to 2D"
            END IF
         END DO
      END IF

  ELSE
    IF  (TRIM(Examples(iExample)%SubExample).EQ.'MeshFile') THEN
       DO iSubExample = 1, MAX(1,Examples(iExample)%SubExampleNumber) 
          IndNum=INDEX(TRIM(Examples(iExample)%SubExampleOption(iSubExample)),'2D')
          IF(IndNum.GT.0)THEN
            Examples(iExample)%SubExampleOption(iSubExample)(IndNum:IndNum)='3'
            SWRITE(UNIT_stdOut,'(A)')"Setting SubExampleOption for MeshFile from 2D to 3D"
          END IF
       END DO
    END IF
  END IF


END IF

! Check for 2D version of the code: use different file names for 2D or 3D code version
SWRITE(UNIT_stdOut,'(A15,L1,A2)',ADVANCE='NO')" Use2D   =[",Use2D,"] "
SWRITE(UNIT_stdOut,'(A)')"Setting option in parameter.ini: MESHFILE (if file name contains '2D' it will be exchanged with '3D)"
! read MeshFile from parameter_ini and search for "3D", then substitute with 2D
CALL GetParameterFromFile(TRIM(Examples(iExample)%PATH)//parameter_ini,'MeshFile',meshFile)
IF(Use2D)THEN
  IndNum=INDEX(meshFile,'3D')
  IF(IndNum.GT.0)THEN
    meshFile(IndNum:IndNum)='2'
    CALL SetSubExample(iExample,-1,parameter_ini,'MeshFile',meshFile)
  END IF
ELSE
  IndNum=INDEX(meshFile,'2D')
  IF(IndNum.GT.0)THEN
    meshFile(IndNum:IndNum)='3'
    CALL SetSubExample(iExample,-1,parameter_ini,'MeshFile',meshFile)
  END IF
END IF

END SUBROUTINE SetParameters




!===================================================================================================================================
!> Set the SubExample in the parameter.ini file
!> Search for "Examples(iExample)%SubExample" + "=" within the parameter.ini file and use "sed" to exchange the line by
!>            "Examples(iExample)%SubExample" + "=" + "Examples(iExample)%SubExampleOption(iSubExample)"
!> Example: TimeDiscMethod = carpenterrk4-5
!===================================================================================================================================
SUBROUTINE SetSubExample(iExample,iSubExample,parameter_ini,ChangeOption,ChangeParameter)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples
USE MOD_RegressionCheck_tools,   ONLY: CheckFileForString
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iExample,iSubExample
CHARACTER(LEN=*),INTENT(IN)          :: parameter_ini
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: ChangeOption,ChangeParameter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iSTATUS           !> status
CHARACTER(LEN=500)                   :: SYSCOMMAND        !> string to fit the system command
LOGICAL                              :: ExistStringInFile !> search a file for the existence of a string
INTEGER                              :: MODE              !> 1: change SubExample option parameter in parameter.ini
                                                          !> 2: change specific parameter in parameter.ini
!===================================================================================================================================
MODE=0 ! default
iSTATUS=0 ! nullify
IF(iSubExample.GT.0)MODE=1
IF(PRESENT(ChangeOption).AND.PRESENT(ChangeParameter))MODE=2

IF((MODE.EQ.1).AND.(Examples(iExample)%SubExampleNumber.GT.0))THEN ! SubExample has been specified
  SWRITE(UNIT_stdOut,'(A)')" "
  SWRITE(UNIT_stdOut,'(A,I2,A,A)')" SubExampleOption(",iSubExample,")=",TRIM(Examples(iExample)%SubExampleOption(iSubExample))
  ! check if SubExampleOption is possible in destination file, e.g. in 'parameter_flexi_navierstokes.ini'
  CALL CheckFileForString(TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini),&
                          TRIM(Examples(iExample)%SubExample)//'=',ExistStringInFile)
ELSEIF(MODE.EQ.2)THEN
  CALL CheckFileForString(TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini),&
                          TRIM(ChangeOption)//'=',ExistStringInFile)
ELSE
  RETURN ! do nothing
END IF

! use 'sed' to replace a parameter in the *.ini file with specific settings
IF(ExistStringInFile)THEN
  IF(MODE.EQ.1)THEN
    SYSCOMMAND=      'cd '//TRIM(Examples(iExample)%PATH)//& ! write the current SubExampleOption(iSubExample) to parameter_ini
     ' && sed -i -e "s/.*'//TRIM(Examples(iExample)%SubExample)//&
                    '=.*/'//TRIM(Examples(iExample)%SubExample)//&
                       '='//TRIM(Examples(iExample)%SubExampleOption(iSubExample))//&
                     '/" '//TRIM(parameter_ini)
    ! E ! 
    ! X ! SYSCOMMAND= cd ~/Flexi/flexi/build_reggie.dev/../regressioncheck/examples/run_basic/ 
    ! A ! && sed -i -e "s/.*TimeDiscMethod
    ! M !               =.*/TimeDiscMethod
    ! P !               =carpenterrk4-5
    ! L !               /" parameter_flexi_navierstokes.ini 
    ! E !
  ELSEIF(MODE.EQ.2)THEN
    SYSCOMMAND=      'cd '//TRIM(Examples(iExample)%PATH)//& ! write the current SubExampleOption(iSubExample) to parameter_ini
     ' && sed -i -e "s/.*'//TRIM(ChangeOption)//&
                    '=.*/'//TRIM(ChangeOption)//&
                       '='//TRIM(ChangeParameter)//&
                     '/" '//TRIM(parameter_ini)
    ! E ! 
    ! X ! SYSCOMMAND= cd ~/Flexi/flexi/build_reggie.dev/../regressioncheck/examples/run_basic/ 
    ! A ! && sed -i -e "s/.*IniExactFunc
    ! M !               =.*/IniExactFunc
    ! P !               =4
    ! L !               /" parameter_flexi_navierstokes.ini 
    ! E !
    SWRITE(UNIT_stdOut,'(A)')"Setting option in parameter.ini: "//TRIM(ChangeOption)//"="//TRIM(ChangeParameter)
  END IF
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
ELSE ! could not change the SubExampleOption parameter, because the the string [Examples(iExample)%SubExample + '='] was not 
     ! found in the destination file ( e.g. "parameter_flexi_navierstokes.ini"). The equation mark "=" must be placed directly
     ! behind the string that is to be changed!
  IF(MODE.EQ.1)THEN
    Examples(iExample)%SubExampleOption(iSubExample)='**failed**'
    SWRITE(UNIT_stdOut,'(A,I2,A,A)')" SubExampleOption(",iSubExample,")=",TRIM(Examples(iExample)%SubExampleOption(iSubExample))
    SWRITE(UNIT_stdOut,'(A)')" The subexample has failed: Could not find the required string in destination parameter file."
    SWRITE(UNIT_stdOut,'(A,A,A)')"   Required String            : [",TRIM(Examples(iExample)%SubExample),"=]"
    SWRITE(UNIT_stdOut,'(A,A,A)')"   Destination parameter file : [",TRIM(parameter_ini),"]"
    SWRITE(UNIT_stdOut,'(A,A,A)')" Maybe the equation mark '=' is not placed directly after the parameter that is to be changed"
  ELSEIF(MODE.EQ.2)THEN
    !Examples(iExample)%SubExampleOption(iSubExample)='**failed**'
    SWRITE(UNIT_stdOut,'(A,A1,A)')TRIM(ChangeOption),"=",TRIM(ChangeParameter)
    SWRITE(UNIT_stdOut,'(A)')" The Parameter change has failed: Could not find the required string in destination parameter file."
    SWRITE(UNIT_stdOut,'(A,A,A)')"   Required String            : [",TRIM(ChangeOption),"=]"
    SWRITE(UNIT_stdOut,'(A,A,A)')"   Destination parameter file : [",TRIM(parameter_ini),"]"
    SWRITE(UNIT_stdOut,'(A,A,A)')" Maybe the equation mark '=' is not placed directly after the parameter that is to be changed"
  END IF
END IF
END SUBROUTINE SetSubExample


!==================================================================================================================================
!> Remove State-files, std.out and err.out files.
!> MODE=0: INITIAL -> delete pre-existing files and folders
!> MODE 1: Delete pre-existing data files before running the code
!> MODE 2: If the example is computed successfully clean up afterwards
!==================================================================================================================================
SUBROUTINE CleanFolder(iExample,MODE)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,MODE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=500)             :: SYSCOMMAND
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)             :: tmp
INTEGER                        :: iSTATUS,ioUnit
!==================================================================================================================================
iSTATUS=0 ! nullify
SELECT CASE(MODE)
  CASE(0) ! MODE=0: INITIAL -> delete pre-existing files and folders
    ! delete "std_files_*" folder
    SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm std_files_* -r > /dev/null 2>&1'
    CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  CASE(1) ! delete pre-existing files before computation
    ! delete pre-existing data files before running the code 
    ! 1.) Files created during simulation which are needed by, e.g., "IntegrateLine" comparison
    IF(Examples(iExample)%IntegrateLine)THEN
      SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm '//TRIM(Examples(iExample)%IntegrateLineFile)//' > /dev/null 2>&1'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! delete, e.g., "TGVAnalysis.dat" or "Database.csv"
    END IF
  CASE(2) ! delete existing files after computation
    ! delete all *.out files
    SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm *.out > /dev/null 2>&1'
    CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
    IF(iSTATUS.NE.0)THEN
      SWRITE(UNIT_stdOut,'(A)')' CleanFolder('//TRIM(Examples(iExample)%PATH)//'): Could not remove *.out files!'
    END IF
    ! delete all *State* files except *reference* state files
    IF((Examples(iExample)%H5DIFFReferenceStateFile.EQ.'').AND. &
       (Examples(iExample)%RestartFileName.EQ.'') ) THEN
      SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm *State* > /dev/null 2>&1'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      IF(iSTATUS.NE.0)THEN
        SWRITE(UNIT_stdOut,'(A)')' CleanFolder('//TRIM(Examples(iExample)%PATH)//'): Could not remove *State* files!'
      END IF
    ELSE
      ! create list of all *State* files and loop them: don't delete *reference* files
      SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && ls *State* > tmp.txt'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      IF(iSTATUS.NE.0)THEN
        SWRITE(UNIT_stdOut,'(A)')' CleanFolder('//TRIM(Examples(iExample)%PATH)//'): Could not create tmp.txt!'
      END IF
      ! read tmp.txt | list of directories if regressioncheck/examples
      FileName=TRIM(Examples(iExample)%PATH)//'tmp.txt'
      OPEN(NEWUNIT = ioUnit, FILE = FileName, STATUS ="OLD", IOSTAT = iSTATUS )
      DO 
        READ(ioUnit,FMT='(A)',IOSTAT=iSTATUS) tmp
        IF (iSTATUS.NE.0) EXIT
        IF((Examples(iExample)%H5DIFFReferenceStateFile.NE.TRIM(tmp)).AND. & ! skip H5DIFFReferenceStateFile and RestartFileName
           (Examples(iExample)%RestartFileName.NE.TRIM(tmp)) ) THEN
           SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm '//TRIM(tmp)
           CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
           IF(iSTATUS.NE.0) THEN
             SWRITE(UNIT_stdOut,'(A)')  ' CleanFolder('//TRIM(Examples(iExample)%PATH)//'): Could not remove state file ',TRIM(tmp)
           END IF
        END IF
      END DO
      CLOSE(ioUnit)
      ! clean tmp.txt
      SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm tmp.txt'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      IF(iSTATUS.NE.0) THEN
        SWRITE(UNIT_stdOut,'(A)')  ' CleanFolder('//TRIM(Examples(iExample)%PATH)//'): Could not remove tmp.txt'
      END IF
    END IF
    ! delete "userblock.tmp"
    SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm userblock.tmp > /dev/null 2>&1'
    CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  CASE DEFAULT
    SWRITE(UNIT_stdOut,'(A,I3,A)') ' SUBROUTINE CleanFolder: MODE=',MODE,' does not exist!'
    ERROR STOP 1
END SELECT

END SUBROUTINE CleanFolder


!===================================================================================================================================
!> Execute the binary and check if the attempt was successful
!===================================================================================================================================
SUBROUTINE RunTheCode(iExample,iSubExample,iScaling,iRun,MPIthreadsStr,EXECPATH,parameter_ini,parameter_ini2,SkipComparison)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples,GlobalRunNumber
USE MOD_RegressionCheck_tools,   ONLY: AddError
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iSubExample,iScaling,iRun
CHARACTER(LEN=*),INTENT(IN)    :: parameter_ini,parameter_ini2,EXECPATH
CHARACTER(LEN=*),INTENT(INOUT) :: MPIthreadsStr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)            :: SkipComparison
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: ComputationSTATUS                               !> simulation successful/failed
INTEGER                        :: iSTATUS                                         !> status
CHARACTER(LEN=1000)            :: SYSCOMMAND                                      !> string to fit the system command
CHARACTER(LEN=255)             :: FileSuffix,FolderSuffix,tempStr                 !> auxiliary vars for file and folder names
CHARACTER(LEN=255)             :: FileOutFolderName,StdOutFileName,ErrOutFileName !> new std.out and err.out files
INTEGER                        :: MPIthreadsInteger,PolynomialDegree              !> character to integer auxiliaray vars
CHARACTER(LEN=255)             :: ComputationResult                               !> auxiliary var to determine a failed computation
                                                                                  !>  even if EXITSTATUS=0
!===================================================================================================================================
SkipComparison=.FALSE.
MPIthreadsInteger=1 ! default: single run
iSTATUS=0 ! nullify
ComputationSTATUS=0 ! nullify
! -----------------------------------------------------------------------------------------------------------------------
! Run the Code
! -----------------------------------------------------------------------------------------------------------------------
IF(Examples(iExample)%MPIrun)THEN ! use "mpirun"
  !IF(Examples(iExample)%MPIthreads.GT.1)THEN
  MPIthreadsStr=ADJUSTL(TRIM(Examples(iExample)%MPIthreadsStr(iScaling))) ! copy string containing the number of MPI threads to temp
  CALL str2int(ADJUSTL(TRIM(MPIthreadsStr)),MPIthreadsInteger,iSTATUS) ! sanity check if the number of threads is correct
  IF((MPIthreadsInteger.LE.0).OR.(iSTATUS.NE.0))CALL abort(&
      __STAMP__&
      ,'RunTheCode(): Number of MPI threads is corrupt = '//ADJUSTL(TRIM(MPIthreadsStr)))
  IF((iScaling.GT.1).AND.(iRun.EQ.1))THEN
    SWRITE(UNIT_stdOut,'(A,A)')" Examples(iExample)%MPIthreads=",Examples(iExample)%MPIthreadsStr(iScaling)
  END IF
  tempStr='' ! default
  IF(Examples(iExample)%MPIcommand.EQ.'mpirun')THEN
    Examples(iExample)%MPIcommand='mpirun -np'
  ELSEIF(Examples(iExample)%MPIcommand.EQ.'aprun')THEN
    Examples(iExample)%MPIcommand='aprun -n' !-N $CoresPerNode'
    IF(MPIthreadsInteger.GT.24)THEN
      tempStr='-N 24'
    ELSE
      tempStr='-N '//ADJUSTL(TRIM(MPIthreadsStr))
    END IF
  END IF


  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && '//TRIM(Examples(iExample)%MPIcommand)//' '//&
                       ADJUSTL(TRIM(MPIthreadsStr))//' '//ADJUSTL(TRIM(tempStr))//' '//TRIM(EXECPATH)//' '//&
           TRIM(parameter_ini)//' '//TRIM(parameter_ini2)//' '//TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
ELSE ! single run
  MPIthreadsStr='-' ! no MPI threads to be used -> single run
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && '//TRIM(EXECPATH)//' '//&
           TRIM(parameter_ini)//' '//TRIM(parameter_ini2)//' '//TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
END IF
GlobalRunNumber=GlobalRunNumber+1
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=ComputationSTATUS) ! run the code
CALL GetParameterFromFile(TRIM(Examples(iExample)%PATH)//'std.out',&
                          'Program abort caused',ComputationResult,DoDisplayInfo=.FALSE.)
IF((TRIM(ComputationResult).NE.'ParameterName does not exist').AND.(TRIM(ComputationResult).NE.'file does not exist'))THEN
  ComputationSTATUS=1 ! an abort was caused (works for MPI=ON/OFF)
  SWRITE(UNIT_stdOut, '(A)')' Program abort caused '//TRIM(ComputationResult)
END IF
! -----------------------------------------------------------------------------------------------------------------------
! copy the std.out and err.out files to sub folder (std_filed_.....)
! -----------------------------------------------------------------------------------------------------------------------
IF(Examples(iExample)%SubExample.EQ.'N')THEN ! when polynomial degree "N" is the SubExample, set special file name
  CALL str2int(Examples(iExample)%SubExampleOption(iSubExample),PolynomialDegree,iSTATUS)
  WRITE(FileSuffix, '(A10,I8.8,A2,I4.4,A5,I4.4)') 'MPIthreads',MPIthreadsInteger,'_N',PolynomialDegree,'_iRun',iRun
  WRITE(FolderSuffix,'(A2,I4.4)')'_N',PolynomialDegree
ELSE
  IF(Examples(iExample)%MPIrun)THEN
    WRITE(FileSuffix, '(A10,I8.8,A6,I4.4,A5,I4.4)') 'MPIthreads',MPIthreadsInteger,'_SubEx',iSubExample,'_iRun',iRun
  ELSE
    WRITE(FileSuffix, '(I4.4,A6,I4.4,A5,I4.4)') iScaling,'_SubEx',iSubExample,'_iRun',iRun
  END IF
  WRITE(FolderSuffix,'(A6,I4.4)')'_SubEx',iSubExample
END IF
FileOutFolderName='std_files'//TRIM(FolderSuffix) ! folder for storing all std.out files
StdOutFileName='std-'//TRIM(FileSuffix)//'.out'
ErrOutFileName='err-'//TRIM(FileSuffix)//'.out'
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && cp std.out '//TRIM(StdOutFileName)
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! copy the std.out file
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && cp err.out '//TRIM(ErrOutFileName)
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! copy the err.out file

SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && mkdir '//TRIM(FileOutFolderName)//' > /dev/null 2>&1'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)

SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && mv '//TRIM(StdOutFileName)//' '//TRIM(FileOutFolderName)//'/'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! move new std.out file
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && mv '//TRIM(ErrOutFileName)//' '//TRIM(FileOutFolderName)//'/'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! move new err.out file
! -----------------------------------------------------------------------------------------------------------------------
! was the run successful? (iSTATUS=0)
! -----------------------------------------------------------------------------------------------------------------------
IF(ComputationSTATUS.EQ.0)THEN ! Computation successful
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='no')  ' Successful Computation ...'
  CALL AddError(MPIthreadsStr,'Successful Computation',iExample,iSubExample,ErrorStatus=0,ErrorCode=0)
ELSE ! Computation failed
  SWRITE(UNIT_stdOut,'(A)')   ' Failed Computation ...'
  SWRITE(UNIT_stdOut,'(A,A)') ' Out-file : ', TRIM(Examples(iExample)%PATH)//TRIM(FileOutFolderName)//'/'//TRIM(StdOutFileName)
  SWRITE(UNIT_stdOut,'(A,A)') ' Errorfile: ', TRIM(Examples(iExample)%PATH)//TRIM(FileOutFolderName)//'/'//TRIM(ErrOutFileName)
  CALL AddError(MPIthreadsStr,'Failed Computation',iExample,iSubExample,ErrorStatus=1,ErrorCode=2)
  SkipComparison=.TRUE. ! when a computation fails no useful output is created for comparison, continue the subexample cycle
END IF
END SUBROUTINE RunTheCode




END MODULE MOD_RegressionCheck_Run
