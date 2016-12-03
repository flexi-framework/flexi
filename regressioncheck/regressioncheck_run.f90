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

PUBLIC::PerformRegressionCheck
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Routine which performs the actual regressioncheck. It triggers the builds and execute commands. Additionally, it performs
!> the checks for L2-error norms, h5-diff and runtime
!==================================================================================================================================
SUBROUTINE PerformRegressionCheck()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Compare, ONLY: CompareResults
USE MOD_RegressionCheck_Tools,   ONLY: InitExample
USE MOD_RegressionCheck_Vars,    ONLY: nExamples,ExampleNames,Examples,EXECPATH,RuntimeOptionType
USE MOD_RegressionCheck_Tools,   ONLY: CleanExample
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
INTEGER                        :: iExample                          !> loop index for example
INTEGER                        :: N_compile_flags                   !> number of compile-flags
INTEGER                        :: iReggieBuild,nReggieBuilds ! field handler unit and ??
INTEGER                        :: iSubExample,iScaling,iRun
LOGICAL                        :: SkipExample,SkipBuild,ExitBuild,SkipFolder,SkipComparison
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Performing tests ...'
ReggieBuildExe=''
!==================================================================================================================================
DO iExample = 1, nExamples ! loop level 1 of 5
!==================================================================================================================================
  CALL CheckExampleName(ExampleNames(iExample),RuntimeOptionType,SkipExample)
  IF(SkipExample)CYCLE ! ignore the example folder and continue with the next

  ! set the build configuration environment when BuildSolver=.TRUE.
  CALL GetnReggieBuilds(iExample,ReggieBuildExe,N_compile_flags,nReggieBuilds)
!==================================================================================================================================
  DO iReggieBuild = 1, nReggieBuilds ! loop level 2 of 5: cycle the number of build configurations (no configuration = only 1 run)
!==================================================================================================================================
    ! Get code binary (build or find it)
    CALL GetCodeBinary(iExample,iReggieBuild,nReggieBuilds,N_compile_flags,ReggieBuildExe,SkipBuild,ExitBuild)
    IF(SkipBuild)CYCLE ! invalid reggie build but not last reggie build
    IF(ExitBuild)EXIT  ! last reggie build -> exit ("cycle" would start an infinite loop)
  
    ! read the parameters for the current example (parameter_reggie.ini)
    CALL InitExample(TRIM(Examples(iExample)%PATH),LEN(TRIM(Examples(iExample)%PATH)),Examples(iExample))

    ! depending on the equation system -> get different Nvar 
    CALL GetNvar(iExample,iReggieBuild)

    ! check if executable is compiled with correct TESTCASE (e.g. for tylorgreenvortex)
    CALL CheckCompilerFlags(iReggieBuild,TESTCASE,TIMEDISCMETHOD)

    ! remove subexample (before printing the case overview) for certain configurations: e.g. Preconditioner when running explicitly
    CALL CheckSubExample(iExample,TIMEDISCMETHOD)

    ! check folder name and decide whether it can be executed with the current binary (e.g. testcases ...)
    CALL CheckFolderName(iExample,TESTCASE,SkipFolder)
    IF(SkipFolder)CYCLE ! e.g. TESTCASE folder and non-TESTCASE binary or vice versa

    ! get list of parameter files for running the simulation
    CALL GetParameterFiles(iExample,TIMEDISCMETHOD,parameter_ini,parameter_ini2)

    ! Output settings (before going into subexamples)
    CALL PrintExampleInfo(iExample,EXECPATH,parameter_ini,parameter_ini2)

!==================================================================================================================================
    DO iSubExample = 1, MAX(1,Examples(iExample)%SubExampleNumber) ! loop level 3 of 5: SubExamples (e.g. different TimeDiscMethods)
!==================================================================================================================================
      ! Set the SubExample in the parameter.ini file
      CALL SetSubExample(iExample,iSubExample,parameter_ini)

      ! delete pre-existing data files before running the code (e.g. "TGVAnalysis.dat" or "Database.csv")
      CALL CleanFolder(iExample)

!==================================================================================================================================
      DO iScaling = 1, Examples(iExample)%MPIthreadsN ! loop level 4 of 5: multiple MPI runs with different MPI threads
!==================================================================================================================================
!==================================================================================================================================
        DO iRun = 1, Examples(iExample)%nRuns ! loop level 5 of 5: repeat the same run multiple times
!==================================================================================================================================
          CALL RunTheCode(iExample,iSubExample,iScaling,iRun,EXECPATH,parameter_ini,parameter_ini2,SkipComparison)
          IF(SkipComparison)CYCLE ! the execution has failed, no comparisons are needed

          ! compare the results and write error messages for the current case
          CALL CompareResults(iExample,iSubExample)

          ! IF all comparisons are successful the error status is 0 -> delete created files in CleanExample(iExample)
          IF(Examples(iExample)%ErrorStatus.EQ.0) CALL CleanExample(iExample)
        END DO ! iScalingRuns = 1, Examples(iExample)%nRuns
      END DO ! iScaling = 1, Examples(iExample)%MPIthreadsN
    END DO ! iSubExample = 1, MAX(1,SubExampleNumber) (for cases without specified SubExamples: SubExampleNumber=0)
  END DO ! iReggieBuild = 1, nReggieBuilds
END DO ! iExample=1,nExamples

END SUBROUTINE PerformRegressionCheck


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
USE MOD_RegressionCheck_Vars,    ONLY: Examples,EXECPATH,BuildEQNSYS,BuildSolver,CodeNameUppCase,CodeNameLowCase
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
CHARACTER(LEN=255)             :: FileName                          ! path to a file or its name
CHARACTER(LEN=255)             :: temp,temp2                        ! temp variables for read in of file lines
LOGICAL                        :: ExistFile                         ! file exists=.true., file does not exist=.false.
INTEGER                        :: iSTATUS                           ! status
INTEGER                        :: ioUnit                            ! IO channel
!===================================================================================================================================
ioUnit=GETFREEUNIT()
IndNum=INDEX(EXECPATH, '/')
IF(IndNum.GT.0)THEN
  IndNum=INDEX(EXECPATH,'/',BACK = .TRUE.) ! get path without binary
  FileName=EXECPATH(1:IndNum)//'configuration.cmake'
  INQUIRE(File=FileName,EXIST=ExistFile)
  IF(ExistFile) THEN
    OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
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
    IF(BuildSolver)THEN ! get EQNSYSNAME from cmake build configuration settings
      Examples(iExample)%EQNSYSNAME=BuildEQNSYS(iReggieBuild)
    ELSE ! stop 
      SWRITE(UNIT_stdOut,'(A12,A)')     ' ERROR: ','no "configuration.cmake" found at the location of the '&
                                                                                                   //CodeNameLowCase//' binary.'
      SWRITE(UNIT_stdOut,'(A12,A)')  ' FileName: ', TRIM(FileName)
      SWRITE(UNIT_stdOut,'(A12,L)') ' ExistFile: ', ExistFile
      ERROR STOP '-1'
    END IF
  END IF
END IF
SELECT CASE (TRIM(Examples(iExample)%EQNSYSNAME))
  CASE ('navierstokes')  
    Examples(iExample)%Nvar=5
  CASE ('linearscalaradvection')  
    Examples(iExample)%Nvar=1
  CASE ('maxwell')  
    Examples(iExample)%Nvar=8
  CASE DEFAULT
    Examples(iExample)%Nvar=-1
    SWRITE(UNIT_stdOut,'(A)')   ' ERROR: missing case select for this '&
                                                                   //CodeNameUppCase//'_EQNSYSNAME with appropriate Nvar. Fix it by'
    SWRITE(UNIT_stdOut,'(A)')   '        adding the correct line of code to ../regressioncheck/regressioncheck_run.f90'
    SWRITE(UNIT_stdOut,'(A,A)') '        Examples(iExample)%EQNSYSNAME=',Examples(iExample)%EQNSYSNAME
    ERROR STOP '77'
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
  SWRITE(UNIT_stdOut,*) ''
END IF
SWRITE(UNIT_stdOut,'(A,2x,A50)',ADVANCE='no') ' Example-Name: ',  TRIM(ExampleName)
IF(dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType)))).NE.RuntimeOptionType)THEN
  SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...skipping'
  SkipExample=.TRUE.
ELSE
  SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...running'
  SkipExample=.FALSE.
END IF
END SUBROUTINE CheckExampleName


!===================================================================================================================================
!> IF build-mode: read the configurations.reggie and determine the number of builds
!> ELSE         : nReggieBuilds=1
!===================================================================================================================================
SUBROUTINE GetnReggieBuilds(iExample,ReggieBuildExe,N_compile_flags,nReggieBuilds)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Build,   ONLY: ReadConfiguration
USE MOD_RegressionCheck_Vars,    ONLY: BuildDir
USE MOD_RegressionCheck_Vars,    ONLY: BuildCounter,BuildSolver
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
! if "BuildSolver" is true, the complete (valid) compiler-flag parameter combination
! is tested (specified in "configuration.reggie", default example is "run_basic")
IF(BuildSolver)THEN
  IF(ReggieBuildExe.EQ.'')THEN
    CALL ReadConfiguration(iExample,nReggieBuilds,N_compile_flags)
    BuildCounter=1 ! reset the counter between read and build (used for selecting the build configuration for compilation)
    SYSCOMMAND='rm -rf '//TRIM(BuildDir)//'build_reggie_bin > /dev/null 2>&1'
    CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
    SYSCOMMAND= 'mkdir '//TRIM(BuildDir)//'build_reggie_bin'
    CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
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
SUBROUTINE GetCodeBinary(iExample,iReggieBuild,nReggieBuilds,N_compile_flags,ReggieBuildExe,SkipBuild,ExitBuild)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: BuildDir,CodeNameLowCase,EXECPATH
USE MOD_RegressionCheck_Vars,    ONLY: BuildValid,BuildSolver
USE MOD_RegressionCheck_Tools,   ONLY: CheckForExecutable
USE MOD_RegressionCheck_Build,   ONLY: BuildConfiguration
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iReggieBuild,nReggieBuilds,N_compile_flags
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
! Get code binary (build or find it)
IF(BuildSolver)THEN
  ! if build is not valid no binary has been built and the lopp can cycle here
  IF(.NOT.BuildValid(iReggieBuild))THEN ! invalid reggie build 
    WRITE (ReggieBuildExe, '(a, i4.4)') "invalid"
  ELSE
    WRITE (ReggieBuildExe, '(a, i4.4)') CodeNameLowCase, COUNT(BuildValid(1:iReggieBuild)) ! e.g. XXXXX0001
  END IF
  ! check if build exists -> if it does, don't build a new executable with cmake
  FileName=TRIM(BuildDir)//'build_reggie_bin/'//ReggieBuildExe
  INQUIRE(File=FileName,EXIST=ExistFile)
  IF(ExistFile) THEN ! 1. build already exists (e.g. XXXX0001 located in ../build_reggie_bin/)
    EXECPATH=TRIM(FileName)
  ELSE ! 2. build does not exists -> create it
    CALL BuildConfiguration(iExample,iReggieBuild,nReggieBuilds,N_compile_flags)
    IF(BuildValid(iReggieBuild))THEN ! only move binary if it has been created (only for valid builds)
      SYSCOMMAND='cd '//TRIM(BuildDir)//' && mv build_reggie/bin/'//CodeNameLowCase//' build_reggie_bin/'//ReggieBuildExe
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
  CALL CheckForExecutable(Mode=2) ! check if executable was created correctly
END IF
END SUBROUTINE GetCodeBinary



!===================================================================================================================================
!> check if executable is compiled with correct TESTCASE (e.g. for tylorgreenvortex)
!===================================================================================================================================
SUBROUTINE CheckCompilerFlags(iReggieBuild,TESTCASE,TIMEDISCMETHOD)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: CodeNameLowCase,EXECPATH
USE MOD_RegressionCheck_Vars,    ONLY: BuildSolver
USE MOD_RegressionCheck_Build,   ONLY: BuildConfiguration
USE MOD_RegressionCheck_Vars,    ONLY: BuildTESTCASE,BuildTIMEDISCMETHOD,CodeNameUppCase
USE MOD_RegressionCheck_Build,   ONLY: GetFlagFromFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iReggieBuild
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: TESTCASE,TIMEDISCMETHOD
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: ExistFile                         !> file exists=.true., file does not exist=.false.
CHARACTER(LEN=255)             :: FileName                          !> path to a file or its name
!===================================================================================================================================
IF(BuildSolver)THEN
  TESTCASE=BuildTESTCASE(iReggieBuild)
  TIMEDISCMETHOD=BuildTIMEDISCMETHOD(iReggieBuild)
ELSE
  FileName=EXECPATH(1:INDEX(EXECPATH,'/',BACK = .TRUE.))//'configuration.cmake'
  INQUIRE(File=FileName,EXIST=ExistFile)
  IF(ExistFile) THEN
    CALL  GetFlagFromFile(FileName,CodeNameUppCase//'_TESTCASE',TESTCASE)
    IF(CodeNameLowCase.EQ.'boltzplatz')TESTCASE='default'! set default (currently no testcases are implemented)
    IF(TRIM(TESTCASE).EQ.'flag does not exist')CALL abort(&
      __STAMP__&
      ,CodeNameUppCase//'_TESTCASE flag not found in configuration.cmake!',999,999.)

    CALL GetFlagFromFile(FileName,CodeNameUppCase//'_TIMEDISCMETHOD',TIMEDISCMETHOD)
    IF(CodeNameLowCase.EQ.'flexi')TIMEDISCMETHOD='default'! set default (TIMEDISCMETHOD is not a compile flag)
    IF(TRIM(TIMEDISCMETHOD).EQ.'flag does not exist')CALL abort(&
      __STAMP__&
      ,CodeNameUppCase//'_TIMEDISCMETHOD flag not found in configuration.cmake!',999,999.)
  ELSE
    SWRITE(UNIT_stdOut,'(A12,A)')     ' ERROR: ','no "configuration.cmake" found at the location of the'//CodeNameLowCase//&
                                                                                                                      ' binary.'
    SWRITE(UNIT_stdOut,'(A12,A)')  ' FileName: ', TRIM(FileName)
    SWRITE(UNIT_stdOut,'(A12,L)') ' ExistFile: ', ExistFile
    ERROR STOP '-1'
  END IF
END IF
END SUBROUTINE CheckCompilerFlags


!===================================================================================================================================
!> Exclude subexamples set in parameter_reggie.ini when the binary is not suited for its execution
!===================================================================================================================================
SUBROUTINE CheckSubExample(iExample,TIMEDISCMETHOD)
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
CHARACTER(LEN=*),INTENT(IN)    :: TIMEDISCMETHOD
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! remove subexample (before printing the case overview) for certain configurations: e.g. Preconditioner when running explicitly
IF((TRIM(TIMEDISCMETHOD).NE.'ImplicitO3').AND.(TRIM(Examples(iExample)%SubExample).EQ.'PrecondType'))THEN
  Examples(iExample)%SubExample       = ''
  Examples(iExample)%SubExampleNumber = 0
  Examples(iExample)%SubExampleOption(1:20) = '-' ! default option is nothing
END IF
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
  FolderName=ExampleNames(iExample) ! e.g. run_TESTCASE-taylorgreenvortex/
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
!===================================================================================================================================
! get list of parameter files for running the simulation
parameter_ini2=''
IF(TRIM(TIMEDISCMETHOD).EQ.'DSMC')THEN
  ! set main parameter.ini file: e.g. parameter_XXX_EQNSYSNAME.ini
  parameter_ini='parameter_'//CodeNameLowCase//'_'//TRIM(ADJUSTL(Examples(iExample)%EQNSYSNAME))//'_DSMC.ini'
  parameter_ini2='parameter_DSMC.ini'
ELSE
  ! set main parameter.ini file: e.g. parameter_XXX_EQNSYSNAME.ini
  parameter_ini='parameter_'//CodeNameLowCase//'_'//TRIM(ADJUSTL(Examples(iExample)%EQNSYSNAME))//'.ini'
END IF
INQUIRE(File=TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini),EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File found under ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' parameter_ini:      ',TRIM(parameter_ini)
  ERROR STOP '-1'
END IF
IF(parameter_ini2.NE.'')THEN
  INQUIRE(File=TRIM(Examples(iExample)%PATH)//TRIM(parameter_ini2),EXIST=ExistFile)
  IF(.NOT.ExistFile) THEN
    SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File found under ',TRIM(Examples(iExample)%PATH)
    SWRITE(UNIT_stdOut,'(A,A)') ' parameter_ini2:     ',TRIM(parameter_ini2)
    ERROR STOP '-1'
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
  SWRITE(UNIT_stdOut,'(A)')    ' MPIthreads (number of threads): ['//TRIM(Examples(iExample)%MPIthreads(1          ))//&
                                                    ' to '//TRIM(Examples(iExample)%MPIthreads(Examples(iExample)%MPIthreadsN))//']'
ELSE
  SWRITE(UNIT_stdOut,'(A)')    ' MPIthreads (number of threads): ['//TRIM(Examples(iExample)%MPIthreads(1))//']'
END IF
SWRITE(UNIT_stdOut,'(A)')      ' Reference:                      ['//TRIM(Examples(iExample)%ReferenceFile)//']'
SWRITE(UNIT_stdOut,'(A)')      ' State:                          ['//TRIM(Examples(iExample)%ReferenceStateFile)//']'
SWRITE(UNIT_stdOut,'(A)')      ' HDF5 dataset:                   ['//TRIM(Examples(iExample)%ReferenceDataSetName)//']'
SWRITE(UNIT_stdOut,'(A)')      ' Restart:                        ['//TRIM(Examples(iExample)%RestartFileName)//']'
SWRITE(UNIT_stdOut,'(A)')      ' Example%SubExample:             ['//TRIM(Examples(iExample)%SubExample)//']'
SWRITE(UNIT_stdOut,'(A,I6,A1)')' Example%SubExampleNumber:       [',      Examples(iExample)%SubExampleNumber,']'
SWRITE(UNIT_stdOut,'(A)')      ' parameter files:                ['//TRIM(parameter_ini)//' '//TRIM(parameter_ini2)//']'
END SUBROUTINE PrintExampleInfo


!===================================================================================================================================
!> Set the SubExample in the parameter.ini file
!===================================================================================================================================
SUBROUTINE SetSubExample(iExample,iSubExample,parameter_ini)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iSubExample
CHARACTER(LEN=*),INTENT(IN)    :: parameter_ini
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSTATUS                           !> status
CHARACTER(LEN=500)             :: SYSCOMMAND                        !> string to fit the system command
!===================================================================================================================================
IF(Examples(iExample)%SubExampleNumber.GT.0)THEN ! SubExample has been specified
  SWRITE(UNIT_stdOut,'(A)')" "
  SWRITE(UNIT_stdOut,'(A,I2,A,A)')" SubExampleOption(",iSubExample,")=",TRIM(Examples(iExample)%SubExampleOption(iSubExample))
  SYSCOMMAND=    'cd '//TRIM(Examples(iExample)%PATH)//& ! print the current SubExampleOption(iSubExample) to parameter_ini
 ' && sed -i -e "s/.*'//TRIM(Examples(iExample)%SubExample)//&
                '=.*/'//TRIM(Examples(iExample)%SubExample)//&
                   '='//TRIM(Examples(iExample)%SubExampleOption(iSubExample))//&
                 '/" '//TRIM(parameter_ini)
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
END IF
END SUBROUTINE SetSubExample


!===================================================================================================================================
!> delete pre-existing data files before running the code
!===================================================================================================================================
SUBROUTINE CleanFolder(iExample)
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=500)             :: SYSCOMMAND                        !> string to fit the system command
INTEGER                        :: iSTATUS                           !> status
!===================================================================================================================================
! delete pre-existing data files before running the code 
! 1.) Files needed by "IntegrateLine" comparison
IF(Examples(iExample)%IntegrateLine)THEN
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm '//TRIM(Examples(iExample)%IntegrateLineFile)//' > /dev/null 2>&1'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! delete, e.g., "TGVAnalysis.dat" or "Database.csv"
END IF
END SUBROUTINE CleanFolder

!===================================================================================================================================
!> Execute the binary and check if the attempt was successful
!===================================================================================================================================
SUBROUTINE RunTheCode(iExample,iSubExample,iScaling,iRun,EXECPATH,parameter_ini,parameter_ini2,SkipComparison)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples
USE MOD_RegressionCheck_tools,   ONLY: AddError,str2int
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iSubExample,iScaling,iRun
CHARACTER(LEN=*),INTENT(IN)    :: parameter_ini,parameter_ini2,EXECPATH
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)            :: SkipComparison
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSTATUS                           !> status
CHARACTER(LEN=500)             :: SYSCOMMAND                        !> string to fit the system command
CHARACTER(LEN=15)              :: MPIthreadsStr                     !> string for the number of MPI threads for execution
CHARACTER(LEN=255)             :: FileSuffix,FolderSuffix 
INTEGER                        :: tempINT,PolynomialDegree,MPIthreads
!===================================================================================================================================
SkipComparison=.FALSE.
! -----------------------------------------------------------------------------------------------------------------------
! Run the Code
! -----------------------------------------------------------------------------------------------------------------------
IF(Examples(iExample)%MPIrun)THEN ! use "mpirun"
  !IF(Examples(iExample)%MPIthreads.GT.1)THEN
  MPIthreadsStr=ADJUSTL(TRIM(Examples(iExample)%MPIthreads(iScaling)))
  CALL str2int(ADJUSTL(TRIM(MPIthreadsStr)),tempINT,iSTATUS) ! sanity check if the number of threads is correct
  IF((tempINT.LE.0).OR.(iSTATUS.NE.0))CALL abort(&
      __STAMP__&
      ,'RunTheCode(): Number of MPI threads is corrupt = '//ADJUSTL(TRIM(MPIthreadsStr)))
  IF((iScaling.GT.1).AND.(iRun.EQ.1))THEN
    SWRITE(*,*)"Examples(iExample)%MPIthreads=",Examples(iExample)%MPIthreads(iScaling)
  END IF

  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && '//TRIM(Examples(iExample)%MPIcommand)//' -np '//&
                       ADJUSTL(TRIM(MPIthreadsStr))//' '//TRIM(EXECPATH)//' '//TRIM(parameter_ini)//' ' &
              //TRIM(parameter_ini2)//' '//TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
ELSE
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && '//TRIM(EXECPATH)//' '//TRIM(parameter_ini)//' ' &
              //TRIM(parameter_ini2)//' '//TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
END IF
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! run the code
! -----------------------------------------------------------------------------------------------------------------------
! was the run successful? (iSTATUS=0)
! -----------------------------------------------------------------------------------------------------------------------
IF(iSTATUS.EQ.0)THEN ! Computation successful
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='no')  ' successful computation ...'
  CALL AddError('successful computation',iExample,iSubExample,ErrorStatus=0,ErrorCode=0)

! copy the std.out file
IF(Examples(iExample)%SubExample.EQ.'N')THEN!when polynomial degree "N" is the SubExample
  CALL str2int(Examples(iExample)%SubExampleOption(iSubExample),PolynomialDegree,iSTATUS)
  CALL str2int(Examples(iExample)%MPIthreads(iScaling),MPIthreads,iSTATUS)
  WRITE(FileSuffix, '(A10,I8.8,A2,I4.4,A5,I4.4)') 'MPIthreads',MPIthreads,'_N',PolynomialDegree,'_iRun',iRun
  WRITE(FolderSuffix,'(A1,I4.4)')'N',PolynomialDegree
ELSE
  WRITE(FileSuffix, '(I4.4,I4.4,I4.4)') iScaling,iSubExample,iRun
  WRITE(FolderSuffix,'(I4.4)')iSubExample
END IF
!print*,"FileSuffix  =",FileSuffix
!print*,"FolderSuffix=",FolderSuffix
!read*
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && cp std.out std-'//TRIM(FileSuffix)//'.out'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! copy the std.out file
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && mkdir std_files'//TRIM(FolderSuffix)//' > /dev/null 2>&1'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && mv std-'//TRIM(FileSuffix)//'.out std_files'//TRIM(FolderSuffix)//'/'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)


ELSE ! Computation failed
  SWRITE(UNIT_stdOut,'(A)')   ' Computation of example failed'
  SWRITE(UNIT_stdOut,'(A,A)') ' Out-file: ', TRIM(Examples(iExample)%PATH)//'std.out'
  SWRITE(UNIT_stdOut,'(A,A)') ' Errorfile: ', TRIM(Examples(iExample)%PATH)//'err.out'
  CALL AddError('Error while executing',iExample,iSubExample,ErrorStatus=1,ErrorCode=2)
  SkipComparison=.TRUE. ! when a computation fails no useful output is created for comparison, continue the subexample cycle
END IF
END SUBROUTINE RunTheCode










END MODULE MOD_RegressionCheck_Run
