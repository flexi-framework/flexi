#include "flexi.h"

!==================================================================================================================================
!> Contains the utilize routines of the regressioncheck
!> -GetExampleList extracts the examples which are subfolders in examples 
!> -CleanExample removes the output in a example after a successful run
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
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE GetCommandLineOption
  MODULE PROCEDURE GetCommandLineOption
END INTERFACE

INTERFACE GetExampleList
  MODULE PROCEDURE GetExampleList
END INTERFACE

INTERFACE CleanExample
  MODULE PROCEDURE CleanExample
END INTERFACE

INTERFACE InitExample
  MODULE PROCEDURE InitExample
END INTERFACE

INTERFACE CheckForExecutable
  MODULE PROCEDURE CheckForExecutable
END INTERFACE

INTERFACE str2real
  MODULE PROCEDURE str2real
END INTERFACE

INTERFACE str2int
  MODULE PROCEDURE str2int
END INTERFACE

INTERFACE str2logical
  MODULE PROCEDURE str2logical
END INTERFACE

PUBLIC::GetExampleList,InitExample,CleanExample, CheckForExecutable,GetCommandLineOption
PUBLIC::str2real,str2int,str2logical
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> reads the command line options for the regressioncheck
!> options are:
!> run [default]:   - runs only the regressioncheck
!> build            - builds all valid combinations of flexi (default uses the configuration.flexi from run_freestream) and 
!>                    performs the tests
!>
!> ./regressioncheck [RuntimeOption] [RuntimeOptionType]
!> 
!> ./regressioncheck                -> uses default "run" and runs the current compiler build and all "run_" examples
!> ./regressioncheck
!> ./regressioncheck build          -> runs "run_freestream" for numerous builds
!> ./regressioncheck build convtest -> runs "feature_convtest" for numerous builds defined in "feature_convtest/configuration.flexi"
!> ./regressioncheck build all      -> runs all examples for numerous builds defined in "run_freestream/configuration.flexi"
!==================================================================================================================================
SUBROUTINE GetCommandLineOption()
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars, ONLY: RuntimeOption,RuntimeOptionType,BuildDebug,RuntimeOptionTypeII
USE MOD_RegressionCheck_Vars, ONLY: RuntimeOptionTypeIII,BuildContinue,BuildSolver
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nArgs                             ! Number of supplied command line arguments
!===================================================================================================================================
RuntimeOption='run'     ! default
RuntimeOptionType='run' ! default
RuntimeOptionTypeII=''  ! default
RuntimeOptionTypeIII='' ! default
! Get number of command line arguments and read in runtime option of regressioncheck
nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs.EQ.0)THEN
  BuildSolver=.FALSE.
ELSE
  CALL GET_COMMAND_ARGUMENT(1,RuntimeOption)
  IF(nArgs.GE.2)CALL GET_COMMAND_ARGUMENT(2,RuntimeOptionType)
  IF(nArgs.GE.3)CALL GET_COMMAND_ARGUMENT(3,RuntimeOptionTypeII)
  IF(nArgs.GE.4)CALL GET_COMMAND_ARGUMENT(4,RuntimeOptionTypeIII)

  ! get number of threads/procs for parallel building
  CALL GetNumberOfProcs(nArgs)
  
  IF(TRIM(RuntimeOption).EQ.'run') THEN
    BuildSolver=.FALSE.
  ELSE IF(TRIM(RuntimeOption(1:5)).EQ.'build') THEN
    IF(TRIM(RuntimeOption).EQ.'build-continue')BuildContinue=.TRUE.
    BuildSolver=.TRUE.
    IF(TRIM(RuntimeOptionType).EQ.'debug')THEN
      BuildDebug=.TRUE.
      RuntimeOptionType='run_freestream' ! debug uses "configuration.flexi" from "run_freestream" and displays the complete 
                                         ! compilation process for debugging
    END IF
    IF(TRIM(RuntimeOptionTypeII).EQ.'debug')BuildDebug=.TRUE. ! e.g. ./regressioncheck build feature_convtest debug
    IF(TRIM(RuntimeOptionType).EQ.'run')RuntimeOptionType='run_freestream'
  ELSE IF((TRIM(RuntimeOption).EQ.'--help').OR.(TRIM(RuntimeOption).EQ.'help').OR.(TRIM(RuntimeOption).EQ.'HELP')) THEN
    CALL Print_Help_Information()
    STOP
  ELSE
    SWRITE(UNIT_stdOut,'(A)') ' ERROR: wrong argument for regressioncheck!'
    ERROR STOP '-2'
  END IF
  
  ! [RuntimeOptionType] = all: run all example folders
  IF((TRIM(RuntimeOptionType).EQ.'all').OR.(TRIM(RuntimeOptionType).EQ.'ALL'))RuntimeOptionType=''
END IF

SWRITE(UNIT_stdOut,'(A,A1,A,A1,A1,A,A1,A1,A,A1,A1,A,A1)')' Running with arguments: ',&
'[',TRIM(RuntimeOption)       ,']',&
'[',TRIM(RuntimeOptionType)   ,']',&
'[',TRIM(RuntimeOptionTypeII) ,']',&
'[',TRIM(RuntimeOptionTypeIII),']'
END SUBROUTINE GetCommandLineOption


!==================================================================================================================================
!> Check if examples exist. Next, scan the folder for all examples. The routine returns the number of examples, their name
!> and nullifies the parameter entries for each example
!> If the regression check is run within an example folder, only said example is executed
!==================================================================================================================================
SUBROUTINE GetExampleList()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: nExamples,ExampleNames,Examples,ExamplesDir,BuildDir,RuntimeOptionType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLESzz
CHARACTER(LEN=500)            :: SYSCOMMAND           ! string to fit the system command
CHARACTER(LEN=500)            :: SYSCOMMANDTWO        ! string to fit the system command
CHARACTER(LEN=255)            :: FilePathName,FileName! file and path strings
INTEGER                       :: ioUnit=33            ! io-unit
INTEGER                       :: iSTATUS              ! status
INTEGER                       :: iExample             ! loop index for example
CHARACTER(len=255)            :: cwd                  ! current cworking directory CALL getcwd(cwd)
LOGICAL                       :: ExistFile            ! T=file exists, F=file does not exist
!==================================================================================================================================
! check if regressioncheck is executed von /bin/regressioncheck directory, otherwise use the current working directory
CALL getcwd(cwd)                                      ! get path of current directory
FileName=cwd(INDEX(cwd,'/',BACK = .TRUE.)+1:LEN(cwd)) ! cut the last folder name from current diectory path
FilePathName=TRIM(cwd)//'/parameter_reggie.ini'       ! check if parameter_reggie.ini is located within cwd
INQUIRE(File=FilePathName,EXIST=ExistFile)            ! inquire

IF(ExistFile.EQV..FALSE.)THEN ! use existing example folder
  ExamplesDir=TRIM(BASEDIR(2:LEN(BASEDIR)-1))//'../regressioncheck/examples/'
  SYSCOMMAND='cd '//TRIM(ExamplesDir)//' && ls -d */ > tmp.txt'
ELSE ! run regressioncheck for a single folder located anywhere from which the reggie is executed
  ExamplesDir='./../'
  SYSCOMMAND='cd '//TRIM(ExamplesDir)//' && ls -d '//TRIM(FileName)//'/ > tmp.txt'
  RuntimeOptionType=TRIM(FileName) ! override RuntimeOptionType in order to select only this directory
END IF
BuildDir=TRIM(BASEDIR(2:LEN(BASEDIR)-1))! use basedir because one cannot use: TRIM(cwd)//'/' because the checked out flexi source 
                                        ! code is needed for building new flexi binaries

! sanity check: change directory into the example folder, if it does not exist this check fails
SYSCOMMANDTWO='cd '//TRIM(ExamplesDir)
CALL EXECUTE_COMMAND_LINE(SYSCOMMANDTWO, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Error: Example folder does not work!'
  ERROR STOP '66'
END IF

! get number of examples by complicated fortran hack:
! ls is called and the output is piped into the file tmp.txt. The number of lines is the number of available examples.
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Could not create tmp.txt to get number of examples'
  ERROR STOP '99'
END IF

! read tmp.txt | list of directories if regressioncheck/examples
FileName=TRIM(ExamplesDir)//'tmp.txt'
OPEN(UNIT = ioUnit, FILE = FileName, STATUS ="OLD", IOSTAT = iSTATUS ) 

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
  Examples(iExample)%CheckedStateFile=''
  Examples(iExample)%ReferenceStateFile=''
  Examples(iExample)%ReferenceDataSetName=''
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
  ERROR STOP '99'
END IF

END SUBROUTINE GetExampleList


!==================================================================================================================================
!> Read the parameter_reggie.ini file of an example given by its relative path. It
!> contains information for the computation of the example:
!>  exec - a mpi or serial example
!>  optional reference files for error-norms, reference state file and tested dataset and name of the checked state file
!>  optional a restart filename
!==================================================================================================================================
SUBROUTINE InitExample(FilePath,FilePathLength,Example)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: tExample
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: FilePathLength
CHARACTER(LEN=FilePathLength),INTENT(IN)  :: FilePath
TYPE(tExample),INTENT(INOUT)              :: Example
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ioUnit
INTEGER                                   :: iSTATUS,IndNum,IndNum2,IndNum3,MaxNum
INTEGER                                   :: iSubExample
CHARACTER(LEN=255)                        :: FileName,temp1,temp2
LOGICAL                                   :: ExistFile
!==================================================================================================================================
! test if file exists and open
ioUnit=GETFREEUNIT()
FileName=TRIM(FilePath)//'parameter_reggie.ini'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A12,A)')     ' ERROR: ','no parameter_reggie.ini found.'
  SWRITE(UNIT_stdOut,'(A12,A)')  ' FileName: ', TRIM(FileName)
  SWRITE(UNIT_stdOut,'(A12,L)') ' ExistFile: ', ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF

! init
Example%EXEC=.FALSE.
Example%nVar=0
Example%ReferenceTolerance     = -1.
Example%SubExampleNumber       = 0   ! init total number of subexamples
Example%SubExampleOption(1:20) = '-' ! default option is nothing
DO ! extract reggie information
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1 ! get first line assuming it is something like "nVar= 5"
  IF(iSTATUS.EQ.-1) EXIT ! end of file (EOF) reached
  IF(INDEX(temp1,'!').GT.0)temp1=temp1(1:INDEX(temp1,'!')-1) ! if temp1 contains a "!", remove it and the following characters
  temp1=TRIM(ADJUSTL(temp1)) ! remove trailing and leading white spaces
  temp2=''
  IndNum=INDEX(temp1,'=')
  MaxNum=LEN(TRIM(temp1))
  IF(IndNum.GT.0)THEN ! definition found
    ! get size of EQNSYS (deprecated)
    IF(TRIM(temp1(1:IndNum-1)).EQ.'nVar')CALL str2int(temp1(IndNum+1:MaxNum),Example%nVar,iSTATUS)
    ! single or parallel (deprecated)
    IF(TRIM(temp1(1:IndNum-1)).EQ.'MPI')CALL str2logical(temp1(IndNum+1:MaxNum),Example%EXEC,iSTATUS)
    ! Reference Norm/State
    IF(TRIM(temp1(1:IndNum-1)).EQ.'ReferenceTolerance')CALL str2real(temp1(IndNum+1:MaxNum),Example%ReferenceTolerance,iSTATUS)
    IF(TRIM(temp1(1:IndNum-1)).EQ.'ReferenceFile')          Example%ReferenceFile         =TRIM(ADJUSTL(temp1(IndNum+1:MaxNum)))
    IF(TRIM(temp1(1:IndNum-1)).EQ.'ReferenceStateFile')     Example%ReferenceStateFile    =TRIM(ADJUSTL(temp1(IndNum+1:MaxNum)))
    IF(TRIM(temp1(1:IndNum-1)).EQ.'CheckedStateFile')       Example%CheckedStateFile      =TRIM(ADJUSTL(temp1(IndNum+1:MaxNum)))
    IF(TRIM(temp1(1:IndNum-1)).EQ.'ReferenceDataSetName')   Example%ReferenceDataSetName  =TRIM(ADJUSTL(temp1(IndNum+1:MaxNum)))
    IF(TRIM(temp1(1:IndNum-1)).EQ.'RestartFileName')        Example%RestartFileName       =TRIM(ADJUSTL(temp1(IndNum+1:MaxNum)))
    ! SubExamples - currently one subexample class is allowed with multiple options
    IF(TRIM(temp1(1:IndNum-1)).EQ.'SubExample')THEN
      temp1=temp1(IndNum+1:MaxNum)
      MaxNum=LEN(TRIM(temp1))
      IndNum2=INDEX(temp1(1:MaxNum),'{')
      IF(IndNum2.GT.0)THEN
        Example%SubExample     = TRIM(ADJUSTL(temp1(1:IndNum2-1)))
        temp1=TRIM(ADJUSTL(temp1(IndNum2+1:MaxNum)))
        IndNum2=INDEX(temp1(1:MaxNum),'}')
        IF(IndNum2.GT.0)THEN
          temp1=temp1(1:IndNum2-1)
          !IndNum2=1
          DO iSubExample=1,19
            IndNum2=INDEX(temp1(1:LEN(TRIM(temp1))),',')
            Example%SubExampleNumber = Example%SubExampleNumber + 1 ! get total number of subexamples
            IF(IndNum2.GT.0)THEN
              Example%SubExampleOption(iSubExample)=TRIM(ADJUSTL(temp1(1:IndNum2-1)))
              temp1=temp1(IndNum2+1:LEN(TRIM(temp1)))
            ELSE
              Example%SubExampleOption(iSubExample)=TRIM(ADJUSTL(temp1(1:LEN(TRIM(temp1)))))
              EXIT
            END IF
          END DO 
        ELSE
          Example%SubExample       = ''
          Example%SubExampleNumber = 0
        END IF      
      ELSE
        Example%SubExample       = ''
        Example%SubExampleNumber = 0
      END IF      
    END IF
    ! Line integration
    IF(TRIM(temp1(1:IndNum-1)).EQ.'IntegrateLine')THEN
       Example%IntegrateLine=.TRUE.
       Example%IntegrateLineRange=0 ! init
       IndNum2=INDEX(temp1(IndNum+1:MaxNum),',')
       IF(IndNum2.GT.0)THEN
         temp2=temp1(IndNum+1:MaxNum)
         Example%IntegrateLineFile=TRIM(ADJUSTL(temp2(1:IndNum2-1)))
         temp2=temp2(IndNum2+1:LEN(TRIM(temp2)))
         IndNum2=INDEX(temp2,',')
         IF(IndNum2.GT.0)THEN ! get column ranges
           IndNum3=INDEX(temp2(1:IndNum2),':')
           IF(IndNum3.GT.0)THEN ! check range
             CALL str2int(temp2(1        :IndNum3-1),Example%IntegrateLineRange(1),iSTATUS)
             CALL str2int(temp2(IndNum3+1:IndNum2-1),Example%IntegrateLineRange(2),iSTATUS)
             temp2=temp2(IndNum2+1:LEN(TRIM(temp2)))
             IndNum2=LEN(temp2)
             IF(IndNum2.GT.0)THEN ! get integral value
               CALL str2real(temp2,Example%IntegrateLineValue,iSTATUS)
             END IF ! get integral value 
           END IF ! check range
         END IF ! get column ranges
       END IF ! get file name
      IF(ANY(Example%IntegrateLineRange(1:2).EQ.0))Example%IntegrateLine=.FALSE. 
      IF(Example%IntegrateLineFile.EQ.'')Example%IntegrateLine=.FALSE.
    END IF ! 'IntegrateLine'
    ! Next feature
    !IF(TRIM(temp1(1:IndNum-1)).EQ.'NextFeature')
  END IF ! IndNum.GT.0 -> definition found
END DO
CLOSE(ioUnit)
END SUBROUTINE InitExample


!==================================================================================================================================
!> This routine goes into each example folder of the regressioncheck. Inside, the not required State-files and std.out and err.out
!> files are removed. The subroutine is called, if the example is computed successfully.
!==================================================================================================================================
SUBROUTINE CleanExample(iExample)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=500)             :: SYSCOMMAND
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)             :: tmp
INTEGER                        :: iSTATUS,ioUnit=22
!==================================================================================================================================
! delete all *.out files
SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm *.out > /dev/null 2>&1'
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0)THEN
  SWRITE(UNIT_stdOut,'(A)')' CleanExample(',Examples(iExample)%PATH,'): Could not remove *.out files!'
END IF

! delete all *State* files except *reference* state files
IF((Examples(iExample)%ReferenceStateFile.EQ.'').AND. &
   (Examples(iExample)%RestartFileName.EQ.'') ) THEN
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm *State* > /dev/null 2>&1'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  IF(iSTATUS.NE.0)THEN
    SWRITE(UNIT_stdOut,'(A)')' CleanExample(',Examples(iExample)%PATH,'): Could not remove *State* files!'
  END IF
ELSE
  ! create list of all *State* files and loop them: don't delete *reference* files
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && ls *State* > tmp.txt'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  IF(iSTATUS.NE.0)THEN
    SWRITE(UNIT_stdOut,'(A)')' CleanExample(',Examples(iExample)%PATH,'): Could not remove tmp.txt!'
  END IF
  ! read tmp.txt | list of directories if regressioncheck/examples
  FileName=TRIM(Examples(iExample)%PATH)//'tmp.txt'
  OPEN(UNIT = ioUnit, FILE = FileName, STATUS ="OLD", IOSTAT = iSTATUS ) 
  DO 
    READ(ioUnit,FMT='(A)',IOSTAT=iSTATUS) tmp
    IF (iSTATUS.NE.0) EXIT
    IF((Examples(iExample)%ReferenceStateFile.NE.TRIM(tmp)).AND. &
       (Examples(iExample)%RestartFileName.NE.TRIM(tmp)) ) THEN
       SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm '//TRIM(tmp)
       CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
       IF(iSTATUS.NE.0) THEN
         SWRITE(UNIT_stdOut,'(A)')  ' CleanExample(',Examples(iExample)%PATH,'): Could not remove state file ',TRIM(tmp)
       END IF
    END IF
  END DO
  CLOSE(ioUnit)
  ! clean tmp.txt
  SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm tmp.txt'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  IF(iSTATUS.NE.0) THEN
    SWRITE(UNIT_stdOut,'(A)')  ' CleanExample(',Examples(iExample)%PATH,'): Could not remove tmp.txt'
  END IF
END IF

END SUBROUTINE CleanExample


!==================================================================================================================================
!>  Check if executable exists
!==================================================================================================================================
SUBROUTINE CheckForExecutable(Mode)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: EXECPATH
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: Mode             ! which mode (1 or 2)
                                                  ! 1: pre-compiled flexi executable
                                                  ! 2: flexi compiled via reggie with defined "EXECPATH", e.g., 
                                                  ! "~/Flexi/flexi/build_reggie.dev/build_reggie/bin/flexi"
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: ExistSolver      ! logical to flag solver
CHARACTER(len=255)            :: cwd             ! current cworking directory CALL getcwd(cwd)
!===================================================================================================================================
IF(Mode.EQ.1) THEN
  ! 1.) check if flexi binary is located at the position of the regressioncheck binary (e.g. in a specific case folder)
  CALL getcwd(cwd) ! get path of current directory
  EXECPATH=TRIM(cwd)//'/flexi'
  INQUIRE(File=EXECPATH,EXIST=ExistSolver)
  IF(.NOT.ExistSolver) THEN
    ! 2.) flexi binary not found in location from where the reggie was executed -> check standard directory instead
    EXECPATH=BASEDIR(2:LEN(BASEDIR)-1)//'bin/flexi'
  END IF
END IF
! EXECPATH has been set, inquire the existence of the flexi binary
INQUIRE(File=EXECPATH,EXIST=ExistSolver)
IF(.NOT.ExistSolver) THEN
  SWRITE(UNIT_stdOut,'(A38,I1,A)') ' CALL CheckForExecutable() with Mode=',&
                                       Mode,': no executable found. Error during compilation or linking?'
  SWRITE(UNIT_stdOut,'(A41,A)')                         ' EXECPATH: ', TRIM(EXECPATH)
  SWRITE(UNIT_stdOut,'(A41,L)')                      ' ExistSolver: ', ExistSolver
  ERROR STOP '77'
END IF
END SUBROUTINE CheckForExecutable


!==================================================================================================================================
!> Get the number of threads/procs for a parallel compilation
!==================================================================================================================================
SUBROUTINE GetNumberOfProcs(nArgs)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars, ONLY: RuntimeOptionType,RuntimeOptionTypeII,RuntimeOptionTypeIII
USE MOD_RegressionCheck_Vars, ONLY: NumberOfProcs,NumberOfProcsStr
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nArgs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSTATUS                           ! Error status
!===================================================================================================================================
IF((nArgs.GE.2))THEN
  CALL str2int(RuntimeOptionType,NumberOfProcs,iSTATUS)
  IF(iSTATUS.EQ.0)THEN
    RuntimeOptionType='run' ! return to default -> needed for setting it to 'run_freestream'
  ELSE
    IF(nArgs.GE.3)THEN
      CALL str2int(RuntimeOptionTypeII,NumberOfProcs,iSTATUS)
      IF(iSTATUS.EQ.0)THEN
        RuntimeOptionTypeII='' ! return to default
      ELSE
        IF(nArgs.GE.4)THEN
          CALL str2int(RuntimeOptionTypeIII,NumberOfProcs,iSTATUS)
          IF(iSTATUS.EQ.0)THEN
            RuntimeOptionTypeIII=''
          ELSE
            NumberOfProcs=1
          END IF
        ELSE
          NumberOfProcs=1
        END IF
      END IF
    ELSE
      NumberOfProcs=1
    END IF
  END IF
  IF(iSTATUS.EQ.0)THEN
    SWRITE(UNIT_stdOut,'(A,I3,A)') ' Building regression checks with',NumberOfProcs,' threads/processors'
    IF(NumberOfProcs.GT.1)THEN
      WRITE(UNIT=NumberOfProcsStr,FMT='(I5)') NumberOfProcs
    ELSE
      NumberOfProcsStr='fail'
    END IF
  END IF
ELSE
  NumberOfProcs=1
END IF
END SUBROUTINE GetNumberOfProcs


!==================================================================================================================================
!> Convert a String to an Integer
!==================================================================================================================================
SUBROUTINE str2int(str,int_number,stat)
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
CHARACTER(len=*),INTENT(IN) :: str
INTEGER,INTENT(OUT)         :: int_number
INTEGER,INTENT(OUT)         :: stat
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  int_number
END SUBROUTINE str2int


!==================================================================================================================================
!> Convert a String to a REAL
!==================================================================================================================================
SUBROUTINE str2real(str,real_number,stat)
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
CHARACTER(len=*),INTENT(IN) :: str
REAL,INTENT(OUT)            :: real_number
INTEGER,INTENT(OUT)         :: stat
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  real_number
END SUBROUTINE str2real


!==================================================================================================================================
!> Convert a String to a LOGICAL
!==================================================================================================================================
SUBROUTINE str2logical(str,logical_number,stat)
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
CHARACTER(len=*),INTENT(IN) :: str
LOGICAL,INTENT(OUT)         :: logical_number
INTEGER,INTENT(OUT)         :: stat
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  logical_number
END SUBROUTINE str2logical


!==================================================================================================================================
!> print regression check information
!> 1.) give regressioncheck option parameters:
!>     ./regressioncheck [RuntimeOption] [RuntimeOptionType]
!>     
!>     ./regressioncheck                -> uses default "run" and runs the current compiler build and all "run_" examples
!>     ./regressioncheck
!>     ./regressioncheck build          -> runs "run_freestream" for numerous builds
!>     ./regressioncheck build convtest -> runs "feature_convtest" for numerous builds defined in "feature_convtest/configuration.flexi"
!>     ./regressioncheck build all      -> runs all examples for numerous builds defined in "run_freestream/configuration.flexi"
!> 2.) information on input files, e.g., comfiguration.flexi and parameter_reggie.ini
!> 3.) give information on error codes for builing/compiling the source code and running the code
!==================================================================================================================================
SUBROUTINE Print_Help_Information()
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
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
SWRITE(UNIT_stdOut,'(A)') '                            | "run_freestream". These tests all inlcude very                           '
SWRITE(UNIT_stdOut,'(A)') '                            | short simulations with <1sec execution time                              '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for on-check-in tests                                          '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' build                      | runs the example "run_freestream" on default and                         '
SWRITE(UNIT_stdOut,'(A)') '                            | (requires locally built HDF5 or loaded HDF5 paths)                       '
SWRITE(UNIT_stdOut,'(A)') '                            | compiles all possible compiler flag combinations                         '
SWRITE(UNIT_stdOut,'(A)') '                            | specified in "comfiguration.flexi" and considers                         '
SWRITE(UNIT_stdOut,'(A)') '                            | the specified exclude list for invalid combinations                      '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for nightly tests                                              '
SWRITE(UNIT_stdOut,'(A)') '                            | The cmake output may be shown on-screen by adding "debug" after possible '
SWRITE(UNIT_stdOut,'(A)') '                            | [RuntimeOptionType] commands                                             '
SWRITE(UNIT_stdOut,'(A)') '                            | Multi-processor compilation is supported by adding "XX"  after possible  '
SWRITE(UNIT_stdOut,'(A)') '                            | [RuntimeOptionType] commands, where "XX" is the number of processors     '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' build-continue             | does the same as "build", but starts the process at the last failed      '
SWRITE(UNIT_stdOut,'(A)') '                            | building position of a previous run, hence, all previously successful    '
SWRITE(UNIT_stdOut,'(A)') '                            | builds are skipped (-> debugging purposes)                               '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' conv_test (ToDo!)          | specific feature test: runs the "conv_test" example                      '
SWRITE(UNIT_stdOut,'(A)') '                            | runs two modes: p-convergence and h-convergence                          '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for weakly tests                                               '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' performance (ToDo!)        | specific feature test: runs the "performance" example                    '
SWRITE(UNIT_stdOut,'(A)') '                            | automatically checks out specified flexi version tag                     '
SWRITE(UNIT_stdOut,'(A)') '                            | and run the example to acquire the reference                             '
SWRITE(UNIT_stdOut,'(A)') '                            | performance wall time                                                    '
SWRITE(UNIT_stdOut,'(A)') '                            | e.g. used for weakly tests                                               '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' Second input argument [RuntimeOptionType] depends on the first input argument                         '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' [RuntimeOptionType]        | [RuntimeOption] | mode                                                   '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' run (default)              |                 |                                                        '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' cavity                     | run             | specific feature test: runs all the the examples that  '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | begin with "cavity", e.g., "cavity_1" or "cavity_xyz"  '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | the example "cavity", e.g., tests long time stability, '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | some BC and ExactFunc ( e.g. used for weakly tests)    '
SWRITE(UNIT_stdOut,'(A)') '                            |                 |                                                        '
SWRITE(UNIT_stdOut,'(A)') '                            | build           | uses the "configurations.flexi" within the "cavity"    '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | directory and builds all compiler flag combinations    '
SWRITE(UNIT_stdOut,'(A)') '                            |                 | that are specified there and runs them all on "cavity" '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------------------------'


! needed/optional files used by the regression check
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' 2.) Which input files are needed, e.g., comfiguration.flexi and parameter_reggie.ini'
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' The following parameter files are supported (within each example folder):           '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' parameter file             | Description                                            '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
SWRITE(UNIT_stdOut,'(A)') ' configuration.flexi        | needed for "./regressioncheck build"                   '
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
SWRITE(UNIT_stdOut,'(A)') '     L2/Linf reference file | ReferenceFile= referencenorm.txt                       '
SWRITE(UNIT_stdOut,'(A)') '  ref state file for h5diff | ReferenceStateFile= cavity_reference_State_0.200.h5    '
SWRITE(UNIT_stdOut,'(A)') '      state file for h5diff | CheckedStateFile= cavity_State_0000000.200000000.h5    '
SWRITE(UNIT_stdOut,'(A)') '      array name for h5diff | ReferenceDataSetName= DG_Solution                      '
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
SWRITE(UNIT_stdOut,'(A)') '                         77 | no flexi executable found for option run               '
SWRITE(UNIT_stdOut,'(A)') '                         99 | fail of execute_system_command                         '
SWRITE(UNIT_stdOut,'(A)') ' ------------------------------------------------------------------------------------'
END SUBROUTINE Print_Help_Information


END MODULE MOD_RegressionCheck_Tools
