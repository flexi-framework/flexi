#include "flexi.h"

!==================================================================================================================================
!> Contains the routines to build the binaries
!==================================================================================================================================
MODULE MOD_RegressionCheck_Build
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ReadConfiguration
  MODULE PROCEDURE ReadConfiguration
END INTERFACE

INTERFACE BuildConfiguration
  MODULE PROCEDURE BuildConfiguration
END INTERFACE

INTERFACE GetFlagFromFile
  MODULE PROCEDURE GetFlagFromFile
END INTERFACE

PUBLIC::ReadConfiguration
PUBLIC::BuildConfiguration
PUBLIC::GetFlagFromFile
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> read the file "configuirations.reggie" and creates multiple compiler flag configurations for cmake that are written to
!> "configurationsX.cmake"
!==================================================================================================================================
SUBROUTINE ReadConfiguration(iExample,nReggieBuilds,N_compile_flags)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples,RuntimeOptionType,BuildEQNSYS,BuildTESTCASE,BuildContinue,BuildContinueNumber
USE MOD_RegressionCheck_Vars,    ONLY: BuildTIMEDISCMETHOD
USE MOD_RegressionCheck_Vars,    ONLY: BuildConfigurations,BuildValid,BuildCounter,BuildIndex
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                           :: iExample
INTEGER,INTENT(INOUT)                        :: N_compile_flags
INTEGER,INTENT(INOUT)                        :: nReggieBuilds
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ioUnit,iSTATUS,I,J,K
INTEGER                                   :: io_error,CurrentIndex,NextIndex,IndNum
INTEGER                                   :: N_subinclude,N_exclude
CHARACTER(LEN=255)                        :: FileName,temp,temp2,COMPILE_FLAG,dummystr
CHARACTER(LEN=255)                        :: EXCLUDE_FLAG_A,EXCLUDE_FLAG_B,EXCLUDE_FLAG_C
LOGICAL                                   :: ExistFile,InvalidA,InvalidB,InvalidC
CHARACTER(LEN=255),ALLOCATABLE            :: ExcludeConfigurations(:,:),BuildValidInfo(:)
INTEGER                                   :: MaxBuildConfigurations=400,N_subinclude_max,N_compile_flags_max

CHARACTER(LEN=255)                        :: FilePath
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
"  Regression Check: Read Cmake Configurations"
SWRITE(UNIT_stdOut,'(132("="))')
ioUnit=GETFREEUNIT()
IF(BuildContinue)CALL GetBuildContinue()
IF(RuntimeOptionType.EQ.'')THEN ! [RuntimeOptionType] has been cleared (set to '') as the input by the user was "all", i.e., use all
                                ! examples use fixed configuration file (maximum number of builds?) but (maximum number of builds?)
  FilePath='./../../regressioncheck/examples/run_basic/'
ELSE
  FilePath=TRIM(Examples(iExample)%PATH)
END IF
FileName=TRIM(FilePath)//'configurations.reggie'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File under: ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:             ','configurations.reggie'
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:            ',ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF

DO I=1,2
  N_compile_flags=0
  N_exclude=0
  nReggieBuilds=1
  N_subinclude_max=1
  N_compile_flags_max=0
  DO
    READ(ioUnit,'(A)',iostat=IO_ERROR)temp
    IF(IO_ERROR.EQ.-1)EXIT
    IF(LEN(trim(temp)).GT.1)THEN
      N_subinclude=0 ! reset
      READ(temp,*)temp2
      IF(TRIM(temp2(1:1)).EQ.'!')CYCLE   !print*,"found !!!"
      IF(INDEX(temp,'!').GT.0)temp=temp(1:INDEX(temp,'!')-1) ! remove '!'
      IF(TRIM(temp(1:7)).EQ.'EXCLUDE')THEN
        IF(INDEX(temp,':').GT.0)THEN
          IF(INDEX(temp,',').GT.0)THEN
            EXCLUDE_FLAG_A=TRIM(ADJUSTL(temp(9                :INDEX(temp,',')-1)))
            temp = temp(INDEX(temp,',')+1:LEN(temp))
            !IndNum = INDEX(temp(INDEX(temp,',')+1:LEN(temp)),',')
            IndNum = INDEX(temp,',')
            IF(IndNum.GT.0)THEN
              EXCLUDE_FLAG_B=TRIM(ADJUSTL(temp(1:IndNum-1        )))
              EXCLUDE_FLAG_C=TRIM(ADJUSTL(temp(IndNum+1:LEN(temp))))
              IndNum = INDEX(temp(IndNum+1:LEN(temp)),',')
              IF (IndNum.GT.0) THEN
                SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: Too many EXCLUDE flags (>3): ',TRIM(temp)
                ERROR STOP '-1'
              ENDIF
            ELSE
              EXCLUDE_FLAG_B=TRIM(ADJUSTL(temp(INDEX(temp,',')+1:LEN(temp)        )))
              EXCLUDE_FLAG_C=''
            ENDIF
            N_exclude=N_exclude+1
            IF(I.EQ.1)write(*, '(A,A45,A25,A45,A5,A45)')"exclude: ",TRIM(EXCLUDE_FLAG_A),' in combination with ',TRIM(EXCLUDE_FLAG_B), &
              ' and ',TRIM(EXCLUDE_FLAG_C)
            IF(I.EQ.2)ExcludeConfigurations(N_exclude,1)=TRIM(EXCLUDE_FLAG_A)
            IF(I.EQ.2)ExcludeConfigurations(N_exclude,2)=TRIM(EXCLUDE_FLAG_B)
            IF(I.EQ.2)ExcludeConfigurations(N_exclude,3)=TRIM(EXCLUDE_FLAG_C)
          END IF
        END IF
      ELSE
        IF(INDEX(temp,'=').GT.0)THEN
          COMPILE_FLAG=TRIM(ADJUSTL(temp(1:INDEX(temp,'=')-1)))
          N_compile_flags=N_compile_flags+1
          IF(I.EQ.1)print*,"include: ",TRIM(COMPILE_FLAG)!,TRIM(temp(INDEX(temp,'='):LEN(temp)))
          IF(I.EQ.2)BuildConfigurations(N_compile_flags,1)=TRIM(ADJUSTL(COMPILE_FLAG))
          temp2=TRIM(ADJUSTL(temp(INDEX(temp,'=')+1:LEN(temp))))
          CurrentIndex=INDEX(temp2,',')
          IF(CurrentIndex.GT.0)THEN
            DO
              N_subinclude=N_subinclude+1
              IF(I.EQ.1)print*,N_subinclude,': ',                     TRIM(ADJUSTL(temp2(1:CurrentIndex-1)))
              IF(I.EQ.2)BuildConfigurations(N_compile_flags,N_subinclude+1)=TRIM(ADJUSTL(temp2(1:CurrentIndex-1)))
              temp2=temp2(CurrentIndex+1:LEN(temp2))
              NextIndex=INDEX(temp2(1:LEN(temp2)),',')
              IF(NextIndex.EQ.0)THEN
                N_subinclude=N_subinclude+1
                IF(I.EQ.1)print*,N_subinclude,': ',                     TRIM(ADJUSTL(temp2(1:LEN(temp2))))
                IF(I.EQ.2)BuildConfigurations(N_compile_flags,N_subinclude+1)=TRIM(ADJUSTL(temp2(1:LEN(temp2))))
                EXIT
              ELSE
                CurrentIndex=NextIndex
              END IF
              CurrentIndex=INDEX(temp2(1:LEN(temp2)),',')            
            END DO
          ELSE
            N_subinclude=N_subinclude+1
            IF(I.EQ.1)print*,N_subinclude,': ',                     TRIM(ADJUSTL(temp2(1:LEN(temp2))))
            IF(I.EQ.2)BuildConfigurations(N_compile_flags,N_subinclude+1)=TRIM(ADJUSTL(temp2(1:LEN(temp2))))
          END IF
          IF(I.EQ.2)BuildIndex(N_compile_flags)=N_subinclude
          nReggieBuilds=nReggieBuilds*N_subinclude
          N_subinclude_max=MAX(N_subinclude_max,N_subinclude)
        END IF
      END IF
    END IF
  END DO
  IF(I.EQ.1)print*,'The number of builds created by Reggie is: ',nReggieBuilds
  IF(I.EQ.1)print*,N_compile_flags
  IF(I.EQ.1)print*,N_subinclude_max
  IF(I.EQ.1)REWIND(ioUnit)
  IF((I.EQ.1).AND.(ALLOCATED(BuildConfigurations)))&
    CALL abort(__STAMP__&
    ,'Fortran runtime error: Attempting to allocate already allocated variable "BuildConfigurations"',iError,999.)
  IF((I.EQ.1).AND.(ALLOCATED(BuildConfigurations)))THEN
    SWRITE(UNIT_stdOut,'(A)') ' Fortran runtime error: Attempting to allocate already allocated variable "BuildConfigurations"'
    STOP
  END IF
  IF(I.EQ.1)ALLOCATE(BuildConfigurations(N_compile_flags,N_subinclude_max+1))
  IF(I.EQ.1)BuildConfigurations=''
  IF(I.EQ.1)ALLOCATE(BuildIndex(N_compile_flags))
  IF(I.EQ.1)BuildIndex=1
  IF(I.EQ.1)ALLOCATE(BuildCounter(N_compile_flags))
  IF(I.EQ.1)BuildIndex=1
  IF(I.EQ.1)ALLOCATE(ExcludeConfigurations(N_exclude,3))
  IF(I.EQ.1)ALLOCATE(BuildValid(nReggieBuilds))
  IF(I.EQ.1)BuildValid=.TRUE.
  IF(I.EQ.1)ALLOCATE(BuildValidInfo(nReggieBuilds))
  IF(I.EQ.1)BuildValidInfo=''
  IF(I.EQ.1)ALLOCATE(BuildEQNSYS(nReggieBuilds))
  IF(I.EQ.1)BuildEQNSYS=''
  IF(I.EQ.1)ALLOCATE(BuildTESTCASE(nReggieBuilds))
  IF(I.EQ.1)BuildTESTCASE='default'
  IF(I.EQ.1)ALLOCATE(BuildTIMEDISCMETHOD(nReggieBuilds))
  IF(I.EQ.1)BuildTIMEDISCMETHOD='default'
END DO


print*,"--- include ---"
print*,"BuildConfigurations(I,J) = ... "
DO I=1,N_compile_flags
  DO J=1,N_subinclude_max+1
      write(*, '(A25)', ADVANCE = "NO") TRIM(BuildConfigurations(I,J))
    IF(J.EQ.N_subinclude_max+1)print*,''
  END DO
END DO
print*,"--- exclude ---"
print*,"ExcludeConfigurations(I,J) = ... "
DO I=1,N_exclude
  DO J=1,3
      write(*, '(A40)', ADVANCE = "NO") TRIM(ExcludeConfigurations(I,J))
    IF(J.EQ.3)print*,''
  END DO
END DO
CLOSE(ioUnit)

!DO I=1,N_compile_flags
  !print*,BuildIndex(I)
!END DO



BuildCounter=1
DO I=1,nReggieBuilds
  !print*,BuildCounter
  DO J=1,N_exclude
    !IF
    InvalidA=.FALSE.
    InvalidB=.FALSE.
    InvalidC=.TRUE.
    DO K=1,N_compile_flags
      dummystr=TRIM(ADJUSTL(BuildConfigurations(K,1)))//'='//TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
      !print*,dummystr
      IF(dummystr.EQ.TRIM(ExcludeConfigurations(J,1)))InvalidA=.TRUE.
      IF(dummystr.EQ.TRIM(ExcludeConfigurations(J,2)))InvalidB=.TRUE.
      IF(dummystr.EQ.TRIM(ExcludeConfigurations(J,3)))InvalidC=.TRUE.

      END DO
  !print*,'tested against: ',TRIM(ExcludeConfigurations(J,1)),' and ',TRIM(ExcludeConfigurations(J,2))
  !print*,InvalidA,InvalidB
    IF(InvalidA.AND.InvalidB.AND.InvalidC)THEN
      BuildValidInfo(I)=TRIM(ExcludeConfigurations(J,1))//'+'//TRIM(ExcludeConfigurations(J,2))//'+'//TRIM(ExcludeConfigurations(J,3))
      BuildValid(I)=.FALSE.
    END IF
  END DO
  

! deprecated  ! TODO!!!! exchange this procedure with the actual compilation (if e.g. the needed flag is not set by hand!!)
! deprecated  ! set EQNSYSNAME
! deprecated  DO K=1,N_compile_flags
! deprecated    IF(TRIM(ADJUSTL(BuildConfigurations(K,1))).EQ.CodeNameUppCase//'_EQNSYSNAME')THEN
! deprecated      BuildEQNSYS(I)=TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
! deprecated    END IF
! deprecated  END DO
! deprecated  ! set BuildTESTCASE
! deprecated  DO K=1,N_compile_flags
! deprecated    IF(TRIM(ADJUSTL(BuildConfigurations(K,1))).EQ.CodeNameUppCase//'_TESTCASE')THEN
! deprecated      BuildTESTCASE(I)=TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
! deprecated    END IF
! deprecated  END DO
  
  ! print cmake compiler flags
  SWRITE(UNIT_stdOut, '(L)', ADVANCE = "NO") BuildValid(I)
  DO K=1,N_compile_flags
    !write(*, '(A)', ADVANCE = "NO") ' '//TRIM(BuildEQNSYS(I))
    SWRITE(UNIT_stdOut, '(A)', ADVANCE = "NO") ' -D'
    SWRITE(UNIT_stdOut, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,1)))
    SWRITE(UNIT_stdOut, '(A)', ADVANCE = "NO") '='
    SWRITE(UNIT_stdOut, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
  END DO
  SWRITE(UNIT_stdOut, '(A)', ADVANCE = "NO") '    '
  SWRITE(UNIT_stdOut, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildValidInfo(I)))
  SWRITE(UNIT_stdOut,*)
  
  
  
  
  ! get next build
  !write(*,*),''
  !read*
  DO J=1,N_compile_flags
    BuildCounter(J)=BuildCounter(J)+1
    IF(BuildCounter(J).GT.BuildIndex(J))THEN
      BuildCounter(J)=1
    ELSE
      EXIT
    END IF
  END DO
END DO

SWRITE(UNIT_stdOut, '(A)')"  "
SWRITE(UNIT_stdOut, '(I5,A4,I5,A12)')COUNT(BuildValid),' of ', nReggieBuilds,' are valid'
IF(BuildContinue)THEN
  SWRITE(UNIT_stdOut, '(A,I5,A1)') 'BuildContinue=.TRUE.    : Skipping builds [1] to [',BuildContinueNumber,']'
  IF(BuildContinueNumber.GT.nReggieBuilds)THEN
    SWRITE(UNIT_stdOut,'(A22,A)')          ' ERROR: ','The number of skipped builds exceeds the maxmum number of allocated builds.'
    SWRITE(UNIT_stdOut,'(A22,I5)') ' BuildContinueNumber: ',BuildContinueNumber 
    SWRITE(UNIT_stdOut,'(A22,I5)') '       nReggieBuilds: ', nReggieBuilds
    ERROR STOP '-1'
  END IF
  BuildValid(1:BuildContinueNumber)=.FALSE.
  SWRITE(UNIT_stdOut, '(I5,A4,I5,A12)')COUNT(BuildValid),' of ', nReggieBuilds,' are valid'
END IF

IF(COUNT(BuildValid).GT.MaxBuildConfigurations)THEN
  SWRITE(UNIT_stdOut,'(A)') ' ERROR: The number of builds exceeds the maxmum number allowed.'
  SWRITE(UNIT_stdOut,'(A,A)') ' COUNT(BuildValid)     :  ', COUNT(BuildValid)
  SWRITE(UNIT_stdOut,'(A,L)') ' MaxBuildConfigurations: ', MaxBuildConfigurations
  ERROR STOP '-1'
END IF

SWRITE(UNIT_stdOut,'(132("="))')

END SUBROUTINE ReadConfiguration

!==================================================================================================================================
!> reads the file "configurationsX.cmake" and creates a binary
!==================================================================================================================================
SUBROUTINE BuildConfiguration(iExample,iReggieBuild,nReggieBuilds,N_compile_flags)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: BuildDebug,BuildNoDebug,BuildEQNSYS,BuildTESTCASE,NumberOfProcs,NumberOfProcsStr
USE MOD_RegressionCheck_Vars,  ONLY: BuildContinue,BuildContinueNumber,BuildDir,BuildTIMEDISCMETHOD
USE MOD_RegressionCheck_Vars,  ONLY: CodeNameLowCase,CodeNameUppCase
USE MOD_RegressionCheck_tools, ONLY: SummaryOfErrors,AddError
USE MOD_RegressionCheck_Vars,  ONLY: BuildConfigurations,BuildValid,BuildCounter,BuildIndex,BuildSolver
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iExample,iReggieBuild,N_compile_flags,nReggieBuilds
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ioUnit,iSTATUS,J,K
!INTEGER                                   :: N_compile_flags,N_subinclude,N_exclude
CHARACTER(LEN=255)                        :: FileName!,temp,temp2,COMPILE_FLAG,dummystr
!CHARACTER(LEN=255)                        :: EXCLUDE_FLAG_A,EXCLUDE_FLAG_B
LOGICAL                                   :: ExistFile!,InvalidA,InvalidB
!CHARACTER(LEN=255),ALLOCATABLE            :: ExcludeConfigurations(:,:),BuildValidInfo(:)
!INTEGER                                   :: MaxBuildConfigurations=400,N_subinclude_max,N_compile_flags_max
!INTEGER,ALLOCATABLE                       :: BuildIndex(:),BuildCounter(:)
CHARACTER(LEN=500)                        :: SYSCOMMAND
CHARACTER(LEN=15)                         :: tempStr
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,I5,A4,I5,A,I5,A,I5,A)') &
"  Regression Check: Build Cmake Configurations",COUNT(BuildValid(1:iReggieBuild)),' of ',COUNT(BuildValid)&
                                            ,'  (',iReggieBuild                     ,'/'   ,nReggieBuilds    ,')'
SWRITE(UNIT_stdOut,'(132("="))')
ioUnit=GETFREEUNIT()
IF(BuildValid(iReggieBuild))THEN
  SYSCOMMAND='rm -rf '//TRIM(BuildDir)//'build_reggie > /dev/null 2>&1' ! clear building folder for next build
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  SYSCOMMAND='cd '//TRIM(BuildDir)//' && mkdir build_reggie'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  SWRITE(tempStr,*)iReggieBuild-1 ! print previously completed build to file for continuation possibility
  SYSCOMMAND='echo '//TRIM(tempStr)//' > '//TRIM(BuildDir)//'build_reggie/BuildContinue.reggie'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  OPEN(UNIT=ioUnit,FILE=TRIM(BuildDir)//'build_reggie/configurationX.cmake',STATUS="NEW",ACTION='WRITE',IOSTAT=iSTATUS)
  DO K=1,N_compile_flags ! print the compiler flag command line to "configurationX.cmake"
    SWRITE(ioUnit, '(A)', ADVANCE = "NO") ' -D'
    SWRITE(ioUnit, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,1)))
    SWRITE(ioUnit, '(A)', ADVANCE = "NO") '='
    SWRITE(ioUnit, '(A)', ADVANCE = "NO") TRIM(ADJUSTL(BuildConfigurations(K,BuildCounter(K)+1)))
  END DO
  CLOSE(ioUnit)
  SYSCOMMAND='cd '//TRIM(BuildDir)//'build_reggie && echo  `cat configurationX.cmake` '
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  SYSCOMMAND='cd '//TRIM(BuildDir)//&
   'build_reggie && cmake `cat configurationX.cmake` ../../ > build_reggie.out  && make '//CodeNameLowCase//' >> build_reggie.out'
  IF(BuildDebug)THEN
    SYSCOMMAND='cd '//TRIM(BuildDir)//'build_reggie && cmake `cat configurationX.cmake` ../../  && make '//CodeNameLowCase
  ELSEIF(BuildNoDebug)THEN
    SYSCOMMAND='cd '//TRIM(BuildDir)//&
   'build_reggie && cmake `cat configurationX.cmake` ../../ > build_reggie.out  && make '//CodeNameLowCase//&
                                                                                                         ' >> build_reggie.out 2>&1'
  END IF
  IF(NumberOfProcs.GT.1)SYSCOMMAND=TRIM(SYSCOMMAND)//' -j '//TRIM(ADJUSTL(NumberOfProcsStr))
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
  ! save compilation flags (even those that are not explicitly selected by the user) for deciding whether a supplied example folder 
  ! can be executed with the compiled executable or not
  IF(iSTATUS.EQ.0)THEN
    CALL GetFlagFromFile(TRIM(BuildDir)//'build_reggie/bin/configuration.cmake',&
                                                               CodeNameUppCase//'_TESTCASE'      ,BuildTESTCASE(iReggieBuild))
    ! set default for, e.g., PICLas code (currently no testcases are implemented)
    IF(BuildTESTCASE(iReggieBuild).EQ.'flag does not exist')BuildTESTCASE(iReggieBuild)='default'
    CALL GetFlagFromFile(TRIM(BuildDir)//'build_reggie/bin/configuration.cmake',&
                                                               CodeNameUppCase//'_TIMEDISCMETHOD',BuildTIMEDISCMETHOD(iReggieBuild))
    ! set default for, e.g., FLEXI (not a compilation flag)
    IF(BuildTESTCASE(iReggieBuild).EQ.'flag does not exist')BuildTESTCASE(iReggieBuild)='default'
    CALL GetFlagFromFile(TRIM(BuildDir)//'build_reggie/bin/configuration.cmake',&
                                                               CodeNameUppCase//'_EQNSYSNAME'     ,BuildEQNSYS(iReggieBuild))
  ELSE ! iSTATUS.NE.0 -> failed to compile cmake configuration build
    CALL AddError('Error while compiling',iExample,1,ErrorStatus=1,ErrorCode=1) !note: iSubExample is set to 1 as argument
    CALL SummaryOfErrors() ! Print the summary or examples and error codes (if they exist) before calling abort
    ! Print the error that caused the compilation abort (if it was written to "build_reggie.out" due to silent compilation)
    FileName=TRIM(BuildDir)//'build_reggie/build_reggie.out'
    INQUIRE(File=TRIM(FileName),EXIST=ExistFile)
    !print*,"ExistFile=",ExistFile
    IF(ExistFile) THEN ! build_reggie.out exists
      SWRITE(UNIT_stdOut, '(A)')' '
      SWRITE(UNIT_stdOut,'(132("="))')
      ! build [no-debug] = T
      ! build [debug]    = F
      ! build [-]        = T (but the error should already be printed to the screen!)
      SWRITE(UNIT_stdOut, '(A)')"Error output from: "//TRIM(FileName)
      SWRITE(UNIT_stdOut, '(A)')' '
      SYSCOMMAND='grep -rin error '//TRIM(FileName)
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      SWRITE(UNIT_stdOut,'(132("="))')
      SWRITE(UNIT_stdOut, '(A)')' '
    !ELSE
      !SWRITE(UNIT_stdOut, '(A)')TRIM(FileName)//" does not exist"
    END IF
    CALL abort(__STAMP__&
    ,'Could not compile '//CodeNameLowCase//'! iSTATUS=',iSTATUS,999.)
  END IF
  SWRITE(UNIT_stdOut,'(A)')"BuildEQNSYS(iReggieBuild)          = ["//TRIM(BuildEQNSYS(iReggieBuild))//"]"
  SWRITE(UNIT_stdOut,'(A)')"BuildTESTCASE(iReggieBuild)        = ["//TRIM(BuildTESTCASE(iReggieBuild))//"]"
  SWRITE(UNIT_stdOut,'(A)')"BuildTIMEDISCMETHOD(iReggieBuild)) = ["//TRIM(BuildTIMEDISCMETHOD(iReggieBuild))//"]"
  SYSCOMMAND='cd '//TRIM(BuildDir)//'build_reggie'
  IF(.NOT.BuildDebug)SYSCOMMAND=TRIM(SYSCOMMAND)//' && tail -n 1 build_reggie.out'
  CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
ELSE
  IF(BuildContinue)THEN
    SWRITE(UNIT_stdOut, '(A,I5,A1)') 'BuildContinue=.TRUE.    : Skipping build [1] to [',BuildContinueNumber,']... skipping...'
  ELSE
    SWRITE(UNIT_stdOut,'(A)')"invalid setup... skipping..."
  END IF
END IF


! get next build
DO J=1,N_compile_flags
  BuildCounter(J)=BuildCounter(J)+1
  IF(BuildCounter(J).GT.BuildIndex(J))THEN
    BuildCounter(J)=1
  ELSE
    EXIT
  END IF
END DO

SWRITE(UNIT_stdOut,'(132("="))')

END SUBROUTINE BuildConfiguration


!==================================================================================================================================
!> read compile flags from a specified file
!==================================================================================================================================
SUBROUTINE GetFlagFromFile(FileName,Flag,output)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals,            ONLY: Getfreeunit
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileName ! e.g. './../build_reggie/bin/configuration.cmake'
CHARACTER(LEN=*),INTENT(IN)  :: Flag     ! e.g. 'XX_EQNSYSNAME'
!INTEGER         :: a
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: output ! e.g. 'navierstokes'
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: ExistFile    ! file exists=.true., file does not exist=.false.
INTEGER                        :: iSTATUS      ! status
CHARACTER(LEN=255)             :: temp,temp2   ! temp variables for read in of file lines
INTEGER                        :: ioUnit       ! field handler unit and ??
INTEGER                        :: IndNum       ! Index Number
!===================================================================================================================================
output=''
INQUIRE(File=TRIM(FileName),EXIST=ExistFile)
IF(ExistFile) THEN
  ioUnit=GETFREEUNIT()
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
  DO
    READ(ioUnit,'(A)',iostat=iSTATUS)temp
    IF(ADJUSTL(temp).EQ.'!') CYCLE  ! complete line is commented out
    IF(iSTATUS.EQ.-1)EXIT           ! end of file is reached
    IF(LEN(trim(temp)).GT.1)THEN    ! exclude empty lines
      IndNum=INDEX(temp,TRIM(Flag)) ! e.g. 'XX_EQNSYSNAME'
      IF(IndNum.GT.0)THEN
        temp2=TRIM(ADJUSTL(temp(IndNum+LEN(TRIM(Flag)):LEN(temp))))
        IndNum=INDEX(temp2, '"')
        IF(IndNum.GT.0)THEN
          temp2=temp2(IndNum+1:LEN(TRIM(temp2)))
          IF(INDEX(temp2(IndNum+1:LEN(TRIM(temp2))), '"')+IndNum.GT.IndNum)THEN ! get binary path up to 2nd '"' in name
            IndNum=INDEX(temp2(IndNum+1:LEN(TRIM(temp2))), '"')+IndNum
          END IF
        END IF
        output=TRIM(ADJUSTL(temp2(1:IndNum-1)))
        EXIT
      END IF
    END IF
  END DO
  CLOSE(ioUnit)
  IF(output.EQ.'')output='flag does not exist'
ELSE 
  output='file does not exist'
END IF
END SUBROUTINE GetFlagFromFile


!==================================================================================================================================
!> Get the number of builds that have been successful in the previous building process in order to skip them
!> If, e.g., the first 10 builds were successful but the 11th failes, the first 10 might want to be skipped (in manual debugging)
!==================================================================================================================================
SUBROUTINE GetBuildContinue()
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_tools, ONLY:str2int
USE MOD_RegressionCheck_vars,  ONLY:BuildContinueNumber,BuildDir
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                                   :: ExistFile
INTEGER                                   :: ioUnit,iSTATUS
CHARACTER(LEN=255)                        :: FileName,temp
!===================================================================================================================================
ioUnit=GETFREEUNIT()
FileName=TRIM(BuildDir)//'build_reggie/BuildContinue.reggie'
INQUIRE(File=TRIM(FileName),EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A)')   ' ERROR: tried to continue building at specific point, file does not exist '
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:             ',TRIM(FileName)
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:            ',ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
END IF
READ(ioUnit,'(A)',iostat=iSTATUS)temp
CLOSE(ioUnit)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A10,A)')   ' ERROR:   ','tried to read BuildContinue.reggie'
  SWRITE(UNIT_stdOut,'(A10,A)') ' temp:    ',temp
  SWRITE(UNIT_stdOut,'(A10,I5)') ' iSTATUS: ',iSTATUS
  ERROR STOP '-1'
END IF
CALL str2int(temp,BuildContinueNumber,iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A22,A)')  ' ERROR:             ','tried to read BuildContinue.reggie, str2int failed'
  SWRITE(UNIT_stdOut,'(A22,I5)') ' BuildContinueNumber: ',BuildContinueNumber
  SWRITE(UNIT_stdOut,'(A22,I5)') ' iSTATUS:             ',iSTATUS
  ERROR STOP '-1'
END IF
IF(BuildContinueNumber.LT.0) THEN
  SWRITE(UNIT_stdOut,'(A22,A)')   ' ERROR:             ','BuildContinueNumber is < 0'
  ERROR STOP '-1'
END IF
END SUBROUTINE GetBuildContinue

END MODULE MOD_RegressionCheck_Build
