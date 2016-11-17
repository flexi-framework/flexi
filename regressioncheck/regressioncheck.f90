#include "flexi.h"

!==================================================================================================================================
!> The regressioncheck tool performs several test to verify the correct behaviour of the solver. The regressioncheck tool is an
!> optional tool which is built by the regressioncheck flag. Only if the tool is built, the regressioncheck examples are checked
!> out.
!> Each example consists of the mesh-file (h5), parameter file, regressiocheck.ini and a reference solution. In order to deal with
!> different compiler, a relative high tolerance is set to 100*epsMach. Please note, that this scaling factor can be modified by
!> the user.
!> Usage: ./regressioncheck run   - uses the already built versions of flexi and only runs the examples with the given executable
!>        ./regressioncheck build - previous to the execution and comparison step, flexi is built with all possible 
!>                                - parameter combinations. each combination is tested with each example
!> error codes are handled by a pointer list and summarized at the end of the program
!> error codes: 0 - no error
!>              1 - failed during build
!>              2 - computation of example failed
!>              3 - mismatch in norms
!>              4 - mismatch in dataset
!>             77 - no flexi executable found for option run
!>             99 - fail of execute_system_command
!==================================================================================================================================
PROGRAM RegressionCheck
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_tools, ONLY: InitExample,CleanExample,GetExampleList,CheckForExecutable,GetCommandLineOption
USE MOD_RegressionCheck_Run,   ONLY: PerformRegressionCheck
USE MOD_RegressionCheck_Vars,  ONLY: ExampleNames,Examples,firstError,aError,BuildSolver
USE MOD_MPI,                   ONLY: InitMPI,DefineParametersMPI
USE MOD_Mesh,                  ONLY: FinalizeMesh
!#ifdef USE_MPI
!USE MOD_MPI_Vars,            ONLY: NbProc,nMPISides_Proc
!#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time                              ! Used to track computation time  
!INTEGER                        :: iExample                       ! Loop counters 
!INTEGER                        :: iSTATUS                           ! system-command status
INTEGER                        :: ioUnit,nReggieBuilds              ! field handler unit and ??
INTEGER                        :: nErrors                           ! number of errors
CHARACTER(LEN=500)             :: SYSCOMMAND                        ! string to fit the system command
CHARACTER(LEN=255)             :: FileName                          ! filename
CHARACTER(LEN=255)             :: tmpstr                           ! tmp variable
!==================================================================================================================================
! errorcodes
ALLOCATE(firstError)
firstError%ErrorCode=-1
NULLIFY(aError)
nReggieBuilds=0
SYSCOMMAND=''
FileName=''
ioUnit=GETFREEUNIT()
CALL InitMPI()
! Define parameters for Converter

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
"  Little ReggressionCheck, add nice ASCII art here"
SWRITE(UNIT_stdOut,'(132("="))')

! get the command line option
CALL GetCommandLineOption()

! set paths for execution
IF(.NOT.BuildSolver) CALL CheckForExecutable(Mode=1)

! Measure regressioncheck runtime 
StartTime=FLEXITIME()

! check if examples are checked out and get list
CALL GetExampleList()

! perform the regressioncheck
CALL PerformRegressionCheck()


!   ! clean all successful tests
!   DO iExample = 1, nExamples
!     ! insert code here
!     IF(Examples(iExample)%ErrorStatus.EQ.0) CALL CleanExample(iExample)
!   END DO ! iExample=1,nExamples

! deallocate example names and example type
DEALLOCATE(ExampleNames)
DEALLOCATE(Examples)

! Measure processing duration
Time=FLEXITIME()
#ifdef USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError,999.)
#endif

SWRITE(UNIT_stdOut,'(132("="))')
nErrors=0
IF(.NOT.ASSOCIATED(aError))THEN
  SWRITE(UNIT_stdOut,'(A)') ' No Examples were executed'
ELSE
  NULLIFY(aError%nextError) ! nullyfy unused next error pointer
  SWRITE(UNIT_stdOut,'(A)') ' Summary of Errors (0=no Error): '
  SWRITE(UNIT_stdOut,'(A)') ' '
  aError=>firstError ! set aError to first error in list
  tmpstr=''
  SWRITE(UNIT_stdOut,'(A45,2x,A20,2x,A10,2x,A10,2x,A65,2x)') 'Example','SubExample','ErrorCode','build','Information'
  DO WHILE (ASSOCIATED(aError))
    IF(TRIM(tmpstr).NE.TRIM(aError%Example))THEN
      SWRITE(UNIT_stdOut,*) ''
    END IF
    tmpstr=TRIM(aError%Example)
    SWRITE(UNIT_stdOut,'(A45,2x)',ADVANCE='no') TRIM(aError%Example)
    SWRITE(UNIT_stdOut,'(A20,2x)',ADVANCE='no') TRIM(aError%SubExampleOption)
    SWRITE(UNIT_stdOut,'(I10,2x)',ADVANCE='no') aError%ErrorCode
    SWRITE(UNIT_stdOut,'(A10,2x)',ADVANCE='no') TRIM(aError%Build)
    SWRITE(UNIT_stdOut,'(A65,2x)',ADVANCE='no') TRIM(aError%Info)
    SWRITE(UNIT_stdOut,*) ''
    IF(aError%ErrorCode.NE.0) nErrors=nErrors+1
    aError=>aError%nextError
  END DO
  SWRITE(UNIT_stdOut,'(A,I4)') ' Number of errors:  ', nErrors
END IF

IF(nErrors.GT.0)THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' RegressionCheck FAILED! [',Time-StartTime,' sec ]'
  SWRITE(UNIT_stdOut,'(132("-"))')
  ERROR STOP '999'
ELSE
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' RegressionCheck SUCCESSFUL! [',Time-StartTime,' sec ]'
  SWRITE(UNIT_stdOut,'(132("-"))')
END IF
SWRITE(UNIT_stdOut,'(132("="))')
END PROGRAM RegressionCheck
