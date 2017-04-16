#include "flexi.h"

!==================================================================================================================================
!> The regressioncheck tool performs several test to verify the correct behaviour of the solver. The regressioncheck tool is an
!> optional tool which is built by the regressioncheck flag. Only if the tool is built, the regressioncheck examples are checked
!> out.
!> Each example consists of the mesh-file (h5), parameter file, regressiocheck.ini and a reference solution. In order to deal with
!> different compiler, a relative high tolerance is set to 100*epsMach. Please note, that this scaling factor can be modified by
!> the user.
!> Usage: ./regressioncheck run   - uses the prev. built binaries and only runs the examples with the given executable
!>        ./regressioncheck build - previous to the execution and comparison step, binaries are built with all possible 
!>                                - parameter combinations. each combination is tested with each example
!> error codes are handled by a pointer list and summarized at the end of the program
!> error codes: 0 - no error
!>              1 - failed during build
!>              2 - computation of example failed
!>              3 - mismatch in norms
!>              4 - mismatch in dataset
!>             77 - no executable found for option run
!>             99 - fail of execute_system_command
!==================================================================================================================================
PROGRAM RegressionCheck
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_tools, ONLY: InitExample,GetExampleList,CheckForExecutable,GetCommandLineOption
USE MOD_RegressionCheck_tools, ONLY: SummaryOfErrors
USE MOD_RegressionCheck_Run,   ONLY: PerformRegressionCheck,PerformFullRegressionCheck
USE MOD_RegressionCheck_Vars,  ONLY: ExampleNames,Examples,firstError,aError,BuildSolver,nErrors
USE MOD_RegressionCheck_Vars,  ONLY: CodeNameUppCase,CodeNameLowCase,DoFullReggie
USE MOD_MPI,                   ONLY: InitMPI
USE MOD_Mesh,                  ONLY: FinalizeMesh
USE MOD_RegressionCheck_tools, ONLY: REGGIETIME
#if USE_MPI
USE MOD_MPI,                     ONLY: FinalizeMPI
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: EndTime             ! time at the end of the reggie execution
INTEGER                        :: nReggieBuilds       ! number of different cmake builds (with different flags)
CHARACTER(LEN=500)             :: SYSCOMMAND          ! string to fit the system command
CHARACTER(LEN=255)             :: FileName            ! filename
!==================================================================================================================================
! errorcodes
ALLOCATE(firstError)
firstError%ErrorCode=-1
NULLIFY(aError)
nReggieBuilds=0
SYSCOMMAND=''
FileName=''
!IF(CodeNameUppCase.EQ.'BOLTZPLATZ')CALL InitGlobals() ! only "boltzplatz"
CALL InitMPI()
! Check Code Names
IF(LEN(CodeNameUppCase).NE.LEN(ADJUSTL(TRIM(CodeNameUppCase))))       CALL abort(&
  __STAMP__&
  ,'CodeNameUppCase=['//CodeNameUppCase//']: the variable has an incorrect length!')
IF((CodeNameUppCase.NE.'FLEXI').AND.(CodeNameUppCase.NE.'BOLTZPLATZ'))CALL abort(&
  __STAMP__&
  ,'CodeNameUppCase=['//CodeNameUppCase//']: the code name is unknown! Add here if it is correct.')
IF(LEN(CodeNameLowCase).NE.LEN(ADJUSTL(TRIM(CodeNameLowCase))))       CALL abort(&
  __STAMP__&
  ,'CodeNameLowCase=['//CodeNameLowCase//']: the variable has an incorrect length!')
IF((CodeNameLowCase.NE.'flexi').AND.(CodeNameLowCase.NE.'boltzplatz'))CALL abort(&
  __STAMP__&
  ,'CodeNameLowCase=['//CodeNameLowCase//']: the code name is unknown! Add here if it is correct.')


SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')"  Little ReggressionCheck, add nice ASCII art here"
SWRITE(UNIT_stdOut,'(132("="))')

! get the command line option
CALL GetCommandLineOption()

! set paths for execution
IF(.NOT.BuildSolver) CALL CheckForExecutable(Mode=1)

! Measure regressioncheck runtime 
StartTime=REGGIETIME()

IF(DoFullReggie)THEN ! call regressioncheck recursivly using the commands from gitlab-ci.yml
  CALL PerformFullRegressionCheck()
ELSE
  ! check if examples are checked out and get list
  CALL GetExampleList()
  
  ! perform the regressioncheck
  CALL PerformRegressionCheck()
  
  ! deallocate example names and example type
  DEALLOCATE(ExampleNames)
  DEALLOCATE(Examples)
END IF

! Measure processing duration
EndTime=REGGIETIME()

!   #if USE_MPI
!   CALL MPI_FINALIZE(iError)
!   IF(iError .NE. 0) CALL abort(&
!     __STAMP__&
!     ,'MPI finalize error',iError,999.)
!   #endif

! Print the summary or examples and error codes (if they exist)
CALL SummaryOfErrors(EndTime)

#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
CALL FinalizeMPI()
#endif

IF(nErrors.GT.0) ERROR STOP '999'
END PROGRAM RegressionCheck
