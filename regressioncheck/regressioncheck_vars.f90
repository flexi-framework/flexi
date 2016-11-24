!==================================================================================================================================
!> Contains global variables required by the regressioncheck 
!==================================================================================================================================
MODULE MOD_RegressionCheck_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=5),PARAMETER    :: CodeNameUppCase='FLEXI'             !> Code name in upper case letters. IMPORTANT: set its length!
CHARACTER(LEN=5),PARAMETER    :: CodeNameLowCase='flexi'             !> Code name in lower case letters. IMPORTANT: set its length!
INTEGER                        :: nErrors                            !> number of errors encountered during reggie execution
INTEGER                        :: NumberOfProcs                      !> number of processors for parallel build
CHARACTER(LEN=20)              :: NumberOfProcsStr                   !> number of processors for parallel build as string
INTEGER                        :: nExamples                          !> number of regressioncheck examples
CHARACTER(LEN=255),ALLOCATABLE :: ExampleNames(:)                    !> name of each example
CHARACTER(LEN=255)             :: RuntimeOption                      !> option for the regressioncheck: default (run), run and build
CHARACTER(LEN=255)             :: RuntimeOptionType                  !> specific option for the regressioncheck: default (run)
CHARACTER(LEN=255)             :: RuntimeOptionTypeII                !> specific option for the regressioncheck: default (empty)
CHARACTER(LEN=255)             :: RuntimeOptionTypeIII               !> specific option for the regressioncheck: default (empty)
CHARACTER(LEN=255)             :: EXECPATH                           !> path to solver incl. executable
CHARACTER(LEN=255)             :: ExamplesDir                        !> path to the regression check example folders
CHARACTER(LEN=255)             :: BuildDir                           !> path to the regression check building environment
CHARACTER(LEN=255),ALLOCATABLE :: BuildEQNSYS(:)                     !> EQNSYS for each build
CHARACTER(LEN=255),ALLOCATABLE :: BuildTESTCASE(:)                   !> TESTCASE for each build: only FLEXI
CHARACTER(LEN=255),ALLOCATABLE :: BuildTIMEDISCMETHOD(:)             !> TIMEDISCMETHOD for each build: only PICLas

LOGICAL                        :: BuildSolver                        !> Flag for automatic building of different flexi cmake configs
LOGICAL                        :: BuildDebug                         !> Prints the complete compilation process for debugging when
                                                                     !> BuildSolver is true 
LOGICAL                        :: BuildNoDebug                       !> Don't print any compiler output (if BuildSolver is true) 
LOGICAL                        :: BuildContinue                      !> allow the building sequence to begin at the last failure
INTEGER                        :: BuildContinueNumber                !> start building sequence from this point

TYPE tExample                                                        !> examples for regressioncheck
  INTEGER                                :: ReferenceType            !> Type of reference
                                                                     !> 0 - no reference
                                                                     !> 1 - L2 and Linf
                                                                     !> 2 - state file   
                                                                     !> 3 - state file  and L2 
  CHARACTER(LEN=255)                     :: EQNSYSNAME               !> Name of current EQNSYS (depends on current exe)
  INTEGER                                :: Nvar                     !> Size of EQNSYS 
  CHARACTER(LEN=255)                     :: PATH                     !> Path to example
  LOGICAL                                :: MPIrun                   !> execution information (MPI)
  INTEGER                                :: MPIthreads               !> number of MPI threads for execution
  CHARACTER(LEN=255)                     :: ReferenceFile            !> Name of references L2/LInf
  REAL                                   :: ReferenceTolerance       !> optional tolerance for L2/LInf
  CHARACTER(LEN=255)                     :: ReferenceStateFile       !> Name of reference state file
  CHARACTER(LEN=255)                     :: CheckedStateFile         !> Name of checked state file
  CHARACTER(LEN=255)                     :: ReferenceDataSetName     !> Name of Dataset in hdf5 file for comparision
  CHARACTER(LEN=255)                     :: RestartFileName          !> Name of RestartFile
  INTEGER                                :: ErrorStatus              !> ErrorStatus
                                                                     !> 0 - success
                                                                     !> 1 - failed during execution
                                                                     !> 2 - test failed
  CHARACTER(LEN=255)                     :: IntegrateLineFile        !> File name with ACSI number columns
  INTEGER                                :: IntegrateLineRange(2)    !> the numerbs of two coulumns with data
  REAL                                   :: IntegrateLineValue       !> the reference integral value
  CHARACTER(LEN=255)                     :: IntegrateLineDelimiter   !> delimiter string for reading the data file
  INTEGER                                :: IntegrateLineHeaderLines !> number of header lines to be ignored from data file
  LOGICAL                                :: IntegrateLine            !> read two columns from a file and integrate over line
                                                                     !> e.g. u(t) is integrated over t for comparison of the 
                                                                     !> integral value
  INTEGER                                :: SubExampleNumber         !> Numbers of sub examples, currently fixed to 1
  CHARACTER(LEN=255)                     :: SubExampleOption(20)     !> for each sub example class, currently 10 options are allowed
  CHARACTER(LEN=255)                     :: SubExample               !> sub example class, e.g., TimeDiscMethod can be chosen for 
                                                                     !> testing multiple time integration schemes
END TYPE

TYPE(tExample), ALLOCATABLE              :: Examples(:)              !> container with variables for each reggie example

TYPE tEC                                                             !> Type to simplify error handling
  INTEGER            :: ErrorCode                                    !> interger code of error
  CHARACTER(LEN=255) :: Example                                      !> name of the example
  CHARACTER(LEN=255) :: SubExample                                   !> name of the subexample
  CHARACTER(LEN=255) :: SubExampleOption                             !> name of the subexample option
  CHARACTER(LEN=255) :: Info                                         !> name of the example
  CHARACTER(LEN=255) :: Build                                        !> flexi cmake build name
  Type(tEC),Pointer  :: nextError                                    !> pointer to next error if several errors occure
END TYPE tEC
TYPE(tEC), POINTER  :: firstError, aError                            !> pointer to first error and looping pointer

!==================================================================================================================================
END MODULE MOD_RegressionCheck_Vars
