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
CHARACTER(LEN=5),PARAMETER    :: CodeNameUppCase='FLEXI'            !> Code name in upper case letters. IMPORTANT: set its length!
CHARACTER(LEN=5),PARAMETER    :: CodeNameLowCase='flexi'            !> Code name in lower case letters. IMPORTANT: set its length!
INTEGER                        :: nErrors                           !> number of errors encountered during reggie execution
INTEGER                        :: GlobalRunNumber                   !> count the number of separate runs for listing in summary
INTEGER                        :: NumberOfProcs                     !> number of processors for parallel build
CHARACTER(LEN=20)              :: NumberOfProcsStr                  !> number of processors for parallel build as string
INTEGER                        :: nExamples                         !> number of regressioncheck examples
CHARACTER(LEN=255),ALLOCATABLE :: ExampleNames(:)                   !> name of each example
CHARACTER(LEN=255)             :: RuntimeOption                     !> option for the regressioncheck: default (run), run and build
CHARACTER(LEN=255)             :: RuntimeOptionType                 !> specific option for the regressioncheck: default (run)
CHARACTER(LEN=255)             :: RuntimeOptionTypeII               !> specific option for the regressioncheck: default (empty)
CHARACTER(LEN=255)             :: RuntimeOptionTypeIII              !> specific option for the regressioncheck: default (empty)
CHARACTER(LEN=255)             :: EXECPATH                          !> path to solver incl. executable
CHARACTER(LEN=255)             :: ExamplesDir                       !> path to the regression check example folders
CHARACTER(LEN=255)             :: readRHS(2)                        !> parameter from parameter_reggie.ini: right hand side 
                                                                    !> parameter name (1) and setting (2)
CHARACTER(LEN=255)             :: BuildDir                          !> path to the regression check building environment
CHARACTER(LEN=255),ALLOCATABLE :: BuildEQNSYS(:)                    !> EQNSYS for each build
CHARACTER(LEN=255),ALLOCATABLE :: BuildTESTCASE(:)                  !> TESTCASE for each build: only FLEXI
CHARACTER(LEN=255),ALLOCATABLE :: BuildTIMEDISCMETHOD(:)            !> TIMEDISCMETHOD for each build: only PICLas
CHARACTER(LEN=255),ALLOCATABLE :: BuildMPI(:)                       !> ON/OFF: build is created with/without MPI
CHARACTER(LEN=255),ALLOCATABLE :: BuildFV(:)                        !> ON/OFF: build is created with/without FV (finite volume)
CHARACTER(LEN=255),ALLOCATABLE :: BuildCODE2D(:)                    !> ON/OFF: build is created with/without 2D only
CHARACTER(LEN=255),ALLOCATABLE :: BuildPARABOLIC(:)                 !> ON/OFF: build is created with/without PARABOLIC terms

LOGICAL                        :: BuildSolver                       !> Flag for automatic building of different flexi cmake configs
LOGICAL                        :: BuildDebug                        !> Prints the complete compilation process for debugging when
                                                                    !> BuildSolver is true 
LOGICAL                        :: BuildNoDebug                      !> Don't print any compiler output (if BuildSolver is true) 
LOGICAL                        :: BuildContinue                     !> allow the building sequence to begin at the last failure
INTEGER                        :: BuildContinueNumber               !> start building sequence from this point

CHARACTER(LEN=255),ALLOCATABLE :: BuildConfigurations(:,:)          !> CMAKE complie flag and value
LOGICAL,ALLOCATABLE            :: BuildValid(:)                     !> use the configuration or don't
INTEGER,ALLOCATABLE            :: BuildCounter(:)                   !> register for creaating all possible cmake configurations
INTEGER,ALLOCATABLE            :: BuildIndex(:)                     !> number of different flag settings for each specified cmake 
                                                                    !> compiler flag

TYPE tExample                                                       !> examples for regressioncheck
  INTEGER                          :: ReferenceType                 !> Type of reference
                                                                    !> 0 - no reference
                                                                    !> 1 - L2 and Linf
                                                                    !> 2 - state file   
                                                                    !> 3 - state file  and L2 
  CHARACTER(LEN=255)               :: EQNSYSNAME                    !> Name of current EQNSYS (depends on current exe)
  INTEGER                          :: Nvar                          !> Size of EQNSYS 
  CHARACTER(LEN=255)               :: PATH                          !> Path to example
  LOGICAL                          :: MPIrun                        !> execution information (MPI)
  CHARACTER(LEN=255)               :: MPIcommand                    !> e.g. 'aprun' or 'mpirun' (default)
  CHARACTER(LEN=15)                :: MPIthreadsStr(100)            !> list (array) with the number of MPI threads for execution
  INTEGER                          :: MPIthreadsN                   !> dimension (vector length) of MPIthreadsStr
  INTEGER                          :: nRuns                         !> number of runs with a specific setup (number of MPI threads)
  CHARACTER(LEN=15)                :: NumberOfCellsStr(100)         !> list (array) with number(s) of DG cells in one direction
  INTEGER                          :: NumberOfCellsN                !> dimension of list with number(s) of DG cells in one direction
  CHARACTER(LEN=255)               :: ReferenceFile                 !> Name of references L2/LInf
  CHARACTER(LEN=255)               :: ReferenceNormFile             !> Name of reference file (arbitrary file, e.g., *.csv)
  REAL                             :: ReferenceTolerance            !> optional tolerance for L2/LInf
  CHARACTER(LEN=255)               :: H5DIFFReferenceStateFile      !> Name of reference state file
  CHARACTER(LEN=255)               :: H5DIFFCheckedStateFile        !> Name of checked state file
  CHARACTER(LEN=255)               :: H5DIFFReferenceDataSetName    !> Name of Dataset in hdf5 file for comparision
  CHARACTER(LEN=255)               :: H5diffToleranceType           !> type of tolerance for h5diff: relative or absolute
  REAL                             :: H5diffTolerance               !> value used for the tolerance check in h5diff
  CHARACTER(LEN=255)               :: RestartFileName               !> Name of RestartFile
  INTEGER                          :: ErrorStatus                   !> ErrorStatus
                                                                    !> 0 - success
                                                                    !> 1 - failed during execution
                                                                    !> 2 - test failed
  CHARACTER(LEN=255)               :: IntegrateLineFile             !> File name with ACSI number columns
  INTEGER                          :: IntegrateLineRange(2)         !> the numerbs of two coulumns with data
  REAL                             :: IntegrateLineValue            !> the reference integral value
  CHARACTER(LEN=255)               :: IntegrateLineDelimiter        !> delimiter string for reading the data file
  INTEGER                          :: IntegrateLineHeaderLines      !> number of header lines to be ignored from data file
  LOGICAL                          :: IntegrateLine                 !> read two columns from a file and integrate over line
                                                                    !> e.g. u(t) is integrated over t for comparison of the 
                                                                    !> integral value
  CHARACTER(LEN=255)               :: CompareDatafileRowFile        !> File name with ACSI number columns
  CHARACTER(LEN=255)               :: CompareDatafileRowRefFile     !> File name with ACSI number columns for reference
  CHARACTER(LEN=255)               :: CompareDatafileRowDelimiter   !> delimiter string for reading the data file
  INTEGER                          :: CompareDatafileRowNumber      !> number of the row chosen for comparison
  REAL                             :: CompareDatafileRowTolerance   !> Tolerance value for comparison
  LOGICAL                          :: CompareDatafileRowReadHeader  !> read the first row of the header (column labels)
  INTEGER                          :: CompareDatafileRowHeaderLines !> number of header lines to be ignored from data file
  LOGICAL                          :: CompareDatafileRow            !> read a single row from a file and compare each entry to
                                                                    !> a reference file (each failed comparison will be dispayed)
  REAL                             :: CompareHDF5ArrayBoundsValue(2)!> value ranges for comparison
  INTEGER                          :: CompareHDF5ArrayBoundsRange(2)!> HDF5 array dim ranges
  CHARACTER(LEN=255)               :: CompareHDF5ArrayBoundsName    !> array name in HDF5 file
  CHARACTER(LEN=255)               :: CompareHDF5ArrayBoundsFile    !> name of HDF5 file
  LOGICAL                          :: CompareHDF5ArrayBounds        !> read an array from a HDF5 file and compare certain entry 
                                                                    !> bounds that must be limited to a supplied value range


  CHARACTER(LEN=255)               :: ConvergenceTestType           !> h- or p-convergence test
  REAL                             :: ConvergenceTestDomainSize     !> length of simulation domain, needed for grid step size
  REAL                             :: ConvergenceTestValue          !> single value for comparison
  REAL                             :: ConvergenceTestTolerance      !> relative tolerance when comparing the "ConvergenceTestValue"
  REAL, ALLOCATABLE                :: ConvergenceTestGridSize(:)    !> array for grid step size: cell length / ( p + 1 )
  REAL, ALLOCATABLE                :: ConvergenceTestError(:,:)     !> array for L2 errors over iteration or polynomial degree
                                                                    !> dimension for "array" will be: [SubExampleNumber]x[nVar]x[2]
                                                                    !>                array(:,1)=[N] 
                                                                    !>                array(:,2)=[iVar]
                                                                    !>                array(:,3)=[L2]
  LOGICAL                          :: ConvergenceTest               !> run convergence test in combination with "SubExample" for
                                                                    !> either "N" (p-convergence) or "MeshFile" (h-convergence)
  INTEGER                          :: SubExampleNumber              !> Numbers of sub examples, currently fixed to 1
  CHARACTER(LEN=255)               :: SubExampleOption(100)         !> for each sub example class, currently 10 options are allowed
  CHARACTER(LEN=255)               :: SubExample                    !> sub example class, e.g., TimeDiscMethod can be chosen for 
                                                                    !> testing multiple time integration schemes
END TYPE

TYPE(tExample), ALLOCATABLE        :: Examples(:)                   !> container with variables for each reggie example

TYPE tEC                                                            !> Type to simplify error handling
  INTEGER            :: ErrorCode                                   !> interger code of error
  INTEGER            :: RunNumber                                   !> number of current run
  CHARACTER(LEN=255) :: Example                                     !> name of the example
  CHARACTER(LEN=255) :: SubExample                                  !> name of the subexample
  CHARACTER(LEN=255) :: SubExampleOption                            !> name of the subexample option
  CHARACTER(LEN=255) :: Info                                        !> name of the example
  CHARACTER(LEN=255) :: MPIthreadsStr                               !> number of used MPI threads or '-' for single computation
  CHARACTER(LEN=255) :: Build                                       !> flexi cmake build name
  Type(tEC),Pointer  :: nextError                                   !> pointer to next error if several errors occure
END TYPE tEC
TYPE(tEC), POINTER  :: firstError, aError                           !> pointer to first error and looping pointer

!==================================================================================================================================
END MODULE MOD_RegressionCheck_Vars
