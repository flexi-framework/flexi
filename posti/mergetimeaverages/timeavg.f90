#include "flexi.h"

!===================================================================================================================================
!> This tool will take pre-averaged files (TimeAvg, Flucs) or simple state files 
!> and perform global temporal averaging
!===================================================================================================================================
PROGRAM TimeAvg
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Input, ONLY:ReadAttribute
USE MOD_HDF5_Output,ONLY:WriteTimeAverage
USE MOD_State,      ONLY:InitState,readStateFile,FinalizeState,VarNames
USE MOD_MPI,        ONLY:InitMPI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: AvgTime,TotalAvgTime,TotalAvgTimeGlobal
INTEGER                              :: iArg,nArgs,iExt
INTEGER                              :: nElemsRef
INTEGER                              :: nVarRef,NRef
CHARACTER(LEN=255)                   :: InputFile,LastInputFile 
CHARACTER(LEN=255)                   :: DataSet=''
CHARACTER(LEN=255)                   :: MeshFileRef,NodeTypeRef
CHARACTER(LEN=255)                   :: ProgramNameRef,FileTypeRef
CHARACTER(LEN=255)                   :: Filename_out,Dummy255
REAL,ALLOCATABLE                     :: UAvg(:,:,:,:,:)
REAL,ALLOCATABLE                     :: UAvg_VMS(:,:,:,:,:)
LOGICAL                              :: isTimeAvg
REAL                                 :: AvgStarttime,Time,AvgEndTime
INTEGER                              :: StartArgs,nFiles,iFile
INTEGER                              :: coarsenFac,iCoarse
REAL                                 :: tol=1.E-7
REAL                                 :: dt
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF (nArgs.GT.2) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: timeavg --start=<starttime> --end=<endtime> FILE1 FILE2 .. FILEN or timeavg --help'// &
  '[option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF

! First check if the first three arguments contain -start=, -end= or -coarsen=
StartArgs=1
AvgStarttime=-HUGE(1.)
AvgEndtime=HUGE(1.)
coarsenFac=999999999
DO iArg=1,MIN(nArgs,3)
  tmp = Args(iArg)

  ! check if the -start= flag is used, using Inputfile as a dummystring
  IF (STRICMP(tmp(1:7), "-start=")) THEN
    StartArgs=StartArgs+1
    READ(TRIM(InputFile(8:LEN(tmp))),*) AvgStarttime
    WRITE(UNIT_stdOut,'(132("="))')
    WRITE(UNIT_stdOut,'(A28,F16.6)') ' Start time for Averaging is ',AvgStarttime
  ELSEIF (STRICMP(tmp(1:5), "-end=")) THEN
    StartArgs=StartArgs+1
    READ(TRIM(InputFile(6:LEN(tmp))) ,*) AvgEndtime
    WRITE(UNIT_stdOut,'(132("="))')
    WRITE(UNIT_stdOut,'(A28,F16.6)') ' End time for Averaging is ',AvgEndtime
  ELSEIF (STRICMP(tmp(1:9), "-coarsen=")) THEN
    StartArgs=StartArgs+1
    READ(TRIM(InputFile(10:LEN(tmp))),*) coarsenFac
    WRITE(UNIT_stdOut,'(132("="))')
    WRITE(UNIT_stdOut,'(A28,F16.6)') ' Coarsening factor is ',coarsenFac
  END IF

END DO

! check if at least 2 timeavg files are there
nFiles=nArgs-StartArgs+1
IF(nFiles.LT.2) CALL CollectiveStop(__STAMP__,'At least two files required for averaging!')

! Get default parameters
CALL InitState(InputFile)
nVarRef    = nVar_HDF5 
NRef       = N_HDF5 
nElemsRef  = nGlobalElems
NodeTypeRef=NodeType_HDF5
MeshFileRef=MeshFile_HDF5

SELECT CASE(TRIM(ProgramName))
CASE('Flexi')
  DataSet=''
CASE DEFAULT
  CALL Abort(__STAMP__,'Unknown program.')
END SELECT
ProgramNameRef=TRIM(ProgramName)

! Check whether the file is a normal state file or already contains time averaged data 
SELECT CASE(TRIM(FileType))
CASE('State')
  isTimeAvg=.FALSE.
CASE('TimeAvg')
  isTimeAvg=.TRUE.
CASE('Fluc')
  isTimeAvg=.TRUE.  ! since we have time-averaged correlations
CASE DEFAULT
  CALL Abort(__STAMP__,'Unknown file type.')
END SELECT
FileTypeRef=TRIM(FileType)
CALL FinalizeState()

ALLOCATE(UAvg(nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,nGlobalElems))
UAvg=0.
!IF(nVarVMS_HDF5.GT.0)THEN
!  ALLOCATE(UAvg_VMS(nVarVMS_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,nGlobalElems))
!  UAvg_VMS=0.
!END IF

! Start the averaging
WRITE(UNIT_stdOut,'(132("="))')
TotalAvgTime=0.
TotalAvgTimeGlobal=0.
iCoarse=0

DO iFile=1,nFiles
  InputFile=Args(iFile+StartArgs-1)
  ! check local time 
  CALL OpenDataFile(InputFile,.FALSE.)
  dt=-Time
  CALL ReadAttribute(File_ID,'Time',1,TRIM(DataSet),RealScalar=Time)
  CALL CloseDataFile()
  dt=dt+Time
  ! If the time in the current file is below the starttime, then cycle
  IF ((Time-AvgStartTime.LT.-tol).OR.(Time-AvgEndTime.GT.tol)) THEN 
    WRITE(UNIT_stdOut,'(A,A)') ' SKIPPING FILE ',TRIM(InputFile)
    nSkipped=nSkipped+1
    ! check if the starttime is later than the last filetime
    IF (StartArgs.GT.nArgs) THEN
      CALL CollectiveStop(__STAMP__,'Starttime is greater than latest Filetime!')
    END IF
    IF (nFiles-nSkipped.LE.1) THEN
      CALL CollectiveStop(__STAMP__,'No file has time lower than endtime!')
    END IF
    CYCLE
  END IF

  WRITE(UNIT_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iFile,' of ',nFiles,' FILES.'
  WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputFile),'" )'
  WRITE(UNIT_stdOut,'(132("="))')

  CALL InitState(InputFile)

  IF(TRIM(FileType).NE.TRIM(FileTypeRef))THEN
    CALL Abort(__STAMP__, 'Mixing files of type '//TRIM(FileTypeRef)//' with files of type '&
                           //TRIM(FileType)//' is not allowed!')
  END IF
  
  ! Read nVar_in, N_in and nElements from file
  avgTime=1.
  CALL OpenDataFile(InputFile,.FALSE.)
  IF(isTimeAvg)THEN
    CALL ReadAttribute(File_ID,'AvgTime',1,TRIM(DataSet),RealScalar=avgTime)
    WRITE(UNIT_stdOut,'(A,F10.5)') ' Time averaged file, averaging time is: ',avgTime 
  ELSE
    WRITE(UNIT_stdOut,'(A,A)') ' Normal state file, each file will be weighted identically.'
  END IF
  CALL CloseDataFile()

  ! Check if input data size has changed
  IF((NRef.NE.N_HDF5).OR.(nVarRef.NE.nVar_HDF5).OR.(nElemsRef.NE.nGlobalElems))&
    CALL Abort(__STAMP__,'Change of polynomial degree and variables not supported!')
  IF(TRIM(NodeTypeRef).NE.NodeType_HDF5)&
    CALL Abort(__STAMP__,'Change of node type not supported yet!')
  IF(TRIM(MeshFileRef).NE.TRIM(MeshFile_HDF5)) &
    CALL Abort(__STAMP__,'Change of meshfile not supported yet!')
  !TODO: Cover switch DG <--> FV

  CALL ReadStateFile(InputFile,nGlobalElems)
  !Perform time averaging
  UAvg       =UAvg       +U_HDF5*AvgTime
  !IF(nVarVMS_HDF5.GT.0)&
  !  UAvg_VMS =UAvg_VMS   +VMSData_HDF5*AvgTime

  TotalAvgTime=TotalAvgTime+AvgTime
  CALL FinalizeState()
  lastInputFile=Inputfile ! we need the last input file for variable names, time etc in the write timeavg to hdf5 routine

  iCoarse=iCoarse+1
  IF((iCoarse.EQ.coarsenFac).OR.(iFile.EQ.nFiles).OR.(ABS(Time-AvgEndTime).LT.dt.AND.Time+dt.GT.AvgEndTime))THEN
    ! Compute total average and write file
    UAvg=UAvg/TotalAvgTime
    !IF(nVarVMS_HDF5.GT.0)THEN
    !  UAvg_VMS   =UAvg_VMS/TotalAvgTime
    !END IF
    FileName_Out=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileType),Time))//'_Merged.h5'
    CALL InitState(lastInputFile)
    CALL WriteTimeAverage(MeshFile,t,TotalAvgTime,FV_Elems_loc,FileType,&
                          (/nVarAvg ,PP_N+1,PP_N+1,PP_NZ+1/),VarNamesAvg,UAvg,FileName_In=FileName_Out)
    CALL FinalizeState()
    iCoarse=0
    UAvg=0.
    !IF(nVarVMS_HDF5.GT.0)THEN
    !  UAvg_VMS=0.
    !END IF
    TotalAvgTimeGlobal=TotalAvgTimeGlobal+TotalAvgTime
    TotalAvgTime=0.
  END IF
END DO
 
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,I5,A,I5,A,F10.5)') "Merging DONE: ",nFiles-nSkipped," of ",nFiles, &
                               " files merged over total averaging timespan ",TotalAvgTimeGlobal
SWRITE(UNIT_stdOut,'(132("="))')

SDEALLOCATE(UAvg)
SDEALLOCATE(UAvg_VMS)

END PROGRAM TimeAvg
