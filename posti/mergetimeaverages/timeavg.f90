#include "flexi.h"

!===================================================================================================================================
!> This tool will take pre-averaged files (TimeAvg, Flucs) or simple state files
!> and perform global temporal averaging
!===================================================================================================================================
PROGRAM TimeAvg
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_StringTools,ONLY:STRICMP
USE MOD_IO_HDF5
USE MOD_HDF5_Input, ONLY:ReadAttribute,ReadArray,GetDataSize
USE MOD_MPI,        ONLY:InitMPI
#if USE_MPI
USE MOD_MPI,        ONLY:FinalizeMPI
#endif
USE MOD_IO_HDF5,    ONLY: InitMPIInfo
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE AND PARAMETER DEFINITIONS
REAL,PARAMETER                       :: tol=1.E-7
INTEGER,PARAMETER                    :: maxDim=16  !< maximum number of permitted array dimension
TYPE tFileSet
  INTEGER                            :: nDataSets  !< number of datasets
  INTEGER(KIND=8)                    :: totalsize  !< 1D size of all arrays
  INTEGER,ALLOCATABLE                :: nDims(:)   !< number of dimensions per dataset
  INTEGER,ALLOCATABLE                :: nVal(:,:)  !< number of entries per dataset
  CHARACTER(LEN=255),ALLOCATABLE     :: DatasetNames(:) !< names of the datasets
  CHARACTER(LEN=255)                 :: FileType,MeshFile,NodeType,ProjectName
  REAL                               :: time       !< time
END TYPE
TYPE(tFileSet)                       :: ref,loc

! LOCAL VARIABLES
REAL                                 :: AvgTime,TotalAvgTime,TotalAvgTimeGlobal
INTEGER                              :: iArg
CHARACTER(LEN=255)                   :: InputFile,LastInputFile
CHARACTER(LEN=255)                   :: DataSet=''
CHARACTER(LEN=255)                   :: FilenameOut,FileTypeOut,tmp,arg
REAL,ALLOCATABLE                     :: UAvg(:),Uloc(:)
LOGICAL                              :: isTimeAvg
REAL                                 :: AvgStarttime,Time,TimeStart,AvgEndTime
INTEGER                              :: StartArgs,nFiles,iFile
INTEGER                              :: coarsenFac,iCoarse
REAL                                 :: dt
INTEGER(KIND=8)                      :: startind,endind,offset
INTEGER                              :: locsize(16),i,n
INTEGER                              :: nSkipped
INTEGER                              :: offsetMPI,nElems
INTEGER                              :: nElemsGlob
CHARACTER(LEN=255),ALLOCATABLE       :: tmpDatasetNames(:)
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL InitMPIInfo()
CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF (nArgs.LT.2) THEN
  ! Print out error message containing valid syntax
  SWRITE(UNIT_stdOut,*) "Please use: timeavg --start=<starttime> --end=<endtime> --coarsen=<factor> FILE1 FILE2 .. FILEN"
  SWRITE(UNIT_stdOut,*) "At least two files are required for merging."
  STOP
END IF

! First check if the first three arguments contain -start=, -end= or -coarsen=
StartArgs=1
AvgStarttime=-HUGE(1.)
AvgEndtime=HUGE(1.)
coarsenFac=999999999
DO iArg=1,MIN(nArgs,3)
  arg = Args(iArg)

  ! check if the -start= flag is used, read requires dummystring
  IF (STRICMP(arg(1:8), "--start=")) THEN
    StartArgs=StartArgs+1
    tmp=TRIM(arg(9:LEN(arg)))
    READ(tmp,*) AvgStarttime
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A28,F16.6)') ' Start time for Averaging is ',AvgStarttime
  ELSEIF (STRICMP(arg(1:6), "--end=")) THEN
    StartArgs=StartArgs+1
    tmp=TRIM(arg(7:LEN(arg)))
    READ(tmp,*) AvgEndtime
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A28,F16.6)') ' End time for Averaging is ',AvgEndtime
  ELSEIF (STRICMP(arg(1:10), "--coarsen=")) THEN
    StartArgs=StartArgs+1
    tmp=TRIM(arg(11:LEN(arg)))
    READ(tmp,*) coarsenFac
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A28,F16.6)') ' Coarsening factor is ',coarsenFac
  END IF

END DO

! check if at least 2 timeavg files are there
nFiles=nArgs-StartArgs+1
IF(nFiles.LT.2) CALL CollectiveStop(__STAMP__,'At least two files required for averaging!')

InputFile=Args(StartArgs+1)
! Partitioning - we partition along the last dimension of the arrays and assume that this
! dimension has the same extend for all data sets
CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDatasetNamesInGroup("/",tmpDatasetNames)
CALL GetDataSize(File_ID,TRIM(tmpDatasetNames(1)),nDims,HSize)
nElemsGlob = INT(HSIZE(nDims))
DEALLOCATE(HSize)
#if USE_MPI
nElems = INT(nElemsGlob/nProcessors)
offsetMPI = myRank*nElems
IF (myRank.EQ.nProcessors-1) THEN
  ! This is the last proc, may have other number of elements
  nElems = nElems + MOD(nElemsGlob,nProcessors)
END IF
#else
offsetMPI = 0
nElems = nElemsGlob
#endif
CALL CloseDataFile()
DEALLOCATE(tmpDatasetNames)
! Read in reference values for data size etc.
CALL GetParams(InputFile,ref)

! Check whether the file is a normal state file or already contains time averaged data
SELECT CASE(TRIM(ref%FileType))
CASE('State')
  isTimeAvg=.FALSE.
  FileTypeOut='TimeAvg'
CASE('TimeAvg')
  isTimeAvg=.TRUE.
  FileTypeOut='TimeAvg'
CASE('Fluc')
  isTimeAvg=.TRUE.  ! since we have time-averaged correlations
  FileTypeOut='Fluc'
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Unknown file type: '//TRIM(ref%FileType))
END SELECT

ALLOCATE(UAvg(ref%totalsize))
ALLOCATE(Uloc(ref%totalsize))

! Start the averaging
SWRITE(UNIT_stdOut,'(132("="))')
TotalAvgTime=0.
TotalAvgTimeGlobal=0.
iCoarse=0
Time=ref%time
nSkipped=0
timestart=-999.
DO iFile=1,nFiles
  InputFile=Args(iFile+StartArgs-1)
  ! check local time
  dt=-Time

  CALL GetParams(InputFile,loc)
  ! Check if input data size has changed
  IF(.NOT.STRICMP(ref%FileType,loc%FileType))&
    CALL CollectiveStop(__STAMP__,'Mixing file types is not allowed: '//TRIM(ref%FileType)//' '//TRIM(loc%FileType))
  IF(.NOT.STRICMP(ref%NodeType,loc%NodeType))&
    CALL CollectiveStop(__STAMP__,'Change of node type not supported yet!')
  IF(.NOT.STRICMP(ref%MeshFile,loc%MeshFile))&
    CALL CollectiveStop(__STAMP__,'Change of mesh file not supported yet!')
  IF(ANY(ref%nVal.NE.loc%nVal))&
    CALL CollectiveStop(__STAMP__,'Change of polynomial degree and variables not supported!')
  ! TODO: check change of FV subcells ?!

  Time=loc%time
  dt=dt+Time
  ! If the time in the current file is below the starttime, then cycle
  IF ((Time-AvgStartTime.LT.-tol).OR.(Time-AvgEndTime.GT.tol)) THEN
    SWRITE(UNIT_stdOut,'(A,A)') ' SKIPPING FILE ',TRIM(InputFile)
    nSkipped=nSkipped+1
    IF (nFiles-nSkipped.LE.1) THEN
      SWRITE(UNIT_stdOut,*) "WARNING: All files have been skipped, no output is performed."
      SWRITE(UNIT_stdOut,*) "         Please check start and end time."
      STOP
    END IF
    CYCLE
  END IF

  IF(iCoarse.EQ.0) TimeStart=Time

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iFile,' of ',nFiles,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputFile),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')

  ! Read AvgTime and Data
  avgTime=1.
  CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  IF(isTimeAvg)THEN
    CALL ReadAttribute(File_ID,'AvgTime',1,TRIM(DataSet),RealScalar=avgTime)
    SWRITE(UNIT_stdOut,'(A,F10.5)') ' Time averaged file, averaging time is: ',avgTime
  ELSE
    SWRITE(UNIT_stdOut,'(A,A)')     ' Normal state file, each file will be weighted identically.'
  END IF

  offset=0
  DO i=1,ref%nDataSets
    n=ref%nDims(i)
    locsize(1:n)=ref%nVal(1:n,i)
    locsize(n) = nElems
    startind=offset+1
    endind  =offset+PRODUCT(locsize(1:n))
    CALL ReadArray(TRIM(ref%DatasetNames(i)),n,locsize(1:n),offsetMPI,n,&
                   RealArray=Uloc(startInd:endInd))
    offset=endInd
  END DO
  CALL CloseDataFile()

  !Perform time averaging
  UAvg         = UAvg         + AvgTime*Uloc
  TotalAvgTime = TotalAvgTime + AvgTime

  lastInputFile=Inputfile ! we need the last input file for variable names, time etc in the write timeavg to hdf5 routine

  iCoarse=iCoarse+1
  IF((iCoarse.EQ.coarsenFac).OR.(iFile.EQ.nFiles).OR.(ABS(Time-AvgEndTime).LT.dt.AND.Time+dt.GT.AvgEndTime))THEN
    ! Compute total average and write file
    UAvg=UAvg/TotalAvgTime
    FileNameOut=TRIM(TIMESTAMP(TRIM(ref%ProjectName)//'_'//TRIM(FileTypeOut)//'_Merged',Time,TimeStart))//'.h5'
    CALL WriteTimeAverageByCopy(InputFile,FileNameOut,FileTypeOut,ref,UAvg,TotalAvgTime)
    iCoarse=0
    UAvg=0.

    TotalAvgTimeGlobal=TotalAvgTimeGlobal+TotalAvgTime
    TotalAvgTime=0.
  END IF
END DO

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,I5,A,I5,A,F10.5)') "Merging DONE: ",nFiles-nSkipped," of ",nFiles, &
                               " files merged over total averaging timespan ",TotalAvgTimeGlobal
SWRITE(UNIT_stdOut,'(132("="))')

SDEALLOCATE(UAvg)
SDEALLOCATE(Uloc)

#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
CALL FinalizeMPI()
#endif

CONTAINS



!===================================================================================================================================
!> Retrieves relevant header and dateset parameters from Flexi files and stores them in a type
!===================================================================================================================================
SUBROUTINE GetParams(filename,f)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename !< input filename
TYPE(tFileSet),INTENT(OUT)  :: f        !< type with infos to be filled
!===================================================================================================================================
CALL OpenDataFile(filename,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDatasetNamesInGroup("/",f%DatasetNames)
f%nDataSets=SIZE(f%DatasetNames)
ALLOCATE(f%nDims(f%nDataSets))
ALLOCATE(f%nVal(maxDim,f%nDataSets))
f%nVal=0
f%totalsize=0
DO i=1,f%nDataSets
  CALL GetDataSize(File_ID,TRIM(f%DatasetNames(i)),f%nDims(i),HSize)
  CHECKSAFEINT(MAXVAL(HSize),4)
  CHECKSAFEINT(MINVAL(HSize),4)
  f%nVal(1:f%nDims(i),i)=INT(HSize)
  DEALLOCATE(HSize)
  f%totalsize=f%totalsize+PRODUCT(f%nVal(1:f%nDims(i)-1,i))*nElems
END DO
! Get default parameters
CALL ReadAttribute(File_ID,'File_Type',   1,StrScalar =f%FileType)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =f%MeshFile)
CALL ReadAttribute(File_ID,'NodeType',    1,StrScalar =f%NodeType)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =f%ProjectName)
CALL ReadAttribute(File_ID,'Time'        ,1,RealScalar=f%Time)
CALL CloseDataFile()
END SUBROUTINE


!===================================================================================================================================
!> Copies an existing state or timeavg file as a basis and write the merged data into it.
!> This ensures that all relevant information is contained in the merged file without too
!> much coding overhead.
!===================================================================================================================================
SUBROUTINE WriteTimeAverageByCopy(filename_in,filename_out,filetype_out,f,uavg,avgTime)
! MODULES
USE MOD_HDF5_Output,ONLY:WriteArray,WriteAttribute
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename_in  !< file to be copied
CHARACTER(LEN=*),INTENT(IN) :: filename_out !< output file
CHARACTER(LEN=*),INTENT(IN) :: filetype_out !< filetype of HDF5
TYPE(tFileSet),INTENT(IN)   :: f            !< type with header and dataset information
REAL,INTENT(IN)             :: UAvg(f%totalsize) !< merged data (1D array with data for all datasets)
REAL,INTENT(IN)             :: avgTime      !< averaged time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: globsize(1:maxDim)
INTEGER,ALLOCATABLE :: offset2(:)
CHARACTER(LEN=255)  :: tmp255
!===================================================================================================================================
IF (MPIRoot) CALL EXECUTE_COMMAND_LINE("cp -f "//TRIM(filename_in)//" "//TRIM(filename_out))
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

CALL OpenDataFile(filename_out,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
offset=0
DO i=1,f%nDataSets
  n=f%nDims(i)
  locsize(1:n)=f%nVal(1:n,i)
  globsize = locsize
  locsize(n) = nElems
  startind=offset+1
  endind  =offset+PRODUCT(locsize(1:n))
  ALLOCATE(offset2(n))
  offset2     = 0
  offset2(n)  = offsetMPI
  CALL WriteArray(TRIM(f%DataSetNames(i)),n,globsize(1:n),locsize(1:n),offset2,&
                          collective=.TRUE.,&
                          RealArray=UAvg(startInd:endInd))
  offset=endInd
  DEALLOCATE(offset2)
END DO
CALL WriteAttribute(File_ID,'AvgTime',1,RealScalar=avgTime)
tmp255=TRIM(filetype_out)
CALL WriteAttribute(File_ID,'File_Type',1,StrScalar=(/tmp255/))
CALL CloseDataFile()
END SUBROUTINE WriteTimeAverageByCopy


END PROGRAM TimeAvg
