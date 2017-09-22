#include "flexi.h"

!===================================================================================================================================
!> This tool will take pre-averaged files (TimeAvg, Flucs) or simple state files 
!> and perform global temporal averaging
!===================================================================================================================================
PROGRAM TimeAvg
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_StringTools,ONLY:STRICMP, GetFileExtension
USE MOD_IO_HDF5
USE MOD_HDF5_Input, ONLY:ReadAttribute,ReadArray,GetDataSize
USE MOD_MPI,        ONLY:InitMPI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: AvgTime,TotalAvgTime,TotalAvgTimeGlobal
INTEGER                              :: iArg
CHARACTER(LEN=255)                   :: InputFile,LastInputFile 
CHARACTER(LEN=255)                   :: DataSet=''
CHARACTER(LEN=255)                   :: Filename_out,tmp,arg
REAL,ALLOCATABLE                     :: UAvg(:),Uloc(:)
LOGICAL                              :: isTimeAvg
REAL                                 :: AvgStarttime,Time,AvgEndTime
INTEGER                              :: StartArgs,nFiles,iFile
INTEGER                              :: coarsenFac,iCoarse
REAL                                 :: tol=1.E-7
REAL                                 :: dt
INTEGER                              :: startind,endind,offset
INTEGER                              :: locsize(16),i,n
INTEGER                              :: nSkipped

TYPE tFileSet
  INTEGER                            :: nDataSets,totalsize
  INTEGER,ALLOCATABLE                :: nDims(:),nVal(:,:)
  CHARACTER(LEN=255),ALLOCATABLE     :: DatasetNames(:)
  CHARACTER(LEN=255)                 :: FileType,MeshFile,NodeType,ProjectName
  REAL                               :: time
END TYPE
TYPE(tFileSet)                       :: ref,loc
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF (nArgs.LT.2) THEN
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
  arg = Args(iArg)

  ! check if the -start= flag is used, read requires dummystring
  IF (STRICMP(arg(1:8), "--start=")) THEN
    StartArgs=StartArgs+1
    tmp=TRIM(arg(9:LEN(arg)))
    READ(tmp,*) AvgStarttime
    WRITE(UNIT_stdOut,'(132("="))')
    WRITE(UNIT_stdOut,'(A28,F16.6)') ' Start time for Averaging is ',AvgStarttime
  ELSEIF (STRICMP(arg(1:6), "--end=")) THEN
    StartArgs=StartArgs+1
    tmp=TRIM(arg(7:LEN(arg)))
    READ(tmp,*) AvgEndtime
    WRITE(UNIT_stdOut,'(132("="))')
    WRITE(UNIT_stdOut,'(A28,F16.6)') ' End time for Averaging is ',AvgEndtime
  ELSEIF (STRICMP(arg(1:10), "--coarsen=")) THEN
    StartArgs=StartArgs+1
    tmp=TRIM(arg(11:LEN(arg)))
    READ(tmp,*) coarsenFac
    WRITE(UNIT_stdOut,'(132("="))')
    WRITE(UNIT_stdOut,'(A28,F16.6)') ' Coarsening factor is ',coarsenFac
  END IF

END DO

! check if at least 2 timeavg files are there
nFiles=nArgs-StartArgs+1
IF(nFiles.LT.2) CALL CollectiveStop(__STAMP__,'At least two files required for averaging!')

InputFile=Args(StartArgs+1)
CALL GetParams(InputFile,ref)

! Check whether the file is a normal state file or already contains time averaged data 
SELECT CASE(TRIM(ref%FileType))
CASE('State')
  isTimeAvg=.FALSE.
CASE('TimeAvg')
  isTimeAvg=.TRUE.
CASE('Fluc')
  isTimeAvg=.TRUE.  ! since we have time-averaged correlations
CASE DEFAULT
  CALL Abort(__STAMP__,'Unknown file type.')
END SELECT

ALLOCATE(UAvg(ref%totalsize))
ALLOCATE(Uloc(ref%totalsize))

! Start the averaging
WRITE(UNIT_stdOut,'(132("="))')
TotalAvgTime=0.
TotalAvgTimeGlobal=0.
iCoarse=0
Time=ref%time
nSkipped=0
DO iFile=1,nFiles
  InputFile=Args(iFile+StartArgs-1)
  ! check local time 
  dt=-Time

  CALL GetParams(InputFile,loc)
  ! Check if input data size has changed
  IF(TRIM(ref%FileType).NE.TRIM(loc%FileType))&
    CALL Abort(__STAMP__,'Mixing file types not allowed: '//TRIM(ref%FileType)//' '//TRIM(loc%FileType))
  IF(ANY(ref%nVal.NE.loc%nVal))&
    CALL Abort(__STAMP__,'Change of polynomial degree and variables not supported!')
  IF(TRIM(ref%NodeType).NE.TRIM(loc%NodeType))&
    CALL Abort(__STAMP__,'Change of node type not supported yet!')
  IF(TRIM(ref%MeshFile).NE.TRIM(loc%MeshFile))&
    CALL Abort(__STAMP__,'Change of meshfile not supported yet!')

  ! TODO : time umbasteln
  Time=loc%time
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

  ! Read AvgTime and Data
  avgTime=1.
  CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  IF(isTimeAvg)THEN
    CALL ReadAttribute(File_ID,'AvgTime',1,TRIM(DataSet),RealScalar=avgTime)
    WRITE(UNIT_stdOut,'(A,F10.5)') ' Time averaged file, averaging time is: ',avgTime 
  ELSE
    WRITE(UNIT_stdOut,'(A,A)') ' Normal state file, each file will be weighted identically.'
  END IF

  offset=0
  DO i=1,ref%nDataSets
    n=ref%nDims(i)
    locsize(1:n)=ref%nVal(1:n,i)
    startind=offset+1
    endind  =offset+PRODUCT(locsize(1:n))
    CALL ReadArray(TRIM(ref%DatasetNames(i)),n,locsize(1:n),0,n,&
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
    FileName_Out=TRIM(TIMESTAMP(TRIM(ref%ProjectName)//'_'//TRIM(ref%FileType),Time))//'_Merged.h5'
    CALL WriteCopyTimeAverage(InputFile,FileName_Out,ref,UAvg,TotalAvgTime)
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

CONTAINS

SUBROUTINE GetParams(filename,f)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename
TYPE(tFileSet),INTENT(OUT)  :: f
!===================================================================================================================================
CALL OpenDataFile(filename,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDatasetNamesInGroup("/",f%DatasetNames)
f%nDataSets=SIZE(f%DatasetNames)
ALLOCATE(f%nDims(f%nDataSets))
ALLOCATE(f%nVal(16,f%nDataSets))
f%nVal=0
f%totalsize=0
DO i=1,f%nDataSets
  CALL GetDataSize(File_ID,TRIM(f%DatasetNames(i)),f%nDims(i),HSize)
  ! TODO: use checksafeint
  f%nVal(1:f%nDims(i),i)=INT(HSize)
  DEALLOCATE(HSize)
  f%totalsize=f%totalsize+PRODUCT(f%nVal(1:f%nDims(i),i))
END DO
! Get default parameters
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar =f%FileType)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar =f%MeshFile)
CALL ReadAttribute(File_ID,'NodeType',1,StrScalar =f%NodeType)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =f%ProjectName)
CALL ReadAttribute(File_ID,'Time'    ,1,RealScalar=f%Time)
CALL CloseDataFile()
END SUBROUTINE


SUBROUTINE WriteCopyTimeAverage(filename_in,filename_out,f,uavg,avgTime)
! MODULES
USE MOD_HDF5_Output,ONLY:WriteArray,WriteAttribute
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename_in,filename_out
TYPE(tFileSet),INTENT(IN)   :: f
REAL,INTENT(IN)             :: UAvg(f%totalsize)
REAL,INTENT(IN)             :: avgTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: offset2(16) = 0
!===================================================================================================================================
CALL EXECUTE_COMMAND_LINE("cp -f "//TRIM(filename_in)//" "//TRIM(filename_out))

CALL OpenDataFile(filename_out,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
offset=0
DO i=1,f%nDataSets
  n=f%nDims(i)
  locsize(1:n)=f%nVal(1:n,i)
  startind=offset+1
  endind  =offset+PRODUCT(locsize(1:n))
  CALL WriteArray(TRIM(f%DataSetNames(i)),n,locsize(1:n),locsize(1:n),offset2(1:n),&
                          collective=.TRUE.,&
                          RealArray=UAvg(startInd:endInd))
  offset=endInd
END DO
CALL WriteAttribute(File_ID,'AvgTime',1,RealScalar=avgTime)
CALL CloseDataFile()
END SUBROUTINE WriteCopyTimeAverage


END PROGRAM TimeAvg
