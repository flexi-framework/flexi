!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!===================================================================================================================================
!> This tool will take pre-averaged files (TimeAvg, Flucs) or simple state files
!> and perform global temporal averaging
!===================================================================================================================================
PROGRAM TimeAvg
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_PreProc              ,ONLY: PP_N
USE MOD_AnalyzeEquation_Vars ,ONLY: UAvg,UFluc
USE MOD_HDF5_Input           ,ONLY: ReadAttribute,ReadArray,GetDataSize
USE MOD_HDF5_Output          ,ONLY: WriteTimeAverage,WriteArray,WriteAttribute
USE MOD_IO_HDF5              ,ONLY: HSIZE,File_ID,nDims
USE MOD_IO_HDF5              ,ONLY: OpenDataFile,CloseDataFile
USE MOD_IO_HDF5              ,ONLY: InitMPIInfo,GetDataSetNamesInGroup
USE MOD_MPI                  ,ONLY: InitMPI
USE MOD_Mesh_Vars            ,ONLY: offsetElem,nElems,nGlobalElems
USE MOD_Output_Vars          ,ONLY: ProjectName
USE MOD_StringTools          ,ONLY: STRICMP
#if USE_MPI
USE MOD_MPI                  ,ONLY: FinalizeMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE AND PARAMETER DEFINITIONS
REAL,PARAMETER                       :: tol   =1.E-7
INTEGER,PARAMETER                    :: maxDim=16            !< maximum number of permitted array dimension

TYPE tFileSet
  INTEGER                            :: nDataSets            !< number of datasets
  INTEGER,ALLOCATABLE                :: nDims(:)             !< number of dimensions per dataset
  INTEGER,ALLOCATABLE                :: nVal(:,:)            !< number of entries per dataset
  CHARACTER(LEN=255),ALLOCATABLE     :: DatasetNames(:)      !< names of the datasets
  CHARACTER(LEN=255),ALLOCATABLE     :: DatasetNamesAvg( :)
  CHARACTER(LEN=255),ALLOCATABLE     :: DatasetNamesFluc(:)
  CHARACTER(LEN=255)                 :: FileType
  CHARACTER(LEN=255)                 :: MeshFile
  CHARACTER(LEN=255)                 :: NodeType
  REAL                               :: time                 !< time
END TYPE
TYPE(tFileSet)                       :: ref  !< first file from command line used as reference
TYPE(tFileSet)                       :: loc  !< current file that is processed

TYPE tAvgData
  REAL,ALLOCATABLE                   :: AvgData(:)
  REAL,ALLOCATABLE                   :: AvgTmp( :)
END TYPE
TYPE(tAvgData),ALLOCATABLE           :: avg(:)

! LOCAL VARIABLES
!> Dataset
INTEGER                              :: iArg,i,nDim
CHARACTER(LEN=255)                   :: tmp,arg
CHARACTER(LEN=255),ALLOCATABLE       :: tmpDatasetNames(:)
INTEGER                              :: nVarAvg(5),nVarFluc(5)
INTEGER                              :: locsize(maxDim),globsize(maxDim),offsetsize(maxDim)

!> Time
REAL                                 :: dt,AvgTime,TotalAvgTime,TotalAvgTimeGlobal

!> FV Elems
INTEGER,ALLOCATABLE                  :: FV_Elems_loc(:)

!> Input
CHARACTER(LEN=255)                   :: InputFile
LOGICAL                              :: isTimeAvg
INTEGER                              :: StartArgs,nFiles,iFile
INTEGER                              :: coarsenFac,iCoarse,nSkipped
REAL                                 :: AvgStarttime,Time,TimeStart,AvgEndTime
CHARACTER(LEN=255),ALLOCATABLE       :: varnames_elems(:)
CHARACTER(LEN=255),ALLOCATABLE       :: varnames_field(:)
LOGICAL                              :: VarNamesElemsExist,VarNamesFieldExist

!> Output
CHARACTER(LEN=255)                   :: FilenameOut,FileTypeOut
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL InitMPIInfo()
CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF (nArgs.LT.2) THEN
  ! Print out error message containing valid syntax
  SWRITE(UNIT_stdOut,*) "Please use: timeavg --start=<starttime> --end=<endtime> --coarsen=<factor> FILE1 FILE2 .. FILEN"
  CALL CollectiveStop(__STAMP__,"At least two files are required for merging.")
END IF

! First check if the first three arguments contain -start=, -end= or -coarsen=
StartArgs    = 1
AvgStarttime = -HUGE(1.)
AvgEndtime   =  HUGE(1.)
coarsenFac   =  HUGE(1)

DO iArg = 1,MIN(nArgs,3)
  arg = Args(iArg)

  ! check if the -start= flag is used, read requires dummystring
  IF (STRICMP(arg(1:8), "--start=")) THEN
    StartArgs = StartArgs+1
    tmp=TRIM(arg(9:LEN(arg)))
    READ(tmp,*) AvgStarttime
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A28,F16.6)') ' Start time for Averaging is ',AvgStarttime
  ELSEIF (STRICMP(arg(1:6), "--end=")) THEN
    StartArgs = StartArgs+1
    tmp=TRIM(arg(7:LEN(arg)))
    READ(tmp,*) AvgEndtime
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A28,F16.6)') ' End time for Averaging is ',AvgEndtime
  ELSEIF (STRICMP(arg(1:10), "--coarsen=")) THEN
    StartArgs = StartArgs+1
    tmp=TRIM(arg(11:LEN(arg)))
    READ(tmp,*) coarsenFac
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A28,F16.6)') ' Coarsening factor is ',coarsenFac
  END IF
END DO

! check if at least 2 timeavg files are there
nFiles = nArgs-StartArgs+1
IF (nFiles.LT.2) CALL CollectiveStop(__STAMP__,'At least two files required for averaging!')

InputFile = Args(StartArgs+1)
! Partitioning - we partition along the last dimension of the arrays and assume that this
! dimension has the same extend for all data sets
CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDatasetNamesInGroup("/",tmpDatasetNames)
CALL GetDataSize(File_ID,TRIM(tmpDatasetNames(1)),nDims,HSize)
nGlobalElems = INT(HSIZE(nDims))
DEALLOCATE(HSize)

#if USE_MPI
nElems     = INT(nGlobalElems/nProcessors)
offsetElem = myRank*nElems
IF (myRank.EQ.nProcessors-1) THEN
  ! This is the last proc, may have other number of elements
  nElems = nElems + MOD(nGlobalElems,nProcessors)
END IF
#else
offsetElem = 0
nElems     = nGlobalElems
#endif

CALL CloseDataFile()
DEALLOCATE(tmpDatasetNames)
! Read in reference values for data size etc.
CALL GetParams(InputFile,ref)

! Check whether the file is a normal state file or already contains time averaged data
SELECT CASE(TRIM(ref%FileType))
  CASE('State')
    isTimeAvg   = .FALSE.
    FileTypeOut = 'TimeAvg'
  CASE('TimeAvg')
    isTimeAvg   = .TRUE.
    FileTypeOut = 'TimeAvg'
  CASE('Fluc')
    isTimeAvg   = .TRUE.  ! since we have time-averaged correlations
    FileTypeOut = 'Fluc'
  CASE DEFAULT
    isTimeAvg   = .TRUE.  ! remove compiler warning
END SELECT

! Allocate arrays for the datasets
ALLOCATE(avg(ref%nDataSets))

DO i = 1,ref%nDataSets
  nDim            = ref%nDims(   i)
  locsize(1:nDim) = ref%nVal(1:nDim,i)
  locsize(  nDim) = nElems

  ALLOCATE(avg(i)%AvgData(PRODUCT(ref%nVal(1:nDim,i))) &
          ,avg(i)%AvgTmp( PRODUCT(ref%nVal(1:nDim,i))))

  avg(i)%AvgData = 0
  avg(i)%AvgTmp  = 0

  ! Identify the datasets
  IF (isTimeAvg) THEN
    SELECT CASE(TRIM(ref%DatasetNames(i)))
      CASE('Mean')
        ALLOCATE(UAvg(locsize(1),locsize(2),locsize(3),locsize(4),locsize(5)))
        PP_N = locsize(2)-1

      CASE('MeanSquare')
        ALLOCATE(UFluc(locsize(1),locsize(2),locsize(3),locsize(4),locsize(5)))
        PP_N = locsize(2)-1
    END SELECT
  ELSE
    SELECT CASE(TRIM(ref%DatasetNames(i)))
      CASE('DG_Solution')
        ALLOCATE(UAvg(locsize(1),locsize(2),locsize(3),locsize(4),locsize(5)))
        PP_N = locsize(2)-1
    END SELECT
  END IF
END DO

! Start the averaging
SWRITE(UNIT_stdOut,'(132("="))')
TotalAvgTime       = 0.
TotalAvgTimeGlobal = 0.
iCoarse            = 0
Time               = ref%time
nSkipped           = 0
timestart          = -999.

DO iFile = 1,nFiles
  InputFile = Args(iFile+StartArgs-1)
  ! check local time
  dt = -Time

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

  Time = loc%time
  dt   = dt + Time
  ! If the time in the current file is below the starttime, then cycle
  IF ((Time-AvgStartTime.LT.-tol).OR.(Time-AvgEndTime.GT.tol)) THEN
    SWRITE(UNIT_stdOut,'(A,A)') ' SKIPPING FILE ',TRIM(InputFile)
    nSkipped = nSkipped + 1

    IF (nFiles-nSkipped.LE.1) THEN
      SWRITE(UNIT_stdOut,'(A)') ' WARNING: All files have been skipped, no output is performed.'
      CALL CollectiveStop(__STAMP__,'Please check start and end time.')
    END IF
    CYCLE
  END IF

  IF (iCoarse.EQ.0) TimeStart = Time

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iFile,' of ',nFiles,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)')       ' ( "',TRIM(InputFile),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')

  ! Read AvgTime and Data
  avgTime = 1.
  CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

  IF (isTimeAvg) THEN
    CALL ReadAttribute(File_ID,'AvgTime',1,RealScalar=avgTime)
    SWRITE(UNIT_stdOut,'(A,F10.5)') ' Time averaged file, averaging time is: ',avgTime
  ELSE
    SWRITE(UNIT_stdOut,'(A,A)')     ' Normal state file, each file will be weighted identically.'
  END IF

  DO i = 1,ref%nDataSets
    nDim            = ref%nDims(i)
    locsize(1:nDim) = ref%nVal(1:nDim,i)
    locsize(nDim)   = nElems

    CALL ReadArray(ArrayName  = TRIM(ref%DatasetNames(i)) &
                  ,Rank       = nDim                      &
                  ,nVal       = locsize(1:nDim)           &
                  ,Offset_in  = offsetElem                &
                  ,Offset_dim = nDim                      &
                  ,RealArray  = avg(i)%AvgTmp)

    avg(i)%AvgData = avg(i)%AvgData + AvgTime*avg(i)%AvgTmp
  END DO

  CALL CloseDataFile()

  ! Perform time averaging
  TotalAvgTime   = TotalAvgTime   + AvgTime
  iCoarse        = iCoarse        + 1

  ! Write the output
  IF ((iCoarse.EQ.coarsenFac).OR.(iFile.EQ.nFiles).OR.(ABS(Time-AvgEndTime).LT.dt.AND.Time+dt.GT.AvgEndTime)) THEN
    FileNameOut = TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileTypeOut)//'_Merged',Time,TimeStart))//'.h5'
    nVarAvg        = 0
    nVarFluc       = 0

    ! Compute total average
    DO i = 1,ref%nDataSets
      avg(i)%AvgData = avg(i)%AvgData/TotalAvgTime

      ! Identify the datasets
      IF (isTimeAvg) THEN
        SELECT CASE(TRIM(ref%DatasetNames(i)))
          CASE('Mean')
            nDim     = ref%nDims(   i)
            nVarAvg  = ref%nVal(1:nDim,i)
            UAvg     = RESHAPE(avg(i)%AvgData,(/nVarAvg( 1),nVarAvg( 2),nVarAvg( 3),nVarAvg( 4),nVarAvg( 5)/))

          CASE('MeanSquare')
            nDim     = ref%nDims(   i)
            nVarFluc = ref%nVal(1:nDim,i)
            UFluc    = RESHAPE(avg(i)%AvgData,(/nVarFluc(1),nVarFluc(2),nVarFluc(3),nVarFluc(4),nVarFluc(5)/))
        END SELECT
      ELSE
        SELECT CASE(TRIM(ref%DatasetNames(i)))
          CASE('DG_Solution')
            nDim     = ref%nDims(   i)
            nVarAvg  = ref%nVal(1:nDim,i)
            UAvg     = RESHAPE(avg(i)%AvgData,(/nVarAvg( 1),nVarAvg( 2),nVarAvg( 3),nVarAvg( 4),nVarAvg( 5)/))

            CALL OpenDataFile(FileNameOut,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
            CALL WriteAttribute(File_ID,'VarNames_Mean',nVarAvg(1),StrArray=ref%DatasetNamesAvg)
            CALL CloseDataFile
        END SELECT
      END IF
    END DO

    ! Allocate dummy datasets
    ALLOCATE(FV_Elems_loc(nElems))
    FV_Elems_loc = FV_ENABLED

    ! Write file
    CALL WriteTimeAverage(MeshFileName = ref%MeshFile         &
                         ,OutputTime   = Time                 &
                         ,dtAvg        = TotalAvgTime         &
                         ,nVal         = nVarAvg(2:4)         &
                         ,nVarAvg      = nVarAvg(1)           &
                         ,VarNamesAvg  = ref%DatasetNamesAvg  &
                         ,UAvg         = UAvg                 &
                         ,nVarFluc     = nVarFluc(1)          &
                         ,VarNamesFluc = ref%DatasetNamesFluc &
                         ,UFluc        = UFluc                &
                         ,FileName_In  = FileNameOut          &
                         ,NodeType_In  = ref%NodeType)

    ! Write the remaining attributes
    CALL OpenDataFile(FileNameOut,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
    IF (VarNamesElemsExist) CALL WriteAttribute(File_ID,'VarNamesAdd'     ,SIZE(varnames_elems),StrArray=varnames_elems)
    IF (VarNamesFieldExist) CALL WriteAttribute(File_ID,'VarNamesAddField',SIZE(varnames_field),StrArray=varnames_field)

    ! Write the remaining datasets
    DO i = 1,ref%nDataSets
      SELECT CASE(TRIM(ref%DatasetNames(i)))
        CASE('DG_Solution','Mean','MeanSquare')
          ! Do nothing, already written above
        CASE DEFAULT
          nDim             = ref%nDims(   i)
          locsize( 1:nDim) = ref%nVal(1:nDim,i)
          locsize(   nDim) = nElems
          globsize(1:nDim) = ref%nVal(1:nDim,i)
          globsize(  nDim) = nGlobalElems
          offsetsize       = 0
          offsetsize(nDim) = offsetElem

          CALL WriteArray(DataSetName  = TRIM(ref%DataSetNames(i)) &
                         ,rank         = nDim                      &
                         ,nValGlobal   = globsize(  1:nDim)        &
                         ,nVal         = locsize(   1:nDim)        &
                         ,offset       = offsetsize(1:nDim)        &
                         ,collective   = .TRUE.                    &
                         ,RealArray    = avg(i)%AvgData)
      END SELECT
    END DO
    CALL CloseDataFile()

    iCoarse            = 0
    TotalAvgTimeGlobal = TotalAvgTimeGlobal+TotalAvgTime
    TotalAvgTime       = 0.
  END IF
END DO

SDEALLOCATE(UAvg)
SDEALLOCATE(UFluc)
SDEALLOCATE(FV_Elems_loc)

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,I5,A,I5,A,F10.5)') "Merging DONE: ",nFiles-nSkipped," of ",nFiles, &
                               " files merged over total averaging timespan ",TotalAvgTimeGlobal
SWRITE(UNIT_stdOut,'(132("="))')

#if USE_MPI
CALL FinalizeMPI()
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

CONTAINS


!===================================================================================================================================
!> Retrieves relevant header and dateset parameters from Flexi files and stores them in a type
!===================================================================================================================================
SUBROUTINE GetParams(filename,f)
! MODULES
USE MOD_Output_Vars          ,ONLY: ProjectName
USE MOD_HDF5_Input           ,ONLY: GetVarNames
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename !< input filename
TYPE(tFileSet),INTENT(OUT)  :: f        !< type with infos to be filled
!===================================================================================================================================
CALL OpenDataFile(filename,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDatasetNamesInGroup("/",f%DatasetNames)

f%nDataSets = SIZE(f%DatasetNames)
ALLOCATE(f%nDims(      f%nDataSets))
ALLOCATE(f%nVal(maxDim,f%nDataSets))
f%nVal = 0

DO i = 1,f%nDataSets
  CALL GetDataSize(File_ID,TRIM(f%DatasetNames(i)),f%nDims(i),HSize)
  CHECKSAFEINT(MAXVAL(HSize),4)
  CHECKSAFEINT(MINVAL(HSize),4)
  f%nVal(1:f%nDims(i),i) = INT(HSize)
  DEALLOCATE(HSize)
END DO

! Get default parameters
CALL ReadAttribute(File_ID,'File_Type',   1,StrScalar =f%FileType)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =f%MeshFile)
CALL ReadAttribute(File_ID,'NodeType',    1,StrScalar =f%NodeType)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
CALL ReadAttribute(File_ID,'Time'        ,1,RealScalar=f%Time)

! Get additional varnames
CALL GetVarNames('VarNamesAdd'     ,varnames_elems,VarNamesElemsExist)
CALL GetVarNames('VarNamesAddField',varnames_field,VarNamesFieldExist)

! Get the VarNames for Mean and MeanSquare
SELECT CASE(f%fileType)
  CASE('TimeAvg','Fluc')
    DO i = 1,f%nDataSets
      SELECT CASE(TRIM(f%DatasetNames(i)))
        CASE('Mean')
          IF (ALLOCATED(f%DatasetNamesAvg )) CYCLE

          ALLOCATE(f%DatasetNamesAvg( f%nVal(1,i)))
          CALL ReadAttribute(File_ID,'VarNames_Mean'      ,f%nVal(1,i),StrArray=f%DatasetNamesAvg)

        CASE('MeanSquare')
          IF (ALLOCATED(f%DatasetNamesFluc)) CYCLE

          ALLOCATE(f%DatasetNamesFluc(f%nVal(1,i)))
          CALL ReadAttribute(File_ID,'VarNames_MeanSquare',f%nVal(1,i),StrArray=f%DatasetNamesFluc)
      END SELECT
    END DO

  CASE('State')
    DO i = 1,f%nDataSets
      IF (STRICMP(f%DatasetNames(i),'DG_Solution')) THEN
          IF (ALLOCATED(f%DatasetNamesAvg )) CYCLE

          ALLOCATE(f%DatasetNamesAvg( f%nVal(1,i)))
          CALL ReadAttribute(File_ID,'VarNames'           ,f%nVal(1,i),StrArray=f%DatasetNamesAvg)
      END IF
    END DO
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'Unknown file type: '//TRIM(ref%FileType))
END SELECT

CALL CloseDataFile()

END SUBROUTINE

END PROGRAM TimeAvg
