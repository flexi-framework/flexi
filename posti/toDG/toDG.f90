#include "flexi.h"

!===================================================================================================================================
!> This tool will convert a TimeAvg file to full DG
!===================================================================================================================================
PROGRAM toDG
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_FV_Basis,           ONLY: FV_GetVandermonde
USE MOD_HDF5_Input,         ONLY: ReadAttribute,ReadArray,GetDataSize,DataSetExists,GetArrayAndName,GetDataProps
USE MOD_HDF5_Output,        ONLY: WriteArray,WriteAttribute
USE MOD_IO_HDF5
USE MOD_Mesh_Vars,          ONLY: nGlobalElems,nElems,OffsetElem
USE MOD_MPI,                ONLY: InitMPI
#if USE_MPI
USE MOD_MPI,                ONLY: FinalizeMPI
#endif
USE MOD_ReadInTools
USE MOD_StringTools,        ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE AND PARAMETER DEFINITIONS
INTEGER                        :: iArg,iElem,iVar
CHARACTER(LEN=255)             :: InputFile
INTEGER                        :: nVarMean,nVarMeanSquare,nVarElemData,nVar
INTEGER                        :: iD,FVIndex
INTEGER                        :: nValMean(15),nValMeanSquare(15),nValAdd(15),nVal(15)
INTEGER                        :: nValMean_glob(15),nValMeanSquare_glob(15),nValAdd_glob(15),nVal_glob(15)
REAL,ALLOCATABLE               :: UMean(:,:,:,:,:),UMeanSquare(:,:,:,:,:),ElemData(:,:)
REAL,ALLOCATABLE               :: U(:,:,:,:,:)
REAL,ALLOCATABLE               :: UMeanDG(:,:,:,:,:),UMeanSquareDG(:,:,:,:,:),ElemDataDG(:,:)
REAL,ALLOCATABLE               :: UDG(:,:,:,:,:)
REAL,ALLOCATABLE               :: UMeanTmp(:),UMeanSquareTmp(:),ElemDataTmp(:)
REAL,ALLOCATABLE               :: UTmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesMean(:),VarNamesMeanSquare(:),VarNamesFluc(:),VarNamesAdd(:),VarNames(:)
INTEGER                        :: nVar_HDF5,N_HDF5,nElems_HDF5,N_HDF5Z
CHARACTER(LEN=255)             :: NodeType_HDF5,FileNameOut
REAL,ALLOCATABLE               :: FV_Vdm(:,:),FV_sVdm(:,:)
CHARACTER(LEN=255)             :: FileType
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL InitMPIInfo()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') '=============================================================   toDG   '//&
                          '============================================================='
SWRITE(UNIT_stdOut,'(132("="))')

CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF (nArgs.LT.1) THEN
  SWRITE(UNIT_stdOut,*) "Please use: posti_toDG <FILE1> <FILE2> ..."
  STOP
END IF

InputFile=Args(1)
! Partitioning - we partition along the last dimension of the arrays
CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5)
N_HDF5Z       = ZDIM(N_HDF5)
#if USE_MPI
nGlobalElems  = nElems_HDF5
nElems        = INT(nElems_HDF5/nProcessors)
OffsetElem    = myRank*nElems
IF (myRank.EQ.nProcessors-1) THEN
  ! This is the last proc, may have other number of elements
  nElems = nElems + MOD(nElems_HDF5,nProcessors)
END IF
#else
OffsetElem    = 0
nElems        = nElems_HDF5
nGlobalElems  = nElems_HDF5
#endif
CALL CloseDataFile()

ALLOCATE(FV_Vdm(0:N_HDF5,0:N_HDF5),FV_sVdm(0:N_HDF5,0:N_HDF5))
CALL FV_GetVandermonde(N_HDF5,NodeType_HDF5,FV_Vdm,FV_sVdm)


DO iArg=1,nArgs
  SWRITE(UNIT_stdOut,"(A,I0,A,I0,A)") "Processing File ", iArg, " of ", nArgs, "..."
  SWRITE(UNIT_stdOut,'(132("-"))')

  InputFile=Args(iArg)

  CALL Readin()

  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Convert Data to DG..."

  !First, Check for FV array in ElemData
  FVIndex = 0
  DO iD=1,nVarElemData
    !index mapping for FV_Elems
    IF(STRICMP(TRIM(VarNamesAdd(iD)),'FV_Elems')) FVIndex=iD
  END DO
  IF(FVIndex.EQ.0) CALL Abort(__STAMP__,"FVElems array not found in HDF5 file.")

  ! Switch between Statefile and TimeAvg mode
  SELECT CASE(TRIM(FileType))
  CASE('State')
    ALLOCATE(UDG(                    nVar,0:N_HDF5,0:N_HDF5,0:N_HDF5Z,nElems))
  CASE('TimeAvg')
    ALLOCATE(UMeanDG(            nVarMean,0:N_HDF5,0:N_HDF5,0:N_HDF5Z,nElems))
    ALLOCATE(UMeanSquareDG(nVarMeanSquare,0:N_HDF5,0:N_HDF5,0:N_HDF5Z,nElems))
  END SELECT
  ALLOCATE(ElemDataDG(nVarElemData,nElems))


  ! Main loop over all elements
  ElemDataDG = ElemData
  DO iElem=1,nElems
    IF(ElemData(FVIndex,iElem).EQ.0) THEN
      SELECT CASE(TRIM(FileType))
      CASE('State')
        UDG(:,:,:,:,iElem)            = U(          :,:,:,:,iElem)
      CASE('TimeAvg')
        UMeanDG(:,:,:,:,iElem)        = UMean(      :,:,:,:,iElem)
        UMeanSquareDG(:,:,:,:,iElem)  = UMeanSquare(:,:,:,:,iElem)
      END SELECT
      ElemDataDG(FVIndex,iElem)       = 0
    ELSE
      SELECT CASE(TRIM(FileType))
      CASE('State')
        DO iVar=1,nVar
          CALL ChangeBasisVolume(N_HDF5,N_HDF5Z,FV_sVdm,U(iVar,:,:,:,iElem)          ,UDG(iVar,:,:,:,iElem))
        END DO
      CASE('TimeAvg')
        DO iVar=1,nVarMean
          CALL ChangeBasisVolume(N_HDF5,N_HDF5Z,FV_sVdm,UMean(iVar,:,:,:,iElem)      ,UMeanDG(iVar,:,:,:,iElem))
        END DO
        DO iVar=1,nVarMeanSquare
          CALL ChangeBasisVolume(N_HDF5,N_HDF5Z,FV_sVdm,UMeanSquare(iVar,:,:,:,iElem),UMeanSquareDG(iVar,:,:,:,iElem))
        END DO
      END SELECT
      ElemDataDG(FVIndex,iElem) = 0
    END IF
  END DO

  SWRITE(UNIT_stdOut,'(A)') " DONE!"
  SWRITE(UNIT_stdOut,'(132("-"))')

  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Write data to HDF5 file..."

  SELECT CASE(TRIM(FileType))
  CASE('State')
    nVal_glob              = nVal
    nValAdd_glob           = nValAdd

    nVal_glob(5)           = nGlobalElems
    nValAdd_glob(2)        = nGlobalElems

    ! Output
    FileNameOut = Args(iArg)(:LEN(TRIM(Args(iArg)))-3)//'_DG.h5'
#ifdef FLANG
    IF (MPIRoot) CALL SYSTEM("cp -f "//TRIM(InputFile)//" "//TRIM(FileNameOut))
#else
    IF (MPIRoot) CALL EXECUTE_COMMAND_LINE("cp -f "//TRIM(InputFile)//" "//TRIM(FileNameOut))
#endif
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
    CALL OpenDataFile(FileNameOut,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
    CALL WriteArray('DG_Solution',5,nVal_glob(1:5)   ,nVal(1:5)   ,(/0,0,0,0,OffsetElem/),collective=.TRUE.,RealArray=UDG)
    CALL WriteArray('ElemData'   ,2,nValAdd_glob(1:2),nValAdd(1:2),(/0,OffsetElem/)      ,collective=.TRUE.,RealArray=ElemDataDG)
    CALL CloseDataFile()
  CASE('TimeAvg')
    nValMean_glob           = nValMean
    nValMeanSquare_glob     = nValMeanSquare
    nValAdd_glob            = nValAdd

    nValMean_glob(5)        = nGlobalElems
    nValMeanSquare_glob(5)  = nGlobalElems
    nValAdd_glob(2)         = nGlobalElems

    ! Output
    FileNameOut = Args(iArg)(:LEN(TRIM(Args(1)))-3)//'_DG.h5'
#ifdef FLANG
    IF (MPIRoot) CALL SYSTEM("cp -f "//TRIM(InputFile)//" "//TRIM(FileNameOut))
#else
    IF (MPIRoot) CALL EXECUTE_COMMAND_LINE("cp -f "//TRIM(InputFile)//" "//TRIM(FileNameOut))
#endif
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
    CALL OpenDataFile(FileNameOut,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
    CALL WriteArray('Mean'      ,5,nValMean_glob(1:5)      ,nValMean(1:5)      ,(/0,0,0,0,OffsetElem/),collective=.TRUE.,RealArray=UMeanDG)
    CALL WriteArray('MeanSquare',5,nValMeanSquare_glob(1:5),nValMeanSquare(1:5),(/0,0,0,0,OffsetElem/),collective=.TRUE.,RealArray=UMeanSquareDG)
    CALL WriteArray('ElemData'  ,2,nValAdd_glob(1:2)       ,nValAdd(1:2)       ,(/0,OffsetElem/)      ,collective=.TRUE.,RealArray=ElemDataDG)
    CALL CloseDataFile()
  END SELECT

  SWRITE(UNIT_stdOut,'(A)') " DONE!"

  SELECT CASE(TRIM(FileType))
  CASE('State')
    SDEALLOCATE(VarNames)
    SDEALLOCATE(U)
    SDEALLOCATE(UDG)
  CASE('TimeAvg')
    SDEALLOCATE(VarNamesMean)
    SDEALLOCATE(VarNamesMeanSquare)
    SDEALLOCATE(VarNamesFluc)
    SDEALLOCATE(UMean)
    SDEALLOCATE(UMeanSquare)
    SDEALLOCATE(UMeanDG)
    SDEALLOCATE(UMeanSquareDG)
  END SELECT
  SDEALLOCATE(ElemData)
  SDEALLOCATE(ElemDataDG)

  SWRITE(UNIT_stdOut,'(132("="))')
END DO !iArg

SDEALLOCATE(Args)
SDEALLOCATE(FV_Vdm)
SDEALLOCATE(FV_sVdm)

SWRITE(UNIT_stdOut,'(A)') "===================================================   CALCULATE  DG  FINISHED!   "//&
                          "==================================================="
SWRITE(UNIT_stdOut,'(132("="))')

#if USE_MPI
CALL FinalizeMPI()
CALL MPI_FINALIZE(iError)
IF (iError.NE.MPI_SUCCESS) STOP 'MPI finalize error'
#endif

CONTAINS

!===================================================================================================================================
!> Read in the mean and mean square data sets from the TimeAvg file
!===================================================================================================================================
SUBROUTINE Readin()
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL  :: found
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A,A,A)') "Read from HDF5 file ", TRIM(InputFile), "..."

CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=FileType)
SELECT CASE(TRIM(FileType))
CASE('State')
  ! Safety check if the number of elements did not change
  CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
  IF (INT(HSIZE(nDims)).NE.nGlobalElems)  STOP 'Number of elements in HDF5 file changed during computation!'
  CALL GetArrayAndName('DG_Solution','VarNames'          ,nVal          ,UTmp          ,VarNames)
  !CALL ReadArray('DG_Solution',5,INT(HSize),0,5,RealArray=UTmp)
  CALL GetArrayAndName('ElemData'  ,'VarNamesAdd'        ,nValAdd       ,ElemDataTmp   ,VarNamesAdd)
  nVar           = nVal(1)
  nVarElemData   = nValAdd(1)
  ALLOCATE(U(          nVar,0:nVal(2)-1,0:nVal(3)-1,0:nVal(4)-1,nElems))
  ALLOCATE(ElemData(   nVarElemData,nElems))
  U           = RESHAPE(UTmp          ,(/nVar      ,nVal(2),nVal(3),nVal(4),nElems/))
  ElemData    = RESHAPE(ElemDataTmp   ,(/nVarElemData,  nElems/))
CASE('TimeAvg')
  ! Safety check if the number of elements did not change
  CALL GetDataSize(File_ID,'Mean',nDims,HSize)
  IF (INT(HSIZE(nDims)).NE.nGlobalElems)  STOP 'Number of elements in HDF5 file changed during computation!'
  CALL GetArrayAndName('Mean'      ,'VarNames_Mean'      ,nValMean      ,UMeanTmp      ,VarNamesMean)
  CALL GetArrayAndName('MeanSquare','VarNames_MeanSquare',nValMeanSquare,UMeanSquareTmp,VarNamesMeanSquare)
  CALL GetArrayAndName('ElemData'  ,'VarNamesAdd'        ,nValAdd       ,ElemDataTmp   ,VarNamesAdd)
  nVarMean       = nValMean(1)
  nVarMeanSquare = nValMeanSquare(1)
  nVarElemData   = nValAdd(1)
  ALLOCATE(UMean(      nVarMean,nValMean(2),nValMean(3),nValMean(4),nElems))
  ALLOCATE(UMeanSquare(nVarMeanSquare,nValMeanSquare(2),nValMeanSquare(3),nValMeanSquare(4),nElems))
  ALLOCATE(ElemData(   nVarElemData,nElems))
  UMean       = RESHAPE(UMeanTmp      ,(/nVarMean      ,nValMean(2),nValMean(3),nValMean(4),nElems/))
  UMeanSquare = RESHAPE(UMeanSquareTmp,(/nVarMeanSquare,nValMeanSquare(2),nValMeanSquare(3),nValMeanSquare(4),nElems/))
  ElemData    = RESHAPE(ElemDataTmp   ,(/nVarElemData  ,nElems/))
  CALL DatasetExists(File_ID, TRIM("Fluc"), found)
  IF(found) STOP 'Run toDG before calculating fluctuations!'
END SELECT

CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(A)') "DONE!"
END SUBROUTINE

END PROGRAM toDG
