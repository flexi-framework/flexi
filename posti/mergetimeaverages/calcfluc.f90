#include "flexi.h"

!===================================================================================================================================
!> This tool will take a TimeAvg file and calculate fluctuations from Mean and MeanSquare values
!===================================================================================================================================
PROGRAM CalcFluc
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_StringTools,ONLY:STRICMP
USE MOD_IO_HDF5
USE MOD_Mesh_Vars,  ONLY:nGlobalElems,nElems
USE MOD_HDF5_Input, ONLY:ReadAttribute,ReadArray,GetDataSize,DataSetExists,GetArrayAndName,GetDataProps
USE MOD_HDF5_Output,ONLY:WriteArray,WriteAttribute
USE MOD_MPI,        ONLY:InitMPI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE AND PARAMETER DEFINITIONS
INTEGER                        :: iArg
CHARACTER(LEN=255)             :: InputFile
INTEGER                        :: nVarMean,nVarMeanSquare,nVar_State,NState
INTEGER                        :: nDOF
INTEGER                        :: iM,iMS
INTEGER                        :: iVarFluc,nVarFluc
INTEGER                        :: velIndexM(3),VelIndexMS(7)
INTEGER,ALLOCATABLE            :: MStoM(:)
INTEGER                        :: nValMean(15),nValMeanSquare(15)
REAL,ALLOCATABLE               :: UMean(:,:),UMeanSquare(:,:),UFluc(:,:)
REAL,ALLOCATABLE               :: UMeanTmp(:),UMeanSquareTmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesMean(:),VarNamesMeanSquare(:),VarNamesFluc(:)
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL Abort(__STAMP__,'Parallel execution not yet implemented. Please run on one processor only.')

WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') '===================================================   CALCULATE   FLUCTUATIONS   '//&
                         '==================================================='
WRITE(UNIT_stdOut,'(132("="))')

CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF (nArgs.LT.1) THEN
  WRITE(UNIT_stdOut,*) "Please use: calcfluc <FILE1> <FILE2> ..."
  STOP
END IF

DO iArg=1,nArgs
  WRITE(UNIT_stdOut,"(A,I0,A,I0,A)") "Processing File ", iArg, " of ", nArgs, "..."
  WRITE(UNIT_stdOut,'(132("-"))')

  InputFile=Args(iArg)

  CALL Readin()

  WRITE(UNIT_stdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Search for variables..."

  ALLOCATE(MStoM(nVarMeanSquare)) !index mapping iVar mean to mean suared
  MStoM=0
  velIndexM=0  !index mapping mean         velocities for Reynolds stresses
  velIndexMS=0 !index mapping mean squared velocities for Reynolds stresses
  nVarFluc=0

  !Check which flucs can be calculated

  !First, Check for equal variable names in Mean and MeanSquare
  DO iMS=1,nVarMeanSquare
    DO iM=1,nVarMean
      IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),TRIM(VarNamesMean(iM)))) THEN
        nVarFluc=nVarFluc+1
        MStoM(iMS)=iM
        EXIT
      END IF
    END DO
    !MS index mapping for Reynolds stresses
    IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),'VelocityX')) velIndexMS(1)=iMS
    IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),'VelocityY')) velIndexMS(2)=iMS
    IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),'VelocityZ')) velIndexMS(3)=iMS
    IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),'uv'       )) velIndexMS(4)=iMS
    IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),'uw'       )) velIndexMS(5)=iMS
    IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),'vw'       )) velIndexMS(6)=iMS
    IF(STRICMP(TRIM(VarNamesMeanSquare(iMS)),'TKE'      )) velIndexMS(7)=iMS
  END DO

  !M index mapping for Reynolds stresses
  DO iM=1,nVarMean
    IF(STRICMP(TRIM(VarNamesMean(iM)),'VelocityX')) velIndexM(1)=iM
    IF(STRICMP(TRIM(VarNamesMean(iM)),'VelocityY')) velIndexM(2)=iM
    IF(STRICMP(TRIM(VarNamesMean(iM)),'VelocityZ')) velIndexM(3)=iM
  END DO

  !Add Reynolds stresses to nVar
  IF((VelIndexMS(4).GT.0).AND.(VelIndexM(1).GT.0).AND.(VelIndexM(2).GT.0)) nVarFluc=nVarFluc+1
  IF((VelIndexMS(5).GT.0).AND.(VelIndexM(1).GT.0).AND.(VelIndexM(3).GT.0)) nVarFluc=nVarFluc+1
  IF((VelIndexMS(6).GT.0).AND.(VelIndexM(2).GT.0).AND.(VelIndexM(3).GT.0)) nVarFluc=nVarFluc+1
  IF(((VelIndexMS(7).GT.0).OR.(MINVAL(velIndexMS(1:3)).GT.0)) .AND. (MINVAL(VelIndexM).GT.0)) nVarFluc=nVarFluc+1

  WRITE(UNIT_stdOut,'(A)') " DONE!"
  WRITE(UNIT_stdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A)') "Perform fluctuation calculations. Calculated variables:"

  ALLOCATE(UFluc(nVarFluc,nDOF))
  ALLOCATE(VarNamesFluc(nVarFluc))

  !calculate RMS values (equal variable names M and MS)
  iVarFluc=0
  DO iMS=1,nVarMeanSquare
    IF(MStoM(iMS).GT.0)THEN
      iVarFluc=iVarFluc+1
      iM=MSToM(iMS)
      UFluc(iVarFluc,:) = UMeanSquare(iMS,:) - UMean(iM,:) * UMean(iM,:)
      VarNamesFluc(iVarFluc)=TRIM(VarNamesMeanSquare(iMS))
      WRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
    END IF
  END DO

  ! (u'v')_avg = ((u -U)*(v-V))_avg = (uv)_avg - (uV)_avg - (Uv)_avg + (UV)_avg = (uv)_avg - UV
  IF((VelIndexMS(4).GT.0).AND.(VelIndexM(1).GT.0).AND.(VelIndexM(2).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:) = UMeanSquare(velIndexMS(4),:) - UMean(velIndexM(1),:) * UMean(velIndexM(2),:)
      VarNamesFluc(iVarFluc)='uv'
      WRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  ! (u'w')_avg = ((u -U)*(w-W))_avg = (uw)_avg - (uW)_avg - (Uw)_avg + (UW)_ag = (uw)_avg - UW
  IF((VelIndexMS(5).GT.0).AND.(VelIndexM(1).GT.0).AND.(VelIndexM(3).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:) = UMeanSquare(velIndexMS(5),:) - UMean(velIndexM(1),:) * UMean(velIndexM(3),:)
      VarNamesFluc(iVarFluc)='uw'
      WRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  ! (v'w')_avg = ((v -V)*(w-W))_avg = (vw)_avg - (vW)_avg - (Vw)_avg + (VW)_avg = (vw)_avg - VW
  IF((VelIndexMS(6).GT.0).AND.(VelIndexM(2).GT.0).AND.(VelIndexM(3).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:) = UMeanSquare(velIndexMS(6),:) - UMean(velIndexM(2),:) * UMean(velIndexM(3),:)
      VarNamesFluc(iVarFluc)='vw'
      WRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  ! TKE can be calculated from MeanSquared variable TKE or from uu+vv+ww
  IF((VelIndexMS(7).GT.0) .AND. (MINVAL(VelIndexM).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:) = UMeanSquare(velIndexMS(7),:) - UMean(velIndexM(1),:) * UMean(velIndexM(1),:) &
                                                       - UMean(velIndexM(2),:) * UMean(velIndexM(2),:) &
                                                       - UMean(velIndexM(3),:) * UMean(velIndexM(3),:)
      VarNamesFluc(iVarFluc)='TKE'
      WRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  ELSEIF((MINVAL(velIndexMS(1:3)).GT.0) .AND. (MINVAL(VelIndexM).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:) = UMeanSquare(velIndexMS(1),:) + UMeanSquare(velIndexMS(2),:) + UMeanSquare(velIndexMS(3),:) &
                                                       - UMean(velIndexM(1),:) * UMean(velIndexM(1),:)               &
                                                       - UMean(velIndexM(2),:) * UMean(velIndexM(2),:)               &
                                                       - UMean(velIndexM(3),:) * UMean(velIndexM(3),:)
      VarNamesFluc(iVarFluc)='TKE'
      WRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  WRITE(UNIT_stdOut,'(A)') "Perform calculations DONE!"
  WRITE(UNIT_stdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Write data to HDF5 file..."

  !nValFluc=nValMean
  nValMean(1)=nVarFluc

  !Output
  CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
  CALL WriteArray('Fluc',5,nValMean,nValMean,0*nValMean,collective=.TRUE.,RealArray=UFluc)
  CALL WriteAttribute(File_ID,'VarNames_Fluc',nVarFluc,StrArray=VarNamesFluc)
  CALL CloseDataFile()

  WRITE(UNIT_stdOut,'(A)') " DONE!"

  SDEALLOCATE(VarNamesMean)
  SDEALLOCATE(VarNamesMeanSquare)
  SDEALLOCATE(VarNamesFluc)
  SDEALLOCATE(UMean)
  SDEALLOCATE(UMeanSquare)
  SDEALLOCATE(UFluc)
  SDEALLOCATE(MStoM)

  WRITE(UNIT_stdOut,'(132("="))')
END DO !iArg

SDEALLOCATE(Args)

WRITE(UNIT_stdOut,'(A)') "==============================================   CALCULATE  FLUCTUATIONS  FINISHED!   "//&
                          "=============================================="
WRITE(UNIT_stdOut,'(132("="))')

CONTAINS



!===================================================================================================================================
!> Reads TimeAvg file, checks for Mean and MeanSquare DataSets, gets data size and reads data
!===================================================================================================================================
SUBROUTINE Readin()
IMPLICIT NONE
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,A,A)') "Read from HDF5 file ", TRIM(InputFile), "..."

! Open the first statefile to read necessary attributes
CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataProps(nVar_State,NState,nGlobalElems)
nElems=nGlobalElems
CALL GetArrayAndName('Mean','VarNames_Mean',nValMean,UMeanTmp,VarNamesMean)
CALL GetArrayAndName('MeanSquare','VarNames_MeanSquare',nValMeanSquare,UMeanSquareTmp,VarNamesMeanSquare)
nVarMean=nValMean(1)
nVarMeanSquare=nValMeanSquare(1)
nDOF=PRODUCT(nValMean(2:5))
ALLOCATE(UMean(nVarMean,nDOF))
ALLOCATE(UMeanSquare(nVarMeanSquare,nDOF))
UMean = RESHAPE(UMeanTmp,(/nVarMean,nDOF/))
UMeanSquare = RESHAPE(UMeanSquareTmp,(/nVarMeanSquare,nDOF/))
CALL CloseDataFile()

WRITE(UNIT_stdOut,'(A)') "DONE!"
END SUBROUTINE


END PROGRAM CalcFluc
