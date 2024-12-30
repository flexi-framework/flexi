!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
!> This tool will take a TimeAvg file and calculate fluctuations from Mean and MeanSquare values
!>
!> During the simulation, FLEXI will write two data sets: The mean value of a variable and the mean of the square of the variable.
!> We can then make usage of the following equality to calculate fluctuations (=mean of the square of the fluctuations):
!> _____   _____________   _____   ______   ____         ___
!> (U*U) = (u+u')*(u+u') = (u*u) + 2*u*u' + u'u' = u*u + u'u'                                                      ___       _
!> where we split the total value of a variable U in the mean u and the fluctuating part u'. Thus, with the stored U*U and u=U we
!> then calculate the fluctuations in here as:
!> ____   ___
!> u'u' = U*U - u*u
!===================================================================================================================================
PROGRAM CalcFluc
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_IO_HDF5
USE MOD_StringTools,ONLY:STRICMP
USE MOD_Mesh_Vars,  ONLY:nGlobalElems,nElems,OffsetElem
USE MOD_HDF5_Input, ONLY:ReadAttribute,ReadArray,GetDataSize,DataSetExists,GetArrayAndName,GetDataProps
USE MOD_HDF5_Output,ONLY:WriteArray,WriteAttribute
USE MOD_MPI,        ONLY:InitMPI
#if USE_MPI
USE MOD_MPI,        ONLY:FinalizeMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE AND PARAMETER DEFINITIONS
INTEGER                        :: iArg
CHARACTER(LEN=255)             :: InputFile
INTEGER                        :: nVarMean,nVarMeanSquare
INTEGER                        :: nDOF
INTEGER                        :: iM,iMS
INTEGER                        :: iVarFluc,nVarFluc
INTEGER                        :: velIndexM(3),VelIndexMS(7)
INTEGER,ALLOCATABLE            :: MStoM(:)
INTEGER                        :: nValMean(15),nValMeanSquare(15)
INTEGER                        :: nValMean_glob(15)
REAL,ALLOCATABLE               :: UMean(:,:,:),UMeanSquare(:,:,:),UFluc(:,:,:)
REAL,ALLOCATABLE               :: UMeanTmp(:),UMeanSquareTmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesMean(:),VarNamesMeanSquare(:),VarNamesFluc(:)
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL InitMPIInfo()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') '===================================================   CALCULATE   FLUCTUATIONS   '//&
                         '==================================================='
SWRITE(UNIT_stdOut,'(132("="))')

CALL ParseCommandlineArguments()

InputFile=Args(1)
! Partitioning - we partition along the last dimension of the arrays
CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataSize(File_ID,'Mean',nDims,HSize)
nGlobalElems = INT(HSize(nDims))
DEALLOCATE(HSize)
#if USE_MPI
nElems = INT(nGlobalElems/nProcessors)
OffsetElem = myRank*nElems
IF (myRank.EQ.nProcessors-1) THEN
  ! This is the last proc, may have other number of elements
  nElems = nElems + MOD(nGlobalElems,nProcessors)
END IF
#else
OffsetElem = 0
nElems = nGlobalElems
#endif
CALL CloseDataFile()

! Check if the number of arguments is correct
IF (nArgs.LT.1) THEN
  SWRITE(UNIT_stdOut,*) "Please use: calcfluc <FILE1> <FILE2> ..."
  STOP
END IF

DO iArg=1,nArgs
  SWRITE(UNIT_stdOut,"(A,I0,A,I0,A)") "Processing File ", iArg, " of ", nArgs, "..."
  SWRITE(UNIT_stdOut,'(132("-"))')

  InputFile=Args(iArg)

  CALL Readin()

  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Search for variables..."

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

  SWRITE(UNIT_stdOut,'(A)') " DONE!"
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') "Perform fluctuation calculations. Calculated variables:"

  ALLOCATE(UFluc(nVarFluc,nDOF,nElems))
  ALLOCATE(VarNamesFluc(nVarFluc))

  !calculate RMS values (equal variable names M and MS)
  iVarFluc=0
  DO iMS=1,nVarMeanSquare
    IF(MStoM(iMS).GT.0)THEN
      iVarFluc=iVarFluc+1
      iM=MSToM(iMS)
      UFluc(iVarFluc,:,:) = UMeanSquare(iMS,:,:) - UMean(iM,:,:) * UMean(iM,:,:)
      VarNamesFluc(iVarFluc)=TRIM(VarNamesMeanSquare(iMS))
      SWRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
    END IF
  END DO

  ! (u'v')_avg = ((u -U)*(v-V))_avg = (uv)_avg - (uV)_avg - (Uv)_avg + (UV)_avg = (uv)_avg - UV
  IF((VelIndexMS(4).GT.0).AND.(VelIndexM(1).GT.0).AND.(VelIndexM(2).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:,:) = UMeanSquare(velIndexMS(4),:,:) - UMean(velIndexM(1),:,:) * UMean(velIndexM(2),:,:)
      VarNamesFluc(iVarFluc)='uv'
      SWRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  ! (u'w')_avg = ((u -U)*(w-W))_avg = (uw)_avg - (uW)_avg - (Uw)_avg + (UW)_ag = (uw)_avg - UW
  IF((VelIndexMS(5).GT.0).AND.(VelIndexM(1).GT.0).AND.(VelIndexM(3).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:,:) = UMeanSquare(velIndexMS(5),:,:) - UMean(velIndexM(1),:,:) * UMean(velIndexM(3),:,:)
      VarNamesFluc(iVarFluc)='uw'
      SWRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  ! (v'w')_avg = ((v -V)*(w-W))_avg = (vw)_avg - (vW)_avg - (Vw)_avg + (VW)_avg = (vw)_avg - VW
  IF((VelIndexMS(6).GT.0).AND.(VelIndexM(2).GT.0).AND.(VelIndexM(3).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:,:) = UMeanSquare(velIndexMS(6),:,:) - UMean(velIndexM(2),:,:) * UMean(velIndexM(3),:,:)
      VarNamesFluc(iVarFluc)='vw'
      SWRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  ! TKE can be calculated from MeanSquared variable TKE or from uu+vv+ww
  IF((VelIndexMS(7).GT.0) .AND. (MINVAL(VelIndexM).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:,:) = UMeanSquare(velIndexMS(7),:,:) - UMean(velIndexM(1),:,:) * UMean(velIndexM(1),:,:) &
                                                           - UMean(velIndexM(2),:,:) * UMean(velIndexM(2),:,:) &
                                                           - UMean(velIndexM(3),:,:) * UMean(velIndexM(3),:,:)
      VarNamesFluc(iVarFluc)='TKE'
      SWRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  ELSEIF((MINVAL(velIndexMS(1:3)).GT.0) .AND. (MINVAL(VelIndexM).GT.0)) THEN
      iVarFluc=iVarFluc+1
      UFluc(iVarFluc,:,:) = UMeanSquare(velIndexMS(1),:,:) + UMeanSquare(velIndexMS(2),:,:) + UMeanSquare(velIndexMS(3),:,:) &
                                                           - UMean(velIndexM(1),:,:) * UMean(velIndexM(1),:,:)               &
                                                           - UMean(velIndexM(2),:,:) * UMean(velIndexM(2),:,:)               &
                                                           - UMean(velIndexM(3),:,:) * UMean(velIndexM(3),:,:)
      VarNamesFluc(iVarFluc)='TKE'
      SWRITE(UNIT_stdOut,*) "  "//TRIM(VarNamesFluc(iVarFluc))
  END IF

  SWRITE(UNIT_stdOut,'(A)') "Perform calculations DONE!"
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Write data to HDF5 file..."

  nValMean(1)      = nVarFluc
  nValMean_glob    = nValMean
  nValMean_glob(5) = nGlobalElems

  !Output
  CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
  CALL WriteArray('Fluc',5,nValMean_glob(1:5),nValMean(1:5),(/0,0,0,0,OffsetElem/),collective=.TRUE.,RealArray=UFluc)
  CALL WriteAttribute(File_ID,'VarNames_Fluc',nVarFluc,StrArray=VarNamesFluc)
  CALL CloseDataFile()

  SWRITE(UNIT_stdOut,'(A)') " DONE!"

  SDEALLOCATE(VarNamesMean)
  SDEALLOCATE(VarNamesMeanSquare)
  SDEALLOCATE(VarNamesFluc)
  SDEALLOCATE(UMean)
  SDEALLOCATE(UMeanSquare)
  SDEALLOCATE(UFluc)
  SDEALLOCATE(MStoM)

  SWRITE(UNIT_stdOut,'(132("="))')
END DO !iArg

SDEALLOCATE(Args)

SWRITE(UNIT_stdOut,'(A)') "==============================================   CALCULATE  FLUCTUATIONS  FINISHED!   "//&
                          "=============================================="
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
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A,A,A)') "Read from HDF5 file ", TRIM(InputFile), "..."

CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
! Safety check if the number of elements did not change
CALL GetDataSize(File_ID,'Mean',nDims,HSize)
IF (INT(HSize(nDims)).NE.nGlobalElems)  STOP 'Number of elements in HDF5 file changed during computation!'
DEALLOCATE(HSize)
CALL GetArrayAndName('Mean'      ,'VarNames_Mean'      ,nValMean      ,UMeanTmp      ,VarNamesMean)
CALL GetArrayAndName('MeanSquare','VarNames_MeanSquare',nValMeanSquare,UMeanSquareTmp,VarNamesMeanSquare)
nVarMean       = nValMean(1)
nVarMeanSquare = nValMeanSquare(1)
nDOF           = PRODUCT(nValMean(2:4))
ALLOCATE(UMean(       nVarMean,     nDOF,nElems))
ALLOCATE(UMeanSquare(nVarMeanSquare,nDOF,nElems))
UMean       = RESHAPE(UMeanTmp,      (/nVarMean,      nDOF,nElems/))
UMeanSquare = RESHAPE(UMeanSquareTmp,(/nVarMeanSquare,nDOF,nElems/))
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(A)') "DONE!"
END SUBROUTINE

END PROGRAM CalcFluc
