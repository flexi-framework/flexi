!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
!> Module containing the main procedures for the POSTI tool: visu_requestInformation is called by ParaView to create a
!> list of available variables and visu is the main routine of POSTI.
!===================================================================================================================================
MODULE MOD_Visu_Avg2D
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitAverage2D
  MODULE PROCEDURE InitAverage2D
END INTERFACE

INTERFACE BuildVandermonds_Avg2D
  MODULE PROCEDURE BuildVandermonds_Avg2D
END INTERFACE

INTERFACE Average2D
  MODULE PROCEDURE Average2D
END INTERFACE

PUBLIC:: InitAverage2D
PUBLIC:: BuildVandermonds_Avg2D
PUBLIC:: Average2D

CONTAINS

SUBROUTINE InitAverage2D()
USE MOD_PreProc
USE MOD_Globals
USE MOD_Visu_Vars
USE MOD_HDF5_Input         ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input         ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile,ReadArray,DatasetExists
USE MOD_Mesh_Vars          ,ONLY: nElems,offsetElem
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: exists
INTEGER                          :: ii,jj,iElem
!===================================================================================================================================
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,'nElems_IJK',exists)
IF (exists) THEN
  SDEALLOCATE(Elem_IJK)
  ALLOCATE(Elem_IJK(3,nElems))
  CALL ReadArray('nElems_IJK',1,(/3/),0,1,IntegerArray=nElems_IJK)
  CALL ReadArray('Elem_IJK',2,(/3,nElems/),offsetElem,2,IntegerArray=Elem_IJK)
ELSE 
  CALL CollectiveStop(__STAMP__,&
      "No Elem_IJK sorting found in mesh file. Required for Avg2D!")
END IF
CALL CloseDataFile()

SDEALLOCATE(FVAmountAvg2D)
ALLOCATE(FVAmountAvg2D(nElems_IJK(1),nElems_IJK(2)))
FVAmountAvg2D = 0.
DO iElem=1,nElems
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  FVAmountAvg2D(ii,jj) = FVAmountAvg2D(ii,jj) + FV_Elems_loc(iElem)
END DO
FVAmountAvg2D = FVAmountAvg2D / REAL(nElems_IJK(3))

nElemsAvg2D_DG = 0
nElemsAvg2D_FV = 0
SDEALLOCATE(mapElemIJToDGElemAvg2D)
SDEALLOCATE(mapElemIJToFVElemAvg2D)
ALLOCATE(mapElemIJToDGElemAvg2D(nElems_IJK(1),nElems_IJK(2)))
ALLOCATE(mapElemIJToFVElemAvg2D(nElems_IJK(1),nElems_IJK(2)))
mapElemIJToDGElemAvg2D = 0
mapElemIJToDGElemAvg2D = 0
DO iElem=1,nElems
  IF (Elem_IJK(3,iElem).EQ.1) THEN
    ii = Elem_IJK(1,iElem)
    jj = Elem_IJK(2,iElem)
    IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN
      nElemsAvg2D_DG = nElemsAvg2D_DG + 1
      mapElemIJToDGElemAvg2D(ii,jj) = nElemsAvg2D_DG
    ELSE
      nElemsAvg2D_FV = nElemsAvg2D_FV + 1
      mapElemIJToFVElemAvg2D(ii,jj) = nElemsAvg2D_FV
    END IF
  END IF
END DO
END SUBROUTINE InitAverage2D

SUBROUTINE BuildVandermonds_Avg2D(NCalc_DG,NCalc_FV) 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars          ,ONLY: Vdm_DGToFV,Vdm_FVToDG,Vdm_DGToVisu,Vdm_FVToVisu,NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu,NVisu_FV
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
#if FV_ENABLED
USE MOD_FV_Basis           ,ONLY: FV_GetVandermonde
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN) :: NCalc_DG,NCalc_FV
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE   :: FV_Vdm(:,:),FV_sVdm(:,:)
REAL,ALLOCATABLE   :: FVdouble(:,:)
INTEGER            :: i
!===================================================================================================================================
SDEALLOCATE(Vdm_DGToFV)
SDEALLOCATE(Vdm_FVToDG)
SDEALLOCATE(Vdm_DGToVisu)
SDEALLOCATE(Vdm_FVToVisu)
#if FV_ENABLED
ALLOCATE(Vdm_DGToFV  (0:NCalc_FV,0:NCalc_DG))
ALLOCATE(Vdm_FVToDG  (0:NCalc_DG,0:NCalc_FV))
ALLOCATE(Vdm_DGToVisu(0:NVisu   ,0:NCalc_DG))
ALLOCATE(Vdm_FVToVisu(0:NVisu_FV,0:NCalc_FV))
ALLOCATE(FVdouble(0:NVisu_FV,0:PP_N))
FVdouble = 0.
DO i = 0, PP_N
  FVdouble(i*2  ,i) = 1. 
  FVdouble(i*2+1,i) = 1.
END DO ! i = 0, PP_N

IF (NCalc_FV.EQ.NCalc_DG) THEN
  CALL FV_GetVandermonde(PP_N,NodeType, Vdm_DGToFV, Vdm_FVToDG)
  Vdm_FVToVisu = FVdouble
ELSE IF (NCalc_FV.EQ.(PP_N+1)*2-1) THEN
  ALLOCATE(FV_Vdm  (0:PP_N    ,0:PP_N))
  ALLOCATE(FV_sVdm (0:PP_N    ,0:PP_N))
  CALL FV_GetVandermonde(PP_N,NodeType, FV_Vdm, FV_sVdm)
  Vdm_DGToFV = MATMUL(FVdouble,FV_Vdm)
  Vdm_FVToDG = MATMUL(FV_sVdm,TRANSPOSE(FVdouble))
  Vdm_FVToVisu = 0.
  DO i = 0, NVisu_FV
    Vdm_FVToVisu(i,i) = 1.
  END DO
  DEALLOCATE(FV_Vdm)
  DEALLOCATE(FV_sVdm)
ELSE
  CALL CollectiveStop(__STAMP__,&
      "BuildVandermonds_Avg2D called with wrong NCalc_DG and NCalc_FV")
END IF
DEALLOCATE(FVdouble)
#else
ALLOCATE(Vdm_DGToFV  (0:0       ,0:NCalc_DG))
ALLOCATE(Vdm_FVToDG  (0:NCalc_DG,0:0       ))
ALLOCATE(Vdm_DGToVisu(0:NVisu   ,0:NCalc_DG))
ALLOCATE(Vdm_FVToVisu(0:NVisu_FV,0:0       ))
#endif
CALL GetVandermonde(NCalc_DG,NodeType,NVisu,NodeTypeVisuPosti,Vdm_DGToVisu,modal=.FALSE.)
END SUBROUTINE BuildVandermonds_Avg2D

SUBROUTINE Average2D(nVarCalc_DG,nVarCalc_FV,NCalc_DG,NCalc_FV,nElems_DG,nElems_FV,&
    NodeTypeCalc_DG,UCalc_DG,UCalc_FV,&
    Vdm_DGToFV,Vdm_FVToDG,Vdm_DGToVisu,Vdm_FVToVisu, &
    startIndexMapVarCalc,endIndexMapVarCalc,mapVarCalc, &
    UVisu_DG,UVisu_FV) 
USE MOD_PreProc
USE MOD_Globals
USE MOD_Visu_Vars          ,ONLY: Elem_IJK,nElems_IJK,FVAmountAvg2D
USE MOD_Visu_Vars          ,ONLY: mapDGElemsToAllElems,mapFVElemsToAllElems
USE MOD_Visu_Vars          ,ONLY: mapElemIJToDGElemAvg2D,mapElemIJToFVElemAvg2D,mapAllVarsToVisuVars
USE MOD_Visu_Vars          ,ONLY: nVarVisu,NVisu,NVisu_FV,nElemsAvg2D_FV,nElemsAvg2D_DG,NodeTypeVisuPosti
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVISUFVEqui
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)            :: nVarCalc_DG,nVarCalc_FV,NCalc_DG,NCalc_FV,nElems_DG,nElems_FV
CHARACTER(LEN=255),INTENT(IN) :: NodeTypeCalc_DG
REAL,INTENT(IN)               :: Vdm_DGToFV(0:NCalc_FV,0:NCalc_DG)
REAL,INTENT(IN)               :: Vdm_FVToDG(0:NCalc_DG,0:NCalc_FV)
REAL,INTENT(IN)               :: Vdm_DGToVisu(0:NVisu   ,0:NCalc_DG)
REAL,INTENT(IN)               :: Vdm_FVToVisu(0:NVisu_FV,0:NCalc_FV)
INTEGER,INTENT(IN)            :: startIndexMapVarCalc,endIndexMapVarCalc,mapVarCalc(startIndexMapVarCalc:endIndexMapVarCalc)
REAL,INTENT(IN)               :: UCalc_DG(0:NCalc_DG,0:NCalc_DG,0:NCalc_DG,1:nElems_DG     ,1:nVarCalc_DG)
REAL,INTENT(IN)               :: UCalc_FV(0:NCalc_FV,0:NCalc_FV,0:NCalc_FV,1:nElems_FV     ,1:nVarCalc_FV)
REAL,INTENT(INOUT)            :: UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,1:nElemsAvg2D_DG,1:nVarVisu)
REAL,INTENT(INOUT)            :: UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,1:nElemsAvg2D_FV,1:nVarVisu)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

REAL,ALLOCATABLE :: Utmp(:,:),Utmp2(:,:)
INTEGER          :: iElem, iElem_DG, iElem_FV, iElemAvg, ii,jj,k,iVar,iVarCalc,iVarVisu
REAL             :: xGP(0:NCalc_DG),wGP(0:NCalc_DG)
REAL,ALLOCATABLE :: UAvg_DG(:,:,:,:)
REAL,ALLOCATABLE :: UAvg_FV(:,:,:,:)
!===================================================================================================================================
ALLOCATE(UAvg_DG(0:NCalc_DG,0:NCalc_DG,nElemsAvg2D_DG,nVarVisu))
ALLOCATE(UAvg_FV(0:NCalc_FV,0:NCalc_FV,nElemsAvg2D_FV,nVarVisu))
UAvg_DG = 0.
UAvg_FV = 0.
CALL GetNodesAndWeights(NCalc_DG,NodeTypeCalc_DG,xGP,wGP)

! average all DG elements first
SWRITE(*,*) " [Avg2D] Average DG elements"
ALLOCATE(Utmp (0:NCalc_DG,0:NCalc_DG))
ALLOCATE(Utmp2(0:NCalc_FV,0:NCalc_FV))
DO iElem_DG = 1,nElems_DG                ! iterate over all DG elements
  iElem = mapDGElemsToAllElems(iElem_DG) ! get global element index
  IF (iElem.EQ.0) CYCLE
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN  ! the averaged ii,jj-th element is rather a DG element and this
                                         ! element is a DG element => just average this element and add to UAvg_DG
    iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        DO k=0,NCalc_DG
          UAvg_DG(:,:,iElemAvg,iVarVisu) = UAvg_DG(:,:,iElemAvg,iVarVisu) + wGP(k)/2. * UCalc_DG(:,:,k,iElem_DG,iVarCalc)
        END DO
      END IF
    END DO
  ELSE ! the averaged ii,jj-th element is rather a FV element, but this element is a DG element
       ! => average this element as DG element and convert the average (Utmp) to FV and add to UAvg_FV
    iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        Utmp = 0.
        DO k=0,NCalc_DG
          Utmp(:,:) = Utmp(:,:) + wGP(k)/2. * UCalc_DG(:,:,k,iElem_DG,iVarCalc)
        END DO
        CALL ChangeBasis2D(NCalc_DG,NCalc_FV,Vdm_DGToFV,Utmp,Utmp2)
        UAvg_FV(:,:,iElemAvg,iVarVisu) = UAvg_FV(:,:,iElemAvg,iVarVisu) + Utmp2
      END IF
    END DO
  END IF
END DO

! all DG elements are averaged and we will now average all FV elements, which requires conversions
! from FV to DG, the other way round then before.
SWRITE(*,*) " [Avg2D] Average FV elements"
DEALLOCATE(Utmp)
DEALLOCATE(Utmp2)
ALLOCATE(Utmp( 0:NCalc_FV,0:NCalc_FV))
ALLOCATE(Utmp2(0:NCalc_DG,0:NCalc_DG))
DO iElem_FV = 1,nElems_FV                ! iterate over all FV elements
  iElem = mapFVElemsToAllElems(iElem_FV) ! get global element index
  IF (iElem.EQ.0) CYCLE
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  IF (FVAmountAvg2D(ii,jj).GT.0.5) THEN  ! the averaged ii,jj-th element is rather a FV element and this 
                                         ! element is a FV element => just average this element and add to UAvg_FV
    iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        DO k=0,NCalc_FV
          UAvg_FV(:,:,iElemAvg,iVarVisu) = UAvg_FV(:,:,iElemAvg,iVarVisu) + 1./(NCalc_FV+1) * UCalc_FV(:,:,k,iElem_FV,iVarCalc)
        END DO
      END IF
    END DO
  ELSE ! the averaged ii,jj-th element is rather a DG element, but this element is a FV element
       ! => average this element as FV element and convert the average (Utmp) to DG and add to UAvg_DG
    iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        Utmp = 0.
        DO k=0,NCalc_FV
          Utmp(:,:) = Utmp(:,:) +  1./(NCalc_FV+1) * UCalc_FV(:,:,k,iElem_FV,iVarCalc)
        END DO
        CALL ChangeBasis2D(NCalc_FV,NCalc_DG,Vdm_FVToDG,Utmp,Utmp2)
        UAvg_DG(:,:,iElemAvg,iVarVisu) = UAvg_DG(:,:,iElemAvg,iVarVisu) + Utmp2
      END IF
    END DO
  END IF
END DO

! divide average solutions by the number of elements in the average direction
UAvg_DG = UAvg_DG / nElems_IJK(3)
UAvg_FV = UAvg_FV / nElems_IJK(3)

DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
  iVarVisu = mapAllVarsToVisuVars(iVar)
  IF (iVarVisu.GT.0) THEN
    DO iElemAvg=1,nElemsAvg2D_DG
      CALL ChangeBasis2D(NCalc_DG,NVisu   ,Vdm_DGToVisu  ,UAvg_DG(:,:,iElemAvg,iVarVisu),UVisu_DG(:,:,0,iElemAvg,iVarVisu))
    END DO
    IF (NCalc_FV.EQ.NVisu_FV) THEN
      UVisu_FV(:,:,0,:,iVarVisu) = UAvg_FV(:,:,:,iVarVisu)
    ELSE
      DO iElemAvg=1,nElemsAvg2D_FV
        CALL ChangeBasis2D(NCalc_FV,NVisu_FV,Vdm_FVToVisu,UAvg_FV(:,:,iElemAvg,iVarVisu),UVisu_FV(:,:,0,iElemAvg,iVarVisu))
      END DO
    END IF
  END IF
END DO


! clear local variables
DEALLOCATE(Utmp)
DEALLOCATE(Utmp2)
DEALLOCATE(UAvg_DG)
DEALLOCATE(UAvg_FV)
END SUBROUTINE Average2D

END MODULE MOD_Visu_Avg2D
