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
!> Contains routines that convert the calculated FV or DG quantities to the visualization grid. There are separate routines
!> to convert the ElemData and FieldData to the visualization grid.
!===================================================================================================================================
MODULE MOD_Posti_ConvertToVisu
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ConvertToVisu_DG
  MODULE PROCEDURE ConvertToVisu_DG
END INTERFACE
PUBLIC:: ConvertToVisu_DG

#if FV_ENABLED
INTERFACE ConvertToVisu_FV
  MODULE PROCEDURE ConvertToVisu_FV
END INTERFACE
PUBLIC:: ConvertToVisu_FV

#if FV_RECONSTRUCT
INTERFACE ConvertToVisu_FV_Reconstruct
  MODULE PROCEDURE ConvertToVisu_FV_Reconstruct
END INTERFACE
PUBLIC:: ConvertToVisu_FV_Reconstruct
#endif /* FV_RECONSTRUCT */
#endif /* FV_ENABLED */

INTERFACE ConvertToVisu_ElemData
  MODULE PROCEDURE ConvertToVisu_ElemData
END INTERFACE

INTERFACE ConvertToVisu_FieldData
  MODULE PROCEDURE ConvertToVisu_FieldData
END INTERFACE

PUBLIC:: ConvertToVisu_ElemData
PUBLIC:: ConvertToVisu_FieldData

CONTAINS

!===================================================================================================================================
!> Perform a ChangeBasis of the calculated DG quantities to the visualization grid.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_DG() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
USE MOD_EOS_Posti_Vars     ,ONLY: nVarTotal
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,iVar,iVarVisu,iVarCalc
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
!===================================================================================================================================

! compute UVisu_DG 
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVisuPosti,Vdm_N_NVisu,modal=.FALSE.)
! convert DG solution to UVisu_DG
SDEALLOCATE(UVisu_DG)
ALLOCATE(UVisu_DG(0:NVisu,0:NVisu,0:NVisu,nElems_DG,nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData))
DO iVar=1,nVarTotal
  IF (mapVisu(iVar).GT.0) THEN
    iVarCalc = mapCalc(iVar) 
    iVarVisu = mapVisu(iVar) 
    DO iElem = 1,nElems_DG
      CALL ChangeBasis3D(PP_N,NVisu,Vdm_N_NVisu,UCalc_DG(:,:,:,iElem,iVarCalc),UVisu_DG(:,:,:,iElem,iVarVisu))
    END DO
  END IF
END DO 
SDEALLOCATE(Vdm_N_NVisu)
END SUBROUTINE ConvertToVisu_DG

#if FV_ENABLED        
!===================================================================================================================================
!> Convert the calculated FV quantities to the visualization grid.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_FV(mapCalc,maskVisu,reallocate)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars         ,ONLY: mapVisu,UVisu_FV,nVarVisu,NVisu_FV,nElems_FV,UCalc_FV
USE MOD_Posti_Vars         ,ONLY: nVarVisu_ElemData, nVarVisu_FieldData
USE MOD_ReadInTools        ,ONLY: GETINT
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
USE MOD_EOS_Posti_Vars     ,ONLY: nVarTotal,DepNames
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)          :: mapCalc(nVarTotal)
INTEGER,INTENT(IN),OPTIONAL :: maskVisu(nVarTotal)
LOGICAL,INTENT(IN),OPTIONAL :: reallocate
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVar,i,j,k,iElem
INTEGER            :: iVarVisu,iVarCalc
LOGICAL            :: reallocate_loc
!===================================================================================================================================
! compute UVisu_FV
reallocate_loc = MERGE(reallocate, .TRUE., PRESENT(reallocate))
IF (reallocate_loc) THEN
  SDEALLOCATE(UVisu_FV)
  ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV,nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData))
END IF
DO iVar=1,nVarTotal
  iVarVisu = mapVisu(iVar) 
  iVarVisu = MERGE(maskVisu(iVar)*iVarVisu, iVarVisu, PRESENT(maskVisu))
  IF (iVarVisu.GT.0) THEN
    SWRITE(*,*) "    ", TRIM(DepNames(iVar))
    iVarCalc = mapCalc(iVar) 
    DO iElem = 1,nElems_FV
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        UVisu_FV(i*2:i*2+1, j*2:j*2+1, k*2:k*2+1,iElem,iVarVisu) = UCalc_FV(i,j,k,iElem,iVarCalc)
      END DO; END DO; END DO
    END DO
  END IF
END DO 

END SUBROUTINE ConvertToVisu_FV


#if FV_RECONSTRUCT
!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE ConvertToVisu_FV_Reconstruct()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_ReadInTools        ,ONLY: GETINT
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
USE MOD_FV_Vars            ,ONLY: gradUxi,gradUeta,gradUzeta
USE MOD_FV_Vars            ,ONLY: FV_dx_XI_L,FV_dx_ETA_L,FV_dx_ZETA_L
USE MOD_FV_Vars            ,ONLY: FV_dx_XI_R,FV_dx_ETA_R,FV_dx_ZETA_R
USE MOD_EOS                ,ONLY: PrimToCons
USE MOD_DG_Vars            ,ONLY: UPrim
USE MOD_EOS_Posti_Vars     ,ONLY: nVarTotal
USE MOD_EOS_Posti          ,ONLY: GetMaskPrim
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar,i,j,k,iElem,iElem_FV
INTEGER             :: iVarCalc
INTEGER             :: nVarPrim,iVarPrim
INTEGER             :: mapUPrim(PP_nVarPrim)
INTEGER             :: mapUCalc(PP_nVarPrim)
INTEGER,ALLOCATABLE :: maskPrim(:)
!===================================================================================================================================
! Build local maps of maximal size PP_nVarPrim:
! - mapUCalc(1:nVarPrim) = indices of the nVarPrim primitive quantities that should be visualized in the UCalc_FV array
! - mapUPrim(1:nVarPrim) = indices of the nVarPrim primitive quantities in the UPrim array
! Example: 
!   If only velocityX and pressure should be visualized then:
!     nVarPrim = 2 
!     mapUPrim(1) = 2     mapUCalc(1) = index of velocityX in UCalc_FV 
!     mapUPrim(2) = 5     mapUCalc(2) = index of pressure  in UCalc_FV 
nVarPrim = 0
iVarPrim = 0
maskPrim = GetMaskPrim()
DO iVar=1,nVarTotal
  IF (maskPrim(iVar).GT.0) THEN
    iVarPrim = iVarPrim + 1
    IF (mapCalc_FV(iVar).GT.0) THEN
      nVarPrim = nVarPrim + 1
      mapUPrim(nVarPrim) = iVarPrim
      mapUCalc(nVarPrim) = mapCalc_FV(iVar)
    END IF
  END IF
END DO
SWRITE(*,*) "  nVarPrim", nVarPrim
SWRITE(*,*) "  mapUPrim", mapUPrim(1:nVarPrim)
SWRITE(*,*) "  mapUCalc", mapUCalc(1:nVarPrim)
DEALLOCATE(maskPrim)


DO iElem_FV=1,nElems_FV
  iElem = mapElems_FV(iElem_FV)  
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO iVar=1,nVarPrim
      iVarPrim = mapUPrim(iVar)
      iVarCalc = mapUCalc(iVar)
      UCalc_FV(i*2  ,j*2  ,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2  ,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2+1,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem)  &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2  ,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem)  &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2+1,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2  ,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2+1,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2+1,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem)  &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
    END DO
  END DO; END DO; END DO
END DO ! iElem_FV
END SUBROUTINE ConvertToVisu_FV_Reconstruct

#endif /* FV_RECONSTRUCT */

#endif /* FV_ENABLED */

!===================================================================================================================================
!> Convert ElemData variables that should be visualized to the visu grid.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_ElemData() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_Restart_Vars       ,ONLY: nVarElemData,ElemData,VarNamesElemData
USE MOD_StringTools        ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVarElemData,iVarVisu,iElem_DG,iElem_FV,iElem
!===================================================================================================================================
SWRITE(*,*) "Convert ElemData to Visu grid"
DO iVarElemData=1,nVarElemData
  DO iVarVisu=1,nVarVisu_ElemData
    IF (STRICMP(VarNamesElemData(iVarElemData),VarNamesVisu_ElemData(iVarVisu))) THEN
      SWRITE(*,*) "  ",TRIM(VarNamesElemData(iVarElemData))
      DO iElem_DG=1,nElems_DG
        iElem = mapElems_DG(iElem_DG)
        UVisu_DG(:,:,:,iElem_DG,nVarVisu+iVarVisu) = ElemData(iVarElemData,iElem)
      END DO
      DO iElem_FV=1,nElems_FV
        iElem = mapElems_FV(iElem_FV)
        UVisu_FV(:,:,:,iElem_FV,nVarVisu+iVarVisu) = ElemData(iVarElemData,iElem)
      END DO
    END IF
  END DO
END DO 
END SUBROUTINE ConvertToVisu_ElemData

!===================================================================================================================================
!> Convert FieldData variables that should be visualized to the visu grid.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_FieldData() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVarFieldData,iVarVisu,iElem_DG,iElem_FV,iElem
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
REAL,ALLOCATABLE   :: Vdm_N_NVisu_FV(:,:)                ! Vandermonde from state to fv visualisation nodes
!===================================================================================================================================
SWRITE(*,*) "Convert FieldData to Visu grid"

ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu_FV,NodeTypeVisuPosti,Vdm_N_NVisu_FV,modal=.FALSE.)
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVisuPosti,Vdm_N_NVisu,modal=.FALSE.)

DO iVarFieldData=1,nVarFieldData
  DO iVarVisu=1,nVarVisu_FieldData
    IF (STRICMP(VarNamesAddField(iVarFieldData),VarNamesVisu_FieldData(iVarVisu))) THEN
      SWRITE(*,*) "  ",TRIM(VarNamesAddField(iVarFieldData))
      DO iElem_DG=1,nElems_DG
        iElem = mapElems_DG(iElem_DG)
        CALL ChangeBasis3D(PP_N,NVisu,Vdm_N_NVisu,FieldData(iVarFieldData,:,:,:,iElem),&
                                                  UVisu_DG(:,:,:,iElem_DG,nVarVisu+nVarVisu_ElemData+iVarVisu))
      END DO
      DO iElem_FV=1,nElems_FV
        iElem = mapElems_FV(iElem_FV)
        CALL ChangeBasis3D(PP_N,NVisu,Vdm_N_NVisu_FV,FieldData(iVarFieldData,:,:,:,iElem),&
                                                     UVisu_FV(:,:,:,iElem_FV,nVarVisu+nVarVisu_ElemData+iVarVisu))
      END DO
    END IF
  END DO
END DO 
END SUBROUTINE ConvertToVisu_FieldData

END MODULE MOD_Posti_ConvertToVisu
