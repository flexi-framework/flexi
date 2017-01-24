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

INTERFACE Average2D
  MODULE PROCEDURE Average2D
END INTERFACE

PUBLIC:: Average2D

CONTAINS

SUBROUTINE Average2D() 
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeFVEqui,wGP
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

REAL,ALLOCATABLE :: Utmp(:,:,:),Utmp2(:,:), Vdm_N_NVisu_FV(:,:), Vdm_NVisu_FV_N(:,:)
INTEGER          :: iElem, iElem_DG, iElem_FV, iElemAvg, ii,jj,k,iVar
!===================================================================================================================================
SDEALLOCATE(UAvg_DG)
SDEALLOCATE(UAvg_FV)
ALLOCATE(UAvg_DG(0:PP_N    ,0:PP_N    ,0:0,nElemsAvg2D_DG,nVarCalc))
ALLOCATE(UAvg_FV(0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV,nVarCalc))
ALLOCATE(Utmp(0:PP_N,0:PP_N,nVarCalc))
ALLOCATE(Utmp2(0:NVisu_FV,0:NVisu_FV))
UAvg_DG = 0.
UAvg_FV = 0.
ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:PP_N))
ALLOCATE(Vdm_NVisu_FV_N(0:PP_N,0:NVisu_FV))
CALL GetVandermonde(PP_N,NodeType,NVisu_FV,NodeTypeFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
CALL GetVandermonde(NVisu_FV,NodeTypeFVEqui,PP_N,NodeType,Vdm_N_NVisu_FV,modal=.FALSE.)
DO iElem_DG = 1,nElems_DG ! iterate over all DG elements
  iElem = mapDGElemsToAllElems(iElem_DG)
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN ! DG
    iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
    DO k=0,PP_N
      UAvg_DG(:,:,0,iElemAvg,:) = UAvg_DG(:,:,0,iElemAvg,:) + wGP(k)/2. * UCalc_DG(:,:,k,iElem_DG,:)
    END DO
  ELSE ! FV
    Utmp = 0.
    DO k=0,PP_N
      Utmp(:,:,:) = Utmp(:,:,:) + wGP(k)/2. * UCalc_DG(:,:,k,iElem_DG,:)
    END DO
    iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
    DO iVar=1,nVarCalc
      CALL ChangeBasis2D(PP_N,NVisu_FV,Vdm_N_NVisu_FV,Utmp(:,:,iVar),Utmp2)
      UAvg_FV(:,:,0,iElemAvg,iVar) = UAvg_FV(:,:,0,iElemAvg,iVar) + Utmp2
    END DO
  END IF
END DO
DEALLOCATE(Utmp)
DEALLOCATE(Utmp2)
ALLOCATE(Utmp(0:NVisu_FV,0:NVisu_FV,nVarCalc))
ALLOCATE(Utmp2(0:PP_N,0:PP_N))
DO iElem_FV = 1,nElems_FV ! iterate over all FV elements
  iElem = mapFVElemsToAllElems(iElem_FV)
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  IF (FVAmountAvg2D(ii,jj).GT.0.5) THEN ! FV
    iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
    DO k=0,NVisu_FV
      UAvg_FV(:,:,0,iElemAvg,:) = UAvg_FV(:,:,0,iElemAvg,:) + 1./(NVisu_FV+1) * UCalc_FV(:,:,k,iElem_FV,:)
    END DO
  ELSE ! DG
    Utmp = 0.
    DO k=0,NVisu_FV
      Utmp(:,:,:) = Utmp(:,:,:) + 1./(NVisu_FV+1) * UCalc_FV(:,:,k,iElem_FV,:)
    END DO
    iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
    DO iVar=1,nVarCalc
      CALL ChangeBasis2D(NVisu_FV,PP_N,Vdm_NVisu_FV_N,Utmp(:,:,iVar),Utmp2)
      UAvg_DG(:,:,0,iElemAvg,iVar) = UAvg_DG(:,:,0,iElemAvg,iVar) + Utmp2
    END DO
  END IF
END DO
UAvg_DG = UAvg_DG / nElems_IJK(3)
UAvg_FV = UAvg_FV / nElems_IJK(3)
DEALLOCATE(Vdm_N_NVisu_FV)
DEALLOCATE(Vdm_NVisu_FV_N)
DEALLOCATE(Utmp)
DEALLOCATE(Utmp2)
END SUBROUTINE Average2D

END MODULE MOD_Visu_Avg2D
