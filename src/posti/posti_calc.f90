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
!> Containes routines that call the equation system specific routines to calculate primitive or derived quantities both for
!> DG and FV.
!===================================================================================================================================
MODULE MOD_Posti_Calc
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE CalcQuantities_DG
  MODULE PROCEDURE CalcQuantities_DG
END INTERFACE
PUBLIC:: CalcQuantities_DG

#if FV_ENABLED 
INTERFACE CalcQuantities_ConvertToVisu_FV
  MODULE PROCEDURE CalcQuantities_ConvertToVisu_FV
END INTERFACE
PUBLIC:: CalcQuantities_ConvertToVisu_FV
#endif

INTERFACE FillCopy
  MODULE PROCEDURE FillCopy
END INTERFACE
PUBLIC:: FillCopy


CONTAINS

!===================================================================================================================================
!> Calc quantities for all DG elements.
!> 1. Allocate UCalc_DG
!> 2. Call CalcQuantities or CalcQuantitiesWithGradients from eos_posti.f90
!===================================================================================================================================
SUBROUTINE CalcQuantities_DG() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_EOS_Posti      ,ONLY: CalcQuantities
USE MOD_Mesh_Vars      ,ONLY: nElems
USE MOD_DG_Vars        ,ONLY: U
USE MOD_StringTools    ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVar,iVar2,iVarCalc
INTEGER            :: maskCalc(nVarTotal)
!===================================================================================================================================
! calc DG solution 
SWRITE(*,*) "[DG] calc quantities"
SDEALLOCATE(UCalc_DG)
ALLOCATE(UCalc_DG(0:PP_N,0:PP_N,0:PP_N,nElems_DG,1:nVarCalc))
maskCalc=1

! Copy exisiting variables from solution array
DO iVar=1,nVarTotal
  iVarCalc = mapCalc(iVar)
  IF (iVarCalc.LT.1) CYCLE
  DO iVar2=1,nVar_State
    IF(STRICMP(VarNamesTotal(iVar),VarNamesHDF5(iVar2)))THEN
      CALL FillCopy(nVar_State,PP_N,nElems,U,nElems_DG,mapElems_DG,UCalc_DG(:,:,:,:,iVarCalc),iVar2)
      maskCalc(iVar)=0
    END IF
  END DO
END DO

IF(TRIM(FileType).EQ.'State')THEN
  CALL CalcQuantities(nVarCalc,PP_N,nElems_DG,mapElems_DG,mapCalc,UCalc_DG,maskCalc) 
END IF
END SUBROUTINE CalcQuantities_DG


#if FV_ENABLED 
!===================================================================================================================================
!> Calc quantities for all FV elements.
!> 
!> If FV_RECONSTRUCT is switched OFF we only have constant cell average values, which can be just used to calc all desired
!> quantities by calling 'CalcQuantities' and afterwards 'ConvertToVisu_FV' to 'interpolate' them to the visu grid.
!>
!> If FV_RECONSTRUCT is switched ON the DGTimeDerivative_weakForm operator is called once to build the limited slopes of the 
!> primitive quantities. This slopes can be used to convert the primitive quantities to the visu grid (call of
!ConvertToVisu_FV_Reconstruct).
!> Since we only have slopes for the primitive quantities this is not possible for the conservative quantities. Therefore the
!> reconstructed primitives on the visu grid are used to calc the conservative quantities directly on the visu grid (call of
!> CalcConsFromPrim). 
!> For all other quantities, that do not depend on gradients the normal 'CalcQuantities' is called on the visu grid.
!> These and the primitive and conservative quantities are then copied from the UCalc_FV array to the UVisu_FV array.
!> Now only quantities depending on gradients remain. They can not be reconstructed and are visualized as cell average values.
!> Therefore the conservative quantities are recomputed from the primitive ones but now on the cell centers and not on the visu
!> grid. Then 'CalcQuantities' is used to calculate gradient quantities only, which are afterwards converted to the FV visu grid.
!===================================================================================================================================
SUBROUTINE CalcQuantities_ConvertToVisu_FV() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV
USE MOD_EOS_Posti           ,ONLY: CalcQuantities
#if FV_RECONSTRUCT
USE MOD_EOS_Posti           ,ONLY: GetMaskCons,GetMaskPrim,GetMaskGrad
USE MOD_EOS_Posti           ,ONLY: AppendNeededPrims
USE MOD_EOS_Posti           ,ONLY: CalcConsFromPrim
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV_Reconstruct
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: maskCalc(nVarTotal)
#if FV_RECONSTRUCT 
INTEGER                      :: maskCons(nVarTotal),maskPrim(nVarTotal),maskGrad(nVarTotal)
INTEGER                      :: iVar
#endif
!===================================================================================================================================
#if FV_RECONSTRUCT
  ! calc FV solution 
  SDEALLOCATE(mapCalc_FV)
  ALLOCATE(mapCalc_FV(1:nVarTotal))
  CALL AppendNeededPrims(mapCalc,mapCalc_FV,nVarCalc_FV)
  SWRITE (*,*) "[FVRE] nVarCalc_FV", nVarCalc_FV
  SWRITE (*,"(A,27I3)") "  mapCalc_FV",  mapCalc_FV
  
  ! convert primitive quantities to the visu grid
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV,1:nVarCalc_FV))
  SWRITE(*,*) "[FVRE] ConvertToVisu_FV_Reconstruct"
  CALL ConvertToVisu_FV_Reconstruct()

  ! calculate all needed conservative variables on the visu grid
  SWRITE(*,*) "[FVRE] CalcConsFromPrim"
  CALL CalcConsFromPrim(mapCalc_FV,nVarCalc_FV,NVisu_FV,nElems_FV,UCalc_FV)

  ! calculate all nonCons, nonPrim, nonGrad quantities on the visu grid
  maskCons=GetMaskCons()
  maskPrim=GetMaskPrim()
  maskGrad=GetMaskGrad()
  maskCalc = 1-MAX(MAX(maskCons,maskPrim),maskGrad)
  SWRITE(*,*) "[FVRE] CalcQuantities (nonCons,nonPrim,nonGrad)"
  CALL CalcQuantities(nVarCalc_FV,NVisu_FV,nElems_FV,mapElems_FV,mapCalc_FV,UCalc_FV,maskCalc)

  ! copy cons, prim and all other non-grad quantities to the UVisu_FV array
  SWRITE(*,*) "[FVRE] copy all non-grad quantities to UVisu_FV"
  SDEALLOCATE(UVisu_FV)
  ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV,nVarVisu+nVarVisu_ElemData))
  DO iVar=1,nVarTotal
    IF (mapVisu(iVar)*(1-maskGrad(iVar)).GT.0) THEN
      SWRITE(*,*) "  ", TRIM(VarNamesTotal(iVar))
      UVisu_FV(:,:,:,:,mapVisu(iVar)) = UCalc_FV(:,:,:,:,mapCalc_FV(iVar))
    END IF
  END DO
  
  ! calc all grad quantities on the normal cell centers and convert them to visu grid
  IF (SUM(maskGrad*mapCalc_FV).GT.0) THEN
    SDEALLOCATE(UCalc_FV)
    ALLOCATE(UCalc_FV(0:PP_N,0:PP_N,0:PP_N,nElems_FV,1:nVarCalc_FV))
    SWRITE(*,*) "[FVRE] CalcQuantitiesWithGradients"
    CALL CalcQuantities(nVarCalc,PP_N,nElems_FV,mapElems_FV,mapCalc_FV,UCalc_FV,maskGrad)
    SWRITE(*,*) "[FVRE] ConvertToVisu_FV"
    CALL ConvertToVisu_FV(mapCalc_FV,maskGrad,reallocate=.FALSE.)
  END IF
#else
  maskCalc=1
  ! calc FV solution 
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:PP_N,0:PP_N,0:PP_N,nElems_FV,1:nVarCalc))
  CALL CalcQuantities(nVarCalc,PP_N,nElems_FV,mapElems_FV,mapCalc,UCalc_FV,withGradients,maskCalc) 
  ! convert FV solution to visu grid
  CALL ConvertToVisu_FV(mapCalc)
#endif /* FV_RECONSTRUCT */
END SUBROUTINE CalcQuantities_ConvertToVisu_FV

#endif /* FV_ENABLED */


!==================================================================================================================================
!> Copies variable for given element range from source array to target array using VTK structure
!==================================================================================================================================
PURE SUBROUTINE FillCopy(nVar,Nloc,nElems,UIn,nElems_calc,indices,UOut,iVar)
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
INTEGER,INTENT(IN) :: nVar
INTEGER,INTENT(IN) :: nElems
INTEGER,INTENT(IN) :: iVar
REAL,INTENT(IN)    :: UIn(nVar,0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(OUT)   :: UOut(    0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER         :: iElem,iElem_calc
!==================================================================================================================================
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
  UOut(:,:,:,iElem_calc) = UIn(iVar,:,:,:,iElem)
END DO ! iElem
END SUBROUTINE FillCopy

END MODULE MOD_Posti_Calc
