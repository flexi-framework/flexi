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
USE MOD_EOS_Posti           ,ONLY: CalcQuantities
IMPLICIT NONE
!===================================================================================================================================
! calc DG solution 
SWRITE(*,*) "[DG] calc quantities"
SDEALLOCATE(UCalc_DG)
ALLOCATE(UCalc_DG(0:PP_N,0:PP_N,0:PP_N,nElems_DG,1:nVarCalc))
CALL CalcQuantities(nVarCalc,PP_N,nElems_DG,mapElems_DG,mapCalc,UCalc_DG,withGradients) 
END SUBROUTINE CalcQuantities_DG


#if FV_ENABLED 
!===================================================================================================================================
!> Calc quantities for all FV elements.
!> TODO
!===================================================================================================================================
SUBROUTINE CalcQuantities_ConvertToVisu_FV() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_EOS_Posti_Vars      ,ONLY: nVarTotal,DepNames
USE MOD_EOS_Posti           ,ONLY: GetMaskCons,GetMaskPrim,GetMaskGrad
USE MOD_EOS_Posti           ,ONLY: CalcQuantities
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV
#if FV_RECONSTRUCT
USE MOD_EOS_Posti           ,ONLY: AppendNeededPrims
USE MOD_EOS_Posti           ,ONLY: CalcConsFromPrim
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV_Reconstruct
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if FV_ENABLED && FV_RECONSTRUCT 
INTEGER                      :: maskCons(nVarTotal),maskPrim(nVarTotal),maskGrad(nVarTotal)
INTEGER                      :: maskNonConsPrim(nVarTotal)
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
  
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV,1:nVarCalc_FV))
  SWRITE(*,*) "[FVRE] ConvertToVisu_FV_Reconstruct"
  CALL ConvertToVisu_FV_Reconstruct() 

  ! calculate all needed conservative variables
  SWRITE(*,*) "[FVRE] CalcConsFromPrim"
  CALL CalcConsFromPrim(mapCalc_FV,nVarCalc_FV,NVisu_FV,nElems_FV,UCalc_FV)

  maskCons=GetMaskCons()
  maskPrim=GetMaskPrim()
  maskNonConsPrim = 1-MAX(maskCons,maskPrim)
  SWRITE(*,*) "[FVRE] CalcQuantities (nonCons,nonPrim)"
  CALL CalcQuantities(nVarCalc_FV,NVisu_FV,nElems_FV,mapElems_FV,mapCalc_FV,UCalc_FV,withGradients,maskNonConsPrim)

  SWRITE(*,*) "[FVRE] copy all non-grad quantities to UVisu_FV"
  maskGrad = GetMaskGrad()
  SDEALLOCATE(UVisu_FV)
  ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV,nVarVisu+nVarVisu_ElemData))
  DO iVar=1,nVarTotal
    IF (mapVisu(iVar)*(1-maskGrad(iVar)).GT.0) THEN
      SWRITE(*,*) "  ", TRIM(DepNames(iVar))
      UVisu_FV(:,:,:,:,mapVisu(iVar)) = UCalc_FV(:,:,:,:,mapCalc_FV(iVar))
    END IF
  END DO
  
  IF (SUM(maskGrad*mapCalc_FV).GT.0) THEN
    SDEALLOCATE(UCalc_FV)
    ALLOCATE(UCalc_FV(0:PP_N,0:PP_N,0:PP_N,nElems_FV,1:nVarCalc_FV))
    SWRITE(*,*) "[FVRE] CalcConsFromPrim"
    CALL CalcConsFromPrim(mapCalc_FV,0,0,0)
    SWRITE(*,*) "[FVRE] CalcQuantitiesWithGradients"
    CALL CalcQuantities(nVarCalc,PP_N,nElems_FV,mapElems_FV,mapCalc_FV,UCalc_FV,withGradients,maskGrad) 
    SWRITE(*,*) "[FVRE] ConvertToVisu_FV"
    CALL ConvertToVisu_FV(mapCalc_FV,maskGrad,reallocate=.FALSE.)
  END IF
#else
  ! calc FV solution 
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:PP_N,0:PP_N,0:PP_N,nElems_FV,1:nVarCalc))
  CALL CalcQuantities(nVarCalc,PP_N,nElems_FV,mapElems_FV,mapCalc,UCalc_FV,withGradients) 
  ! convert FV solution to visu grid
  CALL ConvertToVisu_FV(mapCalc)
#endif /* FV_RECONSTRUCT */
END SUBROUTINE CalcQuantities_ConvertToVisu_FV

#endif /* FV_ENABLED */

END MODULE MOD_Posti_Calc
