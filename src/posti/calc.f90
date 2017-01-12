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
#if PARABOLIC
USE MOD_Lifting_Vars   ,ONLY: gradUx,gradUy,gradUz
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nVal(4)
#if PARABOLIC
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: gradUx_tmp,gradUy_tmp,gradUz_tmp
#endif
!===================================================================================================================================
! calc DG solution 
SWRITE(*,*) "[DG] calc quantities"
SDEALLOCATE(UCalc_DG)
ALLOCATE(UCalc_DG(0:PP_N,0:PP_N,0:PP_N,nElems_DG,1:nVarCalc))
nVal=(/PP_N+1,PP_N+1,PP_N+1,nElems_DG/)

maskCalc = 1
! Copy exisiting variables from solution array
CALL FillCopy(nVar_State,PP_N,nElems,U,nElems_DG,mapElems_DG,UCalc_DG,maskCalc)

IF(TRIM(FileType).EQ.'State')THEN
  IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
    IF(nElems_DG.EQ.nElems)THEN
      CALL CalcQuantities(nVarCalc,nVal,mapElems_DG,mapCalc,UCalc_DG,maskCalc,gradUx,gradUy,gradUz) 
    ELSE
      ALLOCATE(gradUx_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems_DG))
      ALLOCATE(gradUy_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems_DG))
      ALLOCATE(gradUz_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems_DG))
      ! nicer, but only as of gfortran 6+: ALLOCATE(gradUx_tmp,gradUy_tmp,gradUz_tmp,MOLD=gradUx)
      gradUx_tmp=gradUx(:,:,:,:,mapElems_DG)
      gradUy_tmp=gradUy(:,:,:,:,mapElems_DG)
      gradUz_tmp=gradUz(:,:,:,:,mapElems_DG)
      CALL CalcQuantities(nVarCalc,nVal,mapElems_DG,mapCalc,UCalc_DG,maskCalc,gradUx_tmp,gradUy_tmp,gradUz_tmp) 
      DEALLOCATE(gradUx_tmp,gradUy_tmp,gradUz_tmp)
    END IF
#endif
  ELSE
    CALL CalcQuantities(nVarCalc,nVal,mapElems_DG,mapCalc,UCalc_DG,maskCalc) 
  END IF
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
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV_Reconstruct
#else
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_DG_Vars             ,ONLY: U
#endif
#if PARABOLIC
USE MOD_Lifting_Vars        ,ONLY: gradUx,gradUy,gradUz
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: maskCalc(nVarDep)
INTEGER                      :: nVal(4)
#if FV_RECONSTRUCT 
INTEGER                      :: maskPrim(nVarDep),maskGrad(nVarDep)
INTEGER                      :: iVar
#endif
#if PARABOLIC
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems_FV) :: gradUx_tmp,gradUy_tmp,gradUz_tmp
#endif
!===================================================================================================================================
SDEALLOCATE(UVisu_FV)
ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV,nVarVisuTotal))
#if FV_RECONSTRUCT
  ! generate a new mapCalc_FV, which is a copy of the original mapCalc but is extended in the following way.
  ! Since the reconstruction is performed in primitive quantities, the calculation of conservative quantities from them 
  ! introduce for the conservatives dependcies from the primitive ones. Therefore all primitive quantities that
  ! are needed to build the requested conservatives must be added to the mapCalc_FV. 
  nVal=(/NVisu_FV+1,NVisu_FV+1,NVisu_FV+1,nElems_FV/)
  SDEALLOCATE(mapCalc_FV)
  ALLOCATE(mapCalc_FV(1:nVarDep))
  CALL AppendNeededPrims(mapCalc,mapCalc_FV,nVarCalc_FV)
  SWRITE (*,*) "[FVRE] nVarCalc_FV", nVarCalc_FV
  SWRITE (*,"(A,27I3)") "  mapCalc_FV",  mapCalc_FV
  
  ! convert primitive quantities to the visu grid, but store it UCalc_FV, since all dependent calculations including
  ! reconstructed values are performed on the visu grid.
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV,1:nVarCalc_FV))
  SWRITE(*,*) "[FVRE] ConvertToVisu_FV_Reconstruct"
  CALL ConvertToVisu_FV_Reconstruct()

  ! calculate all nonPrim, nonGrad quantities on the visu grid.
  maskPrim=GetMaskPrim()
  maskGrad=GetMaskGrad()
  maskCalc = 1-MAX(maskPrim,maskGrad)
  SWRITE(*,*) "[FVRE] CalcQuantities (nonPrim,nonGrad)"
  CALL CalcQuantities(nVarCalc_FV,nVal,mapElems_FV,mapCalc_FV,UCalc_FV,maskCalc)

  ! copy cons, prim and all other non-grad quantities (namely all quantities that are build with reconstruction)
  ! to the UVisu_FV array
  SWRITE(*,*) "[FVRE] copy all non-grad quantities to UVisu_FV"
  DO iVar=1,nVarDep
    IF (mapVisu(iVar)*(1-maskGrad(iVar)).GT.0) THEN ! visu var, but no gradient var
      SWRITE(*,*) "  ", TRIM(VarNamesTotal(iVar))
      UVisu_FV(:,:,:,:,mapVisu(iVar)) = UCalc_FV(:,:,:,:,mapCalc_FV(iVar))
    END IF
  END DO
  
  ! calc all grad quantities on the normal cell centers and convert them afterwards to visu grid
  IF (SUM(maskGrad*mapCalc_FV).GT.0) THEN
    nVal=(/PP_N+1,PP_N+1,PP_N+1,nElems_FV/)
    SDEALLOCATE(UCalc_FV)
    ALLOCATE(UCalc_FV(0:PP_N,0:PP_N,0:PP_N,nElems_FV,1:nVarCalc_FV))
    SWRITE(*,*) "[FVRE] CalcQuantitiesWithGradients"
#if PARABOLIC
    gradUx_tmp=gradUx(:,:,:,:,mapElems_FV)
    gradUy_tmp=gradUy(:,:,:,:,mapElems_FV)
    gradUz_tmp=gradUz(:,:,:,:,mapElems_FV)
    CALL CalcQuantities(nVarCalc_FV,nVal,mapElems_FV,mapCalc_FV,UCalc_FV,maskGrad,gradUx_tmp,gradUy_tmp,gradUz_tmp)
#else
    CALL CalcQuantities(nVarCalc_FV,nVal,mapElems_FV,mapCalc_FV,UCalc_FV,maskGrad)
#endif
    SWRITE(*,*) "[FVRE] ConvertToVisu_FV"
    CALL ConvertToVisu_FV(mapCalc_FV,maskGrad)
  END IF
#else
  nVal=(/PP_N+1,PP_N+1,PP_N+1,nElems_FV/)
  ! calc FV solution 
  SWRITE(*,*) "[FV] calc quantities"
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:PP_N,0:PP_N,0:PP_N,nElems_FV,1:nVarCalc))

  maskCalc = 1
  ! Copy exisiting variables from solution array
  CALL FillCopy(nVar_State,PP_N,nElems,U,nElems_FV,mapElems_FV,UCalc_FV,maskCalc)

  IF(TRIM(FileType).EQ.'State')THEN
    IF(withDGOperator.AND.PARABOLIC)THEN
#if PARABOLIC
      gradUx_tmp=gradUx(:,:,:,:,mapElems_FV)
      gradUy_tmp=gradUy(:,:,:,:,mapElems_FV)
      gradUz_tmp=gradUz(:,:,:,:,mapElems_FV)
      CALL CalcQuantities(nVarCalc_FV,nVal,nElems_FV,mapElems_FV,mapCalc,UCalc_FV,maskCalc,gradUx_tmp,gradUy_tmp,gradUz_tmp)
#endif
    ELSE
      CALL CalcQuantities(nVarCalc_FV,nVal,nElems_FV,mapElems_FV,mapCalc,UCalc_FV,maskCalc) 
    END IF
  END IF

  ! convert FV solution to visu grid
  CALL ConvertToVisu_FV(mapCalc)
#endif /* FV_RECONSTRUCT */
END SUBROUTINE CalcQuantities_ConvertToVisu_FV

#endif /* FV_ENABLED */


!==================================================================================================================================
!> Copies variable for given element range from source array to target array using VTK structure
!==================================================================================================================================
SUBROUTINE FillCopy(nVar,Nloc,nElems,UIn,nElems_calc,indices,UOut,maskCalc)
USE MOD_Posti_Vars
USE MOD_StringTools ,ONLY: STRICMP
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nVar
INTEGER,INTENT(IN)    :: Nloc
INTEGER,INTENT(IN)    :: nElems
INTEGER,INTENT(IN)    :: nElems_calc
INTEGER,INTENT(IN)    :: indices(nElems_calc)
REAL,INTENT(IN)       :: UIn(nVar,0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(OUT)      :: UOut(0:Nloc,0:Nloc,0:Nloc,nElems_calc,nVarCalc)
INTEGER,INTENT(INOUT) :: maskCalc(nVarDep)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER               :: iVarOut,iVarIn
INTEGER               :: iElem,iElem_calc
!==================================================================================================================================
! Copy exisiting variables from solution array
DO iVarOut=1,nVarDep ! iterate over all out variables
  IF (mapCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarIn=1,nVar_State ! iterate over all in variables
    IF( STRICMP(VarNamesTotal(iVarOut),VarNamesHDF5(iVarIn))) THEN
      DO iElem_calc=1,nElems_calc ! copy variable for all elements
        iElem = indices(iElem_calc)
        UOut(:,:,:,iElem_calc,mapCalc(iVarOut)) = UIn(iVarIn,:,:,:,iElem)
      END DO ! iElem
      maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
    END IF
  END DO
END DO
END SUBROUTINE FillCopy


END MODULE MOD_Posti_Calc
