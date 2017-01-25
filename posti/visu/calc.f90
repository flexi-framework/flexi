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

INTERFACE CalcSurfQuantities_DG
  MODULE PROCEDURE CalcSurfQuantities_DG
END INTERFACE
PUBLIC:: CalcSurfQuantities_DG

#if FV_ENABLED 
INTERFACE CalcQuantities_FV
  MODULE PROCEDURE CalcQuantities_FV
END INTERFACE
PUBLIC:: CalcQuantities_FV

INTERFACE CalcSurfQuantities_FV
  MODULE PROCEDURE CalcSurfQuantities_FV
END INTERFACE
PUBLIC:: CalcSurfQuantities_FV
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
USE MOD_Visu_Vars
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_DG_Vars            ,ONLY: U
USE MOD_StringTools        ,ONLY: STRICMP
#if PARABOLIC
USE MOD_Lifting_Vars       ,ONLY: gradUx,gradUy,gradUz
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
CALL FillCopy(nVar_State,PP_N,nElems,U,nElems_DG,mapDGElemsToAllElems,UCalc_DG,maskCalc)

IF(TRIM(FileType).EQ.'State')THEN
  IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
    IF(nElems_DG.EQ.nElems)THEN
      CALL CalcQuantities(nVarCalc,nVal,mapDGElemsToAllElems,mapDepToCalc,UCalc_DG,maskCalc,gradUx,gradUy,gradUz)
    ELSE
      ALLOCATE(gradUx_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems_DG))
      ALLOCATE(gradUy_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems_DG))
      ALLOCATE(gradUz_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems_DG))
      ! nicer, but only as of gfortran 6+: ALLOCATE(gradUx_tmp,gradUy_tmp,gradUz_tmp,MOLD=gradUx)
      gradUx_tmp=gradUx(:,:,:,:,mapDGElemsToAllElems)
      gradUy_tmp=gradUy(:,:,:,:,mapDGElemsToAllElems)
      gradUz_tmp=gradUz(:,:,:,:,mapDGElemsToAllElems)
      CALL CalcQuantities(nVarCalc,nVal,mapDGElemsToAllElems,mapDepToCalc,UCalc_DG,maskCalc,gradUx_tmp,gradUy_tmp,gradUz_tmp) 
      DEALLOCATE(gradUx_tmp,gradUy_tmp,gradUz_tmp)
    END IF
#endif
  ELSE
    CALL CalcQuantities(nVarCalc,nVal,mapDGElemsToAllElems,mapDepToCalc,UCalc_DG,maskCalc) 
  END IF
END IF
END SUBROUTINE CalcQuantities_DG


!===================================================================================================================================
!> Calc surface quantities for all DG elements.
!> 1. Prolong all independent quantities to the boudary sides that should be visualized, also prolong gradients if they are needed
!> 2. Call CalcQuantities for the surface.
!>
!> This means we only prolong the conservative or additional variables to the boundary. All dependent variables will be calculated
!> on the surface and not prolonged so we don't get interpolation errors.
!===================================================================================================================================
SUBROUTINE CalcSurfQuantities_DG() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_Mesh_Vars          ,ONLY: nBCSides
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2
USE MOD_StringTools        ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nValSide(3)
INTEGER            :: iSide,iSide2
#if PARABOLIC
REAL,ALLOCATABLE   :: gradUxFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUyFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUzFace(:,:,:,:)
#endif
REAL,ALLOCATABLE   :: NormVec_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec1_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec2_loc(:,:,:,:)
!===================================================================================================================================
!------ Surface visualization ----------!
nValSide=(/PP_N+1,PP_N+1,nBCSidesVisu_DG/)

! Allocate array that stores the calculated variables on the visualization boundary.
SDEALLOCATE(USurfCalc_DG)
ALLOCATE(USurfCalc_DG(0:PP_N,0:PP_N,nBCSidesVisu_DG,1:nVarCalc))
#if PARABOLIC
ALLOCATE(gradUxFace(1:PP_nVarPrim,0:PP_N,0:PP_N,nBCSidesVisu_DG))
ALLOCATE(gradUyFace(1:PP_nVarPrim,0:PP_N,0:PP_N,nBCSidesVisu_DG))
ALLOCATE(gradUzFace(1:PP_nVarPrim,0:PP_N,0:PP_N,nBCSidesVisu_DG))
#endif
ALLOCATE(NormVec_loc (1:3,0:PP_N,0:PP_N,nBCSidesVisu_DG))
ALLOCATE(TangVec1_loc(1:3,0:PP_N,0:PP_N,nBCSidesVisu_DG))
ALLOCATE(TangVec2_loc(1:3,0:PP_N,0:PP_N,nBCSidesVisu_DG))

maskCalc=1
CALL ProlongToFace_independent(nVarCalc,nBCSidesVisu_DG,nElems_DG,maskCalc,UCalc_DG,USurfCalc_DG &
#if PARABOLIC    
    ,gradUxFace,gradUyFace,gradUzFace &
#endif    
    )

IF(TRIM(FileType).EQ.'State')THEN
  DO iSide=1,nBCSides
    iSide2 = mapAllBCSidesToDGVisuBCSides(iSide)
    IF (iSide2.GT.0) THEN
      NormVec_loc (:,:,:,iSide2) = NormVec (:,:,:,0,iSide)
      TangVec1_loc(:,:,:,iSide2) = TangVec1(:,:,:,0,iSide)
      TangVec2_loc(:,:,:,iSide2) = TangVec2(:,:,:,0,iSide)
    END IF
  END DO
  IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
    CALL CalcQuantities(nVarCalc,nValSide,mapAllBCSidesToDGVisuBCSides,mapDepToCalc,USurfCalc_DG,maskCalc*(1-DepVolumeOnly),&
        gradUxFace,gradUyFace,gradUzFace,&
        NormVec_loc(:,:,:,:),TangVec1_loc(:,:,:,:),TangVec2_loc(:,:,:,:)) 
#endif
  ELSE
    CALL CalcQuantities(nVarCalc,nValSide,mapAllBCSidesToDGVisuBCSides,mapDepToCalc,USurfCalc_DG,maskCalc*(1-DepVolumeOnly),& 
        NormVec=NormVec_loc(:,:,:,:),TangVec1=TangVec1_loc(:,:,:,:),TangVec2=TangVec2_loc(:,:,:,:)) 
  END IF
END IF

#if PARABOLIC
DEALLOCATE(gradUxFace)
DEALLOCATE(gradUyFace)
DEALLOCATE(gradUzFace)
#endif

END SUBROUTINE CalcSurfQuantities_DG

!===================================================================================================================================
!> Prolong all independent variables (defined as the variables in the state file, e.g. conservative variables for normal
!> Navier Stokes calculations as well as variables that can only be calculated in the volume and must be prolonged anyway)
!> to the faces that should be visualized. If the gradients are needed, they are also prolonged.
!===================================================================================================================================
SUBROUTINE ProlongToFace_independent(nVar,nSides_calc,nElems_calc,maskCalc,UIn,UBoundary&
#if PARABOLIC
    ,gradUxFace,gradUyFace,gradUzFace &
#endif        
    ) 
USE MOD_PreProc
USE MOD_Visu_Vars
#if PARABOLIC
USE MOD_Lifting_Vars       ,ONLY: gradUx,gradUy,gradUz
USE MOD_Mesh_Vars          ,ONLY: SideToElem
#endif
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_Mesh_Vars          ,ONLY: nBCSides,S2V2,ElemToSide
USE MOD_ProlongToFace      ,ONLY: EvalElemFace
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)            :: nVar,nSides_calc,nElems_calc
REAL,INTENT(IN)               :: UIn(0:PP_N,0:PP_N,0:PP_N,nElems_calc,1:nVar)
REAL,INTENT(OUT)              :: UBoundary(0:PP_N,0:PP_N,nSides_calc,1:nVar)
#if PARABOLIC
REAL,INTENT(OUT)              :: gradUxFace(1:PP_nVarPrim,0:PP_N,0:PP_N,nSides_calc)
REAL,INTENT(OUT)              :: gradUyFace(1:PP_nVarPrim,0:PP_N,0:PP_N,nSides_calc)
REAL,INTENT(OUT)              :: gradUzFace(1:PP_nVarPrim,0:PP_N,0:PP_N,nSides_calc)
#endif
INTEGER,INTENT(INOUT)         :: maskCalc(nVarDep)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar,iVarIn,iVarOut,iSide,locSide,iElem,p,q,iElem_DG,iSide_DG
REAL              :: Uface(1,0:PP_N,0:PP_N)
#if PARABOLIC
REAL              :: gradUxFace_tmp(1:PP_nVarPrim,0:PP_N,0:PP_N)
REAL              :: gradUyFace_tmp(1:PP_nVarPrim,0:PP_N,0:PP_N)
REAL              :: gradUzFace_tmp(1:PP_nVarPrim,0:PP_N,0:PP_N)
#endif
!===================================================================================================================================
! Loop over all dependent vairbales
DO iVarOut=1,nVarDep ! iterate over all out variables
  IF (mapDepToCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarIn=1,nVar_State ! iterate over all out variables
    ! Check if this variable is present in the state file, if so define it as independent
    ! Also define variables as independent that can only be computed in the volume
    IF ((STRICMP(VarnamesAll(iVarOut),VarNamesHDF5(iVarIn))).OR.(DepVolumeOnly(iVarOut).EQ.1)) THEN
      WRITE(*,*) "ProlongToFace_independent", TRIM(VarnamesAll(iVarOut))
      iVar=mapDepToCalc(iVarOut)

      DO iElem_DG = 1,nElems_DG                         ! iterate over all DG visu elements
        iElem = mapDGElemsToAllElems(iElem_DG)          ! get global element index
        DO locSide=1,6 
          iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
          IF (iSide.LE.nBCSides) THEN                   ! check if BC side
            iSide_DG = mapAllBCSidesToDGVisuBCSides(iSide)  ! get DG visu side index
            IF (iSide_DG.GT.0) THEN
              IF(PP_NodeType.EQ.1)THEN                  ! prolong solution to face
                CALL EvalElemFace(1,PP_N,UIn(:,:,:,iElem_DG,iVar:iVar),Uface(1:1,:,:),L_Minus,L_Plus,locSide)
              ELSE
                CALL EvalElemFace(1,PP_N,UIn(:,:,:,iElem_DG,iVar:iVar),Uface(1:1,:,:),locSide)
              END IF
              ! Rotate into coordinate system of master side
              DO q=0,PP_N; DO p=0,PP_N
                UBoundary(p,q,iSide_DG,iVar)=Uface(1,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
              END DO; END DO
            END IF   
          END IF
        END DO
      END DO
      maskCalc(iVarOut) = 0
      EXIT
    END IF
  END DO
END DO


! Also prolong the gradients if parabolic terms are needed and the DG operator is called once.
IF(TRIM(FileType).EQ.'State')THEN
  IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
    DO iSide=1,nBCSides
      IF (mapAllBCSidesToDGVisuBCSides(iSide).GT.0) THEN
        iElem = SideToElem(S2E_ELEM_ID,iSide)
        locSide = SideToElem(S2E_LOC_SIDE_ID, iSide)
        IF(PP_NodeType.EQ.1)THEN
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUx(:,:,:,:,iElem),gradUxFace_tmp,L_Minus,L_Plus,locSide)
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUy(:,:,:,:,iElem),gradUyFace_tmp,L_Minus,L_Plus,locSide)
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUz(:,:,:,:,iElem),gradUzFace_tmp,L_Minus,L_Plus,locSide)
        ELSE
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUx(:,:,:,:,iElem),gradUxFace_tmp,locSide)
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUy(:,:,:,:,iElem),gradUyFace_tmp,locSide)
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUz(:,:,:,:,iElem),gradUzFace_tmp,locSide)
        END IF
        iSide_DG = mapAllBCSidesToDGVisuBCSides(iSide)
        DO q=0,PP_N; DO p=0,PP_N
          gradUxFace(:,p,q,iSide_DG)=gradUxFace_tmp(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
          gradUyFace(:,p,q,iSide_DG)=gradUyFace_tmp(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
          gradUzFace(:,p,q,iSide_DG)=gradUzFace_tmp(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
        END DO; END DO
      END IF
    END DO
#endif
  END IF
END IF
END SUBROUTINE ProlongToFace_independent

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
SUBROUTINE CalcQuantities_FV() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV
USE MOD_EOS_Posti           ,ONLY: CalcQuantities
#if FV_RECONSTRUCT
USE MOD_EOS_Posti           ,ONLY: GetMaskCons,GetMaskPrim,GetMaskGrad
USE MOD_EOS_Posti           ,ONLY: AppendNeededPrims
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV_Reconstruct
#else
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_DG_Vars             ,ONLY: U
#if PARABOLIC
USE MOD_Lifting_Vars        ,ONLY: gradUx,gradUy,gradUz
#endif
#endif
USE MOD_StringTools         ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: maskCalc(nVarDep)
INTEGER                      :: nVal(4)
#if PARABOLIC
REAL,ALLOCATABLE             :: gradUx_calc(:,:,:,:,:),gradUy_calc(:,:,:,:,:),gradUz_calc(:,:,:,:,:) 
#endif
!===================================================================================================================================

#if FV_RECONSTRUCT
  ! ================================ WITH RECONSTRUCTION ======================================

  nVal=(/NCalc_FV+1,NCalc_FV+1,NCalc_FV+1,nElems_FV/)
  SWRITE (*,*) "[FVRE] nVarCalc_FV", nVarCalc_FV
  
  ! convert primitive quantities to the visu grid, but store them in UCalc_FV, since all dependent calculations based on
  ! reconstructed values are performed on the visu grid.
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:NCalc_FV,0:NCalc_FV,0:NCalc_FV,nElems_FV,1:nVarCalc_FV))
  SWRITE(*,*) "[FVRE] ConvertToVisu_FV_Reconstruct"

  ! calculate all remaining quantities on the visu grid.
  IF(TRIM(FileType).EQ.'State')THEN
    maskCalc = 1-GetMaskPrim()   ! exclude primitive quantities
    SWRITE(*,*) "[FVRE] CalcQuantities (nonPrim)"
    IF(withDGOperator.AND.(PARABOLIC.EQ.1))THEN
#if PARABOLIC    
      ALLOCATE(gradUx_calc(1:PP_nVarPrim,0:NCalc_FV,0:NCalc_FV,0:NCalc_FV,nElems_FV))
      ALLOCATE(gradUy_calc(1:PP_nVarPrim,0:NCalc_FV,0:NCalc_FV,0:NCalc_FV,nElems_FV))
      ALLOCATE(gradUz_calc(1:PP_nVarPrim,0:NCalc_FV,0:NCalc_FV,0:NCalc_FV,nElems_FV))
      CALL ConvertToVisu_FV_Reconstruct(gradUx_calc,gradUy_calc,gradUz_calc)
      CALL CalcQuantities(nVarCalc_FV,nVal,mapFVElemsToAllElems,mapDepToCalc_FV,UCalc_FV,maskCalc,gradUx_calc,gradUy_calc,gradUz_calc)
#endif
    ELSE
      CALL ConvertToVisu_FV_Reconstruct()
      CALL CalcQuantities(nVarCalc_FV,nVal,mapFVElemsToAllElems,mapDepToCalc_FV,UCalc_FV,maskCalc)
    END IF
  END IF


#else /* FV_RECONSTRUCT */

  ! ================================ WITHOUT RECONSTRUCTION ======================================
  ! Since no reconstruction is involved, we can calculate all dependent quantities from the conservative solution (same
  ! as for DG). Without reconstruction this can be done on the FV cell-centers (PP_N instead of NVisu_FV).
  nVal=(/NCalc_FV+1,NCalc_FV+1,NCalc_FV+1,nElems_FV/)
  SDEALLOCATE(mapDepToCalc_FV)
  ALLOCATE(mapDepToCalc_FV(1:nVarDep))
  mapDepToCalc_FV = mapDepToCalc
  ! calc FV solution 
  SWRITE(*,*) "[FV] calc quantities"
  SDEALLOCATE(UCalc_FV)
  ALLOCATE(UCalc_FV(0:NCalc_FV,0:NCalc_FV,0:NCalc_FV,nElems_FV,1:nVarCalc))

  maskCalc = 1
  ! Copy exisiting variables from solution array
  CALL FillCopy(nVar_State,NCalc_FV,nElems,U,nElems_FV,mapFVElemsToAllElems,UCalc_FV,maskCalc)

  IF(TRIM(FileType).EQ.'State')THEN
    CALL CalcQuantities(nVarCalc_FV,nVal,mapFVElemsToAllElems,mapDepToCalc_FV,UCalc_FV,maskCalc) 
  END IF

#endif /* FV_RECONSTRUCT */
END SUBROUTINE CalcQuantities_FV

!===================================================================================================================================
!> Calc surface quantities for all FV elements.
!> For this, the already calculated quantities as well as the gradients in the volume will simply be copied to to the side. 
!> Also the mesh quantities like normal vectors are prepared and then the CalcQuantities routine is called.
!===================================================================================================================================
SUBROUTINE CalcSurfQuantities_FV() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_Mesh_Vars          ,ONLY: nBCSides,SideToElem,ElemToSide
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2
USE MOD_StringTools        ,ONLY: STRICMP
#if PARABOLIC
USE MOD_Mesh_Vars          ,ONLY: S2V
USE MOD_Lifting_Vars       ,ONLY: gradUx,gradUy,gradUz
#endif
USE MOD_Mappings           ,ONLY: buildMappings
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nValSide(3)
INTEGER            :: iSide,iSide_FV,iElem,iElem_FV,p,q,dir,ijk(3),locSide
INTEGER,ALLOCATABLE          :: S2V_NVisu(:,:,:,:,:,:)
#if PARABOLIC
INTEGER            :: iVar
REAL,ALLOCATABLE   :: gradUxFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUyFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUzFace(:,:,:,:)
#endif
REAL,ALLOCATABLE   :: NormVec_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec1_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec2_loc(:,:,:,:)
!===================================================================================================================================


nValSide=(/NCalc_FV+1,NCalc_FV+1,nBCSidesVisu_FV/)
CALL buildMappings(NCalc_FV,S2V=S2V_NVisu)
SDEALLOCATE(USurfVisu_FV)
ALLOCATE(USurfVisu_FV(0:NCalc_FV,0:NCalc_FV,0:0,nBCSidesVisu_FV,nVarSurfVisuAll))
! ===  Surface visualization ================================
! copy UCalc_FV to USurfCalc_FV
SDEALLOCATE(USurfCalc_FV)
ALLOCATE(USurfCalc_FV(0:NCalc_FV,0:NCalc_FV,nBCSidesVisu_FV,1:nVarCalc_FV))
DO iElem_FV = 1,nElems_FV                         ! iterate over all FV visu elements
  iElem = mapFVElemsToAllElems(iElem_FV)          ! get global element index
  DO locSide=1,6 
    iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
    IF (iSide.LE.nBCSides) THEN                   ! check if BC side
      iSide_FV = mapAllBCSidesToFVVisuBCSides(iSide)  ! get FV visu side index
      IF (iSide_FV.GT.0) THEN
        DO q=0,NCalc_FV; DO p=0,NCalc_FV          ! map volume solution to surface solution
          ijk = S2V_NVisu(:,0,p,q,0,locSide)
          USurfCalc_FV(p,q,iSide_FV,:) = UCalc_FV(ijk(1),ijk(2),ijk(3),iElem_FV,:)
        END DO; END DO
      END IF   
    END IF
  END DO
END DO

ALLOCATE(NormVec_loc (1:3,0:NCalc_FV,0:NCalc_FV,nBCSidesVisu_FV))
ALLOCATE(TangVec1_loc(1:3,0:NCalc_FV,0:NCalc_FV,nBCSidesVisu_FV))
ALLOCATE(TangVec2_loc(1:3,0:NCalc_FV,0:NCalc_FV,nBCSidesVisu_FV))
#if PARABOLIC
ALLOCATE(gradUxFace(1:PP_nVarPrim,0:NCalc_FV,0:NCalc_FV,nBCSidesVisu_FV))
ALLOCATE(gradUyFace(1:PP_nVarPrim,0:NCalc_FV,0:NCalc_FV,nBCSidesVisu_FV))
ALLOCATE(gradUzFace(1:PP_nVarPrim,0:NCalc_FV,0:NCalc_FV,nBCSidesVisu_FV))
#endif
DO iSide=1,nBCSides
  iSide_FV = mapAllBCSidesToFVVisuBCSides(iSide)
  iElem = SideToElem(S2E_ELEM_ID,iSide)
  locSide = SideToElem(S2E_LOC_SIDE_ID, iSide)
  IF (iSide_FV.GT.0) THEN
    DO q=0,PP_N; DO p=0,PP_N
      DO dir=1,3
        NormVec_loc (dir,p*2:p*2+1,q*2:q*2+1,iSide_FV) = NormVec (dir,p,q,1,iSide)
        TangVec1_loc(dir,p*2:p*2+1,q*2:q*2+1,iSide_FV) = TangVec1(dir,p,q,1,iSide)
        TangVec2_loc(dir,p*2:p*2+1,q*2:q*2+1,iSide_FV) = TangVec2(dir,p,q,1,iSide)
      END DO
#if PARABOLIC
      ijk = S2V(:,0,p,q,0,locSide)
      DO iVar=1,PP_nVarPrim
        gradUxFace(iVar,p*2:p*2+1,q*2:q*2+1,iSide_FV) = gradUx(iVar,ijk(1),ijk(2),ijk(3),iElem)
        gradUyFace(iVar,p*2:p*2+1,q*2:q*2+1,iSide_FV) = gradUy(iVar,ijk(1),ijk(2),ijk(3),iElem)
        gradUzFace(iVar,p*2:p*2+1,q*2:q*2+1,iSide_FV) = gradUz(iVar,ijk(1),ijk(2),ijk(3),iElem)
      END DO 
#endif
    END DO; END DO ! p,q=0,PP_N
  END IF
END DO


SWRITE(*,*) "[FVRE] CalcSurfQuantities"
maskCalc = DepSurfaceOnly
IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
  CALL CalcQuantities(nVarCalc_FV,nValSide,mapAllBCSidesToFVVisuBCSides,mapDepToCalc_FV,USurfCalc_FV,maskCalc*(1-DepVolumeOnly),&
      gradUxFace,gradUyFace,gradUzFace,&
      NormVec_loc(:,:,:,:),TangVec1_loc(:,:,:,:),TangVec2_loc(:,:,:,:)) 
#endif
ELSE
  CALL CalcQuantities(nVarCalc_FV,nValSide,mapAllBCSidesToFVVisuBCSides,mapDepToCalc_FV,USurfCalc_FV,maskCalc*(1-DepVolumeOnly),& 
      NormVec=NormVec_loc(:,:,:,:),TangVec1=TangVec1_loc(:,:,:,:),TangVec2=TangVec2_loc(:,:,:,:)) 
END IF

END SUBROUTINE CalcSurfQuantities_FV

#endif /* FV_ENABLED */


!==================================================================================================================================
!> Copies variable for given element range from source array to target array using VTK structure
!==================================================================================================================================
SUBROUTINE FillCopy(nVar,Nloc,nElems,UIn,nElems_calc,indices,UOut,maskCalc)
USE MOD_Visu_Vars
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
  IF (mapDepToCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarIn=1,nVar_State ! iterate over all in variables
    IF( STRICMP(VarnamesAll(iVarOut),VarNamesHDF5(iVarIn))) THEN
      DO iElem_calc=1,nElems_calc ! copy variable for all elements
        iElem = indices(iElem_calc)
        UOut(:,:,:,iElem_calc,mapDepToCalc(iVarOut)) = UIn(iVarIn,:,:,:,iElem)
      END DO ! iElem
      maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
    END IF
  END DO
END DO
END SUBROUTINE FillCopy


END MODULE MOD_Posti_Calc
