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

!===================================================================================================================================
!> Contains global variables provided by the swapmesh routines
!===================================================================================================================================
MODULE MOD_SwapMesh_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: NSuper
REAL                :: maxTol   ! max overshoot in param coords (1+maxTol), used for search and warnings
REAL                :: abortTol ! max tol before code aborts (=inf, if refstate specified)
REAL                :: displacement(3)
LOGICAL             :: InterchangeYZ=.FALSE.
LOGICAL             :: printTroublemakers=.TRUE. ! print warnings for troublemakers 
                                                 ! (very time consuming for many new points outisde the old mesh)

CHARACTER(LEN=255)  :: MeshFileOld,MeshFileNew
INTEGER             :: nElemsOld,  nElemsNew
INTEGER             :: OffsetElemOld,  OffsetElemNew
INTEGER             :: NGeoOld,    NGeoNew
INTEGER             :: NState ! N of state we are interpolating
INTEGER             :: NNew   ! N of new state      (=NState, if not specified otherwise)
INTEGER             :: NInter ! N for interpolation (=NState, if not specified otherwise)
CHARACTER(LEN=255)  :: NodeTypeState
LOGICAL             :: useCurvedsOld,useCurvedsNew
REAL,ALLOCATABLE    :: RefState(:)
INTEGER             :: nVar_State
REAL                :: Time_State


REAL,ALLOCATABLE    :: Vdm_CLNGeo_EquiNSuper(:,:)
REAL,ALLOCATABLE    :: Vdm_CLNInter_GPNNew(:,:)
REAL,ALLOCATABLE    :: Vdm_GPNState_GPNNew(:,:)

REAL,ALLOCATABLE    :: xCLInter(:,:,:,:,:)  ! CL points of interpolation mesh (from new mesh)
REAL,ALLOCATABLE    :: xCLOld(:,:,:,:,:)    ! CL points of old mesh
REAL,ALLOCATABLE    :: xCLNew(:,:,:,:,:)    ! CL points of new mesh

REAL,ALLOCATABLE    :: xiInter(:,:,:,:,:)   ! Parametric coords of new solution xCLInter in old mesh
INTEGER,ALLOCATABLE :: InterToElem(:,:,:,:) ! Mapping from new solution GP to elemID in old mesh
INTEGER,ALLOCATABLE :: equalElem(:)         ! Mapping between equal elements (iElemOld=equalElem(iElemNew)
LOGICAL,ALLOCATABLE :: IPDone(:,:,:,:)      ! Info if IP has been found

REAL,ALLOCATABLE    :: UOld(:,:,:,:,:)      ! Solution from old state

!===================================================================================================================================
END MODULE MOD_SwapMesh_Vars
