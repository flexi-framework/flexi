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
INTEGER             :: NSuper                    !< Polynomial degree used for supersampling of old mesh coordinates
REAL                :: maxTol                    !< max overshoot in param coords (1+maxTol), used for search and warnings
REAL                :: abortTol                  !< max tol before code aborts (=inf, if refstate specified)
REAL                :: displacement(3)           !< optional displacement vector for the old mesh
LOGICAL             :: printTroublemakers=.TRUE. !< print warnings for troublemakers
                                                 !< (very time consuming for many new points outisde the old mesh)
LOGICAL             :: ExtrudeTo3D               !< Switch to perform an extrusion of a one-layer mesh to the 3D version
INTEGER,ALLOCATABLE :: Elem_IJK(:,:)             !< IJK sorting of new mesh
INTEGER             :: nElemsOld_IJK(3)          !< IJK Element numebr of old mesh
INTEGER             :: nElemsNew_IJK(3)          !< IJK Element numebr of new mesh
INTEGER             :: ExtrudeK                  !< Layer which is used in extrusion

LOGICAL             :: ExtrudePeriodic           !< Perform a periodic extrusion of a 3D mesh to a mesh with extended z length

CHARACTER(LEN=255)  :: MeshFileOld               !< Old mesh file (optional, only to overwrite mesh from old state)
CHARACTER(LEN=255)  :: MeshFileNew               !< New mesh file
INTEGER             :: nElemsOld                 !< Number of elements in old mesh
INTEGER             :: nElemsNew                 !< Number of elements in new mesh
INTEGER             :: NGeoOld                   !< Polynomial degree of old mesh
INTEGER             :: NGeoNew                   !< Polynomial degree of new mesh
INTEGER             :: NState                    !< Polynomial degree of old state
INTEGER             :: NNew                      !< Polynomial degree of new state (=NState, if not specified otherwise)
INTEGER             :: NInter                    !< Polynomial degree for interpolation on new mesh (=NState, if not specified otherwise)
CHARACTER(LEN=255)  :: NodeTypeState             !< NodeType of the old state (Gauss/Gauss-Lobatto)
LOGICAL             :: useCurvedsOld             !< Should the old mesh use a curved mesh representation?
LOGICAL             :: useCurvedsNew             !< Should the new mesh use a curved mesh representation?
REAL,ALLOCATABLE    :: RefState(:)               !< Optional reference state used for points that can not be fond
INTEGER             :: nVar_State                !< Number of variables in DG_Solution array
REAL                :: Time_State                !< Output time of the current state file

! Vandermonde matrices
REAL,ALLOCATABLE    :: Vdm_CLNGeo_EquiNSuper(:,:) !< Vandermonde from CL mesh points to equdistant mesh points in supersampling
                                                  !< polynomial degree, used for old mesh points
REAL,ALLOCATABLE    :: Vdm_CLNInter_GPNNew(:,:)   !< Vandermonde from interpolation points on new mesh to solution points
REAL,ALLOCATABLE    :: Vdm_GPNState_GPNNew(:,:)   !< Vandermonde from old solution to new solution (used for equal elements)

REAL,ALLOCATABLE    :: xCLInter(:,:,:,:,:)  !> Mesh coordinates of new mesh represented as CL points on polynoial degree NInter for
                                            !> interpolation
REAL,ALLOCATABLE    :: xCLOld(:,:,:,:,:)    !> CL points of old mesh on NgeoOld
REAL,ALLOCATABLE    :: xCLNew(:,:,:,:,:)    !> CL points of new mesh on NgeoNew

REAL,ALLOCATABLE    :: xiInter(:,:,:,:,:)   !> Parametric coords of interpolation points of new solution in old mesh
INTEGER,ALLOCATABLE :: InterToElem(:,:,:,:) !> Mapping from new interpolation points  to elemID in old mesh
INTEGER,ALLOCATABLE :: equalElem(:)         !> Mapping between equal elements (iElemOld=equalElem(iElemNew))
LOGICAL,ALLOCATABLE :: IPDone(:,:,:,:)      !> Info if IP has been found

REAL,ALLOCATABLE    :: UOld(:,:,:,:,:)      !> Solution from old state

!===================================================================================================================================
END MODULE MOD_SwapMesh_Vars
