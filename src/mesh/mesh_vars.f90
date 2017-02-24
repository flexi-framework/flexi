!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz 
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

!==================================================================================================================================
!> \brief Contains global variables provided by the mesh routines,

!> Contains mesh related variables and data-structures
!> - Gauss-point coordinates (provided by MOD_Metrics)
!> - Surface and volume metrics (provided by MOD_Metrics)
!> - Inner- and inter-element mappings in verious reference frames (provided by MOD_Mappings and MOD_Prepare_Mesh)
!> - support data structures and routines for mesh readin (provided by MOD_MeshReadin)
!==================================================================================================================================
MODULE MOD_Mesh_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! basis
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER           :: NGeo                      !< polynomial degree of geometric transformation
INTEGER           :: NGeoRef                   !< polynomial degree of reference jacobian
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE,TARGET :: NodeCoords(:,:,:,:,:) !< XYZ positions (equidistant,NGeo) of element interpolation points from meshfile
REAL,ALLOCATABLE,TARGET :: Elem_xGP(:,:,:,:,:)          !< XYZ positions (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE :: Elem_xGPO(:,:,:,:,:)         !< XYZ positions (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE :: Face_xGP(:,:,:,:,:)          !< XYZ positions (first index 1:3) of the Face Gauss Point
REAL,ALLOCATABLE :: Face_xGPO(:,:,:,:,:)         !< XYZ positions (first index 1:3) of the Face Gauss Point
!----------------------------------------------------------------------------------------------------------------------------------
! MORTAR DATA FOR NON-CONFORMING MESHES ORIGINATING FROM AN OCTREE BASIS (ONLY ALLOCATED IF isMortarMesh=.TRUE.!!!)
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: isMortarMesh               !< Marker whether non-conforming data is present (false for conforming meshes)
LOGICAL          :: interpolateFromTree        !< Switch whether to build metrics on tree level and interpolate to elements.
                                               !< Only applicable if tree data is present in mesh file
REAL,ALLOCATABLE,TARGET :: TreeCoords(:,:,:,:,:) !< XYZ positions (equidistant,NGeoTree) of tree interpolation points from meshfile
REAL,ALLOCATABLE :: xiMinMax(:,:,:)            !< Position of the 2 bounding nodes of a quadrant in its tree
INTEGER          :: NGeoTree                   !< Polynomial degree of trees geometric transformation
INTEGER          :: nTrees                     !< Local number of trees in mesh
INTEGER          :: nGlobalTrees               !< Global number of trees in mesh
INTEGER          :: offsetTree                 !< Tree offset (for MPI)
INTEGER,ALLOCATABLE :: ElemToTree(:)           !< Index of the tree corresponding to an element
!----------------------------------------------------------------------------------------------------------------------------------
! Metrics on GaussPoints
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: dXCL_N(:,:,:,:,:,:)          !< geometry Jacobian matrix on in reference space
REAL,ALLOCATABLE :: Metrics_fTilde(:,:,:,:,:,:)  !< Metrics for transforming the fluxes f (1:3,0:N,0:N,0:N,nElems,0:FV)
REAL,ALLOCATABLE :: Metrics_gTilde(:,:,:,:,:,:)  !< Metrics for transforming the fluxes g (1:3,0:N,0:N,0:N,nElems,0:FV)
REAL,ALLOCATABLE :: Metrics_hTilde(:,:,:,:,:,:)  !< Metrics for transforming the fluxes h (1:3,0:N,0:N,0:N,nElems,0:FV)
REAL,ALLOCATABLE :: Metrics_fTildeO(:,:,:,:,:)   !< Metrics for transforming the fluxes f for overintegration (1:3,0:M,0:M,0:M,nElems)
REAL,ALLOCATABLE :: Metrics_gTildeO(:,:,:,:,:)   !< Metrics for transforming the fluxes g for overintegration (1:3,0:M,0:M,0:M,nElems)
REAL,ALLOCATABLE :: Metrics_hTildeO(:,:,:,:,:)   !< Metrics for transforming the fluxes h for overintegration (1:3,0:M,0:M,0:M,nElems)
REAL,ALLOCATABLE :: detJac_Ref(:,:,:,:,:)        !< determinant of the mesh Jacobian for each Gauss point at degree 3*NGeo
REAL,ALLOCATABLE :: sJ(:,:,:,:,:)                !< inverse of Jacobian determinent for each Gauss Point at degree N
!----------------------------------------------------------------------------------------------------------------------------------
! surface vectors
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: NormVec(:,:,:,:,:)           !< normal vector for each side       (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec1(:,:,:,:,:)          !< tangential vector 1 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec2(:,:,:,:,:)          !< tangential vector 3 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: SurfElem(:,:,:,:)            !< surface area for each side        (    0:N,0:N,nSides)
REAL,ALLOCATABLE :: Ja_Face(:,:,:,:,:)           !< surface  metrics for each side
REAL,ALLOCATABLE :: NormVecO(:,:,:,:,:)          !< normal vector for each side       for overintegration (1:3,0:M,0:M,nSides)
REAL,ALLOCATABLE :: TangVec1O(:,:,:,:,:)         !< tangential vector 1 for each side for overintegration (1:3,0:M,0:M,nSides)
REAL,ALLOCATABLE :: TangVec2O(:,:,:,:,:)         !< tangential vector 3 for each side for overintegration (1:3,0:M,0:M,nSides)
REAL,ALLOCATABLE :: SurfElemO(:,:,:,:)           !< surface area for each side        for overintegration (    0:M,0:M,nSides)
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: ElemToSide(:,:,:)       !< Array containing element-wise connectivity information to sides
                                               !< SideID    = ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
                                               !< flip      = ElemToSide(E2S_FLIP,ZETA_PLUS,iElem)

INTEGER,ALLOCATABLE :: SideToElem(:,:)         !< Array containing per-side connectivity information to elements and local side ids
                                               !< ElemID      = SideToElem(S2E_ELEM_ID,SideID)
                                               !< NB_ElemID   = SideToElem(S2E_NB_ELEM_ID,SideID)
                                               !< locSideID   = SideToElem(S2E_LOC_SIDE_ID,SideID)
                                               !< nblocSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)

INTEGER,ALLOCATABLE :: BC(:)                   !< BCIndex   = BC(SideID), 1:nBCSides

INTEGER,ALLOCATABLE :: BoundaryType(:,:)       !< List of boundary conditions containing type and state
                                               !< BCType    = BoundaryType(BC(SideID),BC_TYPE)
                                               !< BCState   = BoundaryType(BC(SideID),BC_STATE)

INTEGER,ALLOCATABLE :: AnalyzeSide(:)          !< Marks, wheter a side belongs to a group of analyze sides (e.g. to a BC group)
                                               !< SurfIndex = AnalyzeSide(SideID), 1:nSides

INTEGER,PARAMETER :: NormalDirs(6) = (/ 3 , 2 , 1 , 2 , 1 , 3 /) !< normal vector direction for element local side
INTEGER,PARAMETER :: TangDirs(6)   = (/ 1 , 3 , 2 , 3 , 2 , 1 /) !< first tangential vector direction for element local side
REAL   ,PARAMETER :: NormalSigns(6)= (/-1.,-1., 1., 1.,-1., 1./) !< normal vector sign for element local side

!----------------------------------------------------------------------------------------------------------------------------------
! Volume/Side mappings filled by mappings.f90 - not all available there are currently used!
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: FS2M(:,:,:,:)     !< flip slave side to master and reverse
INTEGER,ALLOCATABLE :: V2S(:,:,:,:,:,:)  !< volume to side mapping
INTEGER,ALLOCATABLE :: S2V(:,:,:,:,:,:)  !< side to volume
INTEGER,ALLOCATABLE :: S2V2(:,:,:,:,:)   !< side to volume 2
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: nGlobalElems=0          !< number of elements in mesh
INTEGER             :: nElems=0                !< number of local elements
INTEGER             :: offsetElem=0            !< for MPI, until now=0 Elems pointer array range: [offsetElem+1:offsetElem+nElems]
INTEGER             :: nSides=0                !< =nInnerSides+nBCSides+nMPISides
INTEGER             :: nSidesMaster=0          !< =sideIDMaster
INTEGER             :: nSidesSlave=0           !< =nInnerSides+nBCSides+nMPISides
INTEGER             :: nInnerSides=0           !< InnerSide index range: sideID [nBCSides+1:nBCSides+nInnerSides]
INTEGER             :: nBCSides=0              !< BCSide index range: sideID [1:nBCSides]
INTEGER             :: nAnalyzeSides=0         !< marker for each side (BC,analyze flag, periodic,...)
INTEGER             :: nMPISides=0             !< number of MPI sides in mesh
INTEGER             :: nMPISides_MINE=0        !< number of MINE MPI sides (on local processor)
INTEGER             :: nMPISides_YOUR=0        !< number of YOUR MPI sides (on neighbour processors)
INTEGER             :: nBCs=0                  !< number of BCs in mesh
INTEGER             :: nUserBCs=0              !< number of BC in inifile
!----------------------------------------------------------------------------------------------------------------------------------
! Define index ranges for all sides in consecutive order for easier access
INTEGER             :: firstBCSide             !< First SideID of BCs (in general 1)
INTEGER             :: firstMortarInnerSide    !< First SideID of Mortars (in general nBCs+1)
INTEGER             :: firstInnerSide          !< First SideID of inner sides
INTEGER             :: firstMPISide_MINE       !< First SideID of MINE MPI sides (on local processor)
INTEGER             :: firstMPISide_YOUR       !< First SideID of YOUR MPI sides (on neighbour processor)
INTEGER             :: firstMortarMPISide      !< First SideID of Mortar MPI sides
INTEGER             :: lastBCSide              !< Last  SideID of BCs (in general nBCs)
INTEGER             :: lastMortarInnerSide     !< Last  SideID of Mortars (in general nBCs+nMortars)
INTEGER             :: lastInnerSide           !< Last  SideID of inner sides
INTEGER             :: lastMPISide_MINE        !< Last  SideID of MINE MPI sides (on local processor)
INTEGER             :: lastMPISide_YOUR        !< Last  SideID of YOUR MPI sides (on neighbour processor)
INTEGER             :: lastMortarMPISide       !< Last  SideID of Mortar MPI sides (in general nSides)
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: nMortarSides=0          !< total number of mortar sides
INTEGER             :: nMortarInnerSides=0     !< number of inner mortar sides
INTEGER             :: nMortarMPISides=0       !< number of mortar MPI sides
INTEGER,ALLOCATABLE :: MortarType(:,:)         !< Type of mortar [1] and position in mortar list [1:nSides]
INTEGER,ALLOCATABLE :: MortarInfo(:,:,:)       !< 1:2,1:4,1:nMortarSides: [1] nbSideID / flip, [2] max 4 mortar sides, [3] sides
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),ALLOCATABLE :: BoundaryName(:) !< names of the boundary conditions read from the mesh file
CHARACTER(LEN=255)             :: MeshFile     !< name of hdf5 meshfile (write with ending .h5!)
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: useCurveds                 !< Marker wheter curved boundaries should be used in the mesh, read from parameter
                                               !< file. If set to false only linear part of curved meshes is used.

LOGICAL          :: CrossProductMetrics        !< Compute metrics in cross-product form instead of curl formulation
                                               !< (not recommended)
!----------------------------------------------------------------------------------------------------------------------------------
! USER DEFINED TYPES

!> Intermediate data type for side pointers (mesh readin only)
TYPE tSidePtr
  TYPE(tSide),POINTER          :: sp              !< side pointer
END TYPE tSidePtr

!> Intermediate data type for element pointers (mesh readin only)
TYPE tElemPtr
  TYPE(tElem),POINTER          :: ep              !< Local element pointer
END TYPE tElemPtr

!> Element data type, containing element mesh information (mesh readin only)
TYPE tElem
  INTEGER                      :: ind             !< global element index
  INTEGER                      :: Type            !< element type (linear/bilinear/curved)
  INTEGER                      :: Zone            !< zone of elements (unused)
  TYPE(tSidePtr),POINTER       :: Side(:)         !< sides connected to element
END TYPE tElem

!> Side data type (mesh readin only)
TYPE tSide
  INTEGER                      :: ind             !< global side ID
  INTEGER                      :: sideID          !< local side ID on Proc
  INTEGER                      :: tmp
  INTEGER                      :: NbProc          !< neighbor processor (if applicable)
  INTEGER                      :: BCindex         !< index in BoundaryType array!
  INTEGER                      :: flip            !< flip of side (0 if master, 1-4 if slave)
  INTEGER                      :: nMortars        !< number of slave mortar sides associated with master mortar
  INTEGER                      :: MortarType      !< type of mortar: Type1 : 1-4 , Type 2: 1-2 in eta, Type 2: 1-2 in xi
  TYPE(tSidePtr),POINTER       :: MortarSide(:)   !< array of side pointers to slave mortar sides
  TYPE(tElem),POINTER          :: Elem            !< pointer to connected element
  TYPE(tSide),POINTER          :: connection      !< pointer to connected neighbour side
END TYPE tSide

!----------------------------------------------------------------------------------------------------------------------------------
TYPE(tElemPtr),POINTER         :: Elems(:)        !< array of mesh elements with geometry and connectivity (only for readin)
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: MeshInitIsDone =.FALSE.       !< marks whether the mesh init routines are finished
!==================================================================================================================================

INTERFACE getNewSide
  MODULE PROCEDURE getNewSide
END INTERFACE

INTERFACE getNewElem
  MODULE PROCEDURE getNewElem
END INTERFACE

INTERFACE deleteMeshPointer
  MODULE PROCEDURE deleteMeshPointer
END INTERFACE

CONTAINS

!==================================================================================================================================
!> Build new side type and initialize values
!==================================================================================================================================
FUNCTION GETNEWSIDE()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tSide),POINTER :: getNewSide !< pointer to new side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
ALLOCATE(getNewSide)
NULLIFY(getNewSide%Elem)
NULLIFY(getNewSide%MortarSide)
NULLIFY(getNewSide%connection)
getNewSide%sideID=0
getNewSide%ind=0
getNewSide%tmp=0
getNewSide%NbProc=-1
getNewSide%BCindex=0
getNewSide%flip=0
getNewSide%nMortars=0
getNewSide%MortarType=0
END FUNCTION GETNEWSIDE

!==================================================================================================================================
!> Build new element type including sides and initialize values
!==================================================================================================================================
FUNCTION GETNEWELEM()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tElem),POINTER :: getNewElem !< pointer to new element
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iLocSide
!==================================================================================================================================
ALLOCATE(getNewElem)
ALLOCATE(getNewElem%Side(6))
DO iLocSide=1,6
  getNewElem%Side(iLocSide)%sp=>getNewSide()
END DO
getNewElem%ind=0
getNewElem%Zone=0
getNewElem%Type=0
END FUNCTION GETNEWELEM



!==================================================================================================================================
!> Deallocates all pointers used for the mesh readin
!==================================================================================================================================
SUBROUTINE deleteMeshPointer()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: FirstElemInd,LastElemInd
INTEGER       :: iElem,iLocSide
INTEGER       :: iMortar
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
!==================================================================================================================================
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    DO iMortar=1,aSide%nMortars
      NULLIFY(aSide%MortarSide(iMortar)%sp)
    END DO
    DEALLOCATE(aSide)
  END DO
  DEALLOCATE(aElem%Side)
  DEALLOCATE(aElem)
END DO
DEALLOCATE(Elems)
END SUBROUTINE deleteMeshPointer


END MODULE MOD_Mesh_Vars
