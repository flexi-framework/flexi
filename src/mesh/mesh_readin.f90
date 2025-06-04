!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
!> \brief Module containing routines to read the mesh and BCs from a HDF5 file
!>
!> This module contains the following routines related to mesh IO
!> - parallel HDF5-based mesh IO
!> - readin of mesh coordinates and connectivity
!> - readin of boundary conditions
!==================================================================================================================================
MODULE MOD_Mesh_Readin
! MODULES
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
!> @defgroup eleminfo ElemInfo parameters
!>  Named parameters for ElemInfo array in mesh file
!> @{
INTEGER,PARAMETER    :: ElemInfoSize=6        !< number of entries in each line of ElemInfo
INTEGER,PARAMETER    :: ELEM_Type=1
INTEGER,PARAMETER    :: ELEM_Zone=2
INTEGER,PARAMETER    :: ELEM_FirstSideInd=3
INTEGER,PARAMETER    :: ELEM_LastSideInd=4
INTEGER,PARAMETER    :: ELEM_FirstNodeInd=5
INTEGER,PARAMETER    :: ELEM_LastNodeInd=6
!> @}

!> @defgroup sideinfo SideInfo parameters
!>  Named parameters for SideInfo array in mesh file
!> @{
INTEGER,PARAMETER    :: SideInfoSize=5        !< number of entries in each line of SideInfo
INTEGER,PARAMETER    :: SIDE_Type=1
INTEGER,PARAMETER    :: SIDE_ID=2
INTEGER,PARAMETER    :: SIDE_nbElemID=3
INTEGER,PARAMETER    :: SIDE_Flip=4
INTEGER,PARAMETER    :: SIDE_BCID=5
!> @}

PUBLIC:: ReadMesh
PUBLIC:: BuildPartition
PUBLIC:: ReadIJKSorting
#if USE_MPI
PUBLIC:: ELEMIPROC
#endif /*USE_MPI*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This subroutine will read boundary conditions from the HDF5 mesh file and from the parameter file.
!> The parameters defined in the mesh file can be overridden by those defined in the parameter file, by specifying the boundary
!> name and a new boundary condition set which consists of a type and a state.
!==================================================================================================================================
SUBROUTINE ReadBCs()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars  ,ONLY:BoundaryName,BoundaryType,nBCs,nUserBCs
USE MOD_ReadInTools,ONLY:GETINTARRAY,CountOption,GETSTR
USE MOD_StringTools,ONLY:LowCase
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: Offset = 0 ! Every process reads all BCs
LOGICAL,ALLOCATABLE            :: UserBCFound(:)
LOGICAL                        :: NameCheck,LengthCheck
CHARACTER(LEN=255), ALLOCATABLE:: BCNames(:)
CHARACTER(LEN=255)             :: currBCName, currBoundaryName
CHARACTER(LEN=255)             :: ErrorString
INTEGER, ALLOCATABLE           :: BCMapping(:),BCType(:,:)
INTEGER                        :: iBC,iUserBC
!==================================================================================================================================
! read in boundary conditions from ini file, will overwrite BCs from meshfile!
nUserBCs = CountOption('BoundaryName')
IF(nUserBCs.GT.0)THEN
  ALLOCATE(BoundaryName(1:nUserBCs))
  ALLOCATE(BoundaryType(1:nUserBCs,2))
  DO iBC=1,nUserBCs
    BoundaryName(iBC)   = GETSTR('BoundaryName')
    BoundaryType(iBC,:) = GETINTARRAY('BoundaryType',2) !(/Type,State/)
  END DO
END IF ! nUserBCs.GT.0

! Read boundary names from data file
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
CHECKSAFEINT(HSize(1),4)
nBCs=INT(HSize(1),4)
DEALLOCATE(HSize)

ALLOCATE(BCNames(nBCs))
ALLOCATE(BCMapping(nBCs))
ALLOCATE(UserBCFound(nUserBCs))
CALL ReadArray('BCNames',1,(/nBCs/),Offset,1,StrArray=BCNames)
! User may have redefined boundaries in the ini file. So we have to create mappings for the boundaries.
BCMapping=0
UserBCFound=.FALSE.
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    DO iUserBC=1,nUserBCs
      ! Check if BoundaryName(iUserBC) is a substring of BCNames(iBC)
      CALL LowCase(BCNames(iBC)           ,currBCName)
      CALL LowCase(BoundaryName(iUserBC)  ,currBoundaryName)
      NameCheck = INDEX(TRIM(currBCName),TRIM(currBoundaryName)).NE.0
      ! Check if both strings have equal length
      LengthCheck = LEN(TRIM(BCNames(iBC))).EQ.LEN(TRIM(BoundaryName(iUserBC)))
      ! Check if both strings are equal (length has to be checked because index checks for substrings!)
      IF(NameCheck.AND.LengthCheck)THEN
        ! Check if the BC was defined multiple times
        IF (BCMapping(iBC).NE.0) THEN
          WRITE(ErrorString,'(A,A,A)') ' Boundary ',TRIM(BCNames(iBC)),' is redefined multiple times in parameter file!'
          CALL CollectiveStop(__STAMP__,ErrorString)
        END IF

        BCMapping(iBC)=iUserBC
        UserBCFound(iUserBC)=.TRUE.
      END IF
    END DO
  END DO
END IF
DO iUserBC=1,nUserBCs
  IF (.NOT.UserBCFound(iUserBC)) CALL Abort(__STAMP__,&
    'Boundary condition specified in parameter file has not been found: '//TRIM(BoundaryName(iUserBC)))
END DO
DEALLOCATE(UserBCFound)

! Read boundary types from data file
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF((HSize(1).NE.4).OR.(HSize(2).NE.nBCs)) CALL CollectiveStop(__STAMP__,'Problem in readBC')
DEALLOCATE(HSize)
ALLOCATE(BCType(4,nBCs))
CALL ReadArray('BCType',2,(/4,nBCs/),Offset,1,IntArray=BCType)
! Now apply boundary mappings
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    IF(BCMapping(iBC) .NE. 0)THEN
      IF((BoundaryType(BCMapping(iBC),1).EQ.1).AND.(BCType(1,iBC).NE.1)) &
        CALL Abort(__STAMP__,&
                   'Remapping non-periodic to periodic BCs is not possible!')
      SWRITE(UNIT_stdOut,'(A,A)')    ' |     Boundary in HDF file found | ',TRIM(BCNames(iBC))
      SWRITE(UNIT_stdOut,'(A,I4,I4)')' |                            was | ',BCType(1,iBC),BCType(3,iBC)
      SWRITE(UNIT_stdOut,'(A,I4,I4)')' |                      is set to | ',BoundaryType(BCMapping(iBC),1:2)
      BCType(1,iBC) = BoundaryType(BCMapping(iBC),BC_TYPE)
      BCType(3,iBC) = BoundaryType(BCMapping(iBC),BC_STATE)
    END IF
  END DO
END IF

IF(ALLOCATED(BoundaryName)) DEALLOCATE(BoundaryName)
IF(ALLOCATED(BoundaryType)) DEALLOCATE(BoundaryType)
ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(nBCs,3))
BoundaryName = BCNames
BoundaryType(:,BC_TYPE)  = BCType(1,:)
BoundaryType(:,BC_STATE) = BCType(3,:)
BoundaryType(:,BC_ALPHA) = BCType(4,:)
SWRITE(UNIT_stdOut,'(132("."))')
SWRITE(UNIT_stdOut,'(A,A15,A20,A10,A10,A10)')' BOUNDARY CONDITIONS','|','Name','Type','State','Alpha'
DO iBC=1,nBCs
  SWRITE(Unit_stdOut,'(A,A33,A20,I10,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:)
END DO

SWRITE(UNIT_stdOut,'(132("."))')

DEALLOCATE(BCNames,BCType,BCMapping)

END SUBROUTINE ReadBCs


!==================================================================================================================================
!> This subroutine reads the mesh from the HDF5 mesh file. The connectivity and further relevant information as flips
!> (i.e. the orientation of sides towards each other) is already contained in the mesh file.
!> For parallel computations the number of elements will be distributed equally onto all processors and each processor only reads
!> its own subset of the mesh.
!> For a documentation of the mesh format see the documentation provided with HOPR (hopr-project.org)
!> The arrays ElemInfo, SideInfo and NodeCoords are read, alongside with the boundary condition data.
!> If the mesh is non-conforming and based on a tree representation, the corresponding tree data (Coords, parameter ranges,
!> connectivity) is also read in.
!==================================================================================================================================
SUBROUTINE ReadMesh(FileString)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,       ONLY:ReadMeshWallTime
USE MOD_Mesh_Vars,          ONLY:tElem,tSide
USE MOD_Mesh_Vars,          ONLY:NGeo,NGeoTree
USE MOD_Mesh_Vars,          ONLY:NodeCoords,TreeCoords
USE MOD_Mesh_Vars,          ONLY:offsetElem,offsetTree,nElems,nTrees,nGlobalTrees
USE MOD_Mesh_Vars,          ONLY:xiMinMax,ElemToTree
USE MOD_Mesh_Vars,          ONLY:nSides,nInnerSides,nBCSides,nMPISides,nAnalyzeSides
USE MOD_Mesh_Vars,          ONLY:nMortarSides,isMortarMesh,meshHasMortars
USE MOD_Mesh_Vars,          ONLY:useCurveds
USE MOD_Mesh_Vars,          ONLY:BoundaryType
USE MOD_Mesh_Vars,          ONLY:MeshInitIsDone
USE MOD_Mesh_Vars,          ONLY:Elems
USE MOD_Mesh_Vars,          ONLY:GETNEWELEM,GETNEWSIDE
USE MOD_Mesh_Vars,          ONLY:ElemInfo,SideInfo
#if USE_MPI
USE MOD_MPI_Vars,           ONLY:nMPISides_Proc,nNbProcs,NbProc
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString !< (IN) mesh filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER            :: aElem
TYPE(tSide),POINTER            :: aSide,bSide
REAL,ALLOCATABLE               :: NodeCoordsTmp(:,:,:,:,:)
#if PP_dim == 2
INTEGER                        :: BCindex
#endif
INTEGER                        :: iElem,nbElemID,nNodes
INTEGER                        :: iLocSide,nbLocSide
INTEGER                        :: iSide
INTEGER                        :: FirstSideInd,LastSideInd,FirstElemInd,LastElemInd
INTEGER                        :: nPeriodicSides,nMPIPeriodics
INTEGER                        :: ReduceData(10)
INTEGER                        :: nSideIDs,offsetSideID
INTEGER                        :: iMortar,jMortar,nMortars
#if USE_MPI
INTEGER                        :: ReduceData_glob(10)
INTEGER                        :: iNbProc
INTEGER                        :: iProc
INTEGER,ALLOCATABLE            :: MPISideCount(:)
#endif
LOGICAL                        :: oriented
LOGICAL                        :: dsExists
REAL                           :: StartT,EndT
!==================================================================================================================================
IF(MESHInitIsDone) RETURN
IF(MPIRoot)THEN
  IF(.NOT.FILEEXISTS(FileString))  CALL CollectiveStop(__STAMP__, &
    'readMesh from data file "'//TRIM(FileString)//'" does not exist')
END IF

SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_stdOut,'(132("-"))')
GETTIME(StartT)

! Open mesh file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL BuildPartition()
FirstElemInd=offsetElem+1
LastElemInd =offsetElem+nElems

CALL readBCs()

!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!
! Build a list of all local elements and set their: global index, type and zone
! but not the pointers to the sides
!----------------------------------------------------------------------------------------------------------------------------
! list of pointers to tElem objects for all local elements
ALLOCATE(Elems(FirstElemInd:LastElemInd))

!read local ElemInfo from mesh file
ALLOCATE(ElemInfo(ElemInfoSize,FirstElemInd:LastElemInd))
CALL ReadArray('ElemInfo',2,(/ElemInfoSize,nElems/),offsetElem,2,IntArray=ElemInfo)

! Fill list of tElem objects by the information stored in ElemInfo
DO iElem=FirstElemInd,LastElemInd
  Elems(iElem)%ep=>GETNEWELEM() ! create a new element object with default values and store pointer of it in the list
  aElem=>Elems(iElem)%ep        ! get the element and set: global index, type and zone
  aElem%Ind    = iElem
  aElem%Type   = ElemInfo(ELEM_Type,iElem)
  aElem%Zone   = ElemInfo(ELEM_Zone,iElem)
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!
! Iterate over all elements and within each element over all sides (6 for hexas).
! For each side set the following informations:
!  - nMortars   (Number of small sides attached to a big Mortar side)
!  - MortarType (zero for non Mortar sides; =1,2 or 3 for big Mortar side; -1 for small sides belonging to a Mortar)
!  - Ind        (global index of the side)
!  - flip
!  - For the big Mortar sides additionally iterate over the small virtual sides and set:
!    - Ind
!    - flip
!----------------------------------------------------------------------------------------------------------------------------

! get offset of local side indices in the global sides
offsetSideID=ElemInfo(ELEM_FirstSideInd,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs    =ElemInfo(ELEM_LastSideInd,LastElemInd)-ElemInfo(ELEM_FirstSideInd,FirstElemInd)
! read local SideInfo from mesh file
FirstSideInd=offsetSideID+1
LastSideInd =offsetSideID+nSideIDs
ALLOCATE(SideInfo(SideInfoSize,FirstSideInd:LastSideInd))
CALL ReadArray('SideInfo',2,(/SideInfoSize,nSideIDs/),offsetSideID,2,IntArray=SideInfo)

! iterate over all local elements and within each element over all sides
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  ! build up sides of the element according to CGNS standard
  ! assign flip
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    nbElemID=SideInfo(SIDE_nbElemID,iSide) ! get neighboring element index of element adjacent to this side

#if PP_dim == 2
    ! In 2D check that there is only one layer of elements in z-direction
    IF ((iLocSide.EQ.1).OR.(iLocSide.EQ.6)) THEN
      BCindex = SideInfo(SIDE_BCID,iSide)
      IF ((iElem.NE.nbElemID).AND.(BCindex.LE.0)) THEN
        CALL Abort(__STAMP__, &
            "Mesh not oriented in z-direction or more than one layer of elements in z-direction! " // &
            "Please set 'orientZ = T' or change number of element in z-direction in HOPR parameter file.")
      END IF
    END IF
#endif

    ! If nbElemID <0, this marks a mortar master side.
    ! The number (-1,-2,-3) is the Type of mortar
    IF(nbElemID.LT.0)THEN ! mortar sides attached!
      aSide%MortarType=ABS(nbElemID)
      SELECT CASE(aSide%MortarType)
      CASE(1)
        aSide%nMortars=4
      CASE(2,3)
        aSide%nMortars=2
      END SELECT

      ALLOCATE(aSide%MortarSide(aSide%nMortars))
      DO iMortar=1,aSide%nMortars
        aSide%MortarSide(iMortar)%sp=>GETNEWSIDE()
      END DO
    ELSE
      aSide%nMortars=0
    END IF

    ! mark side as belonging to a mortar
    IF(SideInfo(SIDE_Type,iSide).LT.0) aSide%MortarType = -1

    ! side is not a big Mortar side
    IF(aSide%MortarType.LE.0)THEN
      aSide%Elem=>aElem
      oriented=(Sideinfo(SIDE_ID,iSide).GT.0)
      aSide%Ind=ABS(SideInfo(SIDE_ID,iSide))
      ! oriented side
      IF(oriented) THEN
        aSide%flip=0
      ! not oriented side
      ELSE
        aSide%flip=MOD(Sideinfo(SIDE_Flip,iSide),10)
        IF((aSide%flip.LT.0).OR.(aSide%flip.GT.4)) CALL Abort(__STAMP__,'NodeID doesnt belong to side')
      END IF

    ! side is a big mortar side
    ELSE
      ! iterate over virtual small mortar sides
      DO iMortar = 1,aSide%nMortars
        iSide=iSide+1
        aSide%mortarSide(iMortar)%sp%Elem=>aElem ! set element pointer to actual element
        IF(SideInfo(SIDE_ID,iSide).LT.0) CALL Abort(__STAMP__,'Problem in Mortar readin,should be flip=0')
        aSide%mortarSide(iMortar)%sp%flip=0
        aSide%mortarSide(iMortar)%sp%Ind =ABS(SideInfo(SIDE_ID,iSide))
      END DO !iMortar
    END IF
  END DO !i=1,locnSides
END DO !iElem


!----------------------------------------------------------------------------------------------------------------------------
!                              CONNECTIONS
!
! Iterate over all elements and within each element over all sides (6 for hexas) and for each big Mortar side over all
! small virtual sides.
! For each side do:
!   - if BC side:
!       Reset side and do not connect
!   - if neighboring element is on local processor:
!       loop over all sides of neighboring element and find side with same index as local side => connect
!   - if neighboring element is on other processor:
!       create a new virtual neighboring side and element
!----------------------------------------------------------------------------------------------------------------------------
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    ! LOOP over mortars, if no mortar, then LOOP is executed once
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0)THEN
        iSide=iSide+1
        aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp ! point to small virtual side
      END IF
      nbElemID      = SideInfo(SIDE_nbElemID,iSide)
      aSide%BCindex = SideInfo(SIDE_BCID,iSide)

      ! BC sides don't need a connection, except for internal (BC_TYPE=0), periodic (BC_TYPE=1) and "dummy" inner BCs (BC_TYPE=100).
      ! For all other BC sides: reset the flip and mortars settings, do not build a connection.
      IF(aSide%BCindex.NE.0)THEN ! BC
        IF((BoundaryType(aSide%BCindex,BC_TYPE).NE.0).AND.&
           (BoundaryType(aSide%BCindex,BC_TYPE).NE.1).AND.&
           (BoundaryType(aSide%BCindex,BC_TYPE).NE.100))THEN
          aSide%flip  =0
          IF(iMortar.EQ.0) aSide%mortarType  = 0
          IF(iMortar.EQ.0) aSide%nMortars    = 0
          CYCLE
        END IF
      END IF

      IF(aSide%mortarType.GT.0) CYCLE        ! no connection for mortar master
      IF(ASSOCIATED(aSide%connection)) CYCLE ! already connected

      IF(nbElemID.NE.0)THEN ! neighboring element exists, since it has a valid global index
        ! check if neighbor on local proc or MPI connection
        IF((nbElemID.LE.LastElemInd).AND.(nbElemID.GE.FirstElemInd))THEN !local
          ! connect sides by looping over the side of the neighboring element and finding the side that
          ! has the same index
          DO nbLocSide=1,6
            bSide=>Elems(nbElemID)%ep%Side(nbLocSide)%sp
            ! LOOP over mortars, if no mortar, then LOOP is executed once
            nMortars=bSide%nMortars
            DO jMortar=0,nMortars
              IF(jMortar.GT.0) bSide=>Elems(nbElemID)%ep%Side(nbLocSide)%sp%mortarSide(jMortar)%sp

              ! check for matching indices of local side and side of the neighboring element
              IF(bSide%ind.EQ.aSide%ind)THEN
                aSide%connection=>bSide
                bSide%connection=>aSide
                EXIT
              END IF
            END DO !jMortar
          END DO !nbLocSide
        ELSE ! MPI connection
#if USE_MPI
          aSide%connection     =>GETNEWSIDE()
          aSide%connection%flip= aSide%flip
          aSide%connection%Elem=>GETNEWELEM()
          aSide%NbProc = ELEMIPROC(nbElemID)
#else
          CALL Abort(__STAMP__, &
            ' ElemID of neighbor not in global Elem list ')
#endif
        END IF
      END IF
    END DO !iMortar
  END DO !iLocSide
END DO !iElem


!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------
! get physical coordinates

IF(useCurveds)THEN
  ALLOCATE(NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),offsetElem*(NGeo+1)**3,2,RealArray=NodeCoords)
ELSE
  ALLOCATE(NodeCoords(   3,0:1,   0:1,   0:1,   nElems))
  ALLOCATE(NodeCoordsTmp(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  ! read all nodes
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),offsetElem*(NGeo+1)**3,2,RealArray=NodeCoordsTmp)
  ! throw away all nodes except the 8 corner nodes of each hexa
  NodeCoords(:,0,0,0,:)=NodeCoordsTmp(:,0,   0,   0,   :)
  NodeCoords(:,1,0,0,:)=NodeCoordsTmp(:,NGeo,0,   0,   :)
  NodeCoords(:,0,1,0,:)=NodeCoordsTmp(:,0,   NGeo,0,   :)
  NodeCoords(:,1,1,0,:)=NodeCoordsTmp(:,NGeo,NGeo,0,   :)
  NodeCoords(:,0,0,1,:)=NodeCoordsTmp(:,0,   0,   NGeo,:)
  NodeCoords(:,1,0,1,:)=NodeCoordsTmp(:,NGeo,0,   NGeo,:)
  NodeCoords(:,0,1,1,:)=NodeCoordsTmp(:,0,   NGeo,NGeo,:)
  NodeCoords(:,1,1,1,:)=NodeCoordsTmp(:,NGeo,NGeo,NGeo,:)
  DEALLOCATE(NodeCoordsTmp)
  NGeo=1 ! linear mesh; set polynomial degree of geometry to 1
ENDIF

! total number of nodes on this processor
nNodes = nElems*(NGeo+1)**3

#if PP_dim == 2
DO iElem=1,nElems
  IF (NodeCoords(3,0,0,0,iElem).GT.NodeCoords(3,0,0,1,iElem))THEN
    CALL Abort(__STAMP__, &
        "Zeta is oriented in negative z-direction, has to be oriented in positive z-direction instead" // &
        "Please set 'orientZ = T' in HOPR parameter file.")
  END IF
END DO
#endif

! Get Mortar specific additional arrays concering Octrees.
! General idea of the octrees:
!   Starting from a conform mesh each element can be halved in each reference direction leading to
!   2, 4 or 8 smaller elements. This bisection can be repeated on the smaller elements. The original
!   coarsest element of the conform mesh is the root of a octree and in the following called a 'tree-element'
!   or just 'tree'. Only the coordinates of the smallest elements are necessary for a valid mesh representation.
!   Nevertheless the mesh format allows to store additional information of the tree-elements (coordinates)
!   and the position of the regular elements inside the trees.
!
! If a mesh is built by refining a conform baseline mesh the refining octree can be
! stored in the mesh file in the following attibutes and arrays:
!   'isMortarMesh' : if present and ==1 indicates a non-conforming Mortar mesh
!   'NgeoTree'     : polynomial degree of the tree-elements
!   'nTrees'       : number of tree-elements
!   'nNodeTrees'   : total number of nodes of the tree-elements (NgeoTree+1)^3 * nTrees
!   'ElemToTree'   : mapping from global element index (ElemID) to octree index (TreeID), which allows to
!                    find to each element the corresponding tree.
!                    size = (1:nElems)
!   'xiMinMax'     : element bounds of an element inside the reference space [-1,1]^3 of its tree-element
!                    size = (1:3, 1:2, 1:nElems)
!                    first index = coordinate
!                    second index =1 for minimum corner node, =2 for maximum corner node
!   'TreeCoords'   : coordinates of all nodes of the tree-elements
dsExists=.FALSE.
iMortar=0
CALL DatasetExists(File_ID,'isMortarMesh',dsExists,.TRUE.)
IF(dsExists)&
  CALL ReadAttribute(File_ID,'isMortarMesh',1,IntScalar=iMortar)
isMortarMesh=(iMortar.EQ.1)
IF(isMortarMesh)THEN
  CALL ReadAttribute(File_ID,'NgeoTree',1,IntScalar=NGeoTree)
  CALL ReadAttribute(File_ID,'nTrees',1,IntScalar=nGlobalTrees)

  ALLOCATE(xiMinMax(3,2,1:nElems))
  xiMinMax=-1.
  CALL ReadArray('xiMinMax',3,(/3,2,nElems/),offsetElem,3,RealArray=xiMinMax)

  ALLOCATE(ElemToTree(1:nElems))
  ElemToTree=0
  CALL ReadArray('ElemToTree',1,(/nElems/),offsetElem,1,IntArray=ElemToTree)

  ! only read trees, connected to a procs elements
  offsetTree=MINVAL(ElemToTree)-1
  ElemToTree=ElemToTree-offsetTree
  nTrees=MAXVAL(ElemToTree)

  ALLOCATE(TreeCoords(3,0:NGeoTree,0:NGeoTree,0:NGeoTree,nTrees))
  TreeCoords=-1.
  CALL ReadArray('TreeCoords',2,(/3,(NGeoTree+1)**3*nTrees/),&
                 (NGeoTree+1)**3*offsetTree,2,RealArray=TreeCoords)
! no mortar mesh
ELSE
  nTrees=0
END IF


CALL CloseDataFile()
! Readin of mesh is now finished


!----------------------------------------------------------------------------------------------------------------------------
!                              COUNT SIDES
!----------------------------------------------------------------------------------------------------------------------------
nBCSides=0
nAnalyzeSides=0
nMortarSides=0
nSides=0
nPeriodicSides=0
nMPIPeriodics=0
nMPISides=0
#if USE_MPI
ALLOCATE(MPISideCount(0:nProcessors-1))
MPISideCount=0
#endif

! Iterate over all elements and within each element over all sides (6 for hexas) and for each big Mortar side over all
! small virtual sides and reset 'tmp' array of each side to zero.
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    ! LOOP over mortars, if no mortar, then LOOP is executed once
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      ! point to small virtual side
      IF(iMortar.GT.0) aSide => aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
      aSide%tmp=0
    END DO !iMortar
  END DO !iLocSide
END DO !iElem

! Iterate over all elements and within each element over all sides (6 for hexas in 3D, 4 for quads in 2D)
! and for each big Mortar side over all small virtual sides
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
#if PP_dim == 3
  DO iLocSide=1,6
#else
  DO iLocSide=2,5
#endif
    aSide=>aElem%Side(iLocSide)%sp
    ! LOOP over mortars, if no mortar, then LOOP is executed once
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      ! point to small virtual side
      IF(iMortar.GT.0) aSide => aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp

      ! if side not counted so far
      IF(aSide%tmp.EQ.0)THEN
        nSides=nSides+1
        aSide%tmp=-1 ! mark side as counted
        IF(ASSOCIATED(aSide%connection))THEN
          ! sanity check that the boundary index is zero
          ! IF(aSide%BCindex.NE.0) CALL Abort(__STAMP__,'Error with BCIndex on connected side')
          aSide%connection%tmp = -1  ! mark connected side as counted
        END IF

        ! side is BC or periodic side.
        IF(aSide%BCindex.NE.0)THEN
          nAnalyzeSides=nAnalyzeSides+1
          IF(ASSOCIATED(aSide%connection))THEN
            SELECT CASE(BoundaryType(aSide%BCindex,BC_TYPE))
              ! internal side
              CASE(0)
                ! do nothing

              ! periodic side
              CASE(1)
                nPeriodicSides=nPeriodicSides+1
#if USE_MPI
                IF(aSide%NbProc.NE.-1) nMPIPeriodics=nMPIPeriodics+1
#endif
              CASE DEFAULT
                CALL Abort(__STAMP__,'Error with BoundaryType on connected side')
            END SELECT
          ! BC side
          ELSE
            IF(aSide%MortarType.EQ.0)     nBCSides      = nBCSides+1
          END IF ! ASSOCIATED(aSide%connection)
        END IF
        IF(aSide%MortarType.GT.0) nMortarSides=nMortarSides+1
#if USE_MPI
        ! count total number of MPI sides and number of MPI sides for each neighboring processor
        IF(aSide%NbProc.NE.-1) THEN
          nMPISides=nMPISides+1
          MPISideCount(aSide%NbProc)=MPISideCount(aSide%NbProc)+1
        END IF
#endif
      END IF
    END DO !iMortar
  END DO !iLocSide
END DO !iElem

! Periodic sides count as inner sides!
nInnerSides = nSides-nBCSides-nMPISides-nMortarSides

LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A22,I8)')'nSides:',      nSides
LOGWRITE(*,'(A22,I8)')'nBCSides:',    nBCSides
LOGWRITE(*,'(A22,I8)')'nMortarSides:',nMortarSides
LOGWRITE(*,'(A22,I8)')'nInnerSides:', nInnerSides
LOGWRITE(*,'(A22,I8)')'nMPISides:',   nMPISides
LOGWRITE(*,*)'-------------------------------------------------------'

#if USE_MPI
! Define variables used for MPI communication:
!   - nBnProcs: count of neighboring processors (sharing a side with this one)
!   - NbProc(1:nNbProcs): global index of each neigboring processor
!   - nMPISides_Proc(1:nNbProcs): count of common sides with each neighboring processor

nNBProcs=0
! iterate over all processors and count processors which are neighbors (having a common side)
DO iProc=0,nProcessors-1
  IF(iProc.EQ.myRank) CYCLE
  IF(MPISideCount(iProc).GT.0) nNBProcs=nNbProcs+1
END DO

! compiled with MPI, but execute on a single processor
IF(nNbProcs.EQ.0)THEN
  ALLOCATE(NbProc(1),nMPISides_Proc(1))
  nNbProcs=1
  NbProc=0
  nMPISides_Proc=0
ELSE
  ! build a mapping 1..nNbProcs to the global number of processors (allows to iterate over all neighboring processors)
  ! and store for each neighboring processor the number of common sides
  ALLOCATE(NbProc(nNbProcs),nMPISides_Proc(1:nNbProcs))
  iNbProc=0
  DO iProc=0,nProcessors-1
    IF(iProc.EQ.myRank) CYCLE
    IF(MPISideCount(iProc).GT.0) THEN
      iNbProc=iNbProc+1
      NbProc(iNbProc)=iProc  ! map neighboring proc number to global proc number
      nMPISides_Proc(iNBProc)=MPISideCount(iProc) ! store number of common sides
    END IF
  END DO
END IF
DEALLOCATE(MPISideCount)

#endif /*USE_MPI*/

! Sum up all element, side and node numbers over all processors (just for stdout)
ReduceData(1)=nElems
ReduceData(2)=nSides
ReduceData(3)=nNodes
ReduceData(4)=nInnerSides
ReduceData(5)=nPeriodicSides
ReduceData(6)=nBCSides
ReduceData(7)=nMPISides
ReduceData(8)=nAnalyzeSides
ReduceData(9)=nMortarSides
ReduceData(10)=nMPIPeriodics

! Store information if mesh has mortars, independent of existence of trees
meshHasMortars = MERGE(.TRUE.,.FALSE.,ReduceData(9).GT.0)

#if USE_MPI
CALL MPI_REDUCE(ReduceData,ReduceData_glob,10,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
IF(MPIRoot) ReduceData=ReduceData_glob
#endif /*USE_MPI*/

IF(MPIRoot)THEN
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nElems               | ',ReduceData(1) !nElems
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nNodes               | ',ReduceData(3) !nNodes
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides               | ',ReduceData(2)-ReduceData(7)/2
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides,    BC        | ',ReduceData(6) !nBCSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides,   MPI        | ',ReduceData(7)/2 !nMPISides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides, Inner        | ',ReduceData(4) !nInnerSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides,Mortar        | ',ReduceData(9) !nMortarSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nPeriodicSides,Total | ',ReduceData(5)-ReduceData(10)/2
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nPeriodicSides,Inner | ',ReduceData(5)-ReduceData(10)
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nPeriodicSides,  MPI | ',ReduceData(10)/2 !nPeriodicSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nAnalyzeSides        | ',ReduceData(8) !nAnalyzeSides
  WRITE(UNIT_stdOut,'(A,A34,L1)')' |','useCurveds           | ',useCurveds
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','Ngeo                 | ',Ngeo
  WRITE(UNIT_stdOut,'(132("."))')
END IF

EndT             = FLEXITIME()
ReadMeshWallTime = EndT-StartT
CALL DisplayMessageAndTime(ReadMeshWallTime,'READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ... DONE',DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE ReadMesh


!==================================================================================================================================
!> Partition the mesh by numbers of processors. Elements are distributed equally to all processors.
!==================================================================================================================================
SUBROUTINE BuildPartition()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY:nElems,nGlobalElems,offsetElem
#if USE_MPI
USE MOD_MPI_Vars,  ONLY:offsetElemMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER           :: iElem
INTEGER           :: iProc
#endif
!===================================================================================================================================
CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
IF(HSize(1).NE.6) CALL Abort(__STAMP__,'ERROR: Wrong size of ElemInfo, should be 6')

CHECKSAFEINT(HSize(2),4)
nGlobalElems=INT(HSize(2),4)
DEALLOCATE(HSize)
#if USE_MPI
IF(nGlobalElems.LT.nProcessors) CALL Abort(__STAMP__,&
  'ERROR: Number of elements (1) is smaller then number of processors (2)!',nGlobalElems,REAL(nProcessors))

! Simple partition: nGlobalelems/nprocs, do this on proc 0
SDEALLOCATE(offsetElemMPI)
ALLOCATE(offsetElemMPI(0:nProcessors))
offsetElemMPI=0
nElems=nGlobalElems/nProcessors
iElem=nGlobalElems-nElems*nProcessors
DO iProc=0,nProcessors-1
  offsetElemMPI(iProc)=nElems*iProc+MIN(iProc,iElem)
END DO
offsetElemMPI(nProcessors)=nGlobalElems

! Local nElems and offset
nElems=offsetElemMPI(myRank+1)-offsetElemMPI(myRank)
offsetElem=offsetElemMPI(myRank)
LOGWRITE(*,*)'offset,nElems',offsetElem,nElems
#else /*USE_MPI*/
nElems=nGlobalElems   !local number of Elements
offsetElem=0          ! offset is the index of first entry, hdf5 array starts at 0-.GT. -1
#endif /*USE_MPI*/

END SUBROUTINE BuildPartition


#if USE_MPI
!==================================================================================================================================
!> Find the id of a processor on which an element with a given ElemID lies, based on the MPI element offsets defined earlier.
!> Use a bisection algorithm for faster search.
!==================================================================================================================================
FUNCTION ELEMIPROC(ElemID)
! MODULES
USE MOD_Globals,   ONLY:nProcessors
USE MOD_MPI_vars,  ONLY:offsetElemMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)                :: ElemID     !< (IN)  NodeID to search for
INTEGER                            :: ELEMIPROC  !< (OUT) processor id
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,high,mid
!==================================================================================================================================
ELEMIPROC=0
maxSteps=INT(LOG(REAL(nProcessors))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=0
high      = nProcessors-1
IF(    ElemID.GT.offsetElemMPI(low)  .AND. ElemID.LE.offsetElemMPI(low +1))THEN
  ELEMIPROC=low
ELSEIF(ElemID.GT.offsetElemMPI(high) .AND. ElemID.LE.offsetElemMPI(high+1))THEN
  ELEMIPROC = high
ELSE
  !bisection
  DO i=1,maxSteps
    mid = (high-low)/2+low
    IF(    ElemID.GT.offsetElemMPI(mid) .AND. ElemID.LE.offsetElemMPI(mid+1))THEN
      ELEMIPROC = mid                            ! index found
      EXIT
    ELSEIF(ElemID .GT. offsetElemMPI(mid+1))THEN ! seek in upper half
      low=mid+1
    ELSE                                         ! seek in lower half
      high = mid
    END IF
  END DO
END IF
END FUNCTION ELEMIPROC
#endif /*USE_MPI*/


!===================================================================================================================================
!> Read arrays nElems_IJK (global number of elements in i,j,k direction) and Elem_IJK (mapping from global element to i,j,k index)
!> for meshes thar are i,j,k sorted.
!===================================================================================================================================
SUBROUTINE ReadIJKSorting(doGlobal)
! MODULES
USE MOD_Mesh_Vars,       ONLY: nElems_IJK,Elem_IJK,offsetElem,nElems,nGlobalElems,MeshFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL      :: doGlobal
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: dsExists
LOGICAL                          :: doGlobal_loc
!===================================================================================================================================
IF (PRESENT(doGlobal)) THEN
  doGlobal_loc = doGlobal
ELSE
  doGlobal_loc = .FALSE.
END IF

CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,'nElems_IJK',dsExists)
IF(dsExists)THEN
  CALL ReadArray('nElems_IJK',1,(/3/),0,1,IntArray=nElems_IJK)
  IF (doGlobal_loc) THEN
    ALLOCATE(Elem_IJK(3,nGlobalElems))
    CALL ReadArray('Elem_IJK',2,(/3,nGlobalElems/),0,2,IntArray=Elem_IJK)
  ELSE
    ALLOCATE(Elem_IJK(3,nElems))
    CALL ReadArray('Elem_IJK',2,(/3,nElems/),offsetElem,2,IntArray=Elem_IJK)
  END IF
END IF
CALL CloseDataFile()
END SUBROUTINE ReadIJKSorting

END MODULE MOD_Mesh_ReadIn
