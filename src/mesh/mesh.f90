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
!> Contains control routines to read high-order meshes, provide mesh data to the solver, build the metrics, partition the domain.
!==================================================================================================================================
MODULE MOD_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE BuildOverintMesh
  MODULE PROCEDURE BuildOverintMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::BuildOverintMesh
PUBLIC::FinalizeMesh
!==================================================================================================================================

PUBLIC::DefineParametersMesh
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateStringOption(  'MeshFile',            "(relative) path to meshfile (mandatory).")
CALL prms%CreateLogicalOption( 'useCurveds',          "Controls usage of high-order information in mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption( 'interpolateFromTree', "For non-conforming meshes, built by refinement from a tree structure, "//&
                                                      "the metrics can be built from the tree geometry if it is contained "//&
                                                      "in the mesh. Can improve free-stream preservation.",&
                                                      '.TRUE.')
CALL prms%CreateRealOption(    'meshScale',           "Scale the mesh by this factor (shrink/enlarge).",&
                                                      '1.0')
CALL prms%CreateLogicalOption( 'meshdeform',          "Apply simple sine-shaped deformation on cartesion mesh (for testing).",&
                                                      '.FALSE.')
CALL prms%CreateLogicalOption( 'crossProductMetrics', "Compute mesh metrics using cross product form. Caution: in this case "//&
                                                      "free-stream preservation is only guaranteed for N=3*NGeo.",&
                                                      '.FALSE.')
CALL prms%CreateIntOption(     'debugmesh',           "Output file with visualization and debug information for the mesh. "//&
                                                      "0: no visualization, 3: Paraview binary",'0')
CALL prms%CreateStringOption(  'BoundaryName',        "Names of boundary conditions to be set (must be present in the mesh!)."//&
                                                      "For each BoundaryName a BoundaryType needs to be specified.",&
                                                      multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType',        "Type of boundary conditions to be set. Format: (BC_TYPE,BC_STATE)",&
                                                      multiple=.TRUE.)
CALL prms%CreateLogicalOption( 'writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
END SUBROUTINE DefineParametersMesh


!==================================================================================================================================
!> Routine controlling the initialization of the mesh.
!> - parameter and mesh reading
!> - domain partitioning
!> - allocation of mesh arrays
!> - build mesh mappings to handle volume/surface operations
!> - compute the mesh metrics
!> - provide mesh metrics for overintegration
!==================================================================================================================================
SUBROUTINE InitMesh(meshMode,MeshFile_IN)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_HDF5_Input
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR,GETREAL,GETINT
USE MOD_Metrics,            ONLY:BuildCoords,CalcMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Mappings,           ONLY:buildMappings
#if USE_MPI
USE MOD_Prepare_Mesh,       ONLY:exchangeFlip
#endif
#if FV_ENABLED
USE MOD_FV_Metrics,         ONLY:InitFV_Metrics
#endif
USE MOD_IO_HDF5,            ONLY:AddToElemData,ElementOut
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode !< 0: only read and build Elem_xGP,
                               !< 1: as 0 + build connectivity, 2: as 1 + calc metrics
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: MeshFile_IN !< file name of mesh to be read
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: x(3),meshScale
REAL,POINTER      :: coords(:,:,:,:,:)
INTEGER           :: iElem,i,j,k,nElemsLoc
LOGICAL           :: validMesh
INTEGER           :: firstMasterSide     ! lower side ID of array U_master/gradUx_master...
INTEGER           :: lastMasterSide      ! upper side ID of array U_master/gradUx_master...
INTEGER           :: firstSlaveSide      ! lower side ID of array U_slave/gradUx_slave...
INTEGER           :: lastSlaveSide       ! upper side ID of array U_slave/gradUx_slave...
INTEGER           :: iSide,LocSideID,SideID
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'InitMesh not ready to be called or already called.')
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,I1,A)') ' INIT MESH IN MODE ',meshMode,'...'

! prepare pointer structure (get nElems, etc.)
IF (PRESENT(MeshFile_IN)) THEN
  MeshFile = MeshFile_IN
ELSE
  MeshFile = GETSTR('MeshFile')
END IF
validMesh = ISVALIDMESHFILE(MeshFile)
IF(.NOT.validMesh) &
    CALL CollectiveStop(__STAMP__,'ERROR - Mesh file not a valid HDF5 mesh.')

useCurveds=GETLOGICAL('useCurveds','.TRUE.')
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo)
CALL CloseDataFile()

IF(useCurveds.AND.(PP_N.LT.NGeo))THEN
  SWRITE(UNIT_stdOut,'(A)') 'WARNING: N<NGeo, for curved hexa normals are only approximated,&
                           & can cause problems on periodic boundaries! Set N>=NGeo'
ENDIF

CALL readMesh(MeshFile) !set nElems

! if trees are available: compute metrics on tree level and interpolate to elements
interpolateFromTree=.FALSE.
IF(isMortarMesh) interpolateFromTree=GETLOGICAL('interpolateFromTree','.TRUE.')
IF(interpolateFromTree)THEN
  coords=>TreeCoords
  NGeo=NGeoTree
  nElemsLoc=nTrees
ELSE
  coords=>NodeCoords
  nElemsLoc=nElems
ENDIF
SWRITE(UNIT_StdOut,'(a3,a30,a3,i0)')' | ','Ngeo',' | ', Ngeo

! scale and deform mesh if desired (warning: no mesh output!)
meshScale=GETREAL('meshScale','1.0')
IF(ABS(meshScale-1.).GT.1e-14)&
  Coords = Coords*meshScale

IF(GETLOGICAL('meshdeform','.FALSE.'))THEN
  DO iElem=1,nElemsLoc
    DO k=0,PP_NGeoZ; DO j=0,NGeo; DO i=0,NGeo
      x(:)=Coords(:,i,j,k,iElem)
#if PP_dim==3
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))*SIN(PP_Pi*x(3))
#else
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))
#endif
    END DO; END DO; END DO;
  END DO
END IF

ALLOCATE(Elem_xGP(3,0:PP_N,0:PP_N,0:PP_N,nElems))
CALL BuildCoords(Elem_xGP)

! Return if no connectivity and metrics are required (e.g. for visualization mode)
IF (meshMode.GT.0) THEN

  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
  CALL setLocalSideIDs()
  
#if USE_MPI
  ! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
  CALL exchangeFlip()
#endif
  
  !RANGES
  !-----------------|-----------------|-------------------|
  !    U_master     | U_slave         |    FLUX           |
  !-----------------|-----------------|-------------------|
  !  BCsides        |                 |    BCSides        |
  !  InnerMortars   |                 |    InnerMortars   |
  !  InnerSides     | InnerSides      |    InnerSides     |
  !  MPI_MINE sides | MPI_MINE sides  |    MPI_MINE sides |
  !                 | MPI_YOUR sides  |    MPI_YOUR sides |
  !  MPIMortars     |                 |    MPIMortars     |
  !-----------------|-----------------|-------------------|
  
  firstBCSide          = 1
  firstMortarInnerSide = firstBCSide         +nBCSides
  firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
  firstMPISide_MINE    = firstInnerSide      +nInnerSides
  firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
  firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR
  
  lastBCSide           = firstMortarInnerSide-1
  lastMortarInnerSide  = firstInnerSide    -1
  lastInnerSide        = firstMPISide_MINE -1
  lastMPISide_MINE     = firstMPISide_YOUR -1
  lastMPISide_YOUR     = firstMortarMPISide-1
  lastMortarMPISide    = nSides
  
  
  firstMasterSide = 1
  lastMasterSide  = nSides
  firstSlaveSide  = firstInnerSide
  lastSlaveSide   = lastMPISide_YOUR
  nSidesMaster    = lastMasterSide-firstMasterSide+1
  nSidesSlave     = lastSlaveSide -firstSlaveSide+1
  
  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8)')   'first/lastMasterSide     ', firstMasterSide,lastMasterSide
  LOGWRITE(*,'(A25,I8)')   'first/lastSlaveSide      ', firstSlaveSide, lastSlaveSide
  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8,I8)')'first/lastBCSide         ', firstBCSide         ,lastBCSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMortarInnerSide', firstMortarInnerSide,lastMortarInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastInnerSide      ', firstInnerSide      ,lastInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_MINE   ', firstMPISide_MINE   ,lastMPISide_MINE
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_YOUR   ', firstMPISide_YOUR   ,lastMPISide_YOUR
  LOGWRITE(*,'(A30,I8,I8)')'first/lastMortarMPISide  ', firstMortarMPISide  ,lastMortarMPISide
  LOGWRITE(*,*)'-------------------------------------------------------'
  
  ! fill ElemToSide, SideToElem,BC
  ALLOCATE(ElemToSide(2,6,nElems))
  ALLOCATE(SideToElem(5,nSides))
  ALLOCATE(BC(1:nBCSides))
  ALLOCATE(AnalyzeSide(1:nSides))
  ElemToSide  = 0
  SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
  BC          = 0
  AnalyzeSide = 0
  
  !NOTE: nMortarSides=nMortarInnerSides+nMortarMPISides
  ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo
  ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides)) ! [1]: 1: Neighbour sides, 2: Flip, [2]: small sides
  MortarType=-1
  MortarInfo=-1
  
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
  CALL fillMeshInfo()
  
  ! dealloacte pointers
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
  CALL deleteMeshPointer()
  
  ! Build necessary mappings 
  CALL buildMappings(PP_N,V2S=V2S,S2V=S2V,S2V2=S2V2,FS2M=FS2M)
END IF

IF (meshMode.GT.1) THEN

  ! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----
  
  ! volume data
  ALLOCATE(      dXCL_N(3,3,0:PP_N,0:PP_N,0:PP_N,nElems)) ! temp
  ALLOCATE(Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_ENABLED))
  ALLOCATE(Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_ENABLED))
  ALLOCATE(Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_ENABLED))
  ALLOCATE(            sJ(  0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_ENABLED))
  NGeoRef=3*NGeo ! build jacobian at higher degree
  ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))
  
  ! surface data
  ALLOCATE(      Face_xGP(3,0:PP_N,0:PP_N,0:FV_ENABLED,1:nSides))
  ALLOCATE(       NormVec(3,0:PP_N,0:PP_N,0:FV_ENABLED,1:nSides))
  ALLOCATE(      TangVec1(3,0:PP_N,0:PP_N,0:FV_ENABLED,1:nSides))
  ALLOCATE(      TangVec2(3,0:PP_N,0:PP_N,0:FV_ENABLED,1:nSides))
  ALLOCATE(      SurfElem(  0:PP_N,0:PP_N,0:FV_ENABLED,1:nSides))
  ALLOCATE(     Ja_Face(3,3,0:PP_N,0:PP_N,             1:nSides)) ! temp
  
  
  ! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
  ! assign 1/detJ (sJ)
  ! assign normal and tangential vectors and surfElems on faces
  
  ! compute metrics using cross product instead of curl form (warning: no free stream preservation!)
  crossProductMetrics=GETLOGICAL('crossProductMetrics','.FALSE.')
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."
  CALL CalcMetrics()     ! DG metrics
#if FV_ENABLED
  CALL InitFV_Metrics()  ! FV metrics
#endif
  ! debugmesh: param specifies format to output, 0: no output, 1: tecplot ascii, 2: tecplot binary, 3: paraview binary
  CALL WriteDebugMesh(GETINT('debugmesh','0'))
END IF

#if PP_dim == 2
  CALL Convert2D(meshMode)
#endif

IF (meshMode.GT.0) THEN
  ALLOCATE(SideToGlobalSide(nSides))
  DO iElem=1,nElems
#if PP_dim == 3
    DO LocSideID=1,6
#else    
    DO LocSideID=2,5
#endif    
      SideID = ElemToSide(E2S_SIDE_ID,LocSideID,iElem)
      iSide = ElemInfo(3,iElem+offsetElem) + LocSideID
      SideToGlobalSide(SideID) = ABS(SideInfo(2,iSide))
    END DO
  END DO ! iElem
END IF

SDEALLOCATE(NodeCoords)
SDEALLOCATE(dXCL_N)
SDEALLOCATE(Ja_Face)
SDEALLOCATE(TreeCoords)
SDEALLOCATE(xiMinMax)
SDEALLOCATE(ElemToTree)

CALL AddToElemData(ElementOut,'myRank',IntScalar=myRank)

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh

!==================================================================================================================================
!> In case the selective version of overintegration is used, all mesh data also needs to be provided on the degree NOver of the
!> higher quadrature rule.
!> This routine interpolates all metric terms from N to NOver, except the surface metrics (normal/tangential vectors, surface area)
!> which have to be recomputed using the Ja interpolated to NOver.
!==================================================================================================================================
SUBROUTINE BuildOverintMesh()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Metrics,             ONLY: CalcSurfMetrics
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_Interpolation_Vars,  ONLY: NodeTypeCL,NodeType
USE MOD_Overintegration_Vars,ONLY: NOver,VdmNToNOver
USE MOD_ChangeBasisByDim,    ONLY: ChangeBasisVolume
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: JaCL_NSurf(3,3,0:NOver,0:NOver,0:PP_NOverZ) !< metric terms P\in NOver
REAL              :: XCL_NSurf(   3,0:NOver,0:NOver,0:PP_NOverZ) !< geometry on P\in NOver
REAL              :: Vdm_N_CLNSurf(0:NOver,0:PP_N)
REAL              :: Vdm_CLNSurf_NSurf(0:NOver,0:NOver)
INTEGER           :: iElem
!==================================================================================================================================

! Build geometry for volume overintegration
ALLOCATE(      Elem_xGPO(3,0:NOver,0:NOver,0:PP_NOverZ,nElems)) ! only needed once by fillini
ALLOCATE(Metrics_fTildeO(3,0:NOver,0:NOver,0:PP_NOverZ,nElems))
ALLOCATE(Metrics_gTildeO(3,0:NOver,0:NOver,0:PP_NOverZ,nElems))
ALLOCATE(Metrics_hTildeO(3,0:NOver,0:NOver,0:PP_NOverZ,nElems))
IF(NOver.GT.PP_N)THEN
  DO iElem=1,nElems
    CALL ChangeBasisVolume(3,PP_N,NOver,VdmNToNOver,Elem_xGP(:,:,:,:,iElem),              Elem_xGPO(:,:,:,:,iElem))
    CALL ChangeBasisVolume(3,PP_N,NOver,VdmNToNOver,Metrics_fTilde(:,:,:,:,iElem,0),Metrics_fTildeO(:,:,:,:,iElem))
    CALL ChangeBasisVolume(3,PP_N,NOver,VdmNToNOver,Metrics_gTilde(:,:,:,:,iElem,0),Metrics_gTildeO(:,:,:,:,iElem))
    CALL ChangeBasisVolume(3,PP_N,NOver,VdmNToNOver,Metrics_hTilde(:,:,:,:,iElem,0),Metrics_hTildeO(:,:,:,:,iElem))
  END DO ! iElem
END IF

! Build geometry for surface overintegration
ALLOCATE(Face_xGPO(3,0:NOver,0:PP_NOverZ,0:0,1:nSides))
ALLOCATE( NormVecO(3,0:NOver,0:PP_NOverZ,0:0,1:nSides))
ALLOCATE(TangVec1O(3,0:NOver,0:PP_NOverZ,0:0,1:nSides))
ALLOCATE(TangVec2O(3,0:NOver,0:PP_NOverZ,0:0,1:nSides))
ALLOCATE(SurfElemO(  0:NOver,0:PP_NOverZ,0:0,1:nSides))

CALL GetVandermonde( PP_N , NodeType  , NOver , NodeTypeCL , Vdm_N_CLNSurf     , modal=.FALSE.)
CALL GetVandermonde( NOver, NodeTypeCL, NOver , NodeType   , Vdm_CLNSurf_NSurf , modal=.FALSE.)
DO iElem=1,nElems
  CALL ChangeBasisVolume(3,PP_N,NOver,Vdm_N_CLNSurf,Metrics_fTilde(:,:,:,:,iElem,0),JaCL_NSurf(1,:,:,:,:))
  CALL ChangeBasisVolume(3,PP_N,NOver,Vdm_N_CLNSurf,Metrics_gTilde(:,:,:,:,iElem,0),JaCL_NSurf(2,:,:,:,:))
  CALL ChangeBasisVolume(3,PP_N,NOver,Vdm_N_CLNSurf,Metrics_hTilde(:,:,:,:,iElem,0),JaCL_NSurf(3,:,:,:,:))
  CALL ChangeBasisVolume(3,PP_N,NOver,Vdm_N_CLNSurf,Elem_xGP(:,:,:,:,iElem),XCL_NSurf)
#if PP_dim==2
  STOP 'Surface metric computation not implemented yet'
#endif
  CALL CalcSurfMetrics(NOver,0,JaCL_NSurf,XCL_NSurf,Vdm_CLNSurf_NSurf,iElem,&
                       NormVecO,TangVec1O,TangVec2O,SurfElemO,Face_xGPO)
END DO
END SUBROUTINE BuildOverintMesh


#if PP_dim == 2
!==================================================================================================================================
!> This routine converts all 3D mesh quantities to a 2D mesh, including mappings, rotations etc.
!> We use CGNS notation, which is for 2D:
!>
!>         ^
!>         |
!>         |  eta
!>
!>       <---
!>     +-------+
!>   | |       | ^
!>   | |       | |   --->  xi
!>   v |       | |
!>     +-------+
!>       --->
!>
!> The arrows along the square indicate the coordinate system on the side. These correspond to the first index of the 3D side 
!> coordinate systems, except for XI_MINUS, where it corresponds to the negative of the second index, so all values need to be 
!> rotated accordingly from 3D on these sides.
!==================================================================================================================================
SUBROUTINE Convert2D(meshMode)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_2D
#if FV_ENABLED
USE MOD_FV_Vars
#endif
USE MOD_Mappings ,ONLY: buildMappings
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,iSide
REAL               :: tmp(3,0:PP_N,0:PP_N,0:FV_ENABLED),zlength
#if FV_ENABLED
REAL               :: tmpFV(0:PP_N,0:PP_N,3)
#endif
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') " CONVERT TO 2D ... "

! Nullify all components in the third dimension
Elem_xGP(3,:,:,:,:) = 0.
CALL to2D_rank5((/1,0,0,0,1/),  (/3,PP_N,PP_N,PP_N,nElems/),4,Elem_xGP)

IF (meshMode.GT.0) THEN
  ! In 2D, there is only one flip for the slave sides (1)
  SideToElem(S2E_FLIP,:) = MIN(1,SideToElem(S2E_FLIP,:))
  ElemToSide(:,1,:) = -999
  ElemToSide(:,6,:) = -999
  ElemToSide(E2S_FLIP,:,:) = MIN(1,ElemToSide(E2S_FLIP,:,:))
  MortarInfo(MI_FLIP,:,:)  = MIN(1,MortarInfo(MI_FLIP,:,:))

  CALL buildMappings(PP_N,V2S=V2S,S2V=S2V,S2V2=S2V2,FS2M=FS2M, dim=2)
END IF

IF (meshMode.GT.1) THEN
  ! Correctly flip and rotate all values on the XI_MINUS sides.
  ! The tangential vector is always stored in TangVec1.
  DO iSide=1,nSides
    SELECT CASE (SideToElem(S2E_LOC_SIDE_ID,iSide))
    CASE(XI_MINUS)
      ! has to be flipped
      tmp=Face_xGP(:,:,:,:,iSide)
      DO j=0,PP_N; DO i=0,PP_N
        Face_xGP(:,PP_N-j,i,:,iSide)=tmp(:,i,j,:)
      END DO; END DO
      tmp=NormVec(:,:,:,:,iSide)
      DO j=0,PP_N; DO i=0,PP_N
        NormVec(:,PP_N-j,i,:,iSide)=tmp(:,i,j,:)
      END DO; END DO
      tmp=-TangVec1(:,:,:,:,iSide)
      DO j=0,PP_N; DO i=0,PP_N
        TangVec1(:,PP_N-j,i,:,iSide)=tmp(:,i,j,:)
      END DO; END DO
      tmp(1,:,:,:)=SurfElem(:,:,:,iSide)
      DO j=0,PP_N; DO i=0,PP_N
        SurfElem(PP_N-j,i,:,iSide)=tmp(1,i,j,:)
      END DO; END DO
    CASE(ETA_MINUS)
      TangVec1(:,:,:,:,iSide) = -TangVec2(:,:,:,:,iSide)
    CASE(ETA_PLUS)
      TangVec1(:,:,:,:,iSide) = -TangVec2(:,:,:,:,iSide)
    END SELECT
  END DO

  ! Nullify all components in the third dimension
  Metrics_fTilde(3,:,:,:,:,:) = 0.
  Metrics_gTilde(3,:,:,:,:,:) = 0.
  Metrics_hTilde(:,:,:,:,:,:) = 0.

  NormVec (3,:,:,:,:) = 0.
  TangVec1(3,:,:,:,:) = 0.
  TangVec2(:,:,:,:,:) = 0.
  Face_xGP(3,:,:,:,:) = 0.

  ! Reduce the bounds of the arrays to be 0:0 in the third dimension (dimension is still kept for compatibility reasons)
  CALL to2D_rank6((/1,0,0,0,1,0/),(/3,PP_N,PP_N,PP_N,nElems,FV_ENABLED/),4,Metrics_fTilde)
  CALL to2D_rank6((/1,0,0,0,1,0/),(/3,PP_N,PP_N,PP_N,nElems,FV_ENABLED/),4,Metrics_gTilde)
  CALL to2D_rank6((/1,0,0,0,1,0/),(/3,PP_N,PP_N,PP_N,nElems,FV_ENABLED/),4,Metrics_hTilde)
  CALL to2D_rank5((/0,0,0,1,0/),  (/PP_N,PP_N,PP_N,nElems,FV_ENABLED/),3,sJ)
  CALL to2D_rank5((/1,0,0,0,1/),  (/1,NgeoRef,NgeoRef,NgeoRef,nElems/),4,DetJac_Ref)

  CALL to2D_rank5((/1,0,0,0,1/),(/3,PP_N,PP_N,FV_ENABLED,nSides/),3,Face_xGP)
  CALL to2D_rank5((/1,0,0,0,1/),(/3,PP_N,PP_N,FV_ENABLED,nSides/),3,NormVec)
  CALL to2D_rank5((/1,0,0,0,1/),(/3,PP_N,PP_N,FV_ENABLED,nSides/),3,TangVec1)
  CALL to2D_rank5((/1,0,0,0,1/),(/3,PP_N,PP_N,FV_ENABLED,nSides/),3,TangVec2)
  CALL to2D_rank4((/0,0,0,1/),  (/PP_N,PP_N,FV_ENABLED,nSides/),2,SurfElem)

  !computation of z length 
  zlength=abs(nodeCoords(3,0,0,Ngeo,1)-nodeCoords(3,0,0,0,1))
  !normalization of geometric terms by z length (reference elemt has length 2 [-1,1])
  sJ=sJ*(zlength/2.)
  ! attention: SurfElem contains NANs (MPI_YOUR sides and MPIMortars). Arithemtic with those causes runtime errors with INTEL debug
  SurfElem(:,:,:,1:lastMPISide_MINE)=SurfElem(:,:,:,1:lastMPISide_MINE)/(zlength/2.)  
  Metrics_fTilde = Metrics_fTilde/(zlength/2.)
  Metrics_gTilde = Metrics_gTilde/(zlength/2.)


#if FV_ENABLED
  FV_TangVec1Eta=-FV_TangVec2Eta


#if FV_RECONSTRUCT
  ! Correctly flip and rotate all values on the XI_MINUS sides.
  ! The tangential vector is always stored in TangVec1.
  DO iSide=1,nSides
    SELECT CASE (SideToElem(S2E_LOC_SIDE_ID,iSide))
    CASE(XI_MINUS)
      ! has to be flipped
      tmpFV=FV_sdx_Face(:,:,:,iSide)
      DO j=0,PP_N; DO i=0,PP_N
        FV_sdx_Face(PP_N-j,i,:,iSide)=tmpFV(i,j,:)
      END DO; END DO

      tmpFV(:,:,1)=FV_dx_slave(1,:,:,iSide)
      DO j=0,PP_N; DO i=0,PP_N
        FV_dx_slave(1,PP_N-j,i,iSide)=tmpFV(i,j,1)
      END DO; END DO
      tmpFV(:,:,1)=FV_dx_master(1,:,:,iSide)
      DO j=0,PP_N; DO i=0,PP_N
        FV_dx_master(1,PP_N-j,i,iSide)=tmpFV(i,j,1)
      END DO; END DO
    END SELECT
  END DO

  CALL to2D_rank4((/0,0,1,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_sdx_XI  ) 
  CALL to2D_rank4((/0,0,1,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_sdx_ETA ) 
  CALL to2D_rank4((/0,0,1,1/),  (/PP_N,PP_N,PP_N,nElems/),3,FV_sdx_ZETA) 

  CALL to2D_rank4((/0,0,1,1/),  (/PP_N,PP_N,3,lastMPISide_MINE/),2,FV_sdx_Face) 

  CALL to2D_rank4((/0,0,0,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_dx_XI_L  )
  CALL to2D_rank4((/0,0,0,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_dx_XI_R  )
  CALL to2D_rank4((/0,0,0,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_dx_ETA_L )
  CALL to2D_rank4((/0,0,0,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_dx_ETA_R )
  CALL to2D_rank4((/0,0,0,1/),  (/PP_N,PP_N,PP_N,nElems/),3,FV_dx_ZETA_L)
  CALL to2D_rank4((/0,0,0,1/),  (/PP_N,PP_N,PP_N,nElems/),3,FV_dx_ZETA_R)

  CALL to2D_rank4((/1,0,0,1/),  (/1,PP_N,PP_N,nSides/),3,FV_dx_slave )
  CALL to2D_rank4((/1,0,0,1/),  (/1,PP_N,PP_N,nSides/),3,FV_dx_master)
#endif
  ! p,q = general face index
  CALL to2D_rank4((/0,0,1,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_SurfElemXi_sw  ) ! Attention: storage order is (p,q,i,iElem)
  CALL to2D_rank4((/0,0,1,1/),  (/PP_N,PP_N,PP_N,nElems/),2,FV_SurfElemEta_sw ) ! Attention: storage order is (p,q,j,iElem)
  CALL to2D_rank4((/0,0,1,1/),  (/PP_N,PP_N,PP_N,nElems/),3,FV_SurfElemZeta_sw) ! Attention: storage order is (p,q,k,iElem)

  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),3,FV_NormVecXi   )  ! Attention: storage order is (p,q,i,iElem)
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),3,FV_TangVec1Xi  )  !  -"-
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),3,FV_TangVec2Xi  )  !  -"-
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),3,FV_NormVecEta  )  ! Attention: storage order is (p,q,j,iElem)
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),3,FV_TangVec1Eta )  !  -"-
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),3,FV_TangVec2Eta )  !  -"-
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),4,FV_NormVecZeta )  ! Attention: storage order is (p,q,k,iElem)
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),4,FV_TangVec1Zeta)  !  -"-
  CALL to2D_rank5((/1,0,0,1,1/),  (/3,PP_N,PP_N,PP_N,nElems/),4,FV_TangVec2Zeta)  !  -"-

#if PARABOLIC
  CALL to2D_rank5((/1,0,0,0,1/),  (/3,PP_N,PP_N,PP_N,nElems/),4,FV_Metrics_fTilde_sJ)
  CALL to2D_rank5((/1,0,0,0,1/),  (/3,PP_N,PP_N,PP_N,nElems/),4,FV_Metrics_gTilde_sJ)
  CALL to2D_rank5((/1,0,0,0,1/),  (/3,PP_N,PP_N,PP_N,nElems/),4,FV_Metrics_hTilde_sJ)
  FV_Metrics_fTilde_sJ = FV_Metrics_fTilde_sJ/(zlength/2.)
  FV_Metrics_gTilde_sJ = FV_Metrics_gTilde_sJ/(zlength/2.)
#endif

  ! attention: SurfElem contains NANs (MPI_YOUR sides and MPIMortars). Arithemtic with those causes runtime errors with INTEL debug
  FV_SurfElemXi_sw=FV_SurfElemXi_sw/(zlength/2.)  
  FV_SurfElemEta_sw=FV_SurfElemEta_sw/(zlength/2.)  
  FV_SurfElemZeta_sw=FV_SurfElemZeta_sw/(zlength/2.)  

#endif /* FV_ENABLED */
END IF

SWRITE(UNIT_stdOut,'(A)')'   DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE Convert2D
#endif

!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE FinalizeMesh()
! MODULES
USE MOD_Mesh_Vars
USE MOD_Mappings  ,ONLY:FinalizeMappings
#if FV_ENABLED
USE MOD_FV_Vars   ,ONLY:FV_Elems_master
USE MOD_FV_Metrics,ONLY:FinalizeFV_Metrics
#endif
IMPLICIT NONE
!============================================================================================================================
!> Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
SDEALLOCATE(AnalyzeSide)

SDEALLOCATE(MortarType)
SDEALLOCATE(MortarInfo)


! allocated during ReadMesh
SDEALLOCATE(NodeCoords)
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)


!> Volume
SDEALLOCATE(Elem_xGP)
SDEALLOCATE(Metrics_fTilde)
SDEALLOCATE(Metrics_gTilde)
SDEALLOCATE(Metrics_hTilde)
SDEALLOCATE(sJ)
SDEALLOCATE(DetJac_Ref)

SDEALLOCATE(Elem_xGPO)
SDEALLOCATE(Metrics_fTildeO)
SDEALLOCATE(Metrics_gTildeO)
SDEALLOCATE(Metrics_hTildeO)

!> surface
SDEALLOCATE(Face_xGP)
SDEALLOCATE(NormVec)
SDEALLOCATE(TangVec1)
SDEALLOCATE(TangVec2)
SDEALLOCATE(SurfElem)

SDEALLOCATE(Face_xGPO)
SDEALLOCATE(NormVecO)
SDEALLOCATE(TangVec1O)
SDEALLOCATE(TangVec2O)
SDEALLOCATE(SurfElemO)

! ijk sorted mesh
SDEALLOCATE(Elem_IJK)
SDEALLOCATE(ElemInfo)
SDEALLOCATE(SideInfo)

!> mappings
CALL FinalizeMappings()

#if FV_ENABLED
SDEALLOCATE(FV_Elems_master) ! moved here from fv.f90
CALL FinalizeFV_Metrics()
#endif

MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
