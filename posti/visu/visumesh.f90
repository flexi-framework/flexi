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

!=================================================================================================================================
!> Routines to build the mesh for visualization.
!=================================================================================================================================
MODULE MOD_Posti_VisuMesh
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE BuildVisuCoords
  MODULE PROCEDURE BuildVisuCoords
END INTERFACE

INTERFACE BuildSurfVisuCoords
  MODULE PROCEDURE BuildSurfVisuCoords
END INTERFACE

INTERFACE VisualizeMesh
  MODULE PROCEDURE VisualizeMesh
END INTERFACE

PUBLIC:: BuildVisuCoords
PUBLIC:: BuildSurfVisuCoords
PUBLIC:: VisualizeMesh

CONTAINS

!=================================================================================================================================
!> Converts the coordinates of the mesh to the visu-mesh.
!=================================================================================================================================
SUBROUTINE BuildVisuCoords()
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeVISUFVEqui,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Mesh_Vars          ,ONLY: nGlobalElems,NodeCoords,NGeo
USE MOD_Visu_Vars          ,ONLY: CoordsVisu_DG
USE MOD_Visu_Vars          ,ONLY: NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu,nElems_DG,mapDGElemsToAllElems
USE MOD_Visu_Vars          ,ONLY: nElemsAvg2D_DG,Avg2D
USE MOD_Visu_Vars          ,ONLY: Elem_IJK_glob,mapElemIJToDGElemAvg2D
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeVISUFVEqui,NodeType
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY: offsetElemMPI
#endif
#if FV_ENABLED
USE MOD_Visu_Vars          ,ONLY: FVAmountAvg2D,mapElemIJToFVElemAvg2D,nElemsAvg2D_FV
USE MOD_Visu_Vars          ,ONLY: NVisu_FV,nElems_FV,mapFVElemsToAllElems,hasFV_Elems
USE MOD_Visu_Vars          ,ONLY: CoordsVisu_FV,changedMeshFile,changedFV_Elems,changedAvg2D
#endif
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY: offsetElemMPI
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,iElem_DG
REAL,ALLOCATABLE   :: Vdm_NGeo_NVisu(:,:)
#if FV_ENABLED
INTEGER            :: iElem_FV
REAL,ALLOCATABLE   :: Vdm_NGeo_NVisu_FV(:,:)
#endif
INTEGER            :: iElemAvg,ii,jj,kk
#if USE_MPI
INTEGER            :: nNodesPerProc(0:nProcessors-1),offsetNodes(0:nProcessors-1)
INTEGER            :: iProc
REAL,ALLOCATABLE   :: NodeCoords_Glob(:,:,:,:,:)
#else
REAL,POINTER       :: NodeCoords_Glob(:,:,:,:,:)
#endif
!===================================================================================================================================

! Convert coordinates to visu grid
SWRITE (Unit_stdOut,'(A)') ' [MESH] Convert coordinates to visu grid (DG)'

! Convert from NodeCoords to create consistent meshes
ALLOCATE(Vdm_NGeo_NVisu(0:NVisu,0:NGeo))
CALL GetVandermonde(NGeo,NodeTypeVisu,NVisu   ,NodeTypeVisuPosti ,Vdm_NGeo_NVisu   ,modal=.FALSE.)
#if FV_ENABLED
ALLOCATE(Vdm_NGeo_NVisu_FV(0:NVisu_FV,0:NGeo))
CALL GetVandermonde(NGeo,NodeTypeVisu,NVisu_FV,NodeTypeVISUFVEqui,Vdm_NGeo_NVisu_FV,modal=.FALSE.)
#endif

! Standard visualization
IF (.NOT.Avg2D) THEN
  ! convert coords of DG elements
  SDEALLOCATE(CoordsVisu_DG)
  ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems_DG))
  DO iElem_DG = 1,nElems_DG
    iElem = mapDGElemsToAllElems(iElem_DG)
    CALL ChangeBasisVolume(3,NGeo,NVisu,Vdm_NGeo_NVisu,NodeCoords(:,:,:,:,iElem),CoordsVisu_DG(:,:,:,:,iElem_DG))
  END DO

#if FV_ENABLED
  IF (hasFV_Elems) THEN
    SWRITE (Unit_stdOut,'(A)') ' [MESH] Convert coordinates to visu grid (FV)'
    ! only NVisu changed, but NVisu_FV is independent of NVisu
    IF ((.NOT.changedMeshFile).AND.(.NOT.changedFV_Elems).AND.(.NOT.changedAvg2D)) RETURN
    ! convert coords of FV elements
    SDEALLOCATE(CoordsVisu_FV)
    ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),nElems_FV))
    DO iElem_FV = 1,nElems_FV
      iElem = mapFVElemsToAllElems(iElem_FV)
      CALL ChangeBasisVolume(3,NGeo,NVisu_FV,Vdm_NGeo_NVisu_FV,NodeCoords(:,:,:,:,iElem),CoordsVisu_FV(:,:,:,:,iElem_FV))
    END DO
    SDEALLOCATE(Vdm_NGeo_NVisu_FV)
  END IF
#endif

! convert coords of DG elements
ELSE
  ! For parallel averaging, the root will gather the whole mesh and convert the first layer to the  visu grid.
#if USE_MPI
  ! For the gather operation, we need to know the number of DOFs per processor
  DO iProc = 0, nProcessors-1
    nNodesPerProc(iProc) = (offsetElemMPI(iProc+1) - offsetElemMPI(iProc)) * 3*(NGeo+1)**(PP_dim)
  END DO ! iProc = 1, nProcessors

  ! On the root, we need the receive array
  IF (MPIRoot) THEN
    ALLOCATE(NodeCoords_Glob(1:3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nGlobalElems))
    DO iProc = 0, nProcessors-1
      offsetNodes(iProc) = offsetElemMPI(iProc) * 3*(NGeo+1)**(PP_dim)
    END DO ! iProc = 1, nProcessors
  END IF

  CALL MPI_GATHERV(NodeCoords     ,nNodesPerProc(myRank)    ,MPI_DOUBLE_PRECISION                          ,&
                   NodeCoords_Glob,nNodesPerProc,offsetNodes,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
#else
  NodeCoords_Glob => NodeCoords
#endif /* USE_MPI */

  ! Now, the root can build the Avg2D mesh
  IF (MPIRoot) THEN
    SDEALLOCATE(CoordsVisu_DG)
    ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:0,nElemsAvg2D_DG))
#if FV_ENABLED
    SDEALLOCATE(CoordsVisu_FV)
    ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV))
#endif
    DO iElem = 1,nGlobalElems
      ii = Elem_IJK_glob(1,iElem)
      jj = Elem_IJK_glob(2,iElem)
      kk = Elem_IJK_glob(3,iElem)
      IF (kk.EQ.1) THEN
#if FV_ENABLED
        IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN ! DG
#endif
          iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
          CALL ChangeBasis2D(3,NGeo,NVisu   ,Vdm_NGeo_NVisu   ,NodeCoords_Glob(:,:,:,0,iElem),CoordsVisu_DG(:,:,:,0,iElemAvg))
#if FV_ENABLED
        ELSE ! FV
          iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
          CALL ChangeBasis2D(3,NGeo,NVisu_FV,Vdm_NGeo_NVisu_FV,NodeCoords_Glob(:,:,:,0,iElem),CoordsVisu_FV(:,:,:,0,iElemAvg))
        END IF
#endif
      END IF
    END DO
#if USE_MPI
    DEALLOCATE(NodeCoords_Glob)
#endif /*USE_MPI*/

  ! .NOT. MPIRoot
  ELSE
    ! All other procs will have 0 average elements
    SDEALLOCATE(CoordsVisu_DG)
    ALLOCATE(CoordsVisu_DG(3,0:NVisu   ,0:NVisu   ,0:0,nElemsAvg2D_DG))
#if FV_ENABLED
    SDEALLOCATE(CoordsVisu_FV)
    ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV))
#endif
  END IF ! MPIRoot

#if !USE_MPI
  NULLIFY(NodeCoords_Glob )
#endif /*!USE_MPI*/
END IF

SDEALLOCATE(Vdm_NGeo_NVisu)

END SUBROUTINE BuildVisuCoords


!=================================================================================================================================
!> Converts the coordinates of the surface mesh to the surface visu-mesh.
!=================================================================================================================================
SUBROUTINE BuildSurfVisuCoords()
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeVISUFVEqui,NodeType
USE MOD_Mappings           ,ONLY: SideToVol2
USE MOD_Mesh_Vars          ,ONLY: NGeo,NodeCoords,nBCSides
USE MOD_Mesh_Vars          ,ONLY: SideToElem
USE MOD_Visu_Vars          ,ONLY: CoordsSurfVisu_DG,nBCSidesVisu_DG,mapAllBCSidesToDGVisuBCSides
USE MOD_Visu_Vars          ,ONLY: NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu
#if FV_ENABLED
USE MOD_Visu_Vars          ,ONLY: CoordsSurfVisu_FV,nBCSidesVisu_FV,mapAllBCSidesToFVVisuBCSides
USE MOD_Visu_Vars          ,ONLY: NVisu_FV,hasFV_Elems
USE MOD_Visu_Vars          ,ONLY: changedMeshFile,changedFV_Elems,changedBCnames
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,iSide,iLocSide,iSideVisu
INTEGER            :: p,q,pq(2)
REAL,ALLOCATABLE   :: Vdm_NGeo_NVisu(:,:)
REAL               :: SurfCoords(1:3,0:NGeo,0:NGeo,0:0,1:nBCSides)   !< XYZ positions (first index 1:3) of the Face Node Coords
REAL               :: tmp       (1:3,0:NGeo,0:NGeo)
#if FV_ENABLED
REAL,ALLOCATABLE   :: Vdm_NGeo_NVisu_FV(:,:)
#endif
!===================================================================================================================================
! Convert coordinates to visu grid
SWRITE (Unit_stdOut,'(A)') ' [MESH] Convert coordinates to surface visu grid (DG)'

! Build surf coords on node coords
DO iSide = 1,nBCSides
  iElem    = SideToElem(S2E_ELEM_ID    ,iSide)
  iLocSide = SideToElem(S2E_LOC_SIDE_ID,iSide)

  SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=NodeCoords(1:3,0   ,:   ,:   ,iElem)
    CASE(XI_PLUS)
      tmp=NodeCoords(1:3,NGeo,:   ,:   ,iElem)
    CASE(ETA_MINUS)
      tmp=NodeCoords(1:3,:   ,0   ,:   ,iElem)
    CASE(ETA_PLUS)
      tmp=NodeCoords(1:3,:   ,NGeo,:   ,iElem)
    CASE(ZETA_MINUS)
      tmp=NodeCoords(1:3,:   ,:   ,0   ,iElem)
    CASE(ZETA_PLUS)
      tmp=NodeCoords(1:3,:   ,:   ,NGeo,iElem)
  END SELECT
  ! CALL ChangeBasisSurf(3,NGeo,NGeo,Vdm_NGeo_N,tmp,tmp2)
  ! turn into right hand system of side
  DO q = 0,ZDIM(NGeo); DO p = 0,NGeo
    pq = SideToVol2(NGeo,p,q,0,iLocSide,PP_dim)
    ! Compute SurfCoords for sides
    SurfCoords(1:3,p,q,0,iSide) = tmp(:,pq(1),pq(2))
  END DO; END DO ! p,q
END DO ! iElem

ALLOCATE(Vdm_NGeo_NVisu(0:NVisu,0:NGeo))
CALL GetVandermonde(NGeo,NodeTypeVisu,NVisu   ,NodeTypeVisuPosti  ,Vdm_NGeo_NVisu   ,modal=.FALSE.)

! convert coords of DG elements
SDEALLOCATE(CoordsSurfVisu_DG)
ALLOCATE(CoordsSurfVisu_DG(3,0:NVisu,0:ZDIM(NVisu),0:0,nBCSidesVisu_DG))
DO iSide = 1,nBCSides
  iSideVisu = mapAllBCSidesToDGVisuBCSides(iSide)
  IF (iSideVisu.GT.0) THEN
    CALL ChangeBasisSurf(3,NGeo,NVisu,   Vdm_NGeo_NVisu,SurfCoords(:,:,:,0,iSide),CoordsSurfVisu_DG(:,:,:,0,iSideVisu))
  END IF
END DO

SDEALLOCATE(Vdm_NGeo_NVisu)

#if FV_ENABLED
IF (hasFV_Elems) THEN
  SWRITE (Unit_stdOut,'(A)') ' [MESH] Convert coordinates to surface visu grid (FV)'
  IF ((.NOT.changedMeshFile).AND.(.NOT.changedFV_Elems).AND.(.NOT.changedBCnames)) RETURN
  ALLOCATE(Vdm_NGeo_NVisu_FV(0:NVisu_FV,0:NGeo))
  CALL GetVandermonde(NGeo,NodeTypeVisu,NVisu_FV,NodeTypeVISUFVEqui,Vdm_NGeo_NVisu_FV,modal=.FALSE.)
  ! convert coords of FV elements
  SDEALLOCATE(CoordsSurfVisu_FV)
  ALLOCATE(CoordsSurfVisu_FV(3,0:NVisu_FV,0:ZDIM(NVisu_FV),0:0,nBCSidesVisu_FV))
  DO iSide = 1,nBCSides
    iSideVisu = mapAllBCSidesToFVVisuBCSides(iSide)
    IF (iSideVisu.GT.0) THEN
      CALL ChangeBasisSurf(3,NGeo,NVisu_FV,Vdm_NGeo_NVisu_FV,SurfCoords(:,:,:,0,iSide),CoordsSurfVisu_FV(:,:,:,0,iSideVisu))
    END IF
  END DO
  SDEALLOCATE(Vdm_NGeo_NVisu_FV)
END IF
#endif

END SUBROUTINE BuildSurfVisuCoords


!=================================================================================================================================
!> Visualize mesh only
!> 1. read mesh
!> 2. BuildVisuCoords
!> 3. Convert scaled jacobian
!> 4. write mesh to VTK array
!> 5. set length of all other output arrays to zero
!=================================================================================================================================
SUBROUTINE VisualizeMesh(postifile,meshfile_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_ReadInTools   ,ONLY: prms,GETINT,GETSTR,GETLOGICAL,CountOption
USE MOD_ReadInTools   ,ONLY: FinalizeParameters
USE MOD_StringTools   ,ONLY: STRICMP
#if USE_MPI
USE MOD_MPI           ,ONLY: FinalizeMPI
#endif
USE MOD_Interpolation ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh_Vars     ,ONLY: nElems,Ngeo,scaledJac
USE MOD_Mesh          ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_VTK           ,ONLY: WriteCoordsToVTK_array
USE MOD_HDF5_Input    ,ONLY: ReadAttribute,File_ID,OpenDataFile,CloseDataFile
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN):: postifile
CHARACTER(LEN=255),INTENT(IN):: meshfile_in
! LOCAL VARIABLES
INTEGER             :: iElem,iVar,jVar,iVarVisu,meshModeLoc
CHARACTER(LEN=255)  :: VarName
!===================================================================================================================================
#if USE_MPI
CALL FinalizeMPI()
#endif
CALL FinalizeMesh()
CALL FinalizeInterpolation()

CALL OpenDataFile(meshfile_in,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=Ngeo)
CALL CloseDataFile()

IF (LEN_TRIM(postifile).GT.0) THEN
  ! read options from parameter file
  CALL DefineParametersInterpolation()
  CALL DefineParametersMesh()
  CALL prms%SetSection("posti")
  CALL prms%CreateIntOption('NVisu', "Number of points at which solution is sampled for visualization.")
  CALL prms%read_options(postifile)
  NVisu     = GETINT('NVisu','1') ! Degree of visualization basis
  HighOrder = GETLOGICAL('HighOrder')
  ! Get number of variables to be visualized
  nVarIni = CountOption("VarName")
ELSE
  NVisu = 2*NGeo ! TODO: correct?
END IF
NVisu_FV = 1

! read mesh, depending if we should visualize the Jacobian or not different mesh modes are needed (calculate metrics or not)
meshModeLoc = 0
IF (nVarIni.GT.0) meshModeLoc=2
CALL InitInterpolation(Ngeo)
CALL InitMesh(meshMode=meshModeLoc, MeshFile_IN=meshfile_in)

! convert to visu grid
nElems_DG = nElems
nElems_FV = 0
SDEALLOCATE(mapDGElemsToAllElems)
ALLOCATE(mapDGElemsToAllElems(nElems))
DO iElem=1,nElems
  mapDGElemsToAllElems(iElem) = iElem
END DO
CALL BuildVisuCoords()
DEALLOCATE(mapDGElemsToAllElems)

! Do we need to visualize the scaled Jacobian, or the max scaled Jacobian?
IF (nVarIni.GT.0) THEN
  ! A very simple mapping is build: There are two depending variables, either one or both of them can be visualized
  NCalc = PP_N
  nVarVisu = nVarIni
  nVarDep = 2
  nVarAll = 2
  SDEALLOCATE(mapDepToCalc)
  SDEALLOCATE(mapAllVarsToVisuVars)
  SDEALLOCATE(mapAllVarsToSurfVisuVars)
  ALLOCATE(mapDepToCalc(nVarDep))
  mapDepToCalc(1) = 1
  mapDepToCalc(2) = 2
  ALLOCATE(mapAllVarsToVisuVars(nVarAll))
  mapAllVarsToVisuVars = 0
  ALLOCATE(mapAllVarsToSurfVisuVars(1:nVarAll))
  mapAllVarsToSurfVisuVars = 0
  iVarVisu = 1
  DO iVar = 1, nVarIni
    VarName = GETSTR("VarName")
    DO jVar = 1, nVarAll
      IF (STRICMP(VarNamesAll(jVar),VarName)) THEN
        mapAllVarsToVisuVars(jVar) = iVarVisu
        iVarVisu = iVarVisu + 1
      END IF
    END DO ! jVar = 1, nVarAll
  END DO ! iVar = 1, nVarIni
  SDEALLOCATE(UCalc_DG)
  ALLOCATE(UCalc_DG(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems_DG,nVarDep))
  UCalc_DG(:,:,:,:,1) = scaledJac
  DO iElem=1,nElems
    UCalc_DG(:,:,:,iElem,2) = MINVAL(UCalc_DG(:,:,:,iElem,1))
  END DO ! iElem

  CALL ConvertToVisu_DG()
ELSE
  nVarVisu = 0
END IF

CALL FinalizeInterpolation()
CALL FinalizeParameters()
END SUBROUTINE VisualizeMesh

END MODULE MOD_Posti_VisuMesh
