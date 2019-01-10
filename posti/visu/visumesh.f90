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
USE MOD_Visu_Vars          ,ONLY: CoordsVisu_DG
USE MOD_Visu_Vars          ,ONLY: NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu,nElems_DG,mapDGElemsToAllElems
USE MOD_Visu_Vars          ,ONLY: nElemsAvg2D_DG,Avg2D
USE MOD_Visu_Vars          ,ONLY: Elem_IJK,mapElemIJToDGElemAvg2D
#if FV_ENABLED
USE MOD_Visu_Vars          ,ONLY: FVAmountAvg2D,mapElemIJToFVElemAvg2D,nElemsAvg2D_FV
USE MOD_Visu_Vars          ,ONLY: NVisu_FV,nElems_FV,mapFVElemsToAllElems,hasFV_Elems
USE MOD_Visu_Vars          ,ONLY: CoordsVisu_FV,changedMeshFile,changedFV_Elems,changedAvg2D
#endif
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeVISUFVEqui,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_Mesh_Vars          ,ONLY: nElems,Elem_xGP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem, iElem_DG
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)    
#if FV_ENABLED
INTEGER            :: iElem_FV
REAL,ALLOCATABLE   :: Vdm_N_NVisu_FV(:,:)
#endif
INTEGER            :: iElemAvg,ii,jj,kk
!===================================================================================================================================
! Convert coordinates to visu grid
SWRITE (*,*) "[MESH] Convert coordinates to visu grid (DG)"
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu   ,NodeTypeVisuPosti  ,Vdm_N_NVisu   ,modal=.FALSE.)
#if FV_ENABLED
ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu_FV,NodeTypeVISUFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
#endif
! convert coords of DG elements
IF (Avg2D) THEN
  SDEALLOCATE(CoordsVisu_DG)
  ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:0,nElemsAvg2D_DG))
#if FV_ENABLED
  SDEALLOCATE(CoordsVisu_FV)
  ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV))
#endif
  DO iElem = 1,nElems
    ii = Elem_IJK(1,iElem)
    jj = Elem_IJK(2,iElem)
    kk = Elem_IJK(3,iElem)
    IF (kk.EQ.1) THEN
#if FV_ENABLED    
      IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN ! DG
#endif
        iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
        CALL ChangeBasis2D(3,PP_N,NVisu,Vdm_N_NVisu,Elem_xGP(:,:,:,0,iElem),CoordsVisu_DG(:,:,:,0,iElemAvg))
#if FV_ENABLED    
      ELSE ! FV
        iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
        CALL ChangeBasis2D(3,PP_N,NVisu_FV,Vdm_N_NVisu_FV,Elem_xGP(:,:,:,0,iElem),CoordsVisu_FV(:,:,:,0,iElemAvg))
      END IF
#endif
    END IF
  END DO
ELSE
  SDEALLOCATE(CoordsVisu_DG)
  ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems_DG))
  DO iElem_DG = 1,nElems_DG
    iElem = mapDGElemsToAllElems(iElem_DG)
    CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_N_NVisu,Elem_xGP(:,:,:,:,iElem),CoordsVisu_DG(:,:,:,:,iElem_DG))
  END DO

#if FV_ENABLED
  IF (hasFV_Elems) THEN
    SWRITE (*,*) "[MESH] Convert coordinates to visu grid (FV)"
    ! only NVisu changed, but NVisu_FV is independent of NVisu
    IF ((.NOT.changedMeshFile).AND.(.NOT.changedFV_Elems).AND.(.NOT.changedAvg2D)) RETURN 
    !ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:PP_N))
    !CALL GetVandermonde(PP_N,NodeType,NVisu_FV,NodeTypeVISUFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
    ! convert coords of FV elements
    SDEALLOCATE(CoordsVisu_FV)
    ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),nElems_FV))
    DO iElem_FV = 1,nElems_FV
      iElem = mapFVElemsToAllElems(iElem_FV)
      CALL ChangeBasisVolume(3,PP_N,NVisu_FV,Vdm_N_NVisu_FV,Elem_xGP(:,:,:,:,iElem),CoordsVisu_FV(:,:,:,:,iElem_FV))
    END DO
    SDEALLOCATE(Vdm_N_NVisu_FV)
  END IF
#endif

END IF
SDEALLOCATE(Vdm_N_NVisu)

END SUBROUTINE BuildVisuCoords

!=================================================================================================================================
!> Converts the coordinates of the surface mesh to the surface visu-mesh.
!=================================================================================================================================
SUBROUTINE BuildSurfVisuCoords()
! MODULES                                                                   
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars          ,ONLY: CoordsSurfVisu_DG,nBCSidesVisu_DG,mapAllBCSidesToDGVisuBCSides
USE MOD_Visu_Vars          ,ONLY: NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu
#if FV_ENABLED
USE MOD_Visu_Vars          ,ONLY: CoordsSurfVisu_FV,nBCSidesVisu_FV,mapAllBCSidesToFVVisuBCSides
USE MOD_Visu_Vars          ,ONLY: NVisu_FV,hasFV_Elems
USE MOD_Visu_Vars          ,ONLY: changedMeshFile,changedFV_Elems,changedBCnames
#endif
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeVISUFVEqui,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
USE MOD_Mesh_Vars          ,ONLY: Face_xGP,nBCSides
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSide,iSideVisu
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)    
#if FV_ENABLED
REAL,ALLOCATABLE   :: Vdm_N_NVisu_FV(:,:)
#endif
CHARACTER(LEN=255) :: NodeType_loc
INTEGER            :: Nloc
!===================================================================================================================================
Nloc = PP_N
NodeType_loc = NodeType

! Convert coordinates to visu grid
SWRITE (*,*) "[MESH] Convert coordinates to surface visu grid (DG)"
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:Nloc))
CALL GetVandermonde(Nloc,NodeType_loc,NVisu   ,NodeTypeVisuPosti  ,Vdm_N_NVisu   ,modal=.FALSE.)
! convert coords of DG elements
SDEALLOCATE(CoordsSurfVisu_DG)
ALLOCATE(CoordsSurfVisu_DG(3,0:NVisu,0:ZDIM(NVisu),0:0,nBCSidesVisu_DG))
DO iSide=1,nBCSides
  iSideVisu = mapAllBCSidesToDGVisuBCSides(iSide)
  IF (iSideVisu.GT.0)THEN
    CALL ChangeBasisSurf(3,Nloc,NVisu,   Vdm_N_NVisu, Face_xGP(:,:,:,0,iSide),CoordsSurfVisu_DG(:,:,:,0,iSideVisu))
  END IF
END DO
SDEALLOCATE(Vdm_N_NVisu)

#if FV_ENABLED
IF (hasFV_Elems) THEN
  SWRITE (*,*) "[MESH] Convert coordinates to surface visu grid (FV)"
  IF ((.NOT.changedMeshFile).AND.(.NOT.changedFV_Elems).AND.(.NOT.changedBCnames)) RETURN 
  ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:Nloc))
  CALL GetVandermonde(Nloc,NodeType_loc,NVisu_FV,NodeTypeVISUFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
  ! convert coords of FV elements
  SDEALLOCATE(CoordsSurfVisu_FV)
  ALLOCATE(CoordsSurfVisu_FV(3,0:NVisu_FV,0:ZDIM(NVisu_FV),0:0,nBCSidesVisu_FV))
  DO iSide=1,nBCSides
    iSideVisu = mapAllBCSidesToFVVisuBCSides(iSide)
    IF (iSideVisu.GT.0)THEN
      CALL ChangeBasisSurf(3,Nloc,NVisu_FV,Vdm_N_NVisu_FV,Face_xGP(:,:,:,0,iSide),CoordsSurfVisu_FV(:,:,:,0,iSideVisu))
    END IF
  END DO
  SDEALLOCATE(Vdm_N_NVisu_FV)
END IF
#endif

END SUBROUTINE BuildSurfVisuCoords


!=================================================================================================================================
!> Visualize mesh only
!> 1. read mesh 
!> 2. BuildVisuCoords
!> 3. write mesh to VTK array
!> 4. set length of all other output arrays to zero
!=================================================================================================================================
SUBROUTINE VisualizeMesh(postifile,meshfile_in)
! MODULES                                                                   
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_ReadInTools   ,ONLY: prms,GETINT
USE MOD_ReadInTools   ,ONLY: FinalizeParameters
#if USE_MPI
USE MOD_MPI           ,ONLY: FinalizeMPI
#endif
USE MOD_Interpolation ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh_Vars     ,ONLY: nElems,Ngeo
USE MOD_Mesh          ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_VTK           ,ONLY: WriteCoordsToVTK_array
USE MOD_HDF5_Input    ,ONLY: ReadAttribute,File_ID,OpenDataFile,CloseDataFile
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),INTENT(IN):: postifile
CHARACTER(LEN=255),INTENT(IN):: meshfile_in
! LOCAL VARIABLES
INTEGER             :: iElem
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
  NVisu = GETINT('NVisu','1') ! Degree of visualization basis
ELSE
  NVisu = 2*NGeo ! TODO: correct?
END IF
NVisu_FV = 1

! read mesh
CALL InitInterpolation(Ngeo)
CALL InitMesh(meshMode=0, MeshFile_IN=meshfile_in)

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

CALL FinalizeInterpolation()
CALL FinalizeParameters()
END SUBROUTINE VisualizeMesh

END MODULE MOD_Posti_VisuMesh
