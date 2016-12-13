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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE BuildVisuCoords
  MODULE PROCEDURE BuildVisuCoords
END INTERFACE

INTERFACE VisualizeMesh
  MODULE PROCEDURE VisualizeMesh
END INTERFACE

PUBLIC:: BuildVisuCoords
PUBLIC:: VisualizeMesh

CONTAINS

!=================================================================================================================================
!> Converts the coordinates of the mesh to the visu-mesh.
!=================================================================================================================================
SUBROUTINE BuildVisuCoords(withGradients)
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars         ,ONLY: NVisu,nElems_DG,mapElems_DG,NVisu_FV,nElems_FV,mapElems_FV,hasFV_Elems
USE MOD_Posti_Vars         ,ONLY: CoordsVisu_DG
#if FV_ENABLED
USE MOD_Posti_Vars         ,ONLY: CoordsVisu_FV,changedMeshFile,changedFV_Elems
#endif
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeFVEqui,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Mesh_Vars          ,ONLY: nElems,Elem_xGP
IMPLICIT NONE
LOGICAL,INTENT(IN) :: withGradients
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem, iElem_DG
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)    
#if FV_ENABLED
INTEGER            :: iElem_FV
REAL,ALLOCATABLE   :: Vdm_N_NVisu_FV(:,:)
#endif
CHARACTER(LEN=255) :: NodeType_loc
INTEGER            :: Nloc
REAL,POINTER       :: NodeCoords_loc(:,:,:,:,:)
!===================================================================================================================================
! Always use Elem_xGP for visualization in case TreeMappings are used or N<NGeo
Nloc = PP_N
NodeType_loc = NodeType
NodeCoords_loc(1:3,0:Nloc,0:Nloc,0:Nloc,1:nElems) => Elem_xGP

! Convert coordinates to visu grid
SWRITE (*,*) "[MESH] Convert coordinates to visu grid (DG)"
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:Nloc))
CALL GetVandermonde(Nloc,NodeType_loc,NVisu   ,NodeTypeVisu  ,Vdm_N_NVisu   ,modal=.FALSE.)
! convert coords of DG elements
SDEALLOCATE(CoordsVisu_DG)
ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:NVisu,nElems_DG))
DO iElem_DG = 1,nElems_DG
  iElem = mapElems_DG(iElem_DG)
  CALL ChangeBasis3D(3,Nloc,NVisu,   Vdm_N_NVisu,   NodeCoords_loc(:,:,:,:,iElem),CoordsVisu_DG   (:,:,:,:,iElem_DG))
END DO
SDEALLOCATE(Vdm_N_NVisu)

#if FV_ENABLED
IF (hasFV_Elems) THEN
  SWRITE (*,*) "[MESH] Convert coordinates to visu grid (FV)"
  IF ((.NOT.changedMeshFile).AND.(.NOT.changedFV_Elems)) RETURN ! only NVisu changed, but NVisu_FV is independent of NVisu
  ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:Nloc))
  CALL GetVandermonde(Nloc,NodeType_loc,NVisu_FV,NodeTypeFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
  ! convert coords of FV elements
  SDEALLOCATE(CoordsVisu_FV)
  ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:NVisu_FV,nElems_FV))
  DO iElem_FV = 1,nElems_FV
    iElem = mapElems_FV(iElem_FV)
    CALL ChangeBasis3D(3,Nloc,NVisu_FV,Vdm_N_NVisu_FV,NodeCoords_loc(:,:,:,:,iElem),CoordsVisu_FV(:,:,:,:,iElem_FV))
  END DO
  SDEALLOCATE(Vdm_N_NVisu_FV)
END IF
#endif

END SUBROUTINE BuildVisuCoords


!=================================================================================================================================
!> Visualize mesh only
!> 1. read mesh 
!> 2. BuildVisuCoords
!> 3. write mesh to VTK array
!> 4. set length of all other output arrays to zero
!=================================================================================================================================
SUBROUTINE VisualizeMesh(postifile,meshfile_in,coordsDG_out,valuesDG_out,nodeidsDG_out, &
        coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)
! MODULES                                                                   
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_ReadInTools   ,ONLY: prms,GETINT
USE MOD_ReadInTools   ,ONLY: FinalizeParameters
USE MOD_MPI           ,ONLY: FinalizeMPI
USE MOD_Interpolation ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh_Vars     ,ONLY: nElems,Ngeo
USE MOD_Mesh          ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_VTK           ,ONLY: WriteCoordsToVTK_array
USE MOD_HDF5_Input    ,ONLY: ReadAttribute,File_ID,OpenDataFile,CloseDataFile
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),INTENT(IN):: postifile
CHARACTER(LEN=255),INTENT(IN):: meshfile_in
TYPE (CARRAY), INTENT(INOUT) :: coordsDG_out
TYPE (CARRAY), INTENT(INOUT) :: valuesDG_out
TYPE (CARRAY), INTENT(INOUT) :: nodeidsDG_out
TYPE (CARRAY), INTENT(INOUT) :: coordsFV_out
TYPE (CARRAY), INTENT(INOUT) :: valuesFV_out
TYPE (CARRAY), INTENT(INOUT) :: nodeidsFV_out
TYPE (CARRAY), INTENT(INOUT) :: varnames_out
TYPE (CARRAY), INTENT(INOUT) :: components_out
! LOCAL VARIABLES
INTEGER             :: iElem
!===================================================================================================================================
#if USE_MPI
CALL FinalizeMPI()
#endif
CALL FinalizeMesh()
CALL FinalizeInterpolation()

CALL OpenDataFile(meshfile_in,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=Ngeo)
CALL CloseDataFile()

IF (LEN_TRIM(postifile).GT.0) THEN
  ! read options from parameter file
  CALL DefineParametersInterpolation()
  CALL DefineParametersMesh()
  CALL prms%SetSection("posti")
  CALL prms%CreateIntOption('NVisu', "Number of points at which solution is sampled for visualization.")
  CALL prms%read_options(postifile)
  NVisu = GETINT('NVisu') ! Degree of visualization basis
ELSE
  NVisu = 2*NGeo ! TODO: correct?
END IF

! read mesh
CALL InitInterpolation(Ngeo)
CALL InitMesh(meshMode=0, MeshFile_IN=meshfile_in)

! convert to visu grid
nElems_DG = nElems
nElems_FV = 0
SDEALLOCATE(mapElems_DG)
ALLOCATE(mapElems_DG(nElems))
DO iElem=1,nElems
  mapElems_DG(iElem) = iElem
END DO
CALL BuildVisuCoords(.FALSE.)
DEALLOCATE(mapElems_DG)

! write to VTK array
CALL WriteCoordsToVTK_array(NVisu,nElems,coordsDG_out,nodeidsDG_out,CoordsVisu_DG,nodeids_DG,dim=3,DGFV=0)
valuesDG_out%len  =0

! set length of all other output arrays to zero
coordsFV_out%len  =0
nodeidsFV_out%len =0
valuesFV_out%len  =0
varnames_out%len  =0
components_out%len=0

CALL FinalizeInterpolation()
CALL FinalizeParameters()
END SUBROUTINE VisualizeMesh

END MODULE MOD_Posti_VisuMesh
