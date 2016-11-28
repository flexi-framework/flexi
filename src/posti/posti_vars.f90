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
!==================================================================================================================================
!> Contains global variables provided by the posti routines 
!==================================================================================================================================
MODULE MOD_Posti_Vars
USE ISO_C_BINDING
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
CHARACTER(LEN=255)                :: prmfile_old = ""
CHARACTER(LEN=255)                :: statefile_old = ""
CHARACTER(LEN=255)                :: MeshFile
CHARACTER(LEN=255)                :: MeshFile_old = ""
INTEGER                           :: NVisu
INTEGER                           :: NVisu_old = -1
INTEGER                           :: nVar_State
INTEGER                           :: nVar_State_old = -1
INTEGER                           :: nElems_DG
LOGICAL                           :: withGradients
LOGICAL                           :: withGradients_old = .FALSE.
INTEGER                           :: nElems_FV
INTEGER                           :: NVisu_FV
REAL                              :: OutputTime
LOGICAL                           :: hasFV_Elems = .FALSE.
INTEGER                           :: VisuDimension

LOGICAL                           :: changedStateFile
LOGICAL                           :: changedMeshFile
LOGICAL                           :: changedNVisu
LOGICAL                           :: changedVarNames
LOGICAL                           :: changedFV_Elems
LOGICAL                           :: changedWithGradients
LOGICAL                           :: changedPrmFile

INTEGER,ALLOCATABLE,TARGET        :: nodeids_DG(:)           ! visu nodeids
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: CoordsVisu_DG(:,:,:,:,:)! visu coordinates
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: UVisu_DG(:,:,:,:,:)     ! state at visu points
INTEGER,ALLOCATABLE,TARGET        :: nodeids_FV(:)           ! visu nodeids
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: CoordsVisu_FV(:,:,:,:,:)! visu coordinates
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: UVisu_FV(:,:,:,:,:)     ! state at visu points
INTEGER                           :: nVarCalc,nVarVisu
INTEGER,ALLOCATABLE               :: mapCalc(:)
#if FV_ENABLED && FV_RECONSTRUCT
INTEGER                           :: nVarCalc_FV
INTEGER,ALLOCATABLE               :: mapCalc_FV(:)
#endif
INTEGER,ALLOCATABLE               :: mapVisu(:)
INTEGER,ALLOCATABLE               :: mapVisu_old(:)
REAL,ALLOCATABLE                  :: UCalc_DG(:,:,:,:,:)
REAL,ALLOCATABLE                  :: UCalc_FV(:,:,:,:,:)

INTEGER,ALLOCATABLE               :: mapElems_DG(:)
INTEGER,ALLOCATABLE               :: mapElems_FV(:)

INTEGER,ALLOCATABLE               :: FV_Elems_loc(:)
INTEGER,ALLOCATABLE               :: FV_Elems_old(:)

INTEGER                           :: nVarVisu_ElemData, nVarVisu_ElemData_old=0
CHARACTER(LEN=255),ALLOCATABLE    :: VarNamesVisu_ElemData(:)
CHARACTER(LEN=255),ALLOCATABLE    :: VarNamesVisu_ElemData_old(:)

INTEGER,ALLOCATABLE,TARGET        :: nodeids_DG_2D(:)           ! visu nodeids
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: CoordsVisu_DG_2D(:,:,:,:,:)! visu coordinates
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: UVisu_DG_2D(:,:,:,:,:)     ! state at visu points
INTEGER,ALLOCATABLE,TARGET        :: nodeids_FV_2D(:)           ! visu nodeids
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: CoordsVisu_FV_2D(:,:,:,:,:)! visu coordinates
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: UVisu_FV_2D(:,:,:,:,:)     ! state at visu points
END MODULE MOD_Posti_Vars
