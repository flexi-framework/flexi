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
!> Module containing the main procedures for the POSTI tool: visu3d_requestInformation is called by ParaView to create a
!> list of available variables and visu3D is the main routine of POSTI.
!===================================================================================================================================
MODULE MOD_Visu3D_Cwrapper
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE visu3d_requestInformation
  MODULE PROCEDURE visu3d_requestInformation
END INTERFACE

INTERFACE visu3D_CWrapper
  MODULE PROCEDURE visu3D_CWrapper
END INTERFACE

INTERFACE visu3d_dealloc_nodeids
  MODULE PROCEDURE visu3d_dealloc_nodeids
END INTERFACE

PUBLIC:: visu3d_requestInformation
PUBLIC:: visu3D_CWrapper
PUBLIC:: visu3d_dealloc_nodeids

CONTAINS


FUNCTION cstrToChar255(cstr, strlen) 
USE ISO_C_BINDING
! INPUT / OUTPUT VARIABLES 
TYPE(C_PTR),TARGET,INTENT(IN)  :: cstr
INTEGER,INTENT(IN)             :: strlen
CHARACTER(LEN=255)             :: cstrToChar255
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(KIND=C_CHAR),POINTER :: tmp(:)
!===================================================================================================================================
CALL C_F_POINTER(C_LOC(cstr), tmp, [strlen])
cstrToChar255 = TRANSFER(tmp(1:strlen), cstrToChar255)
cstrToChar255(strlen+1:255) = ' ' 
END FUNCTION cstrToChar255

!===================================================================================================================================
!> Wrapper to visu3D_InitFile for Paraview plugin
!===================================================================================================================================
SUBROUTINE visu3d_requestInformation(mpi_comm_IN, strlen_state, statefile_IN, varnames, bcnames)
! MODULES
USE MOD_Globals
USE MOD_MPI        ,ONLY: InitMPI
USE MOD_Posti_Vars ,ONLY: VarNamesTotal,BoundaryNamesTotal
USE MOD_Visu3D     ,ONLY: visu3d_getVarNamesAndFileType
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                    :: mpi_comm_IN    
INTEGER,INTENT(IN)                    :: strlen_state
TYPE(C_PTR),TARGET,INTENT(IN)         :: statefile_IN
TYPE (CARRAY), INTENT(INOUT)          :: varnames
TYPE (CARRAY), INTENT(INOUT)          :: bcnames
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
CHARACTER(LEN=255)                    :: statefile
CHARACTER(LEN=255),POINTER            :: varnames_pointer(:)
CHARACTER(LEN=255),POINTER            :: bcnames_pointer(:)
!===================================================================================================================================
statefile = cstrToChar255(statefile_IN, strlen_state)

CALL InitMPI(mpi_comm_IN) 
CALL visu3d_getVarNamesAndFileType(statefile,VarNamesTotal,BoundaryNamesTotal)
IF (ALLOCATED(VarNamesTotal)) THEN
  varnames_pointer => VarNamesTotal
  varnames%len  = SIZE(varnames_pointer)*255 
  varnames%data = C_LOC(varnames_pointer(1))
ELSE
  varnames%len  = 0
  varnames%data = C_NULL_PTR
END IF
IF (ALLOCATED(BoundaryNamesTotal)) THEN
  bcnames_pointer => BoundaryNamesTotal
  bcnames%len  = SIZE(bcnames_pointer)*255 
  bcnames%data = C_LOC(bcnames_pointer(1))
ELSE
  bcnames%len  = 0
  bcnames%data = C_NULL_PTR
END IF
END SUBROUTINE visu3d_requestInformation

!===================================================================================================================================
!> C wrapper routine for the visu3D call from ParaView.
!===================================================================================================================================
SUBROUTINE visu3D_CWrapper(mpi_comm_IN, &
    strlen_prm, prmfile_IN, strlen_posti, postifile_IN, strlen_state, statefile_IN,&
    coordsDG_out,valuesDG_out,nodeidsDG_out, &
    coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out, &
    coordsSurfDG_out,valuesSurfDG_out,nodeidsSurfDG_out, &
    coordsSurfFV_out,valuesSurfFV_out,nodeidsSurfFV_out,varnamesSurf_out)
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_Visu3D      ,ONLY: visu3D
USE MOD_VTK         ,ONLY: WriteCoordsToVTK_array,WriteDataToVTK_array,WriteVarnamesToVTK_array
USE MOD_Output_Vars ,ONLY: ProjectName
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: mpi_comm_IN    
INTEGER,INTENT(IN)            :: strlen_prm    
INTEGER,INTENT(IN)            :: strlen_posti    
INTEGER,INTENT(IN)            :: strlen_state    
TYPE(C_PTR),TARGET,INTENT(IN) :: prmfile_IN
TYPE(C_PTR),TARGET,INTENT(IN) :: postifile_IN
TYPE(C_PTR),TARGET,INTENT(IN) :: statefile_IN
TYPE (CARRAY), INTENT(INOUT)  :: coordsDG_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesDG_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsDG_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsFV_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesFV_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsFV_out
TYPE (CARRAY), INTENT(INOUT)  :: varnames_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsSurfDG_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesSurfDG_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsSurfDG_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsSurfFV_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesSurfFV_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsSurfFV_out
TYPE (CARRAY), INTENT(INOUT)  :: varnamesSurf_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: prmfile
CHARACTER(LEN=255)            :: postifile
CHARACTER(LEN=255)            :: statefile
!===================================================================================================================================
prmfile   = cstrToChar255(prmfile_IN,   strlen_prm)
postifile = cstrToChar255(postifile_IN, strlen_posti)
statefile = cstrToChar255(statefile_IN, strlen_state)
CALL visu3D(mpi_comm_IN, prmfile, postifile, statefile)

! Map Fortran arrays to C pointer

! write UVisu to VTK 2D / 3D arrays (must be done always!)
! write coords, UVisu to VTK  2D / 3D arrays (must be done always!)
IF (VisuDimension.EQ.3) THEN
  ! Volume
  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG,3)
  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV,3)

  CALL WriteCoordsToVTK_array(NVisu   ,nElems_DG,coordsDG_out,nodeidsDG_out,&
      CoordsVisu_DG,nodeids_DG,dim=3,DGFV=0)
  CALL WriteCoordsToVTK_array(NVisu_FV,nElems_FV,coordsFV_out,nodeidsFV_out,&
      CoordsVisu_FV,nodeids_FV,dim=3,DGFV=1)

  CALL WriteVarnamesToVTK_array(nVarTotal,mapVisu,varnames_out,VarNamesTotal,nVarVisuTotal)

  ! Surface
  CALL WriteDataToVTK_array(nVarSurfVisuTotal,NVisu   ,nBCSidesVisu_DG,valuesSurfDG_out,USurfVisu_DG,2)
  CALL WriteDataToVTK_array(nVarSurfVisuTotal,NVisu_FV,nBCSidesVisu_FV,valuesSurfFV_out,USurfVisu_FV,2)

  CALL WriteCoordsToVTK_array(NVisu   ,nBCSidesVisu_DG,coordsSurfDG_out,nodeidsSurfDG_out,&
      CoordsSurfVisu_DG,nodeidsSurf_DG,dim=2,DGFV=0)
  CALL WriteCoordsToVTK_array(NVisu_FV,nBCSidesVisu_FV,coordsSurfFV_out,nodeidsSurfFV_out,&
      CoordsSurfVisu_FV,nodeidsSurf_FV,dim=2,DGFV=1)

  CALL WriteVarnamesToVTK_array(nVarTotal,mapSurfVisu,varnamesSurf_out,VarNamesTotal,nVarSurfVisuTotal)
ELSE IF (VisuDimension.EQ.2) THEN
  STOP 'implement avg2d'

  ! allocate Visu 2D array and copy from first zeta-slice of 3D array
  SDEALLOCATE(UVisu_DG_2D)
  ALLOCATE(UVisu_DG_2D(0:NVisu,0:NVisu,0:0,1:nElems_DG,1:(nVarVisuTotal)))
  UVisu_DG_2D = UVisu_DG(:,:,0:0,:,:)
  SDEALLOCATE(UVisu_FV_2D)
  ALLOCATE(UVisu_FV_2D(0:NVisu_FV,0:NVisu_FV,0:0,1:nElems_FV,1:(nVarVisuTotal)))
#if FV_ENABLED
  UVisu_FV_2D = UVisu_FV(:,:,0:0,:,:)
#else
  CoordsVisu_FV_2D = 0
#endif

  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG_2D,2)
  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV_2D,2)

  ! allocate Coords 2D array and copy from first zeta-slice of 3D array
  SDEALLOCATE(CoordsVisu_DG_2D)
  ALLOCATE(CoordsVisu_DG_2D(1:3,0:NVisu,0:NVisu,0:0,1:nElems_DG))
  CoordsVisu_DG_2D = CoordsVisu_DG(:,:,:,0:0,:)
  SDEALLOCATE(CoordsVisu_FV_2D)
  ALLOCATE(CoordsVisu_FV_2D(1:3,0:NVisu_FV,0:NVisu_FV,0:0,1:nElems_FV))
#if FV_ENABLED    
  CoordsVisu_FV_2D = CoordsVisu_FV(:,:,:,0:0,:)
#else
  CoordsVisu_FV_2D = 0
#endif

  CALL WriteCoordsToVTK_array(NVisu   ,nElems_DG,coordsDG_out,nodeidsDG_out,&
      CoordsVisu_DG_2D,nodeids_DG_2D,dim=2,DGFV=0)
  CALL WriteCoordsToVTK_array(NVisu_FV,nElems_FV,coordsFV_out,nodeidsFV_out,&
      CoordsVisu_FV_2D,nodeids_FV_2D,dim=2,DGFV=1)

  CALL WriteVarnamesToVTK_array(nVarTotal,mapVisu,varnames_out,VarNamesTotal,nVarVisuTotal)
END IF

END SUBROUTINE visu3D_CWrapper

!===================================================================================================================================
!> Deallocate the different NodeID arrays.
!===================================================================================================================================
SUBROUTINE visu3d_dealloc_nodeids() 
USE MOD_Posti_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(nodeids_DG)
SDEALLOCATE(nodeids_FV)
SDEALLOCATE(nodeids_DG_2D)
SDEALLOCATE(nodeids_FV_2D)
END SUBROUTINE visu3d_dealloc_nodeids

END MODULE MOD_Visu3D_Cwrapper
