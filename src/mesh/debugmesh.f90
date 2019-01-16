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
!> This module contains routines for debugging and visualizing purely mesh related data.
!==================================================================================================================================
MODULE MOD_DebugMesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDebugMesh
  MODULE PROCEDURE WriteDebugMesh
END INTERFACE

PUBLIC::WriteDebugMesh
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine will compute a supersampled version of the mesh, to be used for debug purposes and includes various connectivity
!> information, which is built during the (parallel) mesh preprocessing phase. The supersampled mesh data is then output into
!> a visualization file.
!==================================================================================================================================
SUBROUTINE WriteDebugMesh(debugMesh)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_Mesh_Vars,  ONLY:nElems,Elem_xGP,ElemToSide,BC,nBCSides
USE MOD_ChangeBasisByDim,ONLY:ChangeBasisVolume
USE MOD_VTK,        ONLY:WriteDataToVTK
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN):: debugMesh !< file type to be used: 1-2: Tecplot format (deprecated), 3: Paraview format
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,iLocSide,SideID,bcindex_loc
CHARACTER(LEN=32)        :: VarNames(6)
REAL,ALLOCATABLE,TARGET  :: debugVisu(:,:,:,:,:)
REAL,POINTER             :: debugVisu_p(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET  :: X_NVisu(:,:,:,:,:)
REAL,POINTER             :: X_NVisu_p(:,:,:,:,:)
!==================================================================================================================================
IF(debugMesh.LE.0) RETURN

SWRITE(UNIT_stdOut,'(A)')' WRITE DEBUGMESH...'
! WRITE Debugmesh.vtu
ALLOCATE(X_NVisu(3,0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems))

DO iElem=1,nElems
  CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_GaussN_Nvisu,Elem_xGP(:,:,:,:,iElem),X_NVisu(:,:,:,:,iElem))
END DO

VarNames(1)='ElemID'
VarNames(2)='SideID'
VarNames(3)='FLIP'
VarNames(4)='iLocSide'
VarNames(5)='BCIndex'
VarNames(6)='Rank'
ALLOCATE(debugVisu(6,0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems))
debugVisu=-1.
DO iElem=1,nElems
  debugVisu(1,:,:,:,iElem)=REAL(iElem)
  debugVisu(6,:,:,:,iElem)=REAL(myRank)
#if (PP_dim == 3)
  DO iLocSide=1,6
#else
  DO iLocSide=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    bcindex_loc=0
    IF(SideID.LE.nBCSides) bcindex_loc=BC(SideID)
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      debugVisu(2,0,:,:,iElem)=REAL(SideID)
      debugVisu(3,0,:,:,iElem)=REAL(ElemToSide(E2S_FLIP,XI_MINUS,iElem))
      debugVisu(4,0,:,:,iElem)=REAL(iLocSide)
      debugVisu(5,0,:,:,iElem)=REAL(bcindex_loc)
    CASE(XI_PLUS)
      debugVisu(2,NVisu,:,:,iElem)=REAL(SideID)
      debugVisu(3,NVisu,:,:,iElem)=REAL(ElemToSide(E2S_FLIP,XI_PLUS,iElem))
      debugVisu(4,NVisu,:,:,iElem)=REAL(iLocSide)
      debugVisu(5,NVisu,:,:,iElem)=REAL(bcindex_loc)
    CASE(ETA_MINUS)
      debugVisu(2,:,0,:,iElem)=REAL(SideID)
      debugVisu(3,:,0,:,iElem)=REAL(ElemToSide(E2S_FLIP,ETA_MINUS,iElem))
      debugVisu(4,:,0,:,iElem)=REAL(iLocSide)
      debugVisu(5,:,0,:,iElem)=REAL(bcindex_loc)
    CASE(ETA_PLUS)
      debugVisu(2,:,NVisu,:,iElem)=REAL(SideID)
      debugVisu(3,:,NVisu,:,iElem)=REAL(ElemToSide(E2S_FLIP,ETA_PLUS,iElem))
      debugVisu(4,:,NVisu,:,iElem)=REAL(iLocSide)
      debugVisu(5,:,NVisu,:,iElem)=REAL(bcindex_loc)
    CASE(ZETA_MINUS)
      debugVisu(2,:,:,0,iElem)=REAL(SideID)
      debugVisu(3,:,:,0,iElem)=REAL(ElemToSide(E2S_FLIP,ZETA_MINUS,iElem))
      debugVisu(4,:,:,0,iElem)=REAL(iLocSide)
      debugVisu(5,:,:,0,iElem)=REAL(bcindex_loc)
    CASE(ZETA_PLUS)
      debugVisu(2,:,:,NVisu,iElem)=REAL(SideID)
      debugVisu(3,:,:,NVisu,iElem)=REAL(ElemToSide(E2S_FLIP,ZETA_PLUS,iElem))
      debugVisu(4,:,:,NVisu,iElem)=REAL(iLocSide)
      debugVisu(5,:,:,NVisu,iElem)=REAL(bcindex_loc)
    END SELECT
  END DO
END DO

! Tecplot format is deprecated and will probably removed in future versions
SELECT CASE(debugMesh)
CASE(1)
  STOP 'Tecplot Output removed (license issues)'
CASE(2)
  STOP 'Tecplot Output removed (license issues)'
CASE(3)
  X_NVisu_p => X_NVisu
  debugVisu_p => debugVisu
  CALL WriteDataToVTK(6,NVisu,nElems,VarNames,X_NVisu_p,debugVisu_p,'Debugmesh.vtu',dim=PP_dim)
END SELECT

DEALLOCATE(X_NVisu)
DEALLOCATE(debugVisu)
SWRITE(UNIT_stdOut,'(A)')' WRITE DEBUGMESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WriteDebugMesh

END MODULE MOD_DebugMesh
