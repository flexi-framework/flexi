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
!> This module contains additional mortar operators for FV subcells.
!==================================================================================================================================
MODULE MOD_FV_Mortar
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------
#if FV_ENABLED
INTERFACE FV_Elems_Mortar
  MODULE PROCEDURE FV_Elems_Mortar
END INTERFACE

INTERFACE FV_gradU_mortar
  MODULE PROCEDURE FV_gradU_mortar
END INTERFACE

PUBLIC::FV_Elems_Mortar
PUBLIC::FV_gradU_mortar
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Copy FV_Elems information from big mortar sides to the small sides. Compare to U_mortar subroutine.
!==================================================================================================================================
SUBROUTINE FV_Elems_Mortar(FV_Elems_master,FV_Elems_slave,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars ,ONLY: nSides
USE MOD_Mesh_Vars ,ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars ,ONLY: firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(INOUT) :: FV_Elems_master(1:nSides) !< master side FV_elems
INTEGER,INTENT(INOUT) :: FV_Elems_slave( 1:nSides) !< slave  side FV_elems
LOGICAL,INTENT(IN)    :: doMPISides                !< =.TRUE. only MPI sides are filled, =.FALSE. inner sides 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: iMortar,nMortars
INTEGER :: firstMortarSideID,lastMortarSideID
INTEGER :: MortarSideID,SideID,locSide
!==================================================================================================================================
!                         doMPISides==True   doMPISides==False
firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides) 
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides) 

DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide =MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
    SELECT CASE(MortarInfo(MI_FLIP,iMortar,locSide))
      CASE(0) ! master side
        FV_Elems_master(SideID) = FV_Elems_master(MortarSideID)
      CASE(1:4) ! slave side
        FV_Elems_slave( SideID) = FV_Elems_master(MortarSideID)
    END SELECT !flip(iMortar)
  END DO !iMortar 
END DO !MortarSideID
END SUBROUTINE FV_Elems_Mortar

!==================================================================================================================================
!> Fill master/big mortar sides parts of FV_surf_gradU array.
!> Comparable to Flux_Mortar routine, but not totally the same
!==================================================================================================================================
SUBROUTINE FV_gradU_mortar(FV_surf_gradU,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars      ,ONLY: MortarType,nSides
USE MOD_Mesh_Vars      ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars      ,ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_FV_Vars        ,ONLY: FV_Elems_master
USE MOD_FillMortarPrim ,ONLY: Flux_MortarPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT) :: FV_surf_gradU(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) ! slope over interface
LOGICAL,INTENT(IN) :: doMPISides  !< =.TRUE. only MPI sides are filled, =.FALSE. inner sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: firstMortarSideID,lastMortarSideID
INTEGER :: MortarSideID
!==================================================================================================================================
! Attention: we only have one Flux_surf_gradU (no master/slave) 
!            => input it to Flux_Mortar for both fluxes (master/slave)
CALL Flux_MortarPrim(FV_surf_gradU,FV_surf_gradU,doMPISides,weak=.FALSE.,onlyFV=.TRUE.)

!                         doMPISides==True   doMPISides==False
firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides) 
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides) 

DO MortarSideID=firstMortarSideID,lastMortarSideID
  IF(FV_Elems_master(MortarSideID).EQ.0) CYCLE ! DG

#if (PP_dim == 3)
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    FV_surf_gradU(:,:,:,MortarSideID) = 0.25 * FV_surf_gradU(:,:,:,MortarSideID)
  CASE(2) !1->2 in eta
    FV_surf_gradU(:,:,:,MortarSideID) = 0.5  * FV_surf_gradU(:,:,:,MortarSideID)
  CASE(3) !1->2 in xi
#endif
    FV_surf_gradU(:,:,:,MortarSideID) = 0.5  * FV_surf_gradU(:,:,:,MortarSideID)
#if (PP_dim == 3)
  END SELECT
#endif
END DO

END SUBROUTINE FV_gradU_mortar
#endif

END MODULE MOD_FV_Mortar
