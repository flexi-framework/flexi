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
#if FV_ENABLED
#include "flexi.h"

!==================================================================================================================================
!> Prolongate solution of FV subcells to the DG element boundary using the slopes of the 2nd order reconstruction.
!==================================================================================================================================
MODULE MOD_FV_ProlongToFace
#if FV_RECONSTRUCT
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE FV_ProlongToDGFace
  MODULE PROCEDURE FV_ProlongToDGFace
END INTERFACE

PUBLIC::FV_ProlongToDGFace
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Prolongate solution of FV subcells to the DG element boundary using the slopes of the 2nd order reconstruction.
!> Before the call of this function the arrays UPrim_master and UPrim_slave only contain the cell average of the FV subcells
!> (1st order data copied from volume to face data in the prolongtoface.f90, but without reconstruction)
!> When this routine is called the slopes over DG element interfaces are computed but not limited. Therefore first these
!> slopes are limited and then used to calculate the reconstructed solutions UPrim_master/slave at the DG element interfaces.
!==================================================================================================================================
SUBROUTINE FV_ProlongToDGFace(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,FV_surf_gradU,&
    doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_FV_Vars   ,ONLY: FV_Elems_master,FV_Elems_slave
USE MOD_FV_Vars   ,ONLY: FV_dx_master,FV_dx_slave
USE MOD_FV_Limiter,ONLY: FV_Limiter
USE MOD_Mesh_Vars ,ONLY: firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars ,ONLY: firstInnerSide,lastInnerSide,nSides,firstBCSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  !< =.TRUE. only MPI sides are filled, =.FALSE. inner sides
REAL,INTENT(INOUT) :: UPrim_master   (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive master solution (without reconstruction)
REAL,INTENT(INOUT) :: UPrim_slave    (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive slave solution (without reconstruction)
REAL,INTENT(INOUT) :: FV_multi_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< first inner slope of the master element
REAL,INTENT(INOUT) :: FV_multi_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< first inner slope of the slave element
REAL,INTENT(IN)    :: FV_surf_gradU  (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< slope over the interface
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: SideID,firstSideID,lastSideID
INTEGER :: p,q
REAL    :: gradU(PP_nVarPrim,0:PP_N,0:PP_NZ)
!==================================================================================================================================
! reconstruct UPrim_master/slave for sides ranging between firstSideID and lastSideID
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
  lastSideID  = lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstBCSide
  lastSideID  = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  IF (FV_Elems_master(SideID).GT.0) THEN ! FV element
    ! FV_multi_master contains the first inner slope (from sub-cell 0 to sub-cell 1 in normal direction to the interface)
    ! see 'FV_PrepareSurfGradient' subroutine in fv_reconstruction.f90 for details.
    ! FV_surf_gradU contains the gradient over the DG element interface and is calculated in 'FV_SurfCalcGradients'
    ! subroutine in fv_reconstruction.f90
    DO q=0,PP_NZ; DO p=0,PP_N
      CALL FV_Limiter(FV_multi_master(:,p,q,SideID), FV_surf_gradU(:,p,q,SideID), gradU(:,p,q))
      UPrim_master(:,p,q,SideID) = UPrim_master(:,p,q,SideID) - gradU(:,p,q) * FV_dx_master(1,p,q,SideID)
    END DO; END DO ! p,q=0,PP_N
  END IF
  IF (SideID.GE.firstInnerSide) THEN
    IF (FV_Elems_slave(SideID).GT.0) THEN ! FV element
      DO q=0,PP_NZ; DO p=0,PP_N
        CALL FV_Limiter(FV_multi_slave(:,p,q,SideID), -FV_surf_gradU(:,p,q,SideID), gradU(:,p,q))
        UPrim_slave(:,p,q,SideID) = UPrim_slave(:,p,q,SideID) - gradU(:,p,q) * FV_dx_slave(1,p,q,SideID)
      END DO; END DO ! p,q=0,PP_N
    END IF
  END IF
END DO

END SUBROUTINE FV_ProlongToDGFace

#endif /* FV_RECONSTRUCT */
END MODULE MOD_FV_ProlongToFace
#endif /* FV_ENABLED */
