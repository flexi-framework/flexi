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
#if PARABOLIC
#include "flexi.h"

!==================================================================================================================================
!> \brief Contains the different Surface integral formulations for the BR2 lifting operation
!>
!> Computes the Surface integral for all faces using UPrim and updates gradU
!> Computes only inner surface integrals!
!> Surface integrals are separated for each direction
!==================================================================================================================================
MODULE MOD_Lifting_SurfInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_SurfInt
  MODULE PROCEDURE Lifting_SurfInt
END INTERFACE

PUBLIC::Lifting_SurfInt
!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> \brief Surface integral optimized for performance
!==================================================================================================================================
PPURE SUBROUTINE Lifting_SurfInt(Flux,gradU,gradU_master,gradU_slave,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatMinus
USE MOD_Interpolation_Vars, ONLY: L_Minus
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Mesh_Vars,          ONLY: sJ
USE MOD_Mesh_Vars,          ONLY: S2V
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR, lastMPISide_YOUR
USE MOD_Mesh_Vars,          ONLY: firstInnerSide,    lastInnerSide
USE MOD_Mesh_Vars,          ONLY: firstBCSide,       lastMPISide_MINE
USE MOD_Mesh_Vars,          ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Lifting_Vars,       ONLY: etaBR2, etaBR2_wall
USE MOD_Mesh_Vars,          ONLY: BC,BoundaryType,nBCSides
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems_master,FV_Elems_slave
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)   :: doMPISides  !< = .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE
REAL,INTENT(IN)      :: Flux(1:PP_nVarLifting,0:PP_N,0:PP_NZ,nSides) !< Surface flux contribution
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: gradU(       PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Volume contribution to gradients
REAL,INTENT(INOUT)   :: gradU_master(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides)        !< Gradient on the master sides
REAL,INTENT(INOUT)   :: gradU_slave( PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides)        !< Gradient on the slave sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: F_loc(PP_nVarLifting)
INTEGER              :: ElemID,nbElemID,l,p,q,Flip,SideID,locSideID,nblocSideID,ijk(3)
INTEGER              :: firstSideID,lastSideID
REAL                 :: eta
!==================================================================================================================================

!slave Sides
IF(doMPISides)THEN
  ! only MPI YOUR
  firstSideID = firstMPISide_YOUR
  lastSideID  =  lastMPISide_YOUR
ELSE
  ! only inner sides
  firstSideID = firstInnerSide
  lastSideID  =  lastInnerSide
END IF

eta = etaBR2
DO SideID=firstSideID,lastSideID
  nbElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)

  ! slave sides
  IF(nbElemID.GT.0)THEN
    IF (FV_Elems_slave(SideID).EQ.0) THEN ! DG element
      nblocSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
      flip        = SideToElem(S2E_FLIP,SideID)
      DO q=0,PP_NZ; DO p=0,PP_N
#if (PP_NodeType==1)
        DO l=0,PP_N
#elif (PP_NodeType==2)
        DO l=0,0 ! special mapping, returns 0 or PP_N dep. on locSideID and flip
#endif
        ijk=S2V(1:3,l,p,q,flip,nblocSideID) ! get volume indices belonging to the local side indices p,q,l
        F_loc=sJ(ijk(1),ijk(2),ijk(3),nbElemID,0)*Flux(:,p,q,SideID)*L_hatMinus(l)
        ! contribution to volume gradient
        gradU(:,ijk(1),ijk(2),ijk(3),nbElemID)=gradU(:,ijk(1),ijk(2),ijk(3),nbElemID)+F_loc
        ! local contribution to surface gradient
        gradU_slave(:,p,q,SideID) = gradU_slave(:,p,q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! l,q,p
    END IF
  END IF
END DO ! SideID=1,nSides

!master sides
IF(doMPISides)THEN
  ! only mortar MPI sides
  firstSideID = firstMortarMPISide
   lastSideID =  lastMortarMPISide
ELSE
  ! all sides except YOUR and MortarMPI
  firstSideID = firstBCSide
   lastSideID = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  IF(SideID .LE. nBCSides) THEN
    IF (Boundarytype(BC(SideID),BC_TYPE).EQ.4.OR.Boundarytype(BC(SideID),BC_TYPE).EQ.3) THEN
      eta =etaBR2_wall
    ELSE
      eta = etaBR2
    END IF
  ELSE
    eta = etaBR2
  END IF
  ElemID      = SideToElem(S2E_ELEM_ID,   SideID)
  ! master sides
  IF(ElemID.GT.0)THEN
    IF (FV_Elems_master(SideID).EQ.0) THEN ! DG element
      locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
      flip      = 0
      DO q=0,PP_NZ; DO p=0,PP_N
#if (PP_NodeType==1)
        DO l=0,PP_N
#elif (PP_NodeType==2)
        DO l=0,0 ! special mapping, returns 0 or PP_N dep. on locSideID and flip
#endif
        ijk=S2V(1:3,l,p,q,flip,locSideID) ! get volume indices belonging to the local side indices p,q,l
        F_loc=sJ(ijk(1),ijk(2),ijk(3),ElemID,0)*Flux(:,p,q,SideID)*L_hatMinus(l)
        ! contribution to volume gradient
        gradU(:,ijk(1),ijk(2),ijk(3),ElemID) = gradU(:,ijk(1),ijk(2),ijk(3),ElemID)+F_loc
        ! local contribution to surface gradient
        gradU_master(:,p,q,SideID) = gradU_master(:,p,q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! l,q,p
    END IF
  END IF
END DO ! SideID=1,nSides
END SUBROUTINE Lifting_SurfInt

END MODULE MOD_Lifting_SurfInt
#endif /*PARABOLIC*/
