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
SUBROUTINE Lifting_SurfInt(Flux,gradU,gradU_master,gradU_slave,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
#if (PP_NodeType==1) /* Gauss */
USE MOD_Interpolation_Vars, ONLY: L_Plus,L_Minus
#endif
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR, lastMPISide_YOUR
USE MOD_Mesh_Vars,          ONLY: firstInnerSide,    lastInnerSide
USE MOD_Mesh_Vars,          ONLY: firstBCSide,       lastMPISide_MINE
USE MOD_Mesh_Vars,          ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Mesh_Vars,          ONLY: sJ
USE MOD_Mesh_Vars,          ONLY: BC,BoundaryType,nBCSides
USE MOD_Lifting_Vars,       ONLY: etaBR2
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems_master,FV_Elems_slave
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  !< = .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE
REAL,INTENT(IN)    :: Flux(1:PP_nVarPrim,0:PP_N,0:PP_N,nSides) !< Surface flux contribution
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: gradU(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,1:nElems)                   !< Volume contribution to gradients
REAL,INTENT(INOUT)   :: gradU_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)                   !< Gradient on the master sides
REAL,INTENT(INOUT)   :: gradU_slave( PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< Gradient on the slave sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: F_loc(PP_nVarPrim)
INTEGER            :: ElemID,p,q,Flip,SideID,locSideID
INTEGER            :: firstSideID,lastSideID
#if (PP_NodeType==1) /* Gauss */
INTEGER            :: l
#else
REAL               :: L_HatMinus0,L_HatPlusN
REAL               :: eta
#endif
!==================================================================================================================================
#if (PP_NodeType>1)
L_HatMinus0 = L_HatMinus(0)
L_HatPlusN  = L_HatPlus(PP_N)
#endif

IF(doMPISides)THEN
  ! only MPI YOUR
  firstSideID = firstMPISide_YOUR
  lastSideID  =  lastMPISide_YOUR
ELSE
  ! only inner sides
  firstSideID = firstInnerSide
  lastSideID  =  lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  ! neighbor side
  ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip      = SideToElem(S2E_FLIP,SideID)
  ! aditional penalization of wall boundary to prevent high wall velocities
  !HACK!
  IF(SideID .LE. nBCSides) THEN
    IF (Boundarytype(BC(SideID),BC_TYPE).EQ.4.OR.Boundarytype(BC(SideID),BC_TYPE).EQ.3) THEN
      eta =40.
    END IF
  ELSE
    eta = etaBR2
  END IF

  IF (FV_Elems_slave(SideID).EQ.0) THEN ! DG element
  ! update gradients with corresponding SurfInt contribution
#if (PP_NodeType==1)
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    SELECT CASE(flip)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatMinus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,p,q,SideID)=gradU_slave(:,p,q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,PP_N-q,p,SideID)*L_hatMinus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,p,SideID)=gradU_slave(:,PP_N-q,p,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,PP_N-q,SideID)=gradU_slave(:,PP_N-p,PP_N-q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,q,PP_N-p,SideID)*L_hatMinus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,q,PP_N-p,SideID)=gradU_slave(:,q,PP_N-p,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    END SELECT

  CASE(ETA_MINUS)
    SELECT CASE(flip)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,q,p,SideID)*L_hatMinus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,q,p,SideID)=gradU_slave(:,q,p,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatMinus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,q,SideID)=gradU_slave(:,PP_N-p,q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,PP_N-q,PP_N-p,SideID)*L_hatMinus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,PP_N-p,SideID)=gradU_slave(:,PP_N-q,PP_N-p,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,p,PP_N-q,SideID)*L_hatMinus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,p,PP_N-q,SideID)=gradU_slave(:,p,PP_N-q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    END SELECT

  CASE(ZETA_MINUS)
    SELECT CASE(flip)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,p,q,SideID)*L_hatMinus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,p,q,SideID)=gradU_slave(:,p,q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,PP_N-q,p,SideID)*L_hatMinus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,PP_N-q,p,SideID)=gradU_slave(:,PP_N-q,p,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,PP_N-p,PP_N-q,SideID)=gradU_slave(:,PP_N-p,PP_N-q,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,q,PP_N-p,SideID)*L_hatMinus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,q,PP_N-p,SideID)=gradU_slave(:,q,PP_N-p,SideID)+eta*L_Minus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    END SELECT

  CASE(XI_PLUS)
    SELECT CASE(flip)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,q,p,SideID)*L_hatPlus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,q,p,SideID)=gradU_slave(:,q,p,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatPlus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,q,SideID)=gradU_slave(:,PP_N-p,q,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,PP_N-p,SideID)=gradU_slave(:,PP_N-q,PP_N-p,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        F_loc=sJ(l,p,q,ElemID,0)*Flux(:,p,PP_N-q,SideID)*L_hatPlus(l)
        gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
        gradU_slave(:,p,PP_N-q,SideID)=gradU_slave(:,p,PP_N-q,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! l,p,q
    END SELECT

  CASE(ETA_PLUS)
    SELECT CASE(flip)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,q,PP_N-p,SideID)*L_hatPlus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,q,PP_N-p,SideID)=gradU_slave(:,q,PP_N-p,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatPlus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,p,q,SideID)=gradU_slave(:,p,q,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,PP_N-q,p,SideID)*L_hatPlus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,p,SideID)=gradU_slave(:,PP_N-q,p,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,l,q,ElemID,0)*Flux(:,PP_N-p,PP_N-q,SideID)*L_hatPlus(l)
        gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,PP_N-q,SideID)=gradU_slave(:,PP_N-p,PP_N-q,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,l,q
    END SELECT

  CASE(ZETA_PLUS)
    SELECT CASE(flip)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,q,p,SideID)*L_hatPlus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,q,p,SideID)=gradU_slave(:,q,p,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatPlus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,PP_N-p,q,SideID)=gradU_slave(:,PP_N-p,q,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,PP_N-q,PP_N-p,SideID)=gradU_slave(:,PP_N-q,PP_N-p,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,l,ElemID,0)*Flux(:,p,PP_N-q,SideID)*L_hatPlus(l)
        gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
        gradU_slave(:,p,PP_N-q,SideID)=gradU_slave(:,p,PP_N-q,SideID)+eta*L_Plus(l)*F_loc
      END DO; END DO; END DO ! p,q,l
    END SELECT
  END SELECT !locSideID

#else /*PP_NodeType*/

  !update local grid cell
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    SELECT CASE(flip)
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(0,p,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatMinus0
        gradU(:,0,p,q,ElemID)=gradU(:,0,p,q,ElemID)+F_loc
        gradU_slave(:,p,q,SideID)=gradU_slave(:,p,q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(0,p,q,ElemID,0)*Flux(:,PP_N-q,p,SideID)*L_hatMinus0
        gradU(:,0,p,q,ElemID)=gradU(:,0,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,p,SideID)=gradU_slave(:,PP_N-q,p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(0,p,q,ElemID,0)*Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus0
        gradU(:,0,p,q,ElemID)=gradU(:,0,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,PP_N-q,SideID)=gradU_slave(:,PP_N-p,PP_N-q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(0,p,q,ElemID,0)*Flux(:,q,PP_N-p,SideID)*L_hatMinus0
        gradU(:,0,p,q,ElemID)=gradU(:,0,p,q,ElemID)+F_loc
        gradU_slave(:,q,PP_N-p,SideID)=gradU_slave(:,q,PP_N-p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
  CASE(ETA_MINUS)
    SELECT CASE(flip)
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,0,q,ElemID,0)*Flux(:,q,p,SideID)*L_hatMinus0
        gradU(:,p,0,q,ElemID)=gradU(:,p,0,q,ElemID)+F_loc
        gradU_slave(:,q,p,SideID)=gradU_slave(:,q,p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,0,q,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatMinus0
        gradU(:,p,0,q,ElemID)=gradU(:,p,0,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,q,SideID)=gradU_slave(:,PP_N-p,q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,0,q,ElemID,0)*Flux(:,PP_N-q,PP_N-p,SideID)*L_hatMinus0
        gradU(:,p,0,q,ElemID)=gradU(:,p,0,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,PP_N-p,SideID)=gradU_slave(:,PP_N-q,PP_N-p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,0,q,ElemID,0)*Flux(:,p,PP_N-q,SideID)*L_hatMinus0
        gradU(:,p,0,q,ElemID)=gradU(:,p,0,q,ElemID)+F_loc
        gradU_slave(:,p,PP_N-q,SideID)=gradU_slave(:,p,PP_N-q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
  CASE(ZETA_MINUS)
    SELECT CASE(flip)
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,0,ElemID,0)*Flux(:,p,q,SideID)*L_hatMinus0
        gradU(:,p,q,0,ElemID)=gradU(:,p,q,0,ElemID)+F_loc
        gradU_slave(:,p,q,SideID)=gradU_slave(:,p,q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,0,ElemID,0)*Flux(:,PP_N-q,p,SideID)*L_hatMinus0
        gradU(:,p,q,0,ElemID)=gradU(:,p,q,0,ElemID)+F_loc
        gradU_slave(:,PP_N-q,p,SideID)=gradU_slave(:,PP_N-q,p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,0,ElemID,0)*Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus0
        gradU(:,p,q,0,ElemID)=gradU(:,p,q,0,ElemID)+F_loc
        gradU_slave(:,PP_N-p,PP_N-q,SideID)=gradU_slave(:,PP_N-p,PP_N-q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,0,ElemID,0)*Flux(:,q,PP_N-p,SideID)*L_hatMinus0
        gradU(:,p,q,0,ElemID)=gradU(:,p,q,0,ElemID)+F_loc
        gradU_slave(:,q,PP_N-p,SideID)=gradU_slave(:,q,PP_N-p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    END SELECT

  CASE(XI_PLUS)
    SELECT CASE(flip)
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(PP_N,p,q,ElemID,0)*Flux(:,q,p,SideID)*L_hatPlusN
        gradU(:,PP_N,p,q,ElemID)=gradU(:,PP_N,p,q,ElemID)+F_loc
        gradU_slave(:,q,p,SideID)=gradU_slave(:,q,p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(PP_N,p,q,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatPlusN
        gradU(:,PP_N,p,q,ElemID)=gradU(:,PP_N,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,q,SideID)=gradU_slave(:,PP_N-p,q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(PP_N,p,q,ElemID,0)*Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlusN
        gradU(:,PP_N,p,q,ElemID)=gradU(:,PP_N,p,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,PP_N-p,SideID)=gradU_slave(:,PP_N-q,PP_N-p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        F_Loc=sJ(PP_N,p,q,ElemID,0)*Flux(:,p,PP_N-q,SideID)*L_hatPlusN
        gradU(:,PP_N,p,q,ElemID)=gradU(:,PP_N,p,q,ElemID)+F_loc
        gradU_slave(:,p,PP_N-q,SideID)=gradU_slave(:,p,PP_N-q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
  CASE(ETA_PLUS)
    SELECT CASE(flip)
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,PP_N,q,ElemID,0)*Flux(:,q,PP_N-p,SideID)*L_hatPlusN
        gradU(:,p,PP_N,q,ElemID)=gradU(:,p,PP_N,q,ElemID)+F_loc
        gradU_slave(:,q,PP_N-p,SideID)=gradU_slave(:,q,PP_N-p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,PP_N,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatPlusN
        gradU(:,p,PP_N,q,ElemID)=gradU(:,p,PP_N,q,ElemID)+F_loc
        gradU_slave(:,p,q,SideID)=gradU_slave(:,p,q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,PP_N,q,ElemID,0)*Flux(:,PP_N-q,p,SideID)*L_hatPlusN
        gradU(:,p,PP_N,q,ElemID)=gradU(:,p,PP_N,q,ElemID)+F_loc
        gradU_slave(:,PP_N-q,p,SideID)=gradU_slave(:,PP_N-q,p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,PP_N,q,ElemID,0)*Flux(:,PP_N-p,PP_N-q,SideID)*L_hatPlusN
        gradU(:,p,PP_N,q,ElemID)=gradU(:,p,PP_N,q,ElemID)+F_loc
        gradU_slave(:,PP_N-p,PP_N-q,SideID)=gradU_slave(:,PP_N-p,PP_N-q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
  CASE(ZETA_PLUS)
    SELECT CASE(flip)
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,PP_N,ElemID,0)*Flux(:,q,p,SideID)*L_hatPlusN
        gradU(:,p,q,PP_N,ElemID)=gradU(:,p,q,PP_N,ElemID)+F_loc
        gradU_slave(:,q,p,SideID)=gradU_slave(:,q,p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,PP_N,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatPlusN
        gradU(:,p,q,PP_N,ElemID)=gradU(:,p,q,PP_N,ElemID)+F_loc
        gradU_slave(:,PP_N-p,q,SideID)=gradU_slave(:,PP_N-p,q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,PP_N,ElemID,0)*Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlusN
        gradU(:,p,q,PP_N,ElemID)=gradU(:,p,q,PP_N,ElemID)+F_loc
        gradU_slave(:,PP_N-q,PP_N-p,SideID)=gradU_slave(:,PP_N-q,PP_N-p,SideID)+eta*F_loc
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        F_loc=sJ(p,q,PP_N,ElemID,0)*Flux(:,p,PP_N-q,SideID)*L_hatPlusN
        gradU(:,p,q,PP_N,ElemID)=gradU(:,p,q,PP_N,ElemID)+F_loc
        gradU_slave(:,p,PP_N-q,SideID)=gradU_slave(:,p,PP_N-q,SideID)+eta*F_loc
      END DO; END DO ! p,q
    END SELECT
  END SELECT !locSideID
#endif
  END IF
END DO ! SideID=1,nSides


! do master sides
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
  ! master side, flip=0
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)

  IF (FV_Elems_master(SideID).EQ.0) THEN ! DG element
  ! update gradients with corresponding SurfInt contribution
#if (PP_NodeType==1)
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
      F_loc=sJ(l,p,q,ElemID,0)*Flux(:,q,p,SideID)*L_hatMinus(l)
      gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
      gradU_master(:,q,p,SideID)=gradU_master(:,q,p,SideID)+eta*L_Minus(l)*F_loc
    END DO; END DO; END DO ! l,p,q
  CASE(ETA_MINUS)
    DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,l,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatMinus(l)
      gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
      gradU_master(:,p,q,SideID)=gradU_master(:,p,q,SideID)+eta*L_Minus(l)*F_loc
    END DO; END DO; END DO ! p,l,q
  CASE(ZETA_MINUS)
    DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,q,l,ElemID,0)*Flux(:,q,p,SideID)*L_hatMinus(l)
      gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
      gradU_master(:,q,p,SideID)=gradU_master(:,q,p,SideID)+eta*L_Minus(l)*F_loc
    END DO; END DO; END DO ! p,q,l
  CASE(XI_PLUS)
    DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
      F_loc=sJ(l,p,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatPlus(l)
      gradU(:,l,p,q,ElemID)=gradU(:,l,p,q,ElemID)+F_loc
      gradU_master(:,p,q,SideID)=gradU_master(:,p,q,SideID)+eta*L_Plus(l)*F_loc
    END DO; END DO; END DO ! l,p,q
  CASE(ETA_PLUS)
    DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,l,q,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatPlus(l)
      gradU(:,p,l,q,ElemID)=gradU(:,p,l,q,ElemID)+F_loc
      gradU_master(:,PP_N-p,q,SideID)=gradU_master(:,PP_N-p,q,SideID)+eta*L_Plus(l)*F_loc
    END DO; END DO; END DO ! p,l,q
  CASE(ZETA_PLUS)
    DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,q,l,ElemID,0)*Flux(:,p,q,SideID)*L_hatPlus(l)
      gradU(:,p,q,l,ElemID)=gradU(:,p,q,l,ElemID)+F_loc
      gradU_master(:,p,q,SideID)=gradU_master(:,p,q,SideID)+eta*L_Plus(l)*F_loc
    END DO; END DO; END DO ! p,q,l
  END SELECT !locSideID

#else /*PP_NodeType*/

  !update local grid cell
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(0,p,q,ElemID,0)*Flux(:,q,p,SideID)*L_hatMinus0
      gradU(:,0,p,q,ElemID)=gradU(:,0,p,q,ElemID)+F_loc
      gradU_master(:,q,p,SideID)=gradU_master(:,q,p,SideID)+eta*F_loc
    END DO; END DO ! p,q
  ! switch to right hand system for ETA_PLUS direction
  CASE(ETA_MINUS)
    DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,0,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatMinus0
      gradU(:,p,0,q,ElemID)=gradU(:,p,0,q,ElemID)+F_loc
      gradU_master(:,p,q,SideID)=gradU_master(:,p,q,SideID)+eta*F_loc
    END DO; END DO ! p,q
  ! switch to right hand system for ZETA_MINUS direction
  CASE(ZETA_MINUS)
    DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,q,0,ElemID,0)*Flux(:,q,p,SideID)*L_hatMinus0
      gradU(:,p,q,0,ElemID)=gradU(:,p,q,0,ElemID)+F_loc
      gradU_master(:,q,p,SideID)=gradU_master(:,q,p,SideID)+eta*F_loc
    END DO; END DO ! p,q
  CASE(XI_PLUS)
    DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(PP_N,p,q,ElemID,0)*Flux(:,p,q,SideID)*L_hatPlusN
      gradU(:,PP_N,p,q,ElemID)=gradU(:,PP_N,p,q,ElemID)+F_loc
      gradU_master(:,p,q,SideID)=gradU_master(:,p,q,SideID)+eta*F_loc
    END DO; END DO ! p,q
  ! switch to right hand system for ETA_PLUS direction
  CASE(ETA_PLUS)
    DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,PP_N,q,ElemID,0)*Flux(:,PP_N-p,q,SideID)*L_hatPlusN
      gradU(:,p,PP_N,q,ElemID)=gradU(:,p,PP_N,q,ElemID)+F_loc
      gradU_master(:,PP_N-p,q,SideID)=gradU_master(:,PP_N-p,q,SideID)+eta*F_loc
    END DO; END DO ! p,q
  ! switch to right hand system for ZETA_MINUS direction
  CASE(ZETA_PLUS)
    DO q=0,PP_N; DO p=0,PP_N
      F_loc=sJ(p,q,PP_N,ElemID,0)*Flux(:,p,q,SideID)*L_hatPlusN
      gradU(:,p,q,PP_N,ElemID)=gradU(:,p,q,PP_N,ElemID)+F_loc
      gradU_master(:,p,q,SideID)=gradU_master(:,p,q,SideID)+eta*F_loc
    END DO; END DO ! p,q
  END SELECT !locSideID
#endif
  END IF
END DO ! SideID=1,nSides
END SUBROUTINE Lifting_SurfInt


END MODULE MOD_Lifting_SurfInt
#endif /*PARABOLIC*/
