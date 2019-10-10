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
#include "eos.h"

MODULE MOD_Riemann_Deriv
!===================================================================================================================================
!   ! Contains the computation of the local jacobian numerical flux
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Riemann_FD
  MODULE PROCEDURE Riemann_FD
END INTERFACE

PUBLIC::Riemann_FD
!===================================================================================================================================

CONTAINS
!===================================================================================================================================
!> Contains the computation of directional derivative with finite difference of the numerical flux f*_adv+f*_diff
!===================================================================================================================================
SUBROUTINE Riemann_FD(DFDU,U_L,U_R,UPrim_L,UPrim_R,normal,tangent1,tangent2,surf_loc,jk,FV_Elems_Sum,FVElem,FVSide)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Implicit_Vars           ,ONLY: reps0_O1
USE MOD_Riemann                 ,ONLY: Riemann
USE MOD_EOS                     ,ONLY: ConsToPrim
#if FV_ENABLED
USE MOD_ChangeBasisByDim        ,ONLY: ChangeBasisSurf
USE MOD_FV_Vars                 ,ONLY: FV_sVdm,FV_Vdm
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: U_L(    PP_nVar    ,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)                :: U_R(    PP_nVar    ,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)                :: UPrim_L(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)                :: UPrim_R(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)                :: normal(  1:3       ,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)                :: tangent1(1:3       ,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)                :: tangent2(1:3       ,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)                :: surf_loc(           0:PP_N,0:PP_NZ)
INTEGER, INTENT(IN)             :: jk(      2         ,0:PP_N,0:PP_NZ)
INTEGER,INTENT(IN)              :: FV_Elems_Sum
INTEGER, INTENT(IN)             :: FVElem
INTEGER, INTENT(IN)             :: FVSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: DFDU(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: U_L_Tilde    (    1:PP_nVar,0:PP_N,0:PP_NZ)
REAL                            :: UPrim_L_Tilde(1:PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                            :: F            (    1:PP_nVar,0:PP_N,0:PP_NZ)
REAL                            :: F_Tilde_UL   (    1:PP_nVar,0:PP_N,0:PP_NZ)
INTEGER                         :: iVar,jVar
INTEGER                         :: p,q
REAL                            :: reps0_loc,sreps0_loc
#if FV_ENABLED && FV_RECONSTRUCT
REAL                            :: U_R_Tilde    (    1:PP_nVar,0:PP_N,0:PP_NZ)
REAL                            :: UPrim_R_Tilde(1:PP_nVarPrim,0:PP_N,0:PP_NZ)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
!Derivation only in U_L
CALL Riemann(PP_N,F,U_L,U_R,UPrim_L,UPrim_R,normal,tangent1,tangent2,doBC=.FALSE.)
#if FV_ENABLED
! convert flux on FV points to DG points at mixed interface if current element is a DG element
IF((FVSide.EQ.1).AND.(FVElem.EQ.0)) CALL ChangeBasisSurf(PP_nVar,PP_N,PP_N,FV_sVdm,F)
#endif

U_L_Tilde = U_L

DO jVar=1,PP_nVar
#if FV_ENABLED
! transform values at side to DG points and modify slightly with reps0
  IF((FVSide.EQ.1).AND.(FVElem.EQ.0))THEN
    CALL ChangeBasisSurf(PP_nVar,PP_N,PP_N,FV_sVdm,U_L_Tilde)
  END IF
#endif
  reps0_loc  = reps0_O1*(1.+SQRT(NORM2(U_L_Tilde(jVar,:,:))))
  sreps0_loc = 1./reps0_loc
  U_L_Tilde(jVar,:,:) = U_L_Tilde(jVar,:,:) + reps0_loc
  CALL ConsToPrim(PP_N,UPrim_L_Tilde,U_L_Tilde)
#if FV_ENABLED
  ! transform modified values back to FV points, where the Riemann problem is evaluated
  IF((FVSide.EQ.1).AND.(FVElem.EQ.0))THEN
    CALL ChangeBasisSurf(PP_nVar    ,PP_N,PP_N,FV_Vdm,U_L_Tilde)
    CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_Vdm,UPrim_L_Tilde)
  END IF
#endif
  CALL Riemann(PP_N,F_Tilde_UL,U_L_Tilde,U_R,UPrim_L_Tilde,UPrim_R,normal,tangent1,tangent2,doBC=.FALSE.)
#if FV_ENABLED
  ! convert flux on FV points to DG points at mixed interface if current element is a DG element 
  ! as the surface integral is evaluated on DG nodes
  IF((FVSide.EQ.1).AND.(FVElem.EQ.0)) CALL ChangeBasisSurf(PP_nVar,PP_N,PP_N,FV_sVdm,F_Tilde_UL)
#endif
  DO q=0,PP_NZ
    DO p=0,PP_N
      DO iVar=1,PP_nVar
        dFdU(iVar,jVar,jk(1,p,q),jk(2,p,q),1) = surf_loc(p,q)*(F_Tilde_UL(iVar,p,q)-F(iVar,p,q))*sreps0_loc
      END DO ! iVar
    END DO !q
  END DO !p
  U_L_Tilde(:,:,:) = U_L(:,:,:)
END DO !jVar

!----------------------------------------------------------------------------------------------------------------------------------
#if FV_ENABLED && FV_RECONSTRUCT
IF((FVElem.EQ.1).AND.(FV_Elems_Sum.EQ.3))THEN ! do only if current and neighbouring element are fv elements (FV-FV interface)
  ! do the same derivative as above only with respect to U_R instead of U_L
  U_R_Tilde = U_R
  DO jVar=1,PP_nVar
    reps0_loc  = reps0_O1*(1.+SQRT(NORM2(U_R_Tilde(jVar,:,:))))
    sreps0_loc = 1./reps0_loc
    U_R_Tilde(jVar,:,:) = U_R_Tilde(jVar,:,:) + reps0_loc
    CALL ConsToPrim(PP_N,UPrim_R_Tilde,U_R_Tilde)
    CALL Riemann(PP_N,F_Tilde_UL,U_L,U_R_Tilde,UPrim_L,UPrim_R_Tilde,normal,tangent1,tangent2,doBC=.FALSE.)
    DO iVar=1,PP_nVar
      DO q=0,PP_NZ
        DO p=0,PP_N
          dFDU(iVar,jVar,jk(1,p,q),jk(2,p,q),2) = surf_loc(p,q)*(F_Tilde_UL(iVar,p,q)-F(iVar,p,q))*sreps0_loc
        END DO !p
      END DO !q
    END DO ! iVar
    U_R_Tilde(:,:,:) = U_R(:,:,:)
  END DO !jVar
ELSE
#endif
  dFDU(:,:,:,:,2) = 0.
#if FV_ENABLED && FV_RECONSTRUCT
END IF
#endif
END SUBROUTINE Riemann_FD

END MODULE MOD_Riemann_Deriv

