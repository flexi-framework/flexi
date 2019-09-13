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
SUBROUTINE Riemann_FD(DFDU,U_L,U_R,UPrim_L,UPrim_R,normal,tangent1,tangent2,surf_loc,jk)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Implicit_Vars           ,ONLY: reps0,sreps0
USE MOD_Riemann                 ,ONLY: Riemann
USE MOD_EOS                     ,ONLY: ConsToPrim
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: DFDU(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: U_L_Tilde    (    1:PP_nVar,0:PP_N,0:PP_NZ)
REAL                            :: UPrim_L_Tilde(1:PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                            :: F            (    1:PP_nVar,0:PP_N,0:PP_NZ)
REAL                            :: F_Tilde_UL   (    1:PP_nVar,0:PP_N,0:PP_NZ)
INTEGER                         :: iVar,jVar
INTEGER                         :: p, q
!-----------------------------------------------------------------------------------------------------------------------------------
!Derivation only in U_L
CALL Riemann(PP_N,F,U_L,U_R,UPrim_L,UPrim_R,normal,tangent1,tangent2,doBC=.FALSE.)

U_L_Tilde = U_L
DO jVar=1,PP_nVar
  U_L_Tilde(jVar,:,:) = U_L_Tilde(jVar,:,:) + reps0
  CALL ConsToPrim(PP_N,UPrim_L_Tilde,U_L_Tilde)
  CALL Riemann(PP_N,F_Tilde_UL,U_L_Tilde,U_R,UPrim_L_Tilde,UPrim_R,normal,tangent1,tangent2,doBC=.FALSE.)
  DO q=0,PP_NZ
    DO p=0,PP_N
      DO iVar=1,PP_nVar
        dFdU(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surf_loc(p,q)*(F_Tilde_UL(iVar,p,q)-F(iVar,p,q))*sreps0
      END DO ! iVar
    END DO !q
  END DO !p
  U_L_Tilde(jVar,:,:) = U_L(jVar,:,:)
END DO !jVar

END SUBROUTINE Riemann_FD

END MODULE MOD_Riemann_Deriv

