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

!===================================================================================================================================
!> Fills the flux jacobian of the boundary fluxes
!===================================================================================================================================
MODULE MOD_GetBoundaryFlux_FD
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE GetBoundaryFlux_FD
  MODULE PROCEDURE GetBoundaryFlux_FD
END INTERFACE


! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC::GetBoundaryFlux_FD
#if PARABOLIC
PUBLIC::Lifting_GetBoundaryFlux_FD
#endif /*PARABOLIC*/
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!> Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!>              SUBROUTINE CalcSurfInt
!===================================================================================================================================
SUBROUTINE GetBoundaryFlux_FD(SideID,t,DfDU,U_master,UPrim_master,    &
#if PARABOLIC
                              gradUx_Face_cons,gradUy_Face_cons,gradUz_Face_cons,  &
                              gradUx_Face_prim,gradUy_Face_prim,gradUz_Face_prim,  &
                              dF_dQxInner, &
                              dF_dQyInner, &
#if PP_dim==3
                              dF_dQzInner, &
#endif
#endif /*PARABOLIC*/
                              surfElem,xGP_Face,normal,tangent1,tangent2,jk)
! MODULES
USE MOD_Globals      
USE MOD_PreProc
USE MOD_Implicit_Vars           ,ONLY: reps0,sreps0
USE MOD_GetBoundaryFlux         ,ONLY: GetBoundaryFlux
USE MOD_EOS                     ,ONLY: ConsToPrim
#if PARABOLIC
USE MOD_EOS                     ,ONLY: ConsToPrimLifting_loc
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: SideID  
REAL,INTENT(IN)              :: t       !< current time (provided by time integration scheme)
REAL,INTENT(IN)              :: xGP_Face(3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: normal  (3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: tangent1(3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: tangent2(3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: surfElem(0:PP_N,0:PP_NZ)
REAL, INTENT(IN)             :: U_master(PP_nVar,0:PP_N,0:PP_NZ)
REAL, INTENT(IN)             :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ)
INTEGER, INTENT(IN)          :: jk(2,0:PP_N,0:PP_NZ)
#if PARABOLIC
REAL,INTENT(IN)              :: gradUx_Face_cons(PP_nVar,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: gradUy_Face_cons(PP_nVar,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: gradUz_Face_cons(PP_nVar,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: gradUx_Face_prim(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: gradUy_Face_prim(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)              :: gradUz_Face_prim(PP_nVarPrim,0:PP_N,0:PP_NZ)
#endif /*PARABOLIC*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: DfDU(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ)
#if PARABOLIC
REAL,INTENT(OUT),DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ) :: Df_DQxInner
REAL,INTENT(OUT),DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ) :: Df_DQyInner
#if PP_dim==3
REAL,INTENT(OUT),DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N) :: Df_DQzInner
#endif
#endif /*PARABOLIC*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iVar,jVar
INTEGER                      :: p,q
#if PARABOLIC
REAL                         :: gradUx_Face_Tilde_cons(PP_nVar,0:PP_N,0:PP_NZ)
REAL                         :: gradUy_Face_Tilde_cons(PP_nVar,0:PP_N,0:PP_NZ)
REAL                         :: gradUx_Face_Tilde_prim(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                         :: gradUy_Face_Tilde_prim(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                         :: gradUz_Face_Tilde_prim(PP_nVarPrim,0:PP_N,0:PP_N)
#if PP_dim==3
REAL                         :: gradUz_Face_Tilde_cons(PP_nVar,0:PP_N,0:PP_N)
#endif /*PP_dim*/
#endif /*PARABOLCIC*/
REAL                         :: UPrim_master_Tilde(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                         :: U_master_Tilde    (PP_nVar,0:PP_N,0:PP_NZ)
REAL                         :: F_Face_Tilde      (PP_nVar,0:PP_N,0:PP_NZ)
REAL                         :: F_Face            (PP_nVar,0:PP_N,0:PP_NZ)

!===================================================================================================================================
CALL GetBoundaryFlux(SideID,t,PP_N,F_Face,UPrim_master, &
#if PARABOLIC
                           gradUx_Face_prim,gradUy_Face_prim,gradUz_Face_prim, &
#endif /*PARABOLIC*/
                           normal,tangent1,tangent2,xGP_Face)
U_master_Tilde = U_master
DO jVar=1,PP_nVar
  U_master_Tilde(jVar,:,:) = U_master_Tilde(jVar,:,:) + reps0
  CALL ConsToPrim(PP_N,UPrim_master_Tilde,U_master_Tilde) 
#if PARABOLIC
  DO q=0,PP_NZ; DO p=0,PP_N
    CALL ConsToPrimLifting_loc(gradUx_Face_cons(:,p,q),gradUy_Face_cons(:,p,q),gradUz_Face_cons(:,p,q),UPrim_master_Tilde(:,p,q),&
        gradUx_Face_Tilde_prim(:,p,q),gradUy_Face_Tilde_prim(:,p,q),gradUz_Face_Tilde_prim(:,p,q))
  END DO; END DO! p,q
#endif /*PARABOLIC*/

  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_Tilde,UPrim_master_Tilde, &
#if PARABOLIC
                       gradUx_Face_Tilde_prim,gradUy_Face_Tilde_prim,gradUz_Face_Tilde_prim, &
#endif /*PARABOLIC*/
                       normal,tangent1,tangent2,xGP_Face)
  DO iVar=1,PP_nVar
    DO q=0,PP_NZ
      DO p=0,PP_N
        dFdU(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surfElem(p,q)*(F_Face_Tilde(iVar,p,q)-F_Face(iVar,p,q))*sreps0
      END DO !p
    END DO !q
  END DO ! iVar
  U_master_Tilde(jVar,:,:) = U_master(jVar,:,:) 
END DO !jVar


#if PARABOLIC
! dF_dQxInner
gradUx_Face_Tilde_cons = gradUx_Face_cons
DO jVar=1,PP_nVar
  gradUx_Face_Tilde_cons(jVar,:,:) = gradUx_Face_Tilde_cons(jVar,:,:) + reps0
  DO q=0,PP_NZ; DO p=0,PP_N
    CALL ConsToPrimLifting_loc(gradUx_Face_Tilde_cons(:,p,q),gradUy_Face_cons(:,p,q),gradUz_Face_cons(:,p,q),UPrim_master(:,p,q),&
        gradUx_Face_Tilde_prim(:,p,q),gradUy_Face_Tilde_prim(:,p,q),gradUz_Face_Tilde_prim(:,p,q))
  END DO; END DO! p,q
  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_Tilde,UPrim_master, &
                       gradUx_Face_Tilde_prim,&
                       gradUy_Face_Tilde_prim,&
                       gradUz_Face_Tilde_prim,&
                       normal,tangent1,tangent2,xGP_Face)
  DO iVar=1,PP_nVar
    DO q=0,PP_NZ
      DO p=0,PP_N
        dF_dQxInner(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surfElem(p,q)*(F_Face_Tilde(iVar,p,q)-F_Face(iVar,p,q))*sreps0
      END DO !p
    END DO !q
  END DO ! iVar
  gradUx_Face_Tilde_cons(jVar,:,:) = gradUx_Face_cons(jVar,:,:) 
END DO !jVar

! dF_dQyInner
gradUy_Face_Tilde_cons = gradUy_Face_cons
DO jVar=1,PP_nVar
  gradUy_Face_Tilde_cons(jVar,:,:) = gradUy_Face_Tilde_cons(jVar,:,:) + reps0
  DO q=0,PP_NZ; DO p=0,PP_N
    CALL ConsToPrimLifting_loc(gradUx_Face_cons(:,p,q),gradUy_Face_Tilde_cons(:,p,q),gradUz_Face_cons(:,p,q),UPrim_master(:,p,q),&
        gradUx_Face_Tilde_prim(:,p,q),gradUy_Face_Tilde_prim(:,p,q),gradUz_Face_Tilde_prim(:,p,q))
  END DO; END DO! p,q
  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_Tilde,UPrim_master, &
                       gradUx_Face_Tilde_prim,&
                       gradUy_Face_Tilde_prim,&
                       gradUz_Face_Tilde_prim,&
                       normal,tangent1,tangent2,xGP_Face)
  DO iVar=1,PP_nVar
    DO q=0,PP_NZ
      DO p=0,PP_N
        dF_dQyInner(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surfElem(p,q)*(F_Face_Tilde(iVar,p,q)-F_Face(iVar,p,q))*sreps0
      END DO !p
    END DO !q
  END DO ! iVar
  gradUy_Face_Tilde_cons(jVar,:,:) = gradUy_Face_cons(jVar,:,:) 
END DO !jVar

#if PP_dim==3
! dF_dQzInner
gradUz_Face_Tilde_cons = gradUz_Face_cons
DO jVar=1,PP_nVar
  gradUz_Face_Tilde_cons(jVar,:,:) = gradUz_Face_Tilde_cons(jVar,:,:) + reps0
  DO q=0,PP_NZ; DO p=0,PP_N
    CALL ConsToPrimLifting_loc(gradUx_Face_cons(:,p,q),gradUy_Face_cons(:,p,q),gradUz_Face_Tilde_cons(:,p,q),UPrim_master(:,p,q),&
        gradUx_Face_Tilde_prim(:,p,q),gradUy_Face_Tilde_prim(:,p,q),gradUz_Face_Tilde_prim(:,p,q))
  END DO; END DO! p,q
  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_Tilde,UPrim_master, &
                       gradUx_Face_Tilde_prim,&
                       gradUy_Face_Tilde_prim,&
                       gradUz_Face_Tilde_prim,&
                       normal,tangent1,tangent2,xGP_Face)
  DO iVar=1,PP_nVar
    DO q=0,PP_N
      DO p=0,PP_N
        dF_dQzInner(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surfElem(p,q)*(F_Face_Tilde(iVar,p,q)-F_Face(iVar,p,q))*sreps0
        !dF_dQzInner(iVar,jVar,jk(1,p,q),jk(2,p,q)) = (F_Face_Tilde(iVar,p,q)-F_Face(iVar,p,q))
      END DO !p
    END DO !q
  END DO ! iVar
  gradUz_Face_Tilde_cons(jVar,:,:) = gradUz_Face_cons(jVar,:,:) 
END DO !jVar
#endif
#endif /*PARABOLIC*/

END SUBROUTINE GetBoundaryFlux_FD

#if PARABOLIC
!===================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!> Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!>              SUBROUTINE CalcSurfInt
!> Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux_FD(SideID,t,dFdU,U_master,surfElem,xGP_Face,normal,tangent1,tangent2,jk)
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_GetBoundaryFlux         ,ONLY: Lifting_GetBoundaryFlux
USE MOD_Implicit_Vars           ,ONLY: reps0,sreps0
USE MOD_Mesh_Vars               ,ONLY: nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: SideID  
REAL,INTENT(IN)                      :: t
REAL, INTENT(IN)                     :: U_master( PP_nVar,0:PP_N,0:PP_NZ,1:nSides)
REAL,INTENT(IN)                      :: xGP_Face(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: normal  (3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: tangent1(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: tangent2(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: surfElem(  0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
INTEGER,INTENT(IN)                   :: jk(2,0:PP_N,0:PP_NZ)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: dFdU(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q,iVar,jVar
REAL                                 :: U_master_Tilde(PP_nVar,0:PP_N,0:PP_NZ,1:nSides)
REAL                                 :: F_Face        (PP_nVar,0:PP_N,0:PP_NZ,1:nSides)
REAL                                 :: F_Face_Tilde  (PP_nVar,0:PP_N,0:PP_NZ,1:nSides)
!===================================================================================================================================
CALL Lifting_GetBoundaryFlux(SideID,t,U_master(:,:,:,SideID),F_Face(:,:,:,SideID),&
                             normal(:,:,:,0,SideID),tangent1(:,:,:,0,SideID),tangent2(:,:,:,0,SideID),&
                             xGP_Face(:,:,:,0,SideID),surfElem(:,:,0,SideID))
U_master_Tilde = U_master
DO jVar=1,PP_nVar
  U_master_Tilde(jVar,:,:,SideID) = U_master_Tilde(jVar,:,:,SideID) + reps0
  CALL Lifting_GetBoundaryFlux(SideID,t,U_master_Tilde(:,:,:,SideID),F_Face_Tilde(:,:,:,SideID),&
                               normal(:,:,:,0,SideID),tangent1(:,:,:,0,SideID),tangent2(:,:,:,0,SideID),&
                               xGP_Face(:,:,:,0,SideID),surfElem(:,:,0,SideID))
  DO q=0,PP_NZ
    DO p=0,PP_N
      DO iVar=1,PP_nVar
        dFdU(iVar,jVar,jk(1,p,q),jk(2,p,q)) = (F_Face_Tilde(iVar,p,q,SideID)-F_Face(iVar,p,q,SideID))*sreps0
      END DO ! iVar
    END DO !p
  END DO !q
  U_master_Tilde(jVar,:,:,SideID) = U_master(jVar,:,:,SideID) 
END DO !jVar
END SUBROUTINE Lifting_GetBoundaryFlux_FD
#endif /*PARABOLIC*/

END MODULE MOD_GetBoundaryFlux_FD
