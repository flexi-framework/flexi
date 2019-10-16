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
!> Computes the jacobian of the boundary flux for a given face (defined by SideID). Uses the finite difference approach to be able
!> to compute the derivatives for all the different BC types.
!> We compute both the jacobian of the flux w.r.t. the conservative solution, and the jacobians  w.r.t. the gradients in each 
!> direction. We always pre-multiply the jacobians with the surface element.
!===================================================================================================================================
SUBROUTINE GetBoundaryFlux_FD(SideID,t,DfDU,U_master,UPrim_master,    &
#if PARABOLIC
                              gradUx_Face,gradUy_Face,gradUz_Face,  &
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
USE MOD_Implicit_Vars           ,ONLY: reps0_O1,sreps0_O1
USE MOD_GetBoundaryFlux         ,ONLY: GetBoundaryFlux
USE MOD_EOS                     ,ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: SideID                                   !< ID of current side
REAL,INTENT(IN)      :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)      :: xGP_Face(3,0:PP_N,0:PP_NZ)               !< physical coordinates of side GPs
REAL,INTENT(IN)      :: normal  (3,0:PP_N,0:PP_NZ)               !< normal vector
REAL,INTENT(IN)      :: tangent1(3,0:PP_N,0:PP_NZ)               !< first tangential vector
REAL,INTENT(IN)      :: tangent2(3,0:PP_N,0:PP_NZ)               !< second tangential vector
REAL,INTENT(IN)      :: surfElem(0:PP_N,0:PP_NZ)                 !< surface integration element
REAL, INTENT(IN)     :: U_master(PP_nVar,0:PP_N,0:PP_NZ)         !< conservative inner solution 
REAL, INTENT(IN)     :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive inner solution
INTEGER, INTENT(IN)  :: jk(2,0:PP_N,0:PP_NZ)                     !< S2V2 mapping
#if PARABOLIC
REAL,INTENT(IN)      :: gradUx_Face(PP_nVarPrim,0:PP_N,0:PP_NZ)  !< inner gradients in x direction
REAL,INTENT(IN)      :: gradUy_Face(PP_nVarPrim,0:PP_N,0:PP_NZ)  !< inner gradients in y direction
REAL,INTENT(IN)      :: gradUz_Face(PP_nVarPrim,0:PP_N,0:PP_NZ)  !< inner gradients in z direction
#endif /*PARABOLIC*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT),DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2)   :: DfDU         !< jacobian of the boundary flux w.r.t. U inc SurfInt
#if PARABOLIC
REAL,INTENT(OUT),DIMENSION(PP_nVar,PP_nVarPrim,0:PP_N,0:PP_NZ) :: Df_DQxInner  !< jacobian w.r.t. x gradients inc SurfInt
REAL,INTENT(OUT),DIMENSION(PP_nVar,PP_nVarPrim,0:PP_N,0:PP_NZ) :: Df_DQyInner  !< jacobian w.r.t. y gradients inc SurfInt
#if PP_dim==3                                 
REAL,INTENT(OUT),DIMENSION(PP_nVar,PP_nVarPrim,0:PP_N,0:PP_N)  :: Df_DQzInner  !< jacobian w.r.t. z gradients inc SurfInt
#endif
#endif /*PARABOLIC*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iVar,jVar
INTEGER                      :: p,q
#if PARABOLIC
REAL                         :: gradUx_Face_Tilde(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                         :: gradUy_Face_Tilde(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                         :: F_Face_TildeQx   (PP_nVar    ,0:PP_N,0:PP_NZ)
REAL                         :: F_Face_TildeQy   (PP_nVar    ,0:PP_N,0:PP_NZ)
#if PP_dim==3
REAL                         :: gradUz_Face_Tilde(PP_nVarPrim,0:PP_N,0:PP_N)
REAL                         :: F_Face_TildeQz   (PP_nVar    ,0:PP_N,0:PP_NZ)
#endif
#endif /*PARABOLCIC*/
REAL                         :: UPrim_master_Tilde(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                         :: U_master_Tilde    (PP_nVar    ,0:PP_N,0:PP_NZ)
REAL                         :: F_Face_Tilde      (PP_nVar    ,0:PP_N,0:PP_NZ)
REAL                         :: F_Face            (PP_nVar    ,0:PP_N,0:PP_NZ)
REAL                         :: reps0_loc,sreps0_loc
!===================================================================================================================================
! Compute jacobian w.r.t. conservative solution: Store the original state and flux, then pertubate every single DOF, compute
! pertubated flux and use FD approach to approximate the derivative.
CALL GetBoundaryFlux(SideID,t,PP_N,F_Face,UPrim_master,   &
#if PARABOLIC
                     gradUx_Face,gradUy_Face,gradUz_Face, &
#endif /*PARABOLIC*/
                     normal,tangent1,tangent2,xGP_Face)
U_master_Tilde = U_master
DO jVar=1,PP_nVar
  reps0_loc  = reps0_O1*(1.+SQRT(NORM2(U_master_Tilde(jVar,:,:))))
  sreps0_loc = 1./reps0_loc
  U_master_Tilde(jVar,:,:) = U_master_Tilde(jVar,:,:) + reps0_loc
  ! We pertubate the conservative state, but boundary flux takes the primitive state: Call ConsToPrim
  CALL ConsToPrim(PP_N,UPrim_master_Tilde,U_master_Tilde) 
  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_Tilde,UPrim_master_Tilde, &
#if PARABOLIC
                       gradUx_Face,gradUy_Face,gradUz_Face,           &
#endif /*PARABOLIC*/
                       normal,tangent1,tangent2,xGP_Face)
  DO iVar=1,PP_nVar
    DO q=0,PP_NZ
      DO p=0,PP_N
        dFdU(iVar,jVar,jk(1,p,q),jk(2,p,q),1) = surfElem(p,q)*(F_Face_Tilde(iVar,p,q)-F_Face(iVar,p,q))*sreps0_loc
      END DO !p
    END DO !q
  END DO ! iVar
  U_master_Tilde(jVar,:,:) = U_master(jVar,:,:) 
END DO !jVar

#if FV_ENABLED && FV_RECONSTRUCT
dFdU(:,:,:,:,2) = 0. ! only valid for Dirichlet Type BCs
#endif

#if PARABOLIC
! Compute jacobian of the boundary flux w.r.t. the inner gradients for each direction. Again, use the FD approach.
! dF_dQxyzInner
gradUx_Face_Tilde = gradUx_Face
gradUy_Face_Tilde = gradUy_Face
#if PP_dim==3
gradUz_Face_Tilde = gradUz_Face
#endif
DO jVar=1,PP_nVarPrim
  gradUx_Face_Tilde(jVar,:,:) = gradUx_Face_Tilde(jVar,:,:) + reps0_O1
  gradUy_Face_Tilde(jVar,:,:) = gradUy_Face_Tilde(jVar,:,:) + reps0_O1
  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_TildeQx,UPrim_master,   &
                       gradUx_Face_Tilde,gradUy_Face,gradUz_Face, &
                       normal,tangent1,tangent2,xGP_Face)
  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_TildeQy,UPrim_master,   &
                       gradUx_Face,gradUy_Face_Tilde,gradUz_Face, &
                       normal,tangent1,tangent2,xGP_Face)
#if PP_dim==3
  gradUz_Face_Tilde(jVar,:,:) = gradUz_Face_Tilde(jVar,:,:) + reps0_O1
  CALL GetBoundaryFlux(SideID,t,PP_N,F_Face_TildeQz,UPrim_master,   &
                       gradUx_Face,gradUy_Face,gradUz_Face_Tilde, &
                       normal,tangent1,tangent2,xGP_Face)
#endif
  DO iVar=1,PP_nVar
    DO q=0,PP_NZ
      DO p=0,PP_N
        dF_dQxInner(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surfElem(p,q)*(F_Face_TildeQx(iVar,p,q)-F_Face(iVar,p,q))*sreps0_O1
        dF_dQyInner(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surfElem(p,q)*(F_Face_TildeQy(iVar,p,q)-F_Face(iVar,p,q))*sreps0_O1
#if PP_dim==3
        dF_dQzInner(iVar,jVar,jk(1,p,q),jk(2,p,q)) = surfElem(p,q)*(F_Face_TildeQz(iVar,p,q)-F_Face(iVar,p,q))*sreps0_O1
#endif
      END DO !p
    END DO !q
  END DO ! iVar
  gradUx_Face_Tilde(jVar,:,:) = gradUx_Face(jVar,:,:) 
  gradUy_Face_Tilde(jVar,:,:) = gradUy_Face(jVar,:,:) 
#if PP_dim==3
  gradUz_Face_Tilde(jVar,:,:) = gradUz_Face(jVar,:,:) 
#endif
END DO !jVar
#endif /*PARABOLIC*/

END SUBROUTINE GetBoundaryFlux_FD

#if PARABOLIC
!===================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!> Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!>              SUBROUTINE CalcSurfInt
!===================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux_FD(SideID,t,dFdU,UPrim_master,surfElem,xGP_Face,normal,tangent1,tangent2,jk)
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_GetBoundaryFlux         ,ONLY: Lifting_GetBoundaryFlux
USE MOD_Implicit_Vars           ,ONLY: reps0_O1,sreps0_O1
USE MOD_Mesh_Vars               ,ONLY: nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: SideID  
REAL,INTENT(IN)                      :: t
REAL, INTENT(IN)                     :: UPrim_master( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)
REAL,INTENT(IN)                      :: xGP_Face(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: normal  (3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: tangent1(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: tangent2(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)                      :: surfElem(  0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides)
INTEGER,INTENT(IN)                   :: jk(2,0:PP_N,0:PP_NZ)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: dFdU(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_NZ)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q,iVar,jVar
REAL                                 :: UPrim_master_Tilde(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)
REAL                                 :: F_Face        (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)
REAL                                 :: F_Face_Tilde  (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)
!===================================================================================================================================
CALL Lifting_GetBoundaryFlux(SideID,t,UPrim_master(:,:,:,SideID),F_Face(:,:,:,SideID),&
                             normal(:,:,:,0,SideID),tangent1(:,:,:,0,SideID),tangent2(:,:,:,0,SideID),&
                             xGP_Face(:,:,:,0,SideID),surfElem(:,:,0,SideID))
UPrim_master_Tilde = UPrim_master
DO jVar=1,PP_nVarPrim
  UPrim_master_Tilde(jVar,:,:,SideID) = UPrim_master_Tilde(jVar,:,:,SideID) + reps0_O1
  CALL Lifting_GetBoundaryFlux(SideID,t,UPrim_master_Tilde(:,:,:,SideID),F_Face_Tilde(:,:,:,SideID),&
                               normal(:,:,:,0,SideID),tangent1(:,:,:,0,SideID),tangent2(:,:,:,0,SideID),&
                               xGP_Face(:,:,:,0,SideID),surfElem(:,:,0,SideID))
  DO q=0,PP_NZ
    DO p=0,PP_N
      DO iVar=1,PP_nVarPrim
        dFdU(iVar,jVar,jk(1,p,q),jk(2,p,q)) = (F_Face_Tilde(iVar,p,q,SideID)-F_Face(iVar,p,q,SideID))*sreps0_O1
      END DO ! iVar
    END DO !p
  END DO !q
  UPrim_master_Tilde(jVar,:,:,SideID) = UPrim_master(jVar,:,:,SideID) 
END DO !jVar
END SUBROUTINE Lifting_GetBoundaryFlux_FD
#endif /*PARABOLIC*/

END MODULE MOD_GetBoundaryFlux_FD
