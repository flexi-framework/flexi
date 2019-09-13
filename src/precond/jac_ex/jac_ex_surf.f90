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
!> Contains the computation of the local jacobian of the surface integral
!===================================================================================================================================
MODULE MOD_JacSurfInt
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE DGJacSurfInt
  MODULE PROCEDURE DGJacSurfInt
END INTERFACE

PUBLIC::DGJacSurfInt
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Contains the computation of the local jacobian of the surface integral
!> computation is done for one element!
!===================================================================================================================================
SUBROUTINE DGJacSurfInt(t,BJ,iElem)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars                   ,ONLY: U_master,U_slave,UPrim_master,UPrim_slave
USE MOD_Mesh_Vars                 ,ONLY: nBCSides,ElemToSide,S2V2
USE MOD_Mesh_Vars                 ,ONLY: NormVec,TangVec1,TangVec2,SurfElem
USE MOD_Mesh_Vars                 ,ONLY: Face_xGP
USE MOD_Jac_Ex_Vars               ,ONLY: LL_minus, LL_plus
USE MOD_Implicit_Vars             ,ONLY: nDOFVarElem
USE MOD_Riemann_Deriv             ,ONLY: Riemann_FD
USE MOD_GetBoundaryFlux_FD        ,ONLY: GetBoundaryFlux_FD
#if PARABOLIC
USE MOD_Precond_Vars              ,ONLY: NoFillIn
USE MOD_DG_Vars                   ,ONLY: L_Hatminus,L_Hatplus
USE MOD_Precond_Vars              ,ONLY: EulerPrecond
USE MOD_Lifting_Vars              ,ONLY: gradUx_master,gradUx_slave
USE MOD_Lifting_Vars              ,ONLY: gradUx_master_cons
USE MOD_Lifting_Vars              ,ONLY: gradUy_master_cons
USE MOD_Lifting_Vars              ,ONLY: gradUz_master_cons
USE MOD_Lifting_Vars              ,ONLY: gradUx_master_prim
USE MOD_Lifting_Vars              ,ONLY: gradUy_master_prim
USE MOD_Lifting_Vars              ,ONLY: gradUz_master_prim
USE MOD_Lifting_Vars              ,ONLY: gradUy_master,gradUy_slave
USE MOD_Lifting_Vars              ,ONLY: gradUz_master,gradUz_slave
USE MOD_DG_Vars                   ,ONLY: nDOFFace 
USE MOD_Jacobian                  ,ONLY: EvalDiffFluxJacobian
USE MOD_GradJacobian              ,ONLY: EvalFluxGradJacobian
!USE MOD_Jac_Ex_Vars               ,ONLY: PrimConsJac
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                    :: t                                         !< current simulation time
INTEGER,INTENT(IN)                 :: iElem                                     !< index of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                 :: BJ(nDOFVarElem,nDOFVarElem)               !< block-Jacobian of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                                 :: i,j,mm,nn,oo,r,s,vn1,vn2
INTEGER                                                 :: iLocSide,SideID,Flip
REAL                                                    :: df_DUinner(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,6)
#if PARABOLIC
INTEGER                                                 :: jk(2),k,l,p,q
REAL                                                    :: Df_dQ_minus(1:PP_nVar,1:PP_nVar)
REAL                                                    :: Df_dQ_plus(1:PP_nVar,1:PP_nVar)
!REAL                                                    :: Df_dQ_minus_Tilde(1:PP_nVar,1:PP_nVar)
!REAL                                                    :: Df_dQ_plus_Tilde(1:PP_nVar,1:PP_nVar)
REAL                                                    ::dQxvol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
REAL                                                    ::dQyvol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
#if PP_dim==3
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,6)     ::Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter, &
                                                          Df_dQzinner,Df_dQzOuter
REAL                                                    ::dQxInner_dUvol(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL                                                    ::dQyInner_dUvol(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL                                                    ::dQxOuter_dUvol(0:PP_N,0:PP_NZ,6,0:PP_N)
REAL                                                    ::dQyOuter_dUvol(0:PP_N,0:PP_NZ,6,0:PP_N)
REAL                                                    ::dQzvol_dU(0:PP_N,0:PP_N,0:PP_N,0:PP_N,3)
REAL                                                    ::dQzInner_dUvol(PP_nVar,PP_nVar,0:PP_N,0:PP_N,6,0:PP_N)
REAL                                                    ::dQzOuter_dUvol(0:PP_N,0:PP_N,6,0:PP_N)
#else
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_NZ,2:5)  ::Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter
REAL                                                    ::dQxInner_dUvol(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
REAL                                                    ::dQyInner_dUvol(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
REAL                                                    ::dQxOuter_dUvol(0:PP_N,0:PP_NZ,2:5,0:PP_N)
REAL                                                    ::dQyOuter_dUvol(0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif /*PP_dim*/
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ)          ::fJacQx,fJacQz,fJacQy,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ)          ::fJac,gJac,hJac
#endif /*PARABOLIC*/
!===================================================================================================================================
vn1 = PP_nVar * (PP_N + 1)
vn2 = vn1 * (PP_N +1)

#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
  Flip=ElemToSide(E2S_FLIP,iLocSide,iElem)

  IF (SideID.GT.nBCSides) THEN !InnerSide
    IF(Flip==0) THEN !Master
      ! d(f*_ad+f*_diff)_jk/dU_master_jk
      CALL Riemann_FD(df_DUinner(:,:,:,:,iLocSide),U_master(:,:,:,SideID),U_slave(:,:,:,SideID), &
                                 UPrim_master(:,:,:,SideID),UPrim_slave(:,:,:,SideID),&
                                 NormVec(:,:,:,0,SideID),tangVec1(:,:,:,0,SideID),tangVec2(:,:,:,0,SideID), &
                                 SurfElem(:,:,0,SideID),S2V2(:,:,:,Flip,iLocSide))
#if PARABOLIC
      IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
        ! analytic viscous flux derivative: BR2 Flux is 1/2*(Fvisc(U^L,Q^L)+Fvisc(U^R,Q^R) )* (+nvec ), (master)
        ! Derivative with respect to the primitive gradients
#if EQNSYSNR==2
        !Derivative of the diffusive Riemann Flux with respect to U_L
        CALL EvalDiffFluxJacobian(nDOFFace,U_master(:,:,:,SideID),UPrim_master(:,:,:,SideID), &
                                  gradUx_master(:,:,:,SideID), &
                                  gradUy_master(:,:,:,SideID), &
                                  gradUz_master(:,:,:,SideID), &
                                  fJac,gJac,hJac)
#endif /*navierstokes*/
        !Derivative of the diffusive Riemann Flux with respect to gradU_L
        CALL EvalFluxGradJacobian(nDOFFace,U_master(:,:,:,SideID),UPrim_Master(:,:,:,SideID), &
                                   fJacQx,fJacQy,fJacQz, &
                                   gJacQx,gJacQy,gJacQz, &
                                   hJacQx,hJacQy,hJacQz)
        DO q=0,PP_NZ
          DO p=0,PP_N
            ! Df_dQxInner = d(f,g,h)_diff/dQ_x * NormVec
            jk(:)=S2V2(:,p,q,Flip,iLocSide)
#if EQNSYSNR==2
            Df_dUinner(:,:,jk(1),jk(2),iLocSide)= Df_dUinner(:,:,jk(1),jk(2),iLocSide) +&
                                                  0.5*(fJac(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                      +gJac(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                      +hJac(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                      )*SurfElem(p,q,0,SideID)
#endif /*navierstokes*/
            Df_dQxInner(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQx(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQx(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQx(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
            Df_dQyInner(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQy(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQy(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQy(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
#if PP_dim==3
            Df_dQzInner(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQz(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQz(:,:,p,q)*NormVec(2,p,q,0,SideID) &
                                                         +hJacQz(:,:,p,q)*NormVec(3,p,q,0,SideID) )*SurfElem(p,q,0,SideID)
#endif
          END DO !p
        END DO !q
        CALL  EvalFluxGradJacobian(nDOFFace,U_slave(:,:,:,SideID),UPrim_slave(:,:,:,SideID), &
                                   fJacQx,fJacQy,fJacQz, &
                                   gJacQx,gJacQy,gJacQz, &
                                   hJacQx,hJacQy,hJacQz)
        DO q=0,PP_NZ
          DO p=0,PP_N
            jk(:)=S2V2(:,p,q,Flip,iLocSide)
            Df_dQxOuter(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQx(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQx(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQx(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
            Df_dQyOuter(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQy(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQy(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQy(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
#if PP_dim==3
            Df_dQzOuter(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQz(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQz(:,:,p,q)*NormVec(2,p,q,0,SideID) &
                                                         +hJacQz(:,:,p,q)*NormVec(3,p,q,0,SideID) )*SurfElem(p,q,0,SideID)
#endif
          END DO !p
        END DO !q
      END IF ! EulerPrecond==False
#endif /*PARABOLIC*/
    ELSE !Slave
      ! d(f*_ad+f*_diff)_jk/dU_slave_jk
      CALL Riemann_FD(df_DUinner(:,:,:,:,iLocSide),U_slave(:,:,:,SideID),U_master(:,:,:,SideID), &
                      UPrim_slave(:,:,:,SideID),UPrim_master(:,:,:,SideID),&
                      -NormVec(:,:,:,0,SideID),-tangVec1(:,:,:,0,SideID),tangVec2(:,:,:,0,SideID), &
                       SurfElem(:,:,0,SideID),S2V2(:,:,:,Flip,iLocSide))
#if PARABOLIC
      IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
#if EQNSYSNR==2
        !Derivative of the diffusive Riemann Flux with respect to U_L
        CALL EvalDiffFluxJacobian(nDOFFace,U_slave(:,:,:,SideID),UPrim_slave(:,:,:,SideID), &
                                  gradUx_slave(:,:,:,SideID), &
                                  gradUy_slave(:,:,:,SideID), &
                                  gradUz_slave(:,:,:,SideID), &
                                  fJac,gJac,hJac)
        !Derivative of the diffusive Riemann Flux with respect to gradU_L
#endif /*navierstokes*/
        ! analytic viscous flux derivative: BR2 Flux is 1/2*(Fvisc(U^L,Q^L)+Fvisc(U^R,Q^R) )* (-nvec ), (slave)
        CALL  EvalFluxGradJacobian(nDOFFace,U_slave(:,:,:,SideID),UPrim_slave(:,:,:,SideID), &
                                   fJacQx,fJacQy,fJacQz, &
                                   gJacQx,gJacQy,gJacQz, &
                                   hJacQx,hJacQy,hJacQz)
        DO q=0,PP_NZ
          DO p=0,PP_N
            jk(:)=S2V2(:,p,q,Flip,iLocSide)
#if EQNSYSNR==2
            Df_dUinner(:,:,jk(1),jk(2),iLocSide)= Df_dUinner(:,:,jk(1),jk(2),iLocSide) &
                                                  -0.5*(fJac(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                        +gJac(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                        +hJac(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
#endif /*navierstokes*/
            Df_dQxInner(:,:,jk(1),jk(2),iLocSide)= -0.5*( fJacQx(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQx(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQx(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
            Df_dQyInner(:,:,jk(1),jk(2),iLocSide)= -0.5*( fJacQy(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQy(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQy(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
#if PP_dim==3
            Df_dQzInner(:,:,jk(1),jk(2),iLocSide)= -0.5*( fJacQz(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQz(:,:,p,q)*NormVec(2,p,q,0,SideID) &
                                                         +hJacQz(:,:,p,q)*NormVec(3,p,q,0,SideID) )*SurfElem(p,q,0,SideID)
#endif
          END DO !p
        END DO !q
        CALL  EvalFluxGradJacobian(nDOFFace,U_master(:,:,:,SideID),UPrim_master(:,:,:,SideID), &
                                   fJacQx,fJacQy,fJacQz, &
                                   gJacQx,gJacQy,gJacQz, &
                                   hJacQx,hJacQy,hJacQz)
        DO q=0,PP_NZ
          DO p=0,PP_N
            jk(:)=S2V2(:,p,q,Flip,iLocSide)
            Df_dQxOuter(:,:,jk(1),jk(2),iLocSide)= -0.5*( fJacQx(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQx(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQx(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
            Df_dQyOuter(:,:,jk(1),jk(2),iLocSide)= -0.5*( fJacQy(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQy(:,:,p,q)*NormVec(2,p,q,0,SideID) &
#if PP_dim==3
                                                         +hJacQy(:,:,p,q)*NormVec(3,p,q,0,SideID) &
#endif
                                                        )*SurfElem(p,q,0,SideID)
#if PP_dim==3
            Df_dQzOuter(:,:,jk(1),jk(2),iLocSide)= -0.5*( fJacQz(:,:,p,q)*NormVec(1,p,q,0,SideID) &
                                                         +gJacQz(:,:,p,q)*NormVec(2,p,q,0,SideID) &
                                                         +hJacQz(:,:,p,q)*NormVec(3,p,q,0,SideID) )*SurfElem(p,q,0,SideID)
#endif
          END DO !p
        END DO !q
      END IF ! EulerPrecond==False
#endif /*PARABOLIC*/

    END IF !Flip
  ELSE !Boundary
      
#if PARABOLIC
    Df_DQxOuter(:,:,:,:,iLocSide)=0.
    Df_DQyOuter(:,:,:,:,iLocSide)=0.
#if PP_dim==3
    Df_DQzOuter(:,:,:,:,iLocSide)=0.
#endif /*PP_dim*/
#endif /*parabolic*/
    !df*_Boundary_jk/dU_jk*surfElem
    !df*_Boundary_jk/dQ_jk*surfElem
    CALL GetBoundaryFlux_FD(SideID,t,df_DUinner(:,:,:,:,iLocSide),U_master(:,:,:,SideID),UPrim_master(:,:,:,SideID),            &
#if PARABOLIC
                            gradUx_master_cons(:,:,:,SideID),gradUy_master_cons(:,:,:,SideID),gradUz_master_cons(:,:,:,SideID), &
                            gradUx_master_prim(:,:,:,SideID),gradUy_master_prim(:,:,:,SideID),gradUz_master_prim(:,:,:,SideID), &
                            Df_DQxInner(:,:,:,:,iLocSide),                                                                      &
                            Df_DQyInner(:,:,:,:,iLocSide),                                                                      &
#if PP_dim==3
                            Df_DQzInner(:,:,:,:,iLocSide),                                                                      &
#endif
#endif /*PARABOLIC*/
                            SurfElem(:,:,0,SideID),Face_xGP(:,:,:,0,SideID),                                                    &
                            NormVec(:,:,:,0,SideID),TangVec1(:,:,:,0,SideID),TangVec2(:,:,:,0,SideID),S2V2(:,:,:,Flip,iLocSide))
  END IF!SideID
END DO !iLocSide


!Assembling of the preconditioner
!BJ=d(f*_adv+f*_diff)_jk/dU_mno
  DO oo = 0,PP_NZ
    DO nn = 0,PP_N
      DO mm = 0,PP_N
        s = vn2 * oo + vn1 * nn + PP_nVar * mm
        DO i = 0,PP_N
          r = PP_nVar*i + vn1*nn + vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                            + LL_Minus(i,mm)*df_DUinner(:,:,nn,oo,XI_MINUS  )  &
                                            + LL_Plus (i,mm)*df_DUinner(:,:,nn,oo,XI_PLUS   )
        END DO ! i
        DO j = 0,PP_N
          r = PP_nVar*mm + vn1*j + vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                            + LL_Minus(j,nn)*df_DUinner(:,:,mm,oo,ETA_MINUS )  &
                                            + LL_Plus (j,nn)*df_DUinner(:,:,mm,oo,ETA_PLUS  )
        END DO ! j
#if PP_dim==3
        DO k = 0,PP_N
          r = PP_nVar*mm + vn1*nn + vn2*k
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                            + LL_Minus(k,oo)*df_DUinner(:,:,mm,nn,ZETA_MINUS)  &
                                            + LL_Plus (k,oo)*df_DUinner(:,:,mm,nn,ZETA_PLUS )
        END DO ! k
#endif
      END DO ! nn
    END DO ! mm 
  END DO ! oo

#if PARABOLIC
IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
! Analytical derivation of DQInner/dUvol
CALL dQInner(1,iElem,dQxInner_dUvol,dQxVol_dU)
CALL dQInner(2,iElem,dQyInner_dUvol,dQyVol_dU)
CALL dQOuter(1,iElem,dQxOuter_dUvol)
CALL dQOuter(2,iElem,dQyOuter_dUvol)
#if PP_dim==3
CALL dQInner(3,iElem,dQzInner_dUvol,dQzVol_dU)
CALL dQOuter(3,iElem,dQzOuter_dUvol)
#endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!version ohne deltas
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! DF_DU = Df_dQInner ( DQVol_DU + DQinner_DU)  + (dF_dQOuter * DQouter_DU)
! First:  Df_dQInner * DQVol_DU 
!       = Df_DQInner(XI_Minus) * (DQVol_dU_1 + DQVol_DU_2 +DQVol_DU3) 

  DO oo = 0,PP_NZ
    DO nn = 0,PP_N
      DO mm = 0,PP_N
        s = vn2 * oo + vn1 * nn + PP_nVar * mm
        DO l=0,PP_N
          df_dQ_minus(:,:)= (Df_dQxInner(:,:,nn,oo,XI_MINUS  ) * dQxVol_dU(l,nn,oo,mm,1)  &
                                 +Df_dQyInner(:,:,nn,oo,XI_MINUS  ) * dQyVol_dU(l,nn,oo,mm,1) & 
#if PP_dim==3
                                 +Df_dQzInner(:,:,nn,oo,XI_MINUS  ) * dQzVol_dU(l,nn,oo,mm,1) &
#endif
                             )
          !df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = (Df_dQxInner(:,:,nn,oo,XI_PLUS   ) * dQxvol_dU(l,nn,oo,mm,1)  &
                                 +Df_dQyInner(:,:,nn,oo,XI_PLUS   ) * dQyvol_dU(l,nn,oo,mm,1) & 
#if PP_dim==3
                                 +Df_dQzInner(:,:,nn,oo,XI_PLUS   ) * dQzvol_dU(l,nn,oo,mm,1) &
#endif
                             )        
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO i = 0,PP_N
            r = PP_nVar * i + vn1 * nn+ vn2 * oo 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(i,l)* df_dQ_minus(:,:) &
                                               + LL_Plus(i,l) * df_dQ_plus(:,:)
          END DO !i
        END DO ! l
      IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
        DO j = 0,PP_N
          df_dQ_minus(:,:)= (Df_dQxInner(:,:,j,oo,XI_MINUS  ) * dQxVol_dU(mm,j,oo,nn,2)  &
                                 +Df_dQyInner(:,:,j,oo,XI_MINUS  ) * dQyVol_dU(mm,j,oo,nn,2) & 
#if PP_dim==3
                                +Df_dQzInner(:,:,j,oo,XI_MINUS  ) * dQzVol_dU(mm,j,oo,nn,2) &
#endif
                            )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = (Df_dQxInner(:,:,j,oo,XI_PLUS   ) * dQxvol_dU(mm,j,oo,nn,2)  &
                                 +Df_dQyInner(:,:,j,oo,XI_PLUS   ) * dQyvol_dU(mm,j,oo,nn,2) &
#if PP_dim==3
                             &   +Df_dQzInner(:,:,j,oo,XI_PLUS   ) * dQzvol_dU(mm,j,oo,nn,2) &
#endif
                            )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO i = 0,PP_N
            r = PP_nVar * i + vn1 * j+ vn2 * oo 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(i,mm)* df_dQ_minus(:,:) & 
                                               + LL_Plus(i,mm) * df_dQ_plus(:,:) 
          END DO ! i
        END DO ! j
#if PP_dim==3
        DO k = 0,PP_NZ
          df_dQ_minus(:,:)= (Df_dQxInner(:,:,nn,k,XI_MINUS  ) * dQxVol_dU(mm,nn,k,oo,3)  &
                                +Df_dQyInner(:,:,nn,k,XI_MINUS  ) * dQyVol_dU(mm,nn,k,oo,3)  &
                                +Df_dQzInner(:,:,nn,k,XI_MINUS  ) * dQzVol_dU(mm,nn,k,oo,3) &
                            )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = (Df_dQxInner(:,:,nn,k,XI_PLUS   ) * dQxvol_dU(mm,nn,k,oo,3)  &
                                +Df_dQyInner(:,:,nn,k,XI_PLUS   ) * dQyvol_dU(mm,nn,k,oo,3) &
                                +Df_dQzInner(:,:,nn,k,XI_PLUS   ) * dQzvol_dU(mm,nn,k,oo,3) &
                            )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO i = 0,PP_N
            r = PP_nVar * i + vn1 * nn+ vn2 * k 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(i,mm)*df_dQ_minus(:,:) &
                                               + LL_Plus(i,mm) * df_dQ_plus(:,:)
          END DO ! i
        END DO ! k
#endif
        DO i = 0,PP_N
          df_dQ_minus(:,:)= ( Df_dQxInner(:,:,i,oo,ETA_MINUS ) * dQxvol_dU(i,nn,oo,mm,1)  &
                                  +Df_dQyInner(:,:,i,oo,ETA_MINUS ) * dQyvol_dU(i,nn,oo,mm,1)&
#if PP_dim==3
                                  +Df_dQzInner(:,:,i,oo,ETA_MINUS ) * dQzvol_dU(i,nn,oo,mm,1)&
#endif
                            )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = ( Df_dQxInner(:,:,i,oo,ETA_PLUS  ) * dQxvol_dU(i,nn,oo,mm,1)  &
                                 +Df_dQyInner(:,:,i,oo,ETA_PLUS  ) * dQyvol_dU(i,nn,oo,mm,1)  &
#if PP_dim==3
                                 +Df_dQzInner(:,:,i,oo,ETA_PLUS  ) * dQzvol_dU(i,nn,oo,mm,1) &
#endif
                            )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO j = 0,PP_N
            r = PP_nVar * i + vn1 * j+ vn2 * oo 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(j,nn) * df_dQ_minus(:,:) & 
                                               + LL_Plus(j,nn)  * df_dQ_plus(:,:)
          END DO ! j
        END DO ! i
      END IF !NoFillIn
        DO l = 0,PP_N
          df_dQ_minus(:,:)= ( Df_dQxInner(:,:,mm,oo,ETA_MINUS ) * dQxvol_dU(mm,l,oo,nn,2)  &
                                  +Df_dQyInner(:,:,mm,oo,ETA_MINUS ) * dQyvol_dU(mm,l,oo,nn,2)  &
#if PP_dim==3
                                 +Df_dQzInner(:,:,mm,oo,ETA_MINUS ) * dQzvol_dU(mm,l,oo,nn,2) &
#endif
                            )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = ( Df_dQxInner(:,:,mm,oo,ETA_PLUS  ) * dQxvol_dU(mm,l,oo,nn,2)  &
                                  +Df_dQyInner(:,:,mm,oo,ETA_PLUS  ) * dQyvol_dU(mm,l,oo,nn,2) & 
#if PP_dim==3
                                 +Df_dQzInner(:,:,mm,oo,ETA_PLUS  ) * dQzvol_dU(mm,l,oo,nn,2) &
#endif
                             )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO j=0,PP_N
            r = PP_nVar * mm + vn1 * j+ vn2 * oo
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(j,l) * df_dQ_minus(:,:) &
                                               + LL_Plus(j,l)  * df_dQ_plus(:,:) 
          END DO !j
        END DO ! l
#if PP_dim==3
      IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
        DO k = 0,PP_NZ
          df_dQ_minus(:,:)=  (Df_dQxInner(:,:,mm,k,ETA_MINUS ) * dQxvol_dU(mm,nn,k,oo,3)  &
                                  +Df_dQyInner(:,:,mm,k,ETA_MINUS ) * dQyvol_dU(mm,nn,k,oo,3) & 
                                 +Df_dQzInner(:,:,mm,k,ETA_MINUS ) * dQzvol_dU(mm,nn,k,oo,3) &
                             )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = ( Df_dQxInner(:,:,mm,k,ETA_PLUS  ) * dQxvol_dU(mm,nn,k,oo,3)  &
                                  +Df_dQyInner(:,:,mm,k,ETA_PLUS  ) * dQyvol_dU(mm,nn,k,oo,3)  &
                                  +Df_dQzInner(:,:,mm,k,ETA_PLUS  ) * dQzvol_dU(mm,nn,k,oo,3) &
                           )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO j = 0,PP_N
            r = PP_nVar * mm + vn1 * j+ vn2 * k 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(j,nn) * df_dQ_minus(:,:) &
                                               + LL_Plus(j,nn)  * df_dQ_plus(:,:) 
          END DO ! j
        END DO ! k
        DO i = 0,PP_N
          df_dQ_minus(:,:)=  (Df_dQxInner(:,:,i,nn,ZETA_MINUS) * dQxvol_dU(i,nn,oo,mm,1)  &
                                  +Df_dQyInner(:,:,i,nn,ZETA_MINUS) * dQyvol_dU(i,nn,oo,mm,1)  &
                                  +Df_dQzInner(:,:,i,nn,ZETA_MINUS) * dQzvol_dU(i,nn,oo,mm,1) &
                           )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          
          df_dQ_plus(:,:) =(  Df_dQxInner(:,:,i,nn,ZETA_PLUS ) * dQxvol_dU(i,nn,oo,mm,1)  &
                                  +Df_dQyInner(:,:,i,nn,ZETA_PLUS ) * dQyvol_dU(i,nn,oo,mm,1)&
                                 +Df_dQzInner(:,:,i,nn,ZETA_PLUS ) * dQzvol_dU(i,nn,oo,mm,1) &
                           )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO k = 0,PP_NZ
            r = PP_nVar * i + vn1 * nn+ vn2 * k 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(k,oo) * df_dQ_minus(:,:) & 
                                               + LL_Plus(k,oo)  * df_dQ_plus(:,:) 
          END DO ! k
        END DO ! i
        DO j = 0,PP_N
          df_dQ_minus(:,:)= ( Df_dQxInner(:,:,mm,j,ZETA_MINUS) * dQxvol_dU(mm,j,oo,nn,2)  &
                                  +Df_dQyInner(:,:,mm,j,ZETA_MINUS) * dQyvol_dU(mm,j,oo,nn,2)  &
                                  +Df_dQzInner(:,:,mm,j,ZETA_MINUS) * dQzvol_dU(mm,j,oo,nn,2) &
                            )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = ( Df_dQxInner(:,:,mm,j,ZETA_PLUS ) * dQxvol_dU(mm,j,oo,nn,2)  &
                                  +Df_dQyInner(:,:,mm,j,ZETA_PLUS ) * dQyvol_dU(mm,j,oo,nn,2)  &
                                  +Df_dQzInner(:,:,mm,j,ZETA_PLUS ) * dQzvol_dU(mm,j,oo,nn,2) &
                           )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO k = 0,PP_NZ
            r = PP_nVar * mm + vn1 * j+ vn2 * k 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(k,oo) * df_dQ_minus(:,:) &
                                               + LL_Plus(k,oo)  * df_dQ_plus(:,:)
          END DO ! k
        END DO ! j
      END IF !NoFillIn
        DO l = 0,PP_N
          df_dQ_minus(:,:)= ( Df_dQxInner(:,:,mm,nn,ZETA_MINUS) * dQxvol_dU(mm,nn,l,oo,3)  &
                                  +Df_dQyInner(:,:,mm,nn,ZETA_MINUS) * dQyvol_dU(mm,nn,l,oo,3)  &
                                  +Df_dQzInner(:,:,mm,nn,ZETA_MINUS) * dQzvol_dU(mm,nn,l,oo,3) &
                            )
          !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus(:,:) = ( Df_dQxInner(:,:,mm,nn,ZETA_PLUS ) * dQxvol_dU(mm,nn,l,oo,3)  &
                                  +Df_dQyInner(:,:,mm,nn,ZETA_PLUS ) * dQyvol_dU(mm,nn,l,oo,3)  &
                                  +Df_dQzInner(:,:,mm,nn,ZETA_PLUS ) * dQzvol_dU(mm,nn,l,oo,3) &
                            )
          !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO k=0,PP_NZ
            r = PP_nVar * mm + vn1 * nn+ vn2 * k 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                               + LL_Minus(k,l) * df_dQ_minus(:,:) & 
                                               + LL_Plus(k,l)  * df_dQ_plus(:,:) 
          END DO !k
        END DO ! l
#endif
      END DO ! nn
    END DO ! mm 
  END DO ! oo

! Second :( Df_dQInner *  DQinner_DU + dF_dQOuter * DQouter_DU) * PrimConsJac 
  DO oo = 0,PP_NZ
    DO nn = 0,PP_N
      DO mm = 0,PP_N
        s = vn2 * oo + vn1 * nn + PP_nVar * mm
        df_dQ_minus(:,:)=  ( MATMUL(Df_dQxinner(:,:,nn,oo,XI_MINUS  ) ,dQxinner_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
                                +MATMUL(Df_dQyinner(:,:,nn,oo,XI_MINUS  ) ,dQyinner_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
#if PP_dim==3
                                +MATMUL(Df_dQzinner(:,:,nn,oo,XI_MINUS  ) ,dQzinner_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
#endif
                                + Df_dQxOuter(:,:,nn,oo,XI_MINUS  ) * dQxOuter_dUvol(nn,oo,XI_MINUS  ,mm)  &
                                + Df_dQyOuter(:,:,nn,oo,XI_MINUS  ) * dQyOuter_dUvol(nn,oo,XI_MINUS  ,mm)  & 
#if PP_dim==3
                               +Df_dQzOuter(:,:,nn,oo,XI_MINUS  ) * dQzOuter_dUvol(nn,oo,XI_MINUS  ,mm) &
#endif
                               )
        !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        df_dQ_plus(:,:) =  ( MATMUL(Df_dQxinner(:,:,nn,oo,XI_PLUS   ) , dQxinner_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
                                +MATMUL(Df_dQyinner(:,:,nn,oo,XI_PLUS   ) , dQyinner_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
#if PP_dim==3
                                +MATMUL(Df_dQzinner(:,:,nn,oo,XI_PLUS   ) , dQzinner_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
#endif
                                +Df_dQxOuter(:,:,nn,oo,XI_PLUS   ) * dQxOuter_dUvol(nn,oo,XI_PLUS   ,mm)  &
                                +Df_dQyOuter(:,:,nn,oo,XI_PLUS   ) * dQyOuter_dUvol(nn,oo,XI_PLUS   ,mm)  &
#if PP_dim==3
                                +Df_dQzOuter(:,:,nn,oo,XI_PLUS   ) * dQzOuter_dUvol(nn,oo,XI_PLUS   ,mm)  &
#endif
                           )
        !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO l = 0,PP_N
          r = PP_nVar * l + vn1 * nn+ vn2 * oo 
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                             + L_HatMinus(l) * df_dQ_minus(:,:) &
                                             + L_HatPlus(l)  * df_dQ_plus(:,:) 
        END DO ! l
        df_dQ_minus(:,:)=  ( MATMUL(Df_dQxinner(:,:,mm,oo,ETA_MINUS ) , dQxinner_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) &
                                +MATMUL(Df_dQyinner(:,:,mm,oo,ETA_MINUS ) , dQyinner_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) & 
#if PP_dim==3
                                +MATMUL(Df_dQzinner(:,:,mm,oo,ETA_MINUS ) , dQzinner_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) &
#endif
                                +Df_dQxOuter(:,:,mm,oo,ETA_MINUS ) * dQxOuter_dUvol(mm,oo,ETA_MINUS ,nn)  &
                                +Df_dQyOuter(:,:,mm,oo,ETA_MINUS ) * dQyOuter_dUvol(mm,oo,ETA_MINUS ,nn)  &
#if PP_dim==3
                                +Df_dQzOuter(:,:,mm,oo,ETA_MINUS ) * dQzOuter_dUvol(mm,oo,ETA_MINUS ,nn) &
#endif
                            )
        !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        df_dQ_plus(:,:) =  ( MATMUL(Df_dQxinner(:,:,mm,oo,ETA_PLUS  ) , dQxinner_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) &
                                +MATMUL(Df_dQyinner(:,:,mm,oo,ETA_PLUS  ) , dQyinner_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) & 
#if PP_dim==3
                                +MATMUL(Df_dQzinner(:,:,mm,oo,ETA_PLUS  ) , dQzinner_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) &
#endif
                                +Df_dQxOuter(:,:,mm,oo,ETA_PLUS  ) * dQxOuter_dUvol(mm,oo,ETA_PLUS  ,nn)  &
                                +Df_dQyOuter(:,:,mm,oo,ETA_PLUS  ) * dQyOuter_dUvol(mm,oo,ETA_PLUS  ,nn)  &
#if PP_dim==3
                                +Df_dQzOuter(:,:,mm,oo,ETA_PLUS  ) * dQzOuter_dUvol(mm,oo,ETA_PLUS  ,nn) &
#endif
                           )
        !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO l = 0,PP_N
          r = PP_nVar * mm + vn1 * l+ vn2 * oo 
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                             + L_HatMinus(l) *df_dQ_minus(:,:) &
                                             + L_HatPlus(l)  *df_dQ_plus(:,:)
        END DO ! l
#if PP_dim==3
        df_dQ_minus(:,:)= ( MATMUL(Df_dQxinner(:,:,mm,nn,ZETA_MINUS) , dQxinner_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                               +MATMUL(Df_dQyinner(:,:,mm,nn,ZETA_MINUS) , dQyinner_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                               +MATMUL(Df_dQzinner(:,:,mm,nn,ZETA_MINUS) , dQzinner_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                               +Df_dQxOuter(:,:,mm,nn,ZETA_MINUS) * dQxOuter_dUvol(mm,nn,ZETA_MINUS,oo)  &
                               +Df_dQyOuter(:,:,mm,nn,ZETA_MINUS) * dQyOuter_dUvol(mm,nn,ZETA_MINUS,oo)  &
                               +Df_dQzOuter(:,:,mm,nn,ZETA_MINUS) * dQzOuter_dUvol(mm,nn,ZETA_MINUS,oo)  &
                          )
        !df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        df_dQ_plus(:,:) = ( MATMUL(Df_dQxinner(:,:,mm,nn,ZETA_PLUS ) , dQxinner_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                               +MATMUL(Df_dQyinner(:,:,mm,nn,ZETA_PLUS ) , dQyinner_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                               +MATMUL(Df_dQzinner(:,:,mm,nn,ZETA_PLUS ) , dQzinner_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                               +Df_dQxOuter(:,:,mm,nn,ZETA_PLUS ) * dQxOuter_dUvol(mm,nn,ZETA_PLUS ,oo)  &
                               +Df_dQyOuter(:,:,mm,nn,ZETA_PLUS ) * dQyOuter_dUvol(mm,nn,ZETA_PLUS ,oo)  &
                               +Df_dQzOuter(:,:,mm,nn,ZETA_PLUS ) * dQzOuter_dUvol(mm,nn,ZETA_PLUS ,oo)  &
                          )
        !df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO l = 0,PP_N
          r = PP_nVar * mm + vn1 * nn+ vn2 * l 
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                             + L_HatMinus(l) * df_dQ_minus(:,:) &
                                             + L_HatPlus(l)  * df_dQ_plus(:,:) 
        END DO ! l
#endif
      END DO ! nn
    END DO ! mm 
  END DO ! oo

END IF !EulerPrecond
#endif /*PARABOLIC*/
END SUBROUTINE DGJacSurfInt

#if PARABOLIC
!===================================================================================================================================
!> Contains the dervative of the BR2 scheme in U_vol: dQ_dUVol
!> ONLY THE DERIVATIVE OF Q_INNER !!!!!
!> computation is done for one element!
!===================================================================================================================================
SUBROUTINE dQInner(dir,iElem,dQ_dUVolInner,dQVol_dU)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                 ,ONLY: ElemToSide,S2V2
USE MOD_Jac_Ex_Vars               ,ONLY: l_mp
USE MOD_Jac_Ex_Vars               ,ONLY: R_Minus,R_Plus
USE MOD_Jac_Ex_Vars               ,ONLY: JacLiftingFlux 
USE MOD_Mesh_Vars                 ,ONLY: Metrics_fTilde,Metrics_gTilde,sJ   ! metrics
#if PP_dim==3
USE MOD_Mesh_Vars                 ,ONLY: Metrics_hTilde   ! metrics
#endif
USE MOD_DG_Vars                   ,ONLY: D
USE MOD_Lifting_Vars              ,ONLY: etaBR2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: dir
INTEGER,INTENT(IN)                                  :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#if PP_dim == 3
REAL,INTENT(OUT)                                    :: dQ_dUVolInner(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
#else
REAL,INTENT(OUT)                                    :: dQ_dUVolInner(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif
REAL,INTENT(OUT)                                    :: dQVol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: iLocSide,p,q,i,j,k,mm,nn,ll
#if PP_dim == 3
INTEGER                                             :: oo
#endif
INTEGER                                             :: SideID,Flip,jk(2)
REAL                                                :: r (0:PP_N,0:PP_NZ)
REAL                                                :: dQ_dUVolInner_loc(0:PP_N,0:PP_NZ,0:PP_N)
!REAL                                                :: dQVol_dU(0:PP_N,0:PP_N,0:PP_N,0:PP_N,0:PP_N,0:PP_N)
!===================================================================================================================================


dQ_dUVolInner=0.

DO ll=0,PP_N
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        dQVol_dU(i,j,k,ll,1) =  sJ(i,j,k,iElem,0)  *  D(i,ll)*Metrics_fTilde(dir,ll,j,k,iElem,0)
        dQVol_dU(i,j,k,ll,2) =  sJ(i,j,k,iElem,0)  *  D(j,ll)*Metrics_gTilde(dir,i,ll,k,iElem,0)
#if PP_dim==3
        dQVol_dU(i,j,k,ll,3) =  sJ(i,j,k,iElem,0)  *  D(k,ll)*Metrics_hTilde(dir,i,j,ll,iElem,0)
#endif
      END DO !i
    END DO !j
  END DO !k
END DO !ll

!Computation of the dQ_Side/dU_Vol (inner part)
#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  Flip  =ElemToSide(E2S_FLIP,iLocSide,iElem)
  IF(Flip.EQ.0)THEN !master
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Minus(dir,p,q,SideID)
      END DO !p
    END DO !q
  ELSE !slave
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Plus(dir,p,q,SideID)
      END DO !p
    END DO !q
  END IF !Flip=0

  dQ_dUVolInner_loc=0.
  SELECT CASE(iLocSide)

  CASE(XI_MINUS,XI_PLUS)
    DO mm=0,PP_N
      DO k=0,PP_NZ
        DO j=0,PP_N
          dQ_dUVolInner_loc(j,k,mm)   = dQ_dUVolInner_loc(j,k,mm) + etaBR2 * r(j,k) * l_mp(mm,iLocSide)
        END DO !j
      END DO !k
    END DO !mm
  CASE(ETA_MINUS,ETA_PLUS)
    DO nn=0,PP_N
      DO k=0,PP_NZ
        DO i=0,PP_N
          dQ_dUVolInner_loc(i,k,nn) = dQ_dUVolInner_loc(i,k,nn) + etaBR2 * r(i,k) * l_mp(nn,iLocSide)
        END DO !i
      END DO !k
    END DO !nn
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
    DO oo=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          dQ_dUVolInner_loc(i,j,oo) = dQ_dUVolInner_loc(i,j,oo) + etaBR2 * r(i,j) * l_mp(oo,iLocSide)
        END DO !i
      END DO !j
    END DO !oo
#endif
  END SELECT
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO mm=0,PP_N
          dQ_dUVolInner(:,:,j,k,iLocSide,mm)=JacLiftingFlux(:,:,j,k,iLocSide)*dq_dUVolinner_loc(j,k,mm)
      END DO
    END DO !p
  END DO !q
END DO !iLocSide

END SUBROUTINE dQInner

!===================================================================================================================================
!> Contains the dervative of the BR2 scheme in U_vol: dQ_dUVol
!> ONLY THE DERIVATIVE OF Q_OUTER !!!!!
!> computation is done for one element!
!===================================================================================================================================
SUBROUTINE dQOuter(dir,iElem,dQ_dUVolOuter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                 ,ONLY: ElemToSide,S2V2,nBCSides
USE MOD_Jac_Ex_Vars               ,ONLY: Surf,l_mp
USE MOD_Jac_Ex_Vars               ,ONLY: R_Minus,R_Plus
USE MOD_Lifting_Vars              ,ONLY: etaBR2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: dir
INTEGER,INTENT(IN)                                  :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#if PP_dim == 3
REAL,INTENT(OUT)                                    :: dQ_dUVolOuter(0:PP_N,0:PP_NZ,6,0:PP_N)
#else
REAL,INTENT(OUT)                                    :: dQ_dUVolOuter(0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: iLocSide,p,q,i,j,k,mm,nn
#if PP_dim == 3
INTEGER                                             :: oo
#endif
INTEGER                                             :: SideID,Flip,jk(2)
REAL                                                :: r (0:PP_N,0:PP_NZ)
!===================================================================================================================================
dQ_dUVolOuter=0.

!Computation of the dQ_Side/dU_Vol (outer part)
#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  IF(SideID.LE.nBCSides) CYCLE  !for boundary conditions, dQ_dUVol=0.
  Flip  =ElemToSide(E2S_FLIP,iLocSide,iElem)
  IF(Flip.EQ.0)THEN
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Plus(dir,p,q,SideID)
      END DO !p
    END DO !q
  ELSE
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Minus(dir,p,q,SideID)
      END DO !p
    END DO !q
  END IF !Flip=0

  SELECT CASE(iLocSide)
  CASE(XI_MINUS,XI_PLUS)
    DO mm=0,PP_N
      DO k=0,PP_NZ
        DO j=0,PP_N
          dQ_dUVolOuter(j,k,iLocSide,mm) = dQ_dUVolOuter(j,k,iLocSide,mm) + 0.5*etaBR2 * r(j,k) * &
                                           l_mp(mm,iLocSide)*Surf(j,k,iLocSide,iElem)
        END DO !j
      END DO !k
    END DO !mm
  CASE(ETA_MINUS,ETA_PLUS)
    DO nn=0,PP_N
      DO k=0,PP_NZ
        DO i=0,PP_N
          dQ_dUVolOuter(i,k,iLocSide,nn) = dQ_dUVolOuter(i,k,iLocSide,nn) + 0.5*etaBR2 * r(i,k) * &
                                           l_mp(nn,iLocSide)*Surf(i,k,iLocSide,iElem)

        END DO !i
      END DO !k
    END DO !nn
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
    DO oo=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          dQ_dUVolOuter(i,j,iLocSide,oo) = dQ_dUVolOuter(i,j,iLocSide,oo) + 0.5*etaBR2 * r(i,j) * &
                                           l_mp(oo,iLocSide)*Surf(i,j,iLocSide,iElem)

        END DO !i
      END DO !j
    END DO !oo
#endif
  END SELECT
END DO !iLocSide

END SUBROUTINE dQOuter
#endif /*PARABOLIC*/

END MODULE MOD_JacSurfInt
