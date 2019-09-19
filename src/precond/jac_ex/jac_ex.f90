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

!===================================================================================================================================
!> This module contains the routines required for the exact block-Jacobi preconditioner. (Additionally, the exact block-Jacobi
!> preconditioner requires the module Jac_ex_surf)
!===================================================================================================================================
MODULE MOD_Jac_Ex
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitJac_Ex
  MODULE PROCEDURE InitJac_Ex
END INTERFACE

INTERFACE Jac_Ex
  MODULE PROCEDURE Jac_Ex
END INTERFACE

INTERFACE FinalizeJac_Ex
  MODULE PROCEDURE FinalizeJac_Ex
END INTERFACE


PUBLIC::InitJac_Ex,Jac_Ex,FinalizeJac_Ex
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Inititializes the variables required for the exact block-Jacobi preconditioner 
!===================================================================================================================================
SUBROUTINE InitJac_Ex()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Jac_Ex_Vars
USE MOD_ReadInTools        ,ONLY: GETLOGICAL
USE MOD_Interpolation_Vars ,ONLY: L_minus,L_plus
USE MOD_DG_Vars            ,ONLY: L_Hatminus,L_Hatplus
USE MOD_Mesh_Vars          ,ONLY: nElems
#if PARABOLIC
USE MOD_Precond_Vars       ,ONLY: EulerPrecond
#endif
#if FV_ENABLED && FV_RECONSTRUCT
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                      :: i,j,iLocSide
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT BLOCK JACOBIAN...'
ALLOCATE(LL_minus(0:PP_N,0:PP_N), LL_plus(0:PP_N,0:PP_N))
ALLOCATE(l_mp(0:PP_N,6))
ALLOCATE(LL_mp(0:PP_N,6))
!ALLOCATE(PrimConsJac(1:PP_nVarPrim,1:PP_nVar,0:PP_N,0:PP_N,0:PP_N))
#if FV_ENABLED && FV_RECONSTRUCT
ALLOCATE(FV_sdx_XI_extended  (0:PP_N,0:PP_NZ,0:PP_N+1,nElems))          ! 1. / FV_dx_XI   Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_sdx_ETA_extended (0:PP_N,0:PP_NZ,0:PP_N+1,nElems))          ! 1. / FV_dx_ETA  Attention: storage order is (i,k,j,iElem)
#if PP_dim == 3
ALLOCATE(FV_sdx_ZETA_extended(0:PP_N,0:PP_N ,0:PP_N+1,nElems))          ! 1. / FV_dx_ZETA Attention: storage order is (i,j,k,iElem)
ALLOCATE(UPrim_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,1:nElems))! extended primitive solution vector
#else
ALLOCATE(UPrim_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1, 0:PP_NZ ,1:nElems))! extended primitive solution vector
#endif
#endif

DO j=0,PP_N
  DO i=0,PP_N
    LL_minus(i,j) = L_Hatminus(i)*L_minus(j)
    LL_plus(i,j)  = L_Hatplus(i) *L_plus(j)
  END DO
END DO 

#if PP_dim==3
DO iLocSide=1,6
#else
DO iLocSide=2,5
#endif 
  SELECT CASE(iLocSide)
  CASE(XI_MINUS,ETA_MINUS,ZETA_MINUS)
    DO i=0,PP_N
      l_mp (i,iLocSide) = l_minus(i)
      LL_mp(i,iLocSide) = LL_minus(i,i)
    END DO !i
  CASE(XI_PLUS,ETA_PLUS,ZETA_PLUS)
    DO i=0,PP_N
      l_mp (i,iLocSide) = l_plus(i)
      LL_mp(i,iLocSide) = LL_plus(i,i)
    END DO !i
  END SELECT
END DO !iLocSide


#if PP_dim==3
ALLOCATE(nVec (1:3,0:PP_N,0:PP_N,1:6,nElems))
ALLOCATE(Surf (0:PP_N,0:PP_N,1:6,nElems))
#else
ALLOCATE(nVec (1:3,0:PP_N,0:PP_NZ,2:5,nElems))
ALLOCATE(Surf (0:PP_N,0:PP_NZ,2:5,nElems))
#endif
nVec=0.
Surf=0.

#if PARABOLIC
CALL BuildnVecTangSurf()
IF(EulerPrecond.EQV..FALSE.) THEN
  ALLOCATE(JacLiftingFlux(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,6))
  CALL Build_BR2_SurfTerms()
END IF
#endif /*PARABOLIC*/

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT BLOCK JACOBIAN DONE!'
END SUBROUTINE InitJac_Ex

!SUBROUTINE ConsToPrimJac(U,UPrim,PrimConsJac)
!!===================================================================================================================================
!! Calculates the derivative of the primitive Variables UPrim with respect to the conservative Variables U 
!! in the volume coordinates (i,j,k)
!! dUPrim(i,j,k)/DU(i,j,k)
!!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_EOS_Vars      ,ONLY:KappaM1,R
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT / OUTPUT VARIABLES
!REAL,DIMENSION(PP_nVar),INTENT(IN)              :: U
!REAL,DIMENSION(PP_nVarPrim),INTENT(IN)          :: UPrim
!REAL,DIMENSION(PP_nVarPrim,PP_nVar),INTENT(OUT) :: PrimConsJac  ! Derivative of The Primitive Variables in (i,j,k)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES 
!!===================================================================================================================================
!REAL             :: sRho,vv,p_rho,KappaM1sRhoR 
!!==================================================================================================================================
!#if EQNSYSNR==1
!PrimConsJac=1.
!#endif
!#if EQNSYSNR==2
!sRho=1./UPrim(1)
!vv=SUM(UPrim(2:4)*UPrim(2:4))
!!derivative of p=UPrim(5)
!p_rho=KappaM1*0.5*vv
!KappaM1sRhoR=KappaM1*sRho/R
!!rho
!!prim(1)=cons(1)
!PrimConsJac(1,1:5) = (/               1.,   0., 0., 0., 0. /)
!!velocity
!!prim(2:4)=cons(2:4)*sRho
!PrimConsJac(2,1:5) = (/ -UPrim(2)*sRho, sRho, 0., 0., 0. /)
!PrimConsJac(3,1:5) = (/ -UPrim(3)*sRho, 0., sRho, 0., 0. /)
!PrimConsJac(4,1:5) = (/ -UPrim(4)*sRho, 0., 0., sRho, 0. /)
!!pressure
!!prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
!PrimConsJac(5,1:5) = (/ p_rho, -KappaM1*UPrim(2), -KappaM1*UPrim(3), -KappaM1*UPrim(4), KappaM1 /)
!!temperature
!!prim(6) = prim(5)*sRho / R
!!This row is not needed, since the fluxes are independent on $T$
!PrimConsJac(6,1:5) =(/ KappaM1sRhoR*(vv-U(5)*sRho) , -KappaM1sRhoR*UPrim(2), -KappaM1sRhoR*UPrim(3), &
                                                                                  !-KappaM1sRhoR*UPrim(4), KappaM1sRhoR /)
!#endif

!END SUBROUTINE ConsToPrimJac 

SUBROUTINE Jac_Ex(t,iElem,BJ,doVol,doSurf)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars             ,ONLY:nDOFVarElem
USE MOD_JacSurfInt                ,ONLY:DGJacSurfInt
#if PARABOLIC
USE MOD_Precond_Vars              ,ONLY:EulerPrecond
#endif
#if FV_ENABLED
USE MOD_FV_Vars                   ,ONLY:FV_Elems
#endif
!USE MOD_Jac_Ex_Vars               ,ONLY:PrimConsJac
!USE MOD_DG_Vars                   ,ONLY:U,UPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t                                         !< current simulation time
INTEGER,INTENT(IN) :: iElem                                     !< index of current element
LOGICAL,INTENT(IN) :: dovol, dosurf                             !< logicals indicating to do derivative of volume/surface-integral 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: BJ(1:nDOFVarElem,1:nDOFVarElem)             !< block-Jacobian of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: FVEM
!INTEGER            :: i,j,k
!===================================================================================================================================
#if FV_ENABLED
FVEM = FV_Elems(iElem)
#else
FVEM = 0
#endif
! Nullify BJ
BJ = 0.
#if PARABOLIC
IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
  CALL FillJacLiftingFlux(iElem)
END IF
#endif /*PARABOLIC*/
!DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !CALL ConsToPrimJac(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem),PrimConsJac(:,:,i,j,k))
!END DO; END DO; END DO
IF(doVol)THEN
  IF(FVEM.EQ.0)THEN
    CALL DGVolIntJac(BJ,iElem) !without sJ!      !d(F^a+F^v)/dU partial
#if FV_ENABLED
  ELSE
    CALL FVVolIntJac(BJ,iElem)
#endif
  END IF
#if PARABOLIC
  IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
    CALL DGVolIntGradJac(BJ,iElem)               !d(F^v)/dQ
  END IF
#endif /*PARABOLIC*/
END IF!doVol
IF(doSurf) THEN
  CALL DGJacSurfInt(t,BJ,iElem) 
END IF
CALL Apply_sJ(BJ,iElem)

END SUBROUTINE Jac_Ex

!===================================================================================================================================
!> Volume integral Jacobian of the convective flux for DG elements. It uses the analytical flux Jacobian, applies the metric termes
!> and multiplies D_hat
!===================================================================================================================================
SUBROUTINE  DGVolIntJac(BJ,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars       ,ONLY: U,UPrim
USE MOD_Jacobian      ,ONLY: EvalAdvFluxJacobian
#if PARABOLIC
USE MOD_DG_Vars       ,ONLY: nDOFElem
USE MOD_Precond_Vars  ,ONLy: EulerPrecond
USE MOD_Lifting_Vars  ,ONLY: gradUx,gradUy,gradUz
#if EQNSYSNR==2                                         
USE MOD_Jacobian      ,ONLY: EvalDiffFluxJacobian 
#endif /*EQNSYSNR*/
#endif /*PARABOLIC*/
USE MOD_Mesh_Vars     ,ONLY: Metrics_fTilde,Metrics_gTilde
#if PP_dim==3
USE MOD_Mesh_Vars     ,ONLY: Metrics_hTilde
#endif
USE MOD_Implicit_Vars ,ONLY: nDOFVarElem
#ifdef SPLIT_DG
USE MOD_DG_Vars       ,ONLY: D
USE MOD_Jac_Split     ,ONLY: Jac_Split
USE MOD_Interpolation_Vars,ONLY: wGP
#else
USE MOD_DG_Vars       ,ONLY: D_hat,U,UPrim
#endif
!USE MOD_Jac_Ex_Vars   ,ONLY: PrimConsJac
!USE MOD_Implicit_Vars ,ONLY: reps0,sreps0
USE MOD_EOS           ,ONLY: ConsToPrim
USE MOD_Implicit_Vars, ONLY: Xk,sreps0,reps0
USE MOD_VolInt       , ONLY: VolInt
USE MOD_Mesh_Vars    , ONLY: nElems
USE MOD_DG_Vars,       ONLY: Ut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: BJ(1:nDOFVarElem,1:nDOFVarElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!REAL,DIMENSION(PP_nVar) :: f,g,h,u1,u1_tilde,gx,gy,gz
!REAL,DIMENSION(PP_nVarPrim) :: up,up_tilde,gpx,gpy,gpz
!REAL,DIMENSION(PP_nVar) :: f_Tilde,g_Tilde,h_Tilde
!REAL,DIMENSION(PP_nVar,PP_nVar)   :: dfdu,dgdu,dhdu
INTEGER                                                 :: iVar,jVar,i,j,k,r,tt,r11,r22
INTEGER                                                 :: mm,nn,oo
INTEGER                                                 :: s,r1,r2,ll,vn1,vn2
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJacTilde,gJacTilde
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: fJac,gJac,hJac
#if PP_dim==3
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: hJacTilde
INTEGER                                                 :: r3
#endif
#if PARABOLIC
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: fJac_loc,gJac_loc,hJac_loc
#endif
#ifdef SPLIT_DG
REAL                                                    :: DVolSurf_loc(0:PP_N,0:PP_N)
REAL                                                    :: DVolSurf_diag(0:PP_N,0:PP_N)
#endif
REAL                                                    :: Ut_tilde(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL    :: dRdU(1:nDOFVarElem,1:nDOFVarElem)
!===================================================================================================================================
vn1=PP_nVar*(PP_N+1)
vn2=vn1*(PP_N+1)
CALL EvalAdvFluxJacobian(U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem),fJac,gJac,hJac)
!No CALL of EvalDiffFlux Jacobian for EQNSYSNR==1, since the diffusive flux is not depending on U

#if EQNSYSNR==2                                         
#if PARABOLIC
IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
  CALL EvalDiffFluxJacobian(nDOFElem,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem) &
                            ,gradUx(:,:,:,:,iElem) &
                            ,gradUy(:,:,:,:,iElem) &
                            ,gradUz(:,:,:,:,iElem) &
                            ,fJac_loc,gJac_loc,hJac_loc)
  fJac=fJac+fJac_loc                                         
  gJac=gJac+gJac_loc                                         
#if PP_dim==3
  hJac=hJac+hJac_loc                                         
#endif/*PP_dim*/
END IF !EulerPrecond
#endif /*PARABOLIC*/          
#endif /*EQNSYSNR*/

#ifdef SPLIT_DG
DVolSurf_loc = D
!DVolSurf_loc(0,0) = DVolSurf_loc(0,0) + 1.0/(wGP(0))
!DVolSurf_loc(PP_N,PP_N) = DVolSurf_loc(PP_N,PP_N) - 1.0/(wGP(PP_N))
DVolSurf_diag = D
!DVolSurf_diag(0,0) = DVolSurf_diag(0,0) + 1.0/(wGP(0))
!DVolSurf_diag(PP_N,PP_N) = DVolSurf_diag(PP_N,PP_N) - 1.0/(wGP(PP_N))
#endif /*SPLIT_DG*/

s=0
DO oo=0,PP_NZ
  DO nn=0,PP_N
    DO mm=0,PP_N
#ifdef SPLIT_DG
      ! strong form
      r1=           vn1*nn+vn2*oo
      r2=mm*PP_nVar       +vn2*oo
#if PP_dim==3
      r3=mm*PP_nVar+vn1*nn
#endif
      DO ll=0,PP_N
        !CALL Jac_Split(U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem),U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem), &
                       !Metrics_fTilde(:,ll,nn,oo,iElem,0),Metrics_fTilde(:,mm,nn,oo,iElem,0),fJacTilde(:,:))
        !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,mm)*fJacTilde(:,:) 

        !!CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem), &
                       !!Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,ll,nn,oo,iElem,0),fJacTilde(:,:))
        !IF(ll.EQ.mm)THEN
          !BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,ll)*fJacTilde(:,:) 
        !END IF


        !CALL Jac_Split(U(:,mm,ll,oo,iElem),UPrim(:,mm,ll,oo,iElem),U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem), &
                       !Metrics_gTilde(:,mm,ll,oo,iElem,0),Metrics_gTilde(:,mm,nn,oo,iElem,0),gJacTilde(:,:))
        !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,nn)*gJacTilde(:,:) 

        !!CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,ll,oo,iElem),UPrim(:,mm,ll,oo,iElem), &
                       !!Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,ll,oo,iElem,0),gJacTilde(:,:))
        !IF(ll.EQ.nn)THEN
          !BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,ll)*gJacTilde(:,:) 
        !END IF





        !IF(r1.EQ.s)THEN
          !fJacTilde(:,:) = ( fJac(:,:,mm,nn,oo)*Metrics_fTilde(1,mm,nn,oo,iElem,0)  &
                           !+ gJac(:,:,mm,nn,oo)*Metrics_fTilde(2,mm,nn,oo,iElem,0)  &
!#if PP_dim==3
                           !+ hJac(:,:,mm,nn,oo)*Metrics_fTilde(3,mm,nn,oo,iElem,0)  &
!#endif
                           !)

          !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_diag(ll,mm)*fJacTilde(:,:) 

          !CALL Jac_Split(U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem),U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem), &
                         !Metrics_fTilde(:,ll,nn,oo,iElem,0),Metrics_fTilde(:,mm,nn,oo,iElem,0),fJacTilde(:,:))
          !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_diag(ll,mm)*fJacTilde(:,:) 
        !ELSE
          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem), &
                         !Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,ll,nn,oo,iElem,0),fJacTilde(:,:))
          !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,mm)*fJacTilde(:,:) 
        !END IF
        !IF(r2.EQ.s)THEN
          !gJacTilde(:,:) = (fJac(:,:,mm,nn,oo)*Metrics_gTilde(1,mm,nn,oo,iElem,0)  &
                           !+gJac(:,:,mm,nn,oo)*Metrics_gTilde(2,mm,nn,oo,iElem,0)  &
!#if PP_dim==3
                           !+hJac(:,:,mm,nn,oo)*Metrics_gTilde(3,mm,nn,oo,iElem,0)  &
!#endif
                           !)
          !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_diag(ll,nn)*gJacTilde(:,:) 
        !ELSE
          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,ll,oo,iElem),UPrim(:,mm,ll,oo,iElem), &
                         !Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,ll,oo,iElem,0),gJacTilde(:,:))
          !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,nn)*gJacTilde(:,:) 
        !END IF
        !IF((mm.EQ.ll).AND.(nn.EQ.ll).AND.(oo.EQ.ll))THEN
        !ELSE
          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem), &
                         !Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,ll,nn,oo,iElem,0),fJacTilde(:,:))
          !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,mm)*fJacTilde(:,:) 
          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,ll,oo,iElem),UPrim(:,mm,ll,oo,iElem), &
                         !Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,ll,oo,iElem,0),gJacTilde(:,:))
          !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,nn)*gJacTilde(:,:) 
        !END IF

        !IF(mm.EQ.ll)THEN ! main diagonal
          !fJacTilde(:,:) = ( fJac(:,:,mm,nn,oo)*Metrics_fTilde(1,mm,nn,oo,iElem,0)  &
                           !+ gJac(:,:,mm,nn,oo)*Metrics_fTilde(2,mm,nn,oo,iElem,0)  &
!#if PP_dim==3
                           !+ hJac(:,:,mm,nn,oo)*Metrics_fTilde(3,mm,nn,oo,iElem,0)  &
!#endif
                           !)

          !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_diag(ll,mm)*fJacTilde(:,:) 
        !ELSE
        IF(ll.NE.mm)THEN
          CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem), &
                         Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,ll,nn,oo,iElem,0),fJacTilde(:,:))
          BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,mm)*fJacTilde(:,:) 
        ELSE
          fJacTilde(:,:) = ( fJac(:,:,mm,nn,oo)*Metrics_fTilde(1,mm,nn,oo,iElem,0)  &
                           + gJac(:,:,mm,nn,oo)*Metrics_fTilde(2,mm,nn,oo,iElem,0)  &
#if PP_dim==3
                           + hJac(:,:,mm,nn,oo)*Metrics_fTilde(3,mm,nn,oo,iElem,0)  &
#endif
                           )

          BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,mm)*fJacTilde(:,:) 
        END IF

          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem), &
                         !Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,ll,nn,oo,iElem,0),fJacTilde(:,:))
          !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(mm,ll)*fJacTilde(:,:) 
        !END IF
        

        !IF(nn.EQ.ll)THEN
          !gJacTilde(:,:) = (fJac(:,:,mm,nn,oo)*Metrics_gTilde(1,mm,nn,oo,iElem,0)  &
                           !+gJac(:,:,mm,nn,oo)*Metrics_gTilde(2,mm,nn,oo,iElem,0)  &
!#if PP_dim==3
                           !+hJac(:,:,mm,nn,oo)*Metrics_gTilde(3,mm,nn,oo,iElem,0)  &
!#endif
                           !)
          !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_diag(ll,nn)*gJacTilde(:,:) 
        !ELSE
!=====================================================
        !IF(ll.NE.nn)THEN
          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,ll,oo,iElem),UPrim(:,mm,ll,oo,iElem), &
                         !Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,ll,oo,iElem,0),gJacTilde(:,:))
          !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,nn)*gJacTilde(:,:) 
        !ELSE
          !gJacTilde(:,:) = (fJac(:,:,mm,nn,oo)*Metrics_gTilde(1,mm,nn,oo,iElem,0)  &
                           !+gJac(:,:,mm,nn,oo)*Metrics_gTilde(2,mm,nn,oo,iElem,0)  &
!#if PP_dim==3
                           !+hJac(:,:,mm,nn,oo)*Metrics_gTilde(3,mm,nn,oo,iElem,0)  &
!#endif
                           !)
          !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,nn)*gJacTilde(:,:) 
        !END IF
!=====================================================

          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,ll,oo,iElem),UPrim(:,mm,ll,oo,iElem), &
                         !Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,ll,oo,iElem,0),gJacTilde(:,:))
          !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(nn,ll)*gJacTilde(:,:) 
        !END IF

!#if PP_dim==3
        !IF(oo.EQ.ll)THEN
          !hJacTilde(:,:) = ( fJac(:,:,mm,nn,oo)*Metrics_hTilde(1,mm,nn,oo,iElem,0)  &
                            !+gJac(:,:,mm,nn,oo)*Metrics_hTilde(2,mm,nn,oo,iElem,0)  &
                            !+hJac(:,:,mm,nn,oo)*Metrics_hTilde(3,mm,nn,oo,iElem,0)) 
        !ELSE
          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,nn,ll,iElem),UPrim(:,mm,nn,ll,iElem), &
                         !Metrics_hTilde(:,mm,nn,oo,iElem,0),Metrics_hTilde(:,mm,nn,ll,iElem,0),hJacTilde(:,:))
        !END IF
        !BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,oo)*hJacTilde(:,:) 
!#endif
        r1=r1+PP_nVar
        r2=r2+vn1
#if PP_dim==3
        r3=r3+vn2
#endif
      END DO !ll
      DO tt=0,PP_N
        IF(tt.NE.mm)THEN
          CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,tt,nn,oo,iElem),UPrim(:,tt,nn,oo,iElem), &
              Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,tt,nn,oo,iElem,0),fJacTilde(:,:))
          BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(mm,tt)*fJacTilde(:,:)
        END IF
!=====================================================
        !IF(tt.NE.nn)THEN
          !CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,tt,oo,iElem),UPrim(:,mm,tt,oo,iElem), &
              !Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,tt,oo,iElem,0),gJacTilde(:,:))
          !BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(nn,tt)*gJacTilde(:,:) 
        !END IF
!=====================================================
      END DO
#else
      fJacTilde(:,:) = ( fJac(:,:,mm,nn,oo)*Metrics_fTilde(1,mm,nn,oo,iElem,0)  &
                       + gJac(:,:,mm,nn,oo)*Metrics_fTilde(2,mm,nn,oo,iElem,0)  &
#if PP_dim==3
                       + hJac(:,:,mm,nn,oo)*Metrics_fTilde(3,mm,nn,oo,iElem,0) &
#endif
                       )

      gJacTilde(:,:) = (fJac(:,:,mm,nn,oo)*Metrics_gTilde(1,mm,nn,oo,iElem,0)  &
                       +gJac(:,:,mm,nn,oo)*Metrics_gTilde(2,mm,nn,oo,iElem,0)  &
#if PP_dim==3
                       +hJac(:,:,mm,nn,oo)*Metrics_gTilde(3,mm,nn,oo,iElem,0) &
#endif
                       )

#if PP_dim==3
      hJacTilde(:,:) = fJac(:,:,mm,nn,oo)*Metrics_hTilde(1,mm,nn,oo,iElem,0)  &
                       +gJac(:,:,mm,nn,oo)*Metrics_hTilde(2,mm,nn,oo,iElem,0)  &
                       +hJac(:,:,mm,nn,oo)*Metrics_hTilde(3,mm,nn,oo,iElem,0) 
#endif
      r1=           vn1*nn+vn2*oo
      r2=mm*PP_nVar       +vn2*oo
#if PP_dim==3
      r3=mm*PP_nVar+vn1*nn
#endif
      ! weak formulation
      DO ll=0,PP_N
        BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,mm)*fJacTilde(:,:) 
        BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,nn)*gJacTilde(:,:) 
#if PP_dim==3
        BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,oo)*hJacTilde(:,:) 
#endif
        !BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,mm)*fJacTilde(:,:) 
        !BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,nn)*gJacTilde(:,:) 
#if PP_dim==3
        !BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) + DVolSurf_loc(ll,oo)*hJacTilde(:,:) 
#endif
        r1=r1+PP_nVar
        r2=r2+vn1
#if PP_dim==3
        r3=r3+vn2
#endif
      END DO !ll
#endif
      s=s+PP_nVar
    END DO !mm
  END DO !nn
END DO !oo 

!dRdU=0.
!s=1
!CALL VolInt(Ut)
!DO k=0,PP_NZ
  !DO j=0,PP_N
    !DO i=0,PP_N
      !DO jVar=1,PP_nVar
        !U(jVar,i,j,k,iElem) = Xk(jVar,i,j,k,iElem) + reps0
        !CALL ConsToPrim(PP_N,UPrim,U)
        !CALL VolInt(Ut_tilde)
        !U(jVar,i,j,k,iElem) = Xk(jVar,i,j,k,iElem) 
        !r=1
        !DO oo=0,PP_NZ
          !DO nn=0,PP_N
            !DO mm=0,PP_N
              !DO iVar=1,PP_nVar
                !dRdU(r,s) = dRdU(r,s)+(Ut_tilde(iVar,mm,nn,oo,iElem)-Ut(iVar,mm,nn,oo,iElem))*sreps0
                !r=r+1
              !END DO !iVar
            !END DO !mm
          !END DO !nn
        !END DO !oo
        !s=s+1
      !END DO !PP_nVar
    !END DO !i
  !END DO !j
!END DO !k
!WRITE (*,*) MAXVAL(ABS(dRdU-BJ))
!!BJ = dRdU
END SUBROUTINE DGVolIntJac

!===================================================================================================================================
!> Transformation of the block-Jacobian to physical coordinates (multiplication with sJ)
!===================================================================================================================================
SUBROUTINE Apply_sJ(BJ,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Implicit_Vars ,ONLY:nDOFVarElem
USE MOD_Mesh_Vars     ,ONLY:sJ
#if FV_ENABLED
USE MOD_FV_Vars       ,ONLY:FV_Elems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: BJ(nDOFVarElem,nDOFVarElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: r,s,i,j,k
INTEGER :: FVEM
!===================================================================================================================================
#if FV_ENABLED
FVEM = FV_Elems(iElem)
#else
FVEM = 0
#endif
DO s=0,nDOFVarElem-1,PP_nVar
  r=0
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = -BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)*sJ(i,j,k,iElem,FVEM)
        r=r+PP_nVar
      END DO !i
    END DO !j
  END DO !k
END DO ! s

END SUBROUTINE Apply_sJ

#if PARABOLIC
!===================================================================================================================================
! Derivative of the numerical flux h* of the gradient system with respect to U_jk
! h*= 0.5*(U_xi + U_neighbour_(-xi)) (mean value of the face values)
! => h*-U_xi=-0.5*U_xi + 0.5*U_neighbour_(-xi)
! => d(h*-U_xi)_jk(prim)/dU_jk(prim) = -0.5
!===================================================================================================================================
SUBROUTINE FillJacLiftingFlux(iElem)
! MODULES
USE MOD_Globals
USE MOD_Jac_ex_Vars               ,ONLY:JacLiftingFlux,Surf
USE MOD_GetBoundaryFlux_fd        ,ONLY:Lifting_GetBoundaryFlux_FD
USE MOD_Mesh_Vars                 ,ONLY:nBCSides,ElemToSide,S2V2
USE MOD_Mesh_Vars                 ,ONLY:NormVec,TangVec1,TangVec2,SurfElem
USE MOD_Precond_Vars              ,ONLY:tPrecond
USE MOD_DG_Vars                   ,ONLY:U_master
USE MOD_Mesh_Vars                 ,ONLY:Face_xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iVar,iLocSide,SideID
!===================================================================================================================================
#if PP_dim==3
DO iLocSide=1,6
#else
DO iLocSide=2,5
#endif 
  SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
  IF (SideID.LE.nBCSides) THEN !BCSides
    CALL Lifting_GetBoundaryFlux_FD(SideID,tPrecond,JacLiftingFlux(:,:,:,:,iLocSide),U_master, &
                                    SurfElem,Face_xGP,NormVec,TangVec1,TangVec2,S2V2(:,:,:,0,iLocSide)) !flip=0 for BCSide
  ELSE
    JacLiftingFlux(:,:,:,:,iLocSide)=0.
    DO iVar=1,PP_nVar
      JacLiftingFlux(iVar,iVar,:,:,iLocSide)=-0.5*Surf(:,:,iLocSide,iElem)
    END DO !iVar
  END IF !SideID
END DO!iLocSide
END SUBROUTINE FillJacLiftingFlux

SUBROUTINE JacLifting_VolInt(dir,iElem,JacLifting)
!===================================================================================================================================
! Computes the Volume gradient Jacobian of the BR2 scheme dQprim/dUprim (Q= Grad U)
! Normal vectors are supposed to point outwards!
!===================================================================================================================================
! MODULES
USE MOD_Jac_Ex_Vars        ,ONLY: LL_minus,LL_plus,nVec 
USE MOD_Jac_Ex_Vars        ,ONLY: JacLiftingFlux 
USE MOD_DG_Vars            ,ONLY: D
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,sJ   ! metrics
#if PP_dim==3
USE MOD_Mesh_Vars          ,ONLY: Metrics_hTilde
#endif
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                           :: iElem
INTEGER,INTENT(IN)                           :: dir
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                             :: JacLifting(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: i,j,k,ll
INTEGER                                      :: iVar
REAL                                         :: delta(1:PP_nVar,1:PP_nVar) 
!===================================================================================================================================
delta=0.
DO iVar=1,PP_nVar
  delta(iVar,iVar)=1.
END DO

JacLifting=0.
DO ll=0,PP_N
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        JacLifting(:,:,i,j,k,ll,1) = JacLifting(:,:,i,j,k,ll,1) +                                                  &
                                    sJ(i,j,k,iElem,0)*( D(i,ll)*Metrics_fTilde(dir,ll,j,k,iElem,0)*delta(:,:)      &
                                                     + nVec(dir,j,k,   XI_PLUS,iElem)*LL_plus(i,ll)                &
                                                       *JacLiftingFlux(:,:,j,k,XI_PLUS)                            &
                                                     + nVec(dir,j,k,  XI_MINUS,iElem)*LL_minus(i,ll)               &
                                                       *JacLiftingFlux(:,:,j,k,XI_MINUS) )
        JacLifting(:,:,i,j,k,ll,2) = JacLifting(:,:,i,j,k,ll,2) +                                                  &
                                    sJ(i,j,k,iElem,0)*( D(j,ll)*Metrics_gTilde(dir,i,ll,k,iElem,0)*delta(:,:)      &
                                                     + nVec(dir,i,k,  ETA_PLUS,iElem)*LL_plus(j,ll)                &
                                                       *JacLiftingFlux(:,:,i,k,ETA_PLUS)                           &
                                                     + nVec(dir,i,k, ETA_MINUS,iElem)*LL_minus(j,ll)               &
                                                       *JacLiftingFlux(:,:,i,k,ETA_MINUS) )
#if PP_dim==3
        JacLifting(:,:,i,j,k,ll,3) = JacLifting(:,:,i,j,k,ll,3) +                                                  &
                                    sJ(i,j,k,iElem,0)*( D(k,ll)*Metrics_hTilde(dir,i,j,ll,iElem,0)*delta(:,:)      &
                                                     + nVec(dir,i,j, ZETA_PLUS,iElem)*LL_plus(k,ll)                &
                                                       *JacLiftingFlux(:,:,i,k,ZETA_PLUS)                          &
                                                     + nVec(dir,i,j,ZETA_MINUS,iElem)*LL_minus(k,ll)               &
                                                       *JacLiftingFlux(:,:,i,j,ZETA_MINUS) )
#endif
      END DO !i
    END DO !j
  END DO !k
END DO !ll

    !DO ll=0,PP_N
      !DO k=0,PP_NZ
        !DO j=0,PP_N
          !DO i=0,PP_N
            !temp1 = D(i,ll)*Metrics_fTilde(dir,ll,j,k,iElem,0)
            !temp2 = D(j,ll)*Metrics_gTilde(dir,i,ll,k,iElem,0)
!#if PP_dim==3
            !temp3 = D(k,ll)*Metrics_hTilde(dir,i,j,ll,iElem,0)
!#endif
            !DO s=1,PP_nVar
              !!Compute only diagonal elements, since JacLiftingFlux is diagonal for all lifting_geboundaryflux cases
              !JacLifting(s,s,i,j,k,ll,1) = sJ(i,j,k,iElem,0)*(temp1 &
                                                         !+ nVec(dir,j,k,   XI_PLUS,iElem)*LL_plus(i,ll)                &
                                                           !*JacLiftingFlux(s,s,j,k,XI_PLUS)                            &
                                                         !+ nVec(dir,j,k,  XI_MINUS,iElem)*LL_minus(i,ll)               &
                                                           !*JacLiftingFlux(s,s,j,k,XI_MINUS) )
              !JacLifting(s,s,i,j,k,ll,2) = sJ(i,j,k,iElem,0)*( temp2 &
                                                         !+ nVec(dir,i,k,  ETA_PLUS,iElem)*LL_plus(j,ll)                &
                                                           !*JacLiftingFlux(s,s,i,k,ETA_PLUS)                           &
                                                         !+ nVec(dir,i,k, ETA_MINUS,iElem)*LL_minus(j,ll)               &
                                                           !*JacLiftingFlux(s,s,i,k,ETA_MINUS) )
!#if PP_dim==3
              !JacLifting(s,s,i,j,k,ll,3) = sJ(i,j,k,iElem,0)*( temp3 &
                                                     !+ nVec(dir,i,j, ZETA_PLUS,iElem)*LL_plus(k,ll)                &
                                                       !*JacLiftingFlux(s,s,i,k,ZETA_PLUS)                          &
                                                     !+ nVec(dir,i,j,ZETA_MINUS,iElem)*LL_minus(k,ll)               &
                                                       !*JacLiftingFlux(s,s,i,j,ZETA_MINUS) )
!#endif
            !END DO !s
          !END DO !i
        !END DO !j
      !END DO !k
    !END DO !ll
END SUBROUTINE JacLifting_VolInt

SUBROUTINE  DGVolIntGradJac(BJ,iElem)
!===================================================================================================================================
! volume integral: the total derivative of the viscous flux with resprect to U:
!                    dF^v/Du = dF^v/dQ* DQ/DU, with the gradient Q=Grad U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Precond_Vars  ,ONLY: NoFillIn
USE MOD_DG_Vars       ,ONLY: D_hat
USE MOD_DG_Vars       ,ONLY: U,UPrim,nDOFElem
USE MOD_Mesh_Vars     ,ONLY: Metrics_fTilde,Metrics_gTilde
#if PP_dim==3
USE MOD_Mesh_Vars     ,ONLY: Metrics_hTilde
#endif
USE MOD_Implicit_Vars ,ONLY: nDOFVarElem
USE MOD_GradJacobian  ,ONLY: EvalFluxGradJacobian
!USE MOD_Implicit_Vars ,ONLY: reps0,sreps0
!USE MOD_EOS           ,ONLY: ConsToPrim
!#if PARABOLIC
!USE MOD_EOS           ,ONLY: ConsToPrimLifting_loc
!USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
!USE MOD_Lifting_Vars ,ONLY: gradUx_cons,gradUy_cons,gradUz_cons
!#endif
!USE MOD_Jac_Ex_Vars   ,ONLY: PrimConsJac 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: BJ(1:nDOFVarElem,1:nDOFVarElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 

!REAL,DIMENSION(PP_nVar) :: f,g,h,u1,u1_tilde,gx,gy,gz,gy_tilde
!REAL,DIMENSION(PP_nVarPrim) :: up,up_tilde,gpx,gpy,gpz
!REAL,DIMENSION(PP_nVar) :: f_Tilde,g_Tilde,h_Tilde
!REAL,DIMENSION(PP_nVar,PP_nVar)   :: dfdu,dgdu,dhdu
INTEGER                                                 :: i,j,k,l,mm,nn,oo
!INTEGER                                                 :: iVar,jVar
INTEGER                                                 :: r,s,vn1,vn2
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: fJacQx,gJacQx,hJacQx
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: fJacQy,gJacQy,hJacQy
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: fJacQz,gJacQz,hJacQz
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJacTilde,gJacTilde
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJac,gJac
#if PP_dim==3
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: hJac
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: hJacTilde
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N,3):: JacLifting_Z
#endif
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim):: JacLifting_X,JacLifting_Y
!===================================================================================================================================
vn1=PP_nVar*(PP_N+1)
vn2=vn1*(PP_N+1)
! Calculation of the derivative of the diffion flux with respect to gradU
! dF/dQ
CALL EvalFluxGradJacobian(nDOFElem,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem) &
                          ,fJacQx,fJacQy,fJacQz &
                          ,gJacQx,gJacQy,gJacQz &
                          ,hJacQx,hJacQy,hJacQz ) 

DO k=0,PP_NZ
  DO j=0,PP_N
    DO i=0,PP_N
      fJacTilde=fJacQx(:,:,i,j,k)
      gJacTilde=gJacQx(:,:,i,j,k)
#if PP_dim==3
      hJacTilde=hJacQx(:,:,i,j,k)
#endif
      ! Compute the transformed fluxes with the metric terms
      ! Attention 1: we store the transformed fluxes in f,g,h again
      fJacQx(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_fTilde(1,i,j,k,iElem,0) &
                          + gJacTilde(:,:)*Metrics_fTilde(2,i,j,k,iElem,0) &
#if PP_dim==3
                          + hJacTilde(:,:)*Metrics_fTilde(3,i,j,k,iElem,0) &
#endif
                          )
      gJacQx(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_gTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_gTilde(2,i,j,k,iElem,0) &
#if PP_dim==3
                          +hJacTilde(:,:)*Metrics_gTilde(3,i,j,k,iElem,0) &
#endif
                          )
#if PP_dim==3
      hJacQx(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_hTilde(1,i,j,k,iElem,0)  &
                          +gJacTilde(:,:)*Metrics_hTilde(2,i,j,k,iElem,0) &
                          +hJacTilde(:,:)*Metrics_hTilde(3,i,j,k,iElem,0) &
                          )
      hJacTilde=hJacQy(:,:,i,j,k)
#endif
      fJacTilde=fJacQy(:,:,i,j,k)
      gJacTilde=gJacQy(:,:,i,j,k)
      ! Compute the transformed fluxes with the metric terms
      ! Attention 1: we store the transformed fluxes in f,g,h again
      fJacQy(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_fTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_fTilde(2,i,j,k,iElem,0) &
#if PP_dim==3
                          +hJacTilde(:,:)*Metrics_fTilde(3,i,j,k,iElem,0) &
#endif
                          )
      gJacQy(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_gTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_gTilde(2,i,j,k,iElem,0) &
#if PP_dim==3
                          +hJacTilde(:,:)*Metrics_gTilde(3,i,j,k,iElem,0) &
#endif
                          )
#if PP_dim==3
      hJacQy(:,:,i,j,k) = fJacTilde(:,:)*Metrics_hTilde(1,i,j,k,iElem,0) + &
                          gJacTilde(:,:)*Metrics_hTilde(2,i,j,k,iElem,0) + &
                          hJacTilde(:,:)*Metrics_hTilde(3,i,j,k,iElem,0)
      hJacTilde=hJacQz(:,:,i,j,k)
      fJacTilde=fJacQz(:,:,i,j,k)
      gJacTilde=gJacQz(:,:,i,j,k)
      ! Compute the transformed fluxes with the metric terms
      ! Attention 1: we store the transformed fluxes in f,g,h again
      fJacQz(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_fTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_fTilde(2,i,j,k,iElem,0) &
                          +hJacTilde(:,:)*Metrics_fTilde(3,i,j,k,iElem,0) &
                          )
      gJacQz(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_gTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_gTilde(2,i,j,k,iElem,0) &
                          +hJacTilde(:,:)*Metrics_gTilde(3,i,j,k,iElem,0) &
                          )
      hJacQz(:,:,i,j,k) = fJacTilde(:,:)*Metrics_hTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_hTilde(2,i,j,k,iElem,0) &
                          +hJacTilde(:,:)*Metrics_hTilde(3,i,j,k,iElem,0) 
#endif
    END DO ! i
  END DO ! j
END DO ! k

CALL JacLifting_VolInt(1,iElem,JacLifting_X) !d(Q^1)/dU(1:3) (3 sums)
CALL JacLifting_VolInt(2,iElem,JacLifting_Y) !d(Q^2)/dU(1:3)
#if PP_dim==3
CALL JacLifting_VolInt(3,iElem,JacLifting_Z) !d(Q^3)/dU(1:3)
#endif

s=0
DO oo=0,PP_NZ
  DO nn=0,PP_N
    DO mm=0,PP_N
      IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
      DO j=0,PP_N
        fJac(:,:)= ( MATMUL(fJacQx(:,:,mm,j,oo) , JacLifting_X(:,:,mm,j,oo,nn,2) )  &
                   +MATMUL(fJacQy(:,:,mm,j,oo) , JacLifting_Y(:,:,mm,j,oo,nn,2) )  &
#if PP_dim==3
                   +MATMUL(fJacQz(:,:,mm,j,oo) , JacLifting_Z(:,:,mm,j,oo,nn,2) )  &
#endif
                   )
        !fJac(:,:) = MATMUL(fJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO i=0,PP_N
          r=PP_nVar*i+vn1*j+vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(i,mm)*fJac(:,:)
        END DO !i
      END DO !j
      DO i=0,PP_N
        gJac(:,:)= ( MATMUL(gJacQx(:,:,i,nn,oo) , JacLifting_X(:,:,i,nn,oo,mm,1) )  &
                   + MATMUL(gJacQy(:,:,i,nn,oo) , JacLifting_Y(:,:,i,nn,oo,mm,1) )  &
#if PP_dim==3
                   + MATMUL(gJacQz(:,:,i,nn,oo) , JacLifting_Z(:,:,i,nn,oo,mm,1) )  &
#endif
                   )
        !gJac(:,:) = MATMUL(gJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO j=0,PP_N
          r=PP_nVar*i+vn1*j+vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(j,nn)*gJac(:,:)
        END DO !j
      END DO !i
#if PP_dim==3
      DO k=0,PP_N
        fJac(:,:)= ( MATMUL(fJacQx(:,:,mm,nn,k) , JacLifting_X(:,:,mm,nn,k,oo,3) )  &
                   +MATMUL(fJacQy(:,:,mm,nn,k) , JacLifting_Y(:,:,mm,nn,k,oo,3) )  &
                   +MATMUL(fJacQz(:,:,mm,nn,k) , JacLifting_Z(:,:,mm,nn,k,oo,3) ) &
                   )
        !fJac(:,:) = MATMUL(fJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO i=0,PP_N
          r=PP_nVar*i+vn1*nn+vn2*k
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(i,mm)*fJac(:,:)
        END DO !i
      END DO !k 
      DO i=0,PP_N
        hJac(:,:)=  MATMUL(hJacQx(:,:,i,nn,oo) , JacLifting_X(:,:,i,nn,oo,mm,1) )  &
                   +MATMUL(hJacQy(:,:,i,nn,oo) , JacLifting_Y(:,:,i,nn,oo,mm,1) )  &
                   +MATMUL(hJacQz(:,:,i,nn,oo) , JacLifting_Z(:,:,i,nn,oo,mm,1) ) 
        !hJac(:,:) = MATMUL(hJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO k=0,PP_N
          r=PP_nVar*i+vn1*nn+vn2*k
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(k,oo)*hJac(:,:)
        END DO !k
      END DO !i 
      DO k=0,PP_NZ
        gJac(:,:)= ( MATMUL(gJacQx(:,:,mm,nn,k) , JacLifting_X(:,:,mm,nn,k,oo,3) )  &
                   +MATMUL(gJacQy(:,:,mm,nn,k) , JacLifting_Y(:,:,mm,nn,k,oo,3) )  &
                   +MATMUL(gJacQz(:,:,mm,nn,k) , JacLifting_Z(:,:,mm,nn,k,oo,3) )  &
                   )
        !gJac(:,:) = MATMUL(gJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO j=0,PP_N
          r=PP_nVar*mm+vn1*j+vn2*k
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(j,nn)*gJac(:,:)
        END DO !j
      END DO !k 
      DO j=0,PP_N
        hJac(:,:)=  MATMUL(hJacQx(:,:,mm,j,oo) , JacLifting_X(:,:,mm,j,oo,nn,2) )  &
                   +MATMUL(hJacQy(:,:,mm,j,oo) , JacLifting_Y(:,:,mm,j,oo,nn,2) )  &
                   +MATMUL(hJacQz(:,:,mm,j,oo) , JacLifting_Z(:,:,mm,j,oo,nn,2) )
        !hJac(:,:) = MATMUL(hJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO k=0,PP_NZ
          r=PP_nVar*mm+vn1*j+vn2*k
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(k,oo)*hJac(:,:)
        END DO !k
      END DO !j 
#endif
    END IF !NoFillIn
      DO l=0,PP_N
        fJac(:,:)= ( MATMUL(fJacQx(:,:,l,nn,oo) , JacLifting_X(:,:,l,nn,oo,mm,1) )  &
                   +MATMUL(fJacQy(:,:,l,nn,oo) , JacLifting_Y(:,:,l,nn,oo,mm,1) )  &
#if PP_dim==3
                   +MATMUL(fJacQz(:,:,l,nn,oo) , JacLifting_Z(:,:,l,nn,oo,mm,1) )  &
#endif
                   )
        !TODO:Besser umschreiben, Summe ueber l ausserhalb, so dass nur einmal mit ConToPrimJac multipliziert werden muss
        !fJac(:,:) = MATMUL(fJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO i=0,PP_N
          r=PP_nVar*i+vn1*nn+vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(i,l)*fJac(:,:)
        END DO !i
        gJac(:,:)= ( MATMUL(gJacQx(:,:,mm,l,oo) , JacLifting_X(:,:,mm,l,oo,nn,2) )  &
                   +MATMUL(gJacQy(:,:,mm,l,oo) , JacLifting_Y(:,:,mm,l,oo,nn,2) )  &
#if PP_dim==3
                   +MATMUL(gJacQz(:,:,mm,l,oo) , JacLifting_Z(:,:,mm,l,oo,nn,2) )  &
#endif
                  )
        !gJac(:,:) = MATMUL(gJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO j=0,PP_N
          r=PP_nVar*mm+vn1*j+vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(j,l)*gJac(:,:)
        END DO !j
#if PP_dim==3
        hJac(:,:)=  MATMUL(hJacQx(:,:,mm,nn,l) , JacLifting_X(:,:,mm,nn,l,oo,3) )  &
                   +MATMUL(hJacQy(:,:,mm,nn,l) , JacLifting_Y(:,:,mm,nn,l,oo,3) )  &
                   +MATMUL(hJacQz(:,:,mm,nn,l) , JacLifting_Z(:,:,mm,nn,l,oo,3) )
        !hJac(:,:) = MATMUL(hJacTilde(:,:),PrimConsJac(:,:,mm,nn,oo))
        DO k=0,PP_N
          r=PP_nVar*mm+vn1*nn+vn2*k
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(k,l)*hJac(:,:)
        END DO !k
#endif
      END DO !l 
      s=s+PP_nVar
    END DO !mm
  END DO !nn
END DO !oo 


! Debugging test for derivatives of advective flux (fd and exact)
!U1=U(:,1,1,0,iElem)
!UP=UPrim(:,1,1,0,iElem)
!gx=gradux_cons(:,1,1,0,iElem)
!gy=graduy_cons(:,1,1,0,iElem)
!gz=graduz_cons(:,1,1,0,iElem)
!gpx=gradux(:,1,1,0,iElem)
!gpy=graduy(:,1,1,0,iElem)
!gpz=graduz(:,1,1,0,iElem)

!dfdu=0.
!dgdu=0.
!dhdu=0.
!CALL EvalDiffFlux3D(UP,gpx,gpy,gpz,f,g,h)
!!CALL EvalDiffFlux3D(UP,gx,gy,gz,f,g,h)
!gy_tilde=gy
!DO jVar=1,PP_nVar
  !IF(jVar==4) CYCLE
  !gy_tilde(jVar) = gy(jVar) + reps0
  !CALL ConsToPrimLifting_loc(gx,gy_Tilde,gz,UP,gpx,gpy,gpz)
  !CALL EvalDiffFlux3D(UP,gpx,gpy,gpz,f_tilde,g_tilde,h_tilde)
  !!CALL EvalDiffFlux3D(UP,gx_Tilde,gy,gz,f_tilde,g_tilde,h_tilde)
  !gy_tilde(jVar)=gy(jVar)
      !DO iVar=1,PP_nVar
        !dFdU(iVar,jVar) = (f_tilde(iVar)-f(iVar))*sreps0
        !dgdU(iVar,jVar) = (g_tilde(iVar)-g(iVar))*sreps0
        !dhdU(iVar,jVar) = (h_tilde(iVar)-h(iVar))*sreps0
      !END DO ! iVar
!END DO
!DO jVar=1,PP_nVar
  !WRITE(*,*) jVar
!WRITE(*,*) dfdu(:,jVar)
!WRITE(*,*)fJacQy(:,jVar,1,1,0)
!WRITE(*,*) '--------------------------------------------'
!WRITE(*,*) dgdu(:,jVar)
!WRITE(*,*)gJacQy(:,jVar,1,1,0)
!WRITE(*,*) '--------------------------------------------'
!WRITE(*,*) dhdu(:,jVar)
!WRITE(*,*)hJacQy(:,jVar,1,1,0)
!WRITE(*,*) '--------------------------------------------'
!eND DO
!stop

END SUBROUTINE DGVolIntGradJac

SUBROUTINE BuildnVecTangSurf()
!===================================================================================================================================
! used for BR2: normal vectors, outward pointing and sorted in ijk element fashion!!! 
! The usual NormVec is only outward poiting for the master Sides.
! nVec are also outward pointing for the slave Sides
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: Normvec,SurfElem,ElemToSide,nElems 
USE MOD_Jac_ex_Vars        ,ONLY: nVec,Surf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: iElem
INTEGER                                      :: p,q,SideID,Flip
!===================================================================================================================================
DO iElem=1,nElems
  SideID=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,XI_MINUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      nVec( :,q,p,XI_MINUS,iElem)          =  NormVec( :,p,q,0,SideID)
      Surf(   q,p,XI_MINUS,iElem)          =  SurfElem(  p,q,0,SideID)
#else
      nVec (:,PP_N-p,q,XI_MINUS,iElem)          =  NormVec (:,p,q,0,SideID)
      Surf (  PP_N-p,q,XI_MINUS,iElem)          =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      nVec( :,p,q,XI_MINUS,iElem)           = -NormVec( :,p,q,0,SideID)
      Surf(   p,q,XI_MINUS,iElem)           =  SurfElem(  p,q,0,SideID)
#else
      nVec (:,p,q,XI_MINUS,iElem)      = -NormVec (:,p,q,0,SideID)
      Surf (  p,q,XI_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
#if PP_dim==3
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,PP_N-q,p,XI_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 4!
      nVec( :,q,PP_N-p,XI_MINUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   q,PP_N-p,XI_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-p,PP_N-q,XI_MINUS,iElem) = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-p,PP_N-q,XI_MINUS,iElem) =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,q,PP_N-p,XI_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 2!
      nVec( :,PP_N-q,p,XI_MINUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-q,p,XI_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
#endif
  END SELECT

  
  SideID=ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,XI_PLUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      nVec( :,p,q,XI_PLUS,iElem)           =  NormVec( :,p,q,0,SideID)
      Surf(   p,q,XI_PLUS,iElem)           =  SurfElem(  p,q,0,SideID)
#else
      nVec (:,p,q,XI_PLUS,iElem)          =  NormVec (:,p,q,0,SideID)
      Surf (  p,q,XI_PLUS,iElem)          =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      nVec( :,q,p,XI_PLUS,iElem)           = -NormVec( :,p,q,0,SideID)
      Surf(   q,p,XI_PLUS,iElem)           =  SurfElem(  p,q,0,SideID)
#else
      nVec (:,PP_N-p,q,XI_PLUS,iElem)          = -NormVec (:,p,q,0,SideID)
      Surf (  PP_N-p,q,XI_PLUS,iElem)          =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
#if PP_dim==3
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-p,q,XI_PLUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-p,q,XI_PLUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-q,PP_N-p,XI_PLUS,iElem) = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-q,PP_N-p,XI_PLUS,iElem) =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,p,PP_N-q,XI_PLUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   p,PP_N-q,XI_PLUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
#endif
  END SELECT
  SideID=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ETA_MINUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      nVec( :,p,q,ETA_MINUS,iElem)           =  NormVec( :,p,q,0,SideID)
      Surf(   p,q,ETA_MINUS,iElem)           =  SurfElem(  p,q,0,SideID)
#else
      nVec (:,p,q,ETA_MINUS,iElem)          =  NormVec (:,p,q,0,SideID)
      Surf (  p,q,ETA_MINUS,iElem)          =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      nVec( :,q,p,ETA_MINUS,iElem)           = -NormVec( :,p,q,0,SideID)
      Surf(   q,p,ETA_MINUS,iElem)           =  SurfElem(  p,q,0,SideID)
#else
      nVec (:,PP_N-p,q,ETA_MINUS,iElem)      = -NormVec (:,p,q,0,SideID)
      Surf (  PP_N-p,q,ETA_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
#if PP_dim==3
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-p,q,ETA_MINUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-p,q,ETA_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-q,PP_N-p,ETA_MINUS,iElem) = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-q,PP_N-p,ETA_MINUS,iElem) =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,p,PP_N-q,ETA_MINUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   p,PP_N-q,ETA_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
#endif
  END SELECT

  SideID=ElemToSide(E2S_SIDE_ID,ETA_PLUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ETA_PLUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      nVec( :,PP_N-p,q,ETA_PLUS,iElem)           = NormVec( :,p,q,0,SideID)
      Surf(   PP_N-p,q,ETA_PLUS,iElem)           = SurfElem(  p,q,0,SideID)
#else
      nVec (:,PP_N-p,q,ETA_PLUS,iElem)          =  NormVec (:,p,q,0,SideID)
      Surf (  PP_N-p,q,ETA_PLUS,iElem)          =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_NZ; DO p=0,PP_N
      !nVecSurf(:,q,PP_N-p,ETA_PLUS,iElem)           = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 3!
#if PP_dim==3
      nVec( :,PP_N-q,p,ETA_PLUS,iElem)           = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-q,p,ETA_PLUS,iElem)           =  SurfElem(  p,q,0,SideID)
#else
      nVec (:,p,q,ETA_PLUS,iElem)      = -NormVec (:,p,q,0,SideID)
      Surf (  p,q,ETA_PLUS,iElem)      =  SurfElem(  p,q,0,SideID)
#endif
    END DO; END DO 
#if PP_dim==3
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,p,q,ETA_PLUS,iElem)                = -NormVec (:,p,q,0,SideID)
      Surf(   p,q,ETA_PLUS,iElem)                =  SurfElem  (p,q,0,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,PP_N-q,p,ETA_PLUS,iElem)           = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 1!
      nVec( :,q,PP_N-p,ETA_PLUS,iElem)           = -NormVec( :,p,q,0,SideID)
      Surf(   q,PP_N-p,ETA_PLUS,iElem)           =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-p,PP_N-q,ETA_PLUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-p,PP_N-q,ETA_PLUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
#endif
  END SELECT
#if PP_dim==3
  SideID=ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ZETA_MINUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,q,p,ZETA_MINUS,iElem)           = NormVec( :,p,q,0,SideID)
      Surf(   q,p,ZETA_MINUS,iElem)           = SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,p,q,ZETA_MINUS,iElem)           = -NormVec( :,p,q,0,SideID)
      Surf(   p,q,ZETA_MINUS,iElem)           =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,PP_N-q,p,ZETA_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 4!
      nVec( :,q,PP_N-p,ZETA_MINUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   q,PP_N-p,ZETA_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-p,PP_N-q,ZETA_MINUS,iElem) = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-p,PP_N-q,ZETA_MINUS,iElem) =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,q,PP_N-p,ZETA_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 2!
      nVec( :,PP_N-q,p,ZETA_MINUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-q,p,ZETA_MINUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  END SELECT
  SideID=ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ZETA_PLUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,p,q,ZETA_PLUS,iElem)           = NormVec( :,p,q,0,SideID)
      Surf(   p,q,ZETA_PLUS,iElem)           = SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,q,p,ZETA_PLUS,iElem)           = -NormVec( :,p,q,0,SideID)
      Surf(   q,p,ZETA_PLUS,iElem)           =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-p,q,ZETA_PLUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-p,q,ZETA_PLUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,PP_N-q,PP_N-p,ZETA_PLUS,iElem) = -NormVec( :,p,q,0,SideID)
      Surf(   PP_N-q,PP_N-p,ZETA_PLUS,iElem) =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec( :,p,PP_N-q,ZETA_PLUS,iElem)      = -NormVec( :,p,q,0,SideID)
      Surf(   p,PP_N-q,ZETA_PLUS,iElem)      =  SurfElem(  p,q,0,SideID)
    END DO; END DO 
  END SELECT
#endif
END DO !iElem

END SUBROUTINE BuildnVecTangSurf

SUBROUTINE Build_BR2_SurfTerms()
!===================================================================================================================================
!used for BR2: normal vectors, outward pointing and sorted in ijk element fashion!!! 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: sJ,ElemToSide,nElems
USE MOD_Mesh_Vars          ,ONLY: nSides,NormVec
USE MOD_Jac_ex_Vars        ,ONLY: R_Minus,R_Plus,LL_Minus,LL_plus
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY:nNbProcs
USE MOD_MPI                ,ONLY:StartSendMPIData,StartReceiveMPIData,FinishExchangeMPIData
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: iElem,p,q,l,iLocSide,SideID,Flip
REAL                                         :: RFace(0:PP_N,0:PP_NZ)
#if USE_MPI
!INTEGER :: MPIRequest_R(nNbProcs,2,2)
INTEGER :: MPIRequest_RPlus(nNbProcs,2)
INTEGER :: MPIRequest_RMinus(nNbProcs,2)
#endif 
!===================================================================================================================================

ALLOCATE(R_Minus(3,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(R_Plus(3,0:PP_N,0:PP_NZ,1:nSides))
R_Minus=0.
R_Plus =0.
DO iElem=1,nElems
#if PP_dim==3
  DO iLocSide=1,6
#else
  DO iLocSide=2,5
#endif 
    RFace=0.
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    Flip  =ElemToSide(E2S_FLIP,iLocSide,iElem)
    SELECT CASE(ilocSide)
    CASE(XI_MINUS)
      DO q=0,PP_NZ
        DO p=0,PP_N
          DO l=0,PP_N
            ! switch to right hand system
#if PP_dim==3
            RFace(q,p)=RFace(q,p)+sJ(l,p,q,iElem,0)*LL_Minus(l,l)
#else
            RFace(PP_N-p,0)=RFace(PP_N-p,0)+sJ(l,p,q,iElem,0)*LL_Minus(l,l)
#endif
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_MINUS)
      DO q=0,PP_NZ
        DO p=0,PP_N
          DO l=0,PP_N
#if PP_dim==3
            RFace(p,q)=RFace(p,q)+sJ(p,l,q,iElem,0)*LL_Minus(l,l)
#else
            RFace(p,0)=RFace(p,0)+sJ(p,l,q,iElem,0)*LL_Minus(l,l)
#endif
          END DO ! l
        END DO ! p
      END DO ! q
#if PP_dim==3
    CASE(ZETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          DO l=0,PP_N
            ! switch to right hand system
            RFace(q,p)=RFace(q,p)+sJ(p,q,l,iElem,0)*LL_Minus(l,l)
          END DO ! l
        END DO ! p
      END DO ! q
#endif
    CASE(XI_PLUS)
      DO q=0,PP_NZ
        DO p=0,PP_N
          DO l=0,PP_N
#if PP_dim==3
            RFace(p,q)=RFace(p,q)+sJ(l,p,q,iElem,0)*LL_Plus(l,l)
#else
            RFace(p,0)=RFace(p,0)+sJ(l,p,q,iElem,0)*LL_Plus(l,l)
#endif
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_PLUS)
      DO q=0,PP_NZ
        DO p=0,PP_N
          DO l=0,PP_N
            ! switch to right hand system
#if PP_dim==3
            RFace(PP_N-p,q)=RFace(PP_N-p,q)+sJ(p,l,q,iElem,0)*LL_Plus(l,l)
#else
            RFace(PP_N-p,0)=RFace(PP_N-p,0)+sJ(p,l,q,iElem,0)*LL_Plus(l,l)
#endif
          END DO ! l
        END DO ! p
      END DO ! q
#if PP_dim==3
    CASE(ZETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          DO l=0,PP_N
            RFace(p,q)=RFace(p,q)+sJ(p,q,l,iElem,0)*LL_Plus(l,l)
          END DO ! l
        END DO ! p
      END DO ! q
#endif
    END SELECT
    SELECT CASE(Flip)
      CASE(0) ! master side
        DO q=0,PP_NZ
          DO p=0,PP_N
            R_Minus(:,p,q,SideID)=RFace(p,q)*NormVec(:,p,q,0,SideID)
          END DO ! p
        END DO ! q
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_NZ
          DO p=0,PP_N
#if PP_dim==3
            R_Plus(:,p,q,SideID)=-RFace(q,p)*NormVec(:,p,q,0,SideID)
#else
            R_Plus(:,p,q,SideID)=-RFace(PP_N-p,0)*NormVec(:,p,q,0,SideID)
#endif
          END DO ! p
        END DO ! q
#if PP_dim==3
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N
          DO p=0,PP_N
            R_Plus(:,p,q,SideID)=-RFace(PP_N-p,q)*NormVec(:,p,q,0,SideID)
          END DO ! p
        END DO ! q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N
          DO p=0,PP_N
            R_Plus(:,p,q,SideID)=-RFace(PP_N-q,PP_N-p)*NormVec(:,p,q,0,SideID)
          END DO ! p
        END DO ! q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N
          DO p=0,PP_N
            R_Plus(:,p,q,SideID)=-RFace(p,PP_N-q)*NormVec(:,p,q,0,SideID)
          END DO ! p
        END DO ! q
#endif
    END SELECT
  END DO !iLocSide
END DO !iElem

#if USE_MPI
!EXCHANGE R_Minus and R_Plus vice versa !!!!!
MPIRequest_RPlus=0
MPIRequest_RMinus=0
CALL StartReceiveMPIData(R_Plus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RPlus(:,RECV),SendID=2) ! Receive MINE / Geo: slave -> master
CALL StartSendMPIData(   R_Plus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RPlus(:,SEND),SendID=2) ! SEND YOUR / Geo: slave -> master
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_RPlus) 

CALL StartReceiveMPIData(R_Minus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RMinus(:,RECV),SendID=1) ! Receive YOUR / Geo: master -> slave
CALL StartSendMPIData(R_Minus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RMinus(:,SEND),SendID=1) ! SEND MINE / Geo: master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_RMinus) 
#endif

END SUBROUTINE Build_BR2_SurfTerms
#endif /*PARABOLIC*/

#if FV_ENABLED
!===================================================================================================================================
!> Volume integral Jacobian of the convective flux for FV elements.
!> Each finite volume subcell requires the solution of a Riemann problem. Hence, the derivation is done using a finite difference
!> for the derivation of the Riemann problem (same is always done for the surface fluxes). It also takes into account the additional
!> dependencies caused by the 2nd order reconstruction procedure. Here, a primitive reconstruction is assumed resulting in a back
!> and forth transformation from conservative to primitive. The additional dependencies are stored in dUdUvol_minus/plus and then
!> multiplied with the finite difference of the Reiemann flux. Finally the derivatives are assembled at the right position of
!> the block-Jacobian. Note that the derivation of the volume integral has to be done directionwise.
!===================================================================================================================================
SUBROUTINE  FVVolIntJac(BJ,iElem)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars             ,ONLY: nDOFElem
#if PP_dim == 3
USE MOD_FV_Vars             ,ONLY: FV_NormVecZeta,FV_TangVec1Zeta,FV_TangVec2Zeta
USE MOD_FV_Vars             ,ONLY: FV_SurfElemZeta_sw
#endif
USE MOD_FV_Vars             ,ONLY: FV_NormVecXi,FV_TangVec1Xi,FV_TangVec2Xi
USE MOD_FV_Vars             ,ONLY: FV_NormVecEta,FV_TangVec1Eta,FV_TangVec2Eta
USE MOD_FV_Vars             ,ONLY: FV_SurfElemXi_sw,FV_SurfElemEta_sw
USE MOD_Implicit_Vars       ,ONLY: rEps0
USE MOD_Riemann             ,ONLY: Riemann_Point
USE MOD_EOS                 ,ONLY: ConsToPrim,PrimToCons
USE MOD_DG_Vars             ,ONLY: UPrim
#if FV_RECONSTRUCT
USE MOD_Mesh_Vars           ,ONLY: ElemToSide
USE MOD_Jac_Ex_Vars         ,ONLY: UPrim_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended
USE MOD_Jac_Reconstruction  ,ONLY: FV_Reconstruction_Derivative
USE MOD_FV_Vars             ,ONLY: gradUxi,gradUeta,FV_dx_XI_L,FV_dx_XI_R,FV_dx_ETA_L,FV_dx_ETA_R
#if PP_dim == 3
USE MOD_Jac_Ex_Vars         ,ONLY: FV_sdx_ZETA_extended
USE MOD_FV_Vars             ,ONLY: gradUzeta,FV_dx_ZETA_L,FV_dx_ZETA_R
#endif
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: iElem                                     !< current element index
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: BJ(1:nDOFElem*PP_nVar,1:nDOFElem*PP_nVar) !< block-Jacobian of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ):: dFdU_plus,dFdU_minus
INTEGER                                                  :: i,j,p,q,iVar,jVar
INTEGER                                                  :: s,r,vn1,vn2
REAL,DIMENSION(PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ)        :: F,F_Tilde
REAL,DIMENSION(PP_nVar    ,0:PP_N)                       :: FV_U_plus_Tilde,FV_U_minus_Tilde
REAL,DIMENSION(PP_nVarPrim,0:PP_N)                       :: FV_UPrim_plus_Tilde,FV_UPrim_minus_Tilde
REAL,DIMENSION(PP_nVar    ,0:PP_N)                       :: FV_U_Ximinus_loc,FV_U_Xiplus_loc,FV_U_Etaminus_loc,FV_U_Etaplus_loc
REAL,DIMENSION(PP_nVarPrim,0:PP_N)                       :: FV_UPrim_Ximinus_loc,FV_UPrim_Xiplus_loc
REAL,DIMENSION(PP_nVarPrim,0:PP_N)                       :: FV_UPrim_Etaminus_loc,FV_UPrim_Etaplus_loc
#if PP_dim==3
REAL,DIMENSION(PP_nVar    ,0:PP_N)                       :: FV_U_Zetaminus_loc,FV_U_Zetaplus_loc
REAL,DIMENSION(PP_nVarPrim,0:PP_N)                       :: FV_UPrim_Zetaminus_loc,FV_UPrim_Zetaplus_loc
INTEGER                                                  :: k
#endif
REAL                                                     :: reps0_L,sreps0_L,reps0_R,sreps0_R
#ifdef LEVELSET
REAL                                                     :: LS_L,LS_R
REAL                                                     :: curv 
REAL                                                     :: curv_L,curv_R 
REAL                                                     :: U_dummy1(PP_nVar),U_dummy2(PP_nVar),vel_dummy(3) 
REAL,DIMENSION(3)                                        :: NormVec_L,NormVec_R,NormVec_LS,TangVec1_LS,TangVec2_LS
#endif
#if FV_RECONSTRUCT
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,1:3)           :: dUdUvol_plus,dUdUvol_minus
INTEGER                                                  :: m,l
INTEGER                                                  :: SideID_minus,SideID_plus
#endif
REAL,DIMENSION(PP_nVar,PP_nVar)                          :: matrix !< has to be used due to intel compiler error with matmul
REAL,DIMENSION(PP_nVar,PP_nVar)                          :: vector !< has to be used due to intel compiler error with matmul
!===================================================================================================================================
vn1 = PP_nVar*(PP_N+1)
vn2 = vn1*(PP_N+1)
dfDU_plus  = 0
dFdU_minus = 0


! === Xi-Direction ================================================================================================================
#if FV_RECONSTRUCT
SideID_minus         = ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
SideID_plus          = ElemToSide(E2S_SIDE_ID,XI_PLUS ,iElem)
#endif

DO p=0,PP_N
  DO q=0,PP_NZ
    DO i=0,PP_N
#if FV_RECONSTRUCT
      FV_UPrim_Xiplus_loc( :,i) = UPrim(:,i,p,q,iElem) + gradUxi(:,p,q,i,iElem) * FV_dx_XI_R(p,q,i,iElem)
      FV_UPrim_Ximinus_loc(:,i) = UPrim(:,i,p,q,iElem) - gradUxi(:,p,q,i,iElem) * FV_dx_XI_L(p,q,i,iElem)
#else
      FV_UPrim_Xiplus_loc( :,i) = UPrim(:,i,p,q,iElem)
      FV_UPrim_Ximinus_loc(:,i) = UPrim(:,i,p,q,iElem)
#endif
      CALL PrimToCons(FV_UPrim_Ximinus_loc(:,i),FV_U_Ximinus_loc(:,i))
      CALL PrimToCons(FV_UPrim_Xiplus_loc( :,i),FV_U_Xiplus_loc( :,i))
    END DO
    FV_UPrim_plus_Tilde  = FV_UPrim_Xiplus_loc( :,:)
    FV_UPrim_minus_Tilde = FV_UPrim_Ximinus_loc(:,:)
    FV_U_plus_Tilde      = FV_U_Xiplus_loc(     :,:)
    FV_U_minus_Tilde     = FV_U_Ximinus_loc(    :,:)
#if FV_RECONSTRUCT
    CALL FV_Reconstruction_Derivative(FV_sdx_XI_extended(p,q,:,iElem),FV_dx_XI_L(p,q,:,iElem),FV_dx_XI_R(p,q,:,iElem), &
                                      FV_UPrim_Xiplus_loc(:,:),FV_UPrim_Ximinus_loc(:,:),                              &
                                      UPrim_extended(:,:,p,q,iElem),dUdUvol_plus(:,:,:,:),dUdUvol_minus(:,:,:,:))
#endif
    DO i=1,PP_N
      CALL Riemann_Point(F(:,i-1,p,q),                             &
                   FV_U_Xiplus_loc(     :,i-1),                    &
                   FV_U_Ximinus_loc(    :,i  ),                    &
                   FV_UPrim_Xiplus_loc( :,i-1),                    &
                   FV_UPrim_Ximinus_loc(:,i  ),                    &
                   FV_NormVecXi(        :,p  ,q,i,iElem),          &
                   FV_TangVec1Xi(       :,p  ,q,i,iElem),          &
                   FV_TangVec2Xi(       :,p  ,q,i,iElem),.FALSE.)
      DO jVar=1,PP_nVar
        ! modify U_plus at i-1
        reps0_L  = reps0*(1.+ABS(FV_U_plus_Tilde(jVar,i-1)))
        sreps0_L = 1./reps0_L
        FV_U_plus_Tilde(jVar,i-1) = FV_U_plus_Tilde(jVar,i-1) + reps0_L
        CALL ConsToPrim(FV_UPrim_plus_Tilde(:,i-1),FV_U_plus_Tilde(:,i-1))
        ! modify U_minus at i
        reps0_R  = reps0*(1.+ABS(FV_U_minus_Tilde(jVar,i)))
        sreps0_R = 1./reps0_R
        FV_U_minus_Tilde(jVar,i) = FV_U_minus_Tilde(jVar,i) + reps0_R
        CALL ConsToPrim(FV_UPrim_minus_Tilde(:,i),FV_U_minus_Tilde(:,i))
 
        CALL Riemann_Point(F_Tilde(:,i-1,p,q),                        &
                     FV_U_plus_Tilde(     :,i-1),                     &
                     FV_U_Ximinus_loc(    :,i  ),                     &
                     FV_UPrim_plus_Tilde( :,i-1),                     &
                     FV_UPrim_Ximinus_loc(:,i  ),                     &
                     FV_NormVecXi(        :,p  ,q,i,iElem),           &
                     FV_TangVec1Xi(       :,p  ,q,i,iElem),           &
                     FV_TangVec2Xi(       :,p  ,q,i,iElem),.FALSE.)
        CALL Riemann_Point(F_Tilde(:,i,p,q),                          &
                     FV_U_Xiplus_loc(     :,i-1),                     &
                     FV_U_minus_Tilde(    :,i  ),                     &
                     FV_UPrim_Xiplus_loc( :,i-1),                     &
                     FV_UPrim_minus_Tilde(:,i  ),                     &
                     FV_NormVecXi(        :,p  ,q,i,iElem),           &
                     FV_TangVec1Xi(       :,p  ,q,i,iElem),           &
                     FV_TangVec2Xi(       :,p  ,q,i,iElem),.FALSE.)
        DO iVar=1,PP_nVar
          dFdU_plus( iVar,jVar,i-1,p,q) =  FV_SurfElemXi_sw(p,q,i,iElem)*(F_Tilde(iVar,i-1,p,q) - F(iVar,i-1,p,q))*sreps0_L
          dFdU_plus( iVar,jVar,i  ,p,q) = -FV_SurfElemXi_sw(p,q,i,iElem)*(F_Tilde(iVar,i-1,p,q) - F(iVar,i-1,p,q))*sreps0_L
          dFdU_minus(iVar,jVar,i-1,p,q) =  FV_SurfElemXi_sw(p,q,i,iElem)*(F_Tilde(iVar,i  ,p,q) - F(iVar,i-1,p,q))*sreps0_R
          dFdU_minus(iVar,jVar,i  ,p,q) = -FV_SurfElemXi_sw(p,q,i,iElem)*(F_Tilde(iVar,i  ,p,q) - F(iVar,i-1,p,q))*sreps0_R
        END DO !iVar
        ! reset U_plus
        FV_U_plus_Tilde(jVar,i-1) = FV_U_Xiplus_loc(jVar,i-1)
        ! reset U_minus
        FV_U_minus_Tilde(jVar,i) = FV_U_Ximinus_loc(jVar,i)
      END DO !jVar
      ! assemble preconditioner
      r = vn1*p + vn2*q + PP_nVar*(i-1)
      s = vn1*p + vn2*q + PP_nVar*i
#if FV_RECONSTRUCT
      m = vn1*p + vn2*q + PP_nVar*(i-2)
      l = vn1*p + vn2*q + PP_nVar*(i+1)
      ! derivatives with respect to i-2
      IF(i.GT.1)THEN ! i-2 is not at neighboring element
        matrix=dFdU_plus( :,:,i-1,p,q); vector=dUdUvol_plus( :,:,i-1,1)
        BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) = BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) + MATMUL(matrix,vector)
        matrix=dFdU_plus( :,:,i  ,p,q); vector=dUdUvol_plus( :,:,i-1,1)
        BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) = BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) + MATMUL(matrix,vector)
      END IF
      ! derivatives with respect to i-1
      matrix=dFdU_plus( :,:,i-1,p,q); vector=dUdUvol_plus( :,:,i-1,2)
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,i  ,p,q); vector=dUdUvol_plus( :,:,i-1,2)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,i-1,p,q); vector=dUdUvol_minus(:,:,i  ,1)
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,i  ,p,q); vector=dUdUvol_minus(:,:,i  ,1)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      ! derivatives with respect to i
      matrix=dFdU_minus(:,:,i-1,p,q); vector=dUdUvol_minus(:,:,i  ,2)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,i  ,p,q); vector=dUdUvol_minus(:,:,i  ,2)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,i-1,p,q); vector=dUdUvol_plus( :,:,i-1,3)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,i  ,p,q); vector=dUdUvol_plus( :,:,i-1,3)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      ! derivatives with respect to i+1
      IF(i.LT.PP_N)THEN ! i+1 is not at neighboring element
      matrix=dFdU_minus(:,:,i-1,p,q); vector=dUdUvol_minus(:,:,i  ,3)
        BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) = BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) + MATMUL(matrix,vector)
        matrix=dFdU_minus(:,:,i  ,p,q); vector=dUdUvol_minus(:,:,i  ,3)
        BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) = BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) + MATMUL(matrix,vector)
      END IF
#else
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + dFdU_plus( :,:,i-1,p,q)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + dFdU_minus(:,:,i-1,p,q)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + dFdU_minus(:,:,i  ,p,q)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + dFdU_plus( :,:,i  ,p,q)
#endif
    END DO !i
  END DO !p
END DO !q

! === Eta-Direction ===============================================================================================================
#if FV_RECONSTRUCT
SideID_minus         = ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
SideID_plus          = ElemToSide(E2S_SIDE_ID,ETA_PLUS ,iElem)
#endif

DO p=0,PP_N
  DO q=0,PP_NZ
    DO j=0,PP_N
#if FV_RECONSTRUCT
      FV_UPrim_Etaplus_loc( :,j) = UPrim(:,p,j,q,iElem) + gradUeta(:,p,q,j,iElem) * FV_dx_ETA_R(p,q,j,iElem)
      FV_UPrim_Etaminus_loc(:,j) = UPrim(:,p,j,q,iElem) - gradUeta(:,p,q,j,iElem) * FV_dx_ETA_L(p,q,j,iElem)
#else
      FV_UPrim_Etaplus_loc( :,j) = UPrim(:,p,j,q,iElem)
      FV_UPrim_Etaminus_loc(:,j) = UPrim(:,p,j,q,iElem)
#endif
      CALL PrimToCons(FV_UPrim_Etaminus_loc(:,j),FV_U_Etaminus_loc(:,j))
      CALL PrimToCons(FV_UPrim_Etaplus_loc( :,j),FV_U_Etaplus_loc( :,j))
    END DO
    FV_UPrim_plus_Tilde  = FV_UPrim_Etaplus_loc( :,:)
    FV_UPrim_minus_Tilde = FV_UPrim_Etaminus_loc(:,:)
    FV_U_plus_Tilde      = FV_U_Etaplus_loc(     :,:)
    FV_U_minus_Tilde     = FV_U_Etaminus_loc(    :,:)
#if FV_RECONSTRUCT
    CALL FV_Reconstruction_Derivative(FV_sdx_ETA_extended(p,q,:,iElem),FV_dx_ETA_L(p,q,:,iElem),FV_dx_ETA_R(p,q,:,iElem), &
                                      FV_UPrim_Etaplus_loc(:,:),FV_UPrim_Etaminus_loc(:,:),                               &
                                      UPrim_extended(:,p,:,q,iElem),dUdUvol_plus(:,:,:,:),dUdUvol_minus(:,:,:,:))
#endif
    DO j=1,PP_N
      CALL Riemann_Point(F(:,p,j-1,q),                              &
                   FV_U_Etaplus_loc(     :,j-1),                    &
                   FV_U_Etaminus_loc(    :,j  ),                    &
                   FV_UPrim_Etaplus_loc( :,j-1),                    &
                   FV_UPrim_Etaminus_loc(:,j  ),                    &
                   FV_NormVecEta(        :,p,q  ,j,iElem),          &
                   FV_TangVec1Eta(       :,p,q  ,j,iElem),          &
                   FV_TangVec2Eta(       :,p,q  ,j,iElem),.FALSE.)
 
      DO jVar=1,PP_nVar
        ! modify U_plus at j-1
        reps0_L  = reps0*(1.+ABS(FV_U_plus_Tilde(jVar,j-1)))
        sreps0_L = 1./reps0_L
        FV_U_plus_Tilde(jVar,j-1) = FV_U_plus_Tilde(jVar,j-1) + reps0_L
        CALL ConsToPrim(FV_UPrim_plus_Tilde(:,j-1),FV_U_plus_Tilde(:,j-1))
        ! modify U_minus at j
        reps0_R  = reps0*(1.+ABS(FV_U_minus_Tilde(jVar,j)))
        sreps0_R = 1./reps0_R
        FV_U_minus_Tilde(jVar,j) = FV_U_minus_Tilde(jVar,j) + reps0_R
        CALL ConsToPrim(FV_UPrim_minus_Tilde(:,j),FV_U_minus_Tilde(:,j))
        CALL Riemann_Point(F_Tilde(:,p,j-1,q),                         &
                     FV_U_plus_Tilde(      :,j-1),                     &
                     FV_U_Etaminus_loc(    :,j  ),                     &
                     FV_UPrim_plus_Tilde(  :,j-1),                     &
                     FV_UPrim_Etaminus_loc(:,j  ),                     &
                     FV_NormVecEta(        :,p,q  ,j,iElem),           &
                     FV_TangVec1Eta(       :,p,q  ,j,iElem),           &
                     FV_TangVec2Eta(       :,p,q  ,j,iElem),.FALSE.)
        CALL Riemann_Point(F_Tilde(:,p,j,q),                           &
                     FV_U_Etaplus_loc(    :,j-1),                      &
                     FV_U_minus_Tilde(    :,j  ),                      &
                     FV_UPrim_Etaplus_loc(:,j-1),                      &
                     FV_UPrim_minus_Tilde(:,j  ),                      &
                     FV_NormVecEta(       :,p,q  ,j,iElem),            &
                     FV_TangVec1Eta(      :,p,q  ,j,iElem),            &
                     FV_TangVec2Eta(      :,p,q  ,j,iElem),.FALSE.)
        DO iVar=1,PP_nVar
          dFdU_plus( iVar,jVar,p,j-1,q) =  FV_SurfElemEta_sw(p,q,j,iElem)*(F_Tilde(iVar,p,j-1,q) - F(iVar,p,j-1,q))*sreps0_L
          dFdU_plus( iVar,jVar,p,j  ,q) = -FV_SurfElemEta_sw(p,q,j,iElem)*(F_Tilde(iVar,p,j-1,q) - F(iVar,p,j-1,q))*sreps0_L
          dFdU_minus(iVar,jVar,p,j-1,q) =  FV_SurfElemEta_sw(p,q,j,iElem)*(F_Tilde(iVar,p,j  ,q) - F(iVar,p,j-1,q))*sreps0_R
          dFdU_minus(iVar,jVar,p,j  ,q) = -FV_SurfElemEta_sw(p,q,j,iElem)*(F_Tilde(iVar,p,j  ,q) - F(iVar,p,j-1,q))*sreps0_R
        END DO !iVar
        ! reset U_plus
        FV_U_plus_Tilde(jVar,j-1) = FV_U_Etaplus_loc(jVar,j-1)
        ! reset U_minus
        FV_U_minus_Tilde(jVar,j) = FV_U_Etaminus_loc(jVar,j)
      END DO !jVar
      ! assemble preconditioner
      r = vn1*(j-1) + vn2*q + PP_nVar*p
      s = vn1*j     + vn2*q + PP_nVar*p
#if FV_RECONSTRUCT
      m = vn1*(j-2) + vn2*q + PP_nVar*p
      l = vn1*(j+1) + vn2*q + PP_nVar*p
      ! derivatives with respect to j-2
      IF(j.GT.1)THEN
        matrix=dFdU_plus( :,:,p,j-1,q); vector=dUdUvol_plus( :,:,j-1,1)
        BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) = BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) + MATMUL(matrix,vector)
        matrix=dFdU_plus( :,:,p,j  ,q); vector=dUdUvol_plus( :,:,j-1,1)
        BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) = BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) + MATMUL(matrix,vector)
      END IF
      ! derivatives with respect to j-1
      matrix=dFdU_plus( :,:,p,j-1,q); vector=dUdUvol_plus( :,:,j-1,2)
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,p,j  ,q); vector=dUdUvol_plus( :,:,j-1,2)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,p,j-1,q); vector=dUdUvol_minus(:,:,j  ,1)
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,p,j  ,q); vector=dUdUvol_minus(:,:,j  ,1)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      ! derivatives with respect to j
      matrix=dFdU_minus(:,:,p,j-1,q); vector=dUdUvol_minus(:,:,j  ,2)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,p,j  ,q); vector=dUdUvol_minus(:,:,j  ,2)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,p,j-1,q); vector=dUdUvol_plus( :,:,j-1,3)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,p,j  ,q); vector=dUdUvol_plus( :,:,j-1,3)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      ! derivatives with respect to j+1
      IF(j.LT.PP_N)THEN
        matrix=dFdU_minus(:,:,p,j-1,q); vector=dUdUvol_minus(:,:,j  ,3)
        BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) = BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) + MATMUL(matrix,vector)
        matrix=dFdU_minus(:,:,p,j  ,q); vector=dUdUvol_minus(:,:,j  ,3)
        BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) = BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) + MATMUL(matrix,vector)
      END IF
#else
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + dFdU_plus( :,:,p,j-1,q)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + dFdU_minus(:,:,p,j-1,q)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + dFdU_minus(:,:,p,j  ,q)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + dFdU_plus( :,:,p,j  ,q)
#endif
    END DO !j
  END DO !p
END DO !q

#if PP_dim == 3
! === Zeta-Direction ==============================================================================================================
#if FV_RECONSTRUCT
SideID_minus         = ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
SideID_plus          = ElemToSide(E2S_SIDE_ID,ZETA_PLUS ,iElem)
#endif

DO p=0,PP_N
  DO q=0,PP_N
    DO k=0,PP_N
#if FV_RECONSTRUCT
      FV_UPrim_Zetaplus_loc( :,k) = UPrim(:,p,q,k,iElem) + gradUzeta(:,p,q,k,iElem) * FV_dx_ZETA_R(p,q,k,iElem)
      FV_UPrim_Zetaminus_loc(:,k) = UPrim(:,p,q,k,iElem) - gradUzeta(:,p,q,k,iElem) * FV_dx_ZETA_L(p,q,k,iElem)
#else
      FV_UPrim_Zetaplus_loc( :,k) = UPrim(:,p,q,k,iElem)
      FV_UPrim_Zetaminus_loc(:,k) = UPrim(:,p,q,k,iElem)
#endif
      CALL PrimToCons(FV_UPrim_Zetaminus_loc(:,k),FV_U_Zetaminus_loc(:,k))
      CALL PrimToCons(FV_UPrim_Zetaplus_loc( :,k),FV_U_Zetaplus_loc( :,k))
    END DO
    FV_UPrim_plus_Tilde  = FV_UPrim_Zetaplus_loc( :,:)
    FV_UPrim_minus_Tilde = FV_UPrim_Zetaminus_loc(:,:)
    FV_U_plus_Tilde      = FV_U_Zetaplus_loc(     :,:)
    FV_U_minus_Tilde     = FV_U_Zetaminus_loc(    :,:)
#if FV_RECONSTRUCT
    CALL FV_Reconstruction_Derivative(FV_sdx_ZETA_extended(p,q,:,iElem),FV_dx_ZETA_L(p,q,:,iElem),FV_dx_ZETA_R(p,q,:,iElem), &
                                      FV_UPrim_Zetaplus_loc(:,:),FV_UPrim_Zetaminus_loc(:,:),                                &
                                      UPrim_extended(:,p,q,:,iElem),dUdUvol_plus(:,:,:,:),dUdUvol_minus(:,:,:,:))
#endif
    DO k=1,PP_N
      CALL Riemann_Point(F(:,p,q,k-1),                               &
                   FV_U_Zetaplus_loc(     :,k-1),                    &
                   FV_U_Zetaminus_loc(    :,k  ),                    &
                   FV_UPrim_Zetaplus_loc( :,k-1),                    &
                   FV_UPrim_Zetaminus_loc(:,k  ),                    &
                   FV_NormVecZeta(        :,p,q,k  ,iElem),          &
                   FV_TangVec1Zeta(       :,p,q,k  ,iElem),          &
                   FV_TangVec2Zeta(       :,p,q,k  ,iElem),.FALSE.)
 
      DO jVar=1,PP_nVar
        ! modify U_plus at k-1
        reps0_L  = reps0*(1.+ABS(FV_U_plus_Tilde(jVar,k-1)))
        sreps0_L = 1./reps0_L
        FV_U_plus_Tilde(jVar,k-1) = FV_U_plus_Tilde(jVar,k-1) + reps0_L
        CALL ConsToPrim(FV_UPrim_plus_Tilde(:,k-1),FV_U_plus_Tilde(:,k-1))
        ! modify U_minus at k
        reps0_R  = reps0*(1.+ABS(FV_U_minus_Tilde(jVar,k)))
        sreps0_R = 1./reps0_R
        FV_U_minus_Tilde(jVar,k) = FV_U_minus_Tilde(jVar,k) + reps0_R
        CALL ConsToPrim(FV_UPrim_minus_Tilde(:,k),FV_U_minus_Tilde(:,k))
 
        CALL Riemann_Point(F_Tilde(:,p,q,k-1),                          &
                     FV_U_plus_Tilde(       :,k-1),                     &
                     FV_U_Zetaminus_loc(    :,k),                       &
                     FV_UPrim_plus_Tilde(   :,k-1),                     &
                     FV_UPrim_Zetaminus_loc(:,k),                       &
                     FV_NormVecZeta(        :,p,q,k  ,iElem),           &
                     FV_TangVec1Zeta(       :,p,q,k  ,iElem),           &
                     FV_TangVec2Zeta(       :,p,q,k  ,iElem),.FALSE.)
        CALL Riemann_Point(F_Tilde(:,p,q,k),                            &
                     FV_U_Zetaplus_loc(    :,k-1),                      &
                     FV_U_minus_Tilde(     :,k),                        &
                     FV_UPrim_Zetaplus_loc(:,k-1),                      &
                     FV_UPrim_minus_Tilde( :,k),                        &
                     FV_NormVecZeta(       :,p,q,k  ,iElem),            &
                     FV_TangVec1Zeta(      :,p,q,k  ,iElem),            &
                     FV_TangVec2Zeta(      :,p,q,k  ,iElem),.FALSE.)
        DO iVar=1,PP_nVar
          dFdU_plus( iVar,jVar,p,q,k-1) =  FV_SurfElemZeta_sw(p,q,k,iElem)*(F_Tilde(iVar,p,q,k-1) - F(iVar,p,q,k-1))*sreps0_L
          dFdU_plus( iVar,jVar,p,q,k  ) = -FV_SurfElemZeta_sw(p,q,k,iElem)*(F_Tilde(iVar,p,q,k-1) - F(iVar,p,q,k-1))*sreps0_L
          dFdU_minus(iVar,jVar,p,q,k-1) =  FV_SurfElemZeta_sw(p,q,k,iElem)*(F_Tilde(iVar,p,q,k  ) - F(iVar,p,q,k-1))*sreps0_R
          dFdU_minus(iVar,jVar,p,q,k  ) = -FV_SurfElemZeta_sw(p,q,k,iElem)*(F_Tilde(iVar,p,q,k  ) - F(iVar,p,q,k-1))*sreps0_R
        END DO !iVar
        ! reset U_plus
        FV_U_plus_Tilde(jVar,k-1) = FV_U_Zetaplus_loc(jVar,k-1)
        ! reset U_minus
        FV_U_minus_Tilde(jVar,k) = FV_U_Zetaminus_loc(jVar,k)
      END DO !jVar
      ! assemble preconditioner
      r = vn1*q + vn2*(k-1) + PP_nVar*p
      s = vn1*q + vn2*k     + PP_nVar*p
#if FV_RECONSTRUCT
      m = vn1*q + vn2*(k-2) + PP_nVar*p
      l = vn1*q + vn2*(k+1) + PP_nVar*p
      ! derivatives with respect to k-2
      IF(k.GT.1)THEN
        matrix=dFdU_plus( :,:,p,q,k-1); vector=dUdUvol_plus( :,:,k-1,1)
        BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) = BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) + MATMUL(matrix,vector)
        matrix=dFdU_plus( :,:,p,q,k  ); vector=dUdUvol_plus( :,:,k-1,1)
        BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) = BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) + MATMUL(matrix,vector)
      END IF
      ! derivatives with respect to k-1
      matrix=dFdU_plus( :,:,p,q,k-1); vector=dUdUvol_plus( :,:,k-1,2)
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,p,q,k  ); vector=dUdUvol_plus( :,:,k-1,2)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,p,q,k-1); vector=dUdUvol_minus(:,:,k  ,1)
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,p,q,k  ); vector=dUdUvol_minus(:,:,k  ,1)
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + MATMUL(matrix,vector)
      ! derivatives with respect to k
      matrix=dFdU_minus(:,:,p,q,k-1); vector=dUdUvol_minus(:,:,k  ,2)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_minus(:,:,p,q,k  ); vector=dUdUvol_minus(:,:,k  ,2)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,p,q,k-1); vector=dUdUvol_plus( :,:,k-1,3)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      matrix=dFdU_plus( :,:,p,q,k  ); vector=dUdUvol_plus( :,:,k-1,3)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + MATMUL(matrix,vector)
      ! derivatives with respect to k+1
      IF(k.LT.PP_N)THEN
        matrix=dFdU_minus(:,:,p,q,k-1); vector=dUdUvol_minus(:,:,k  ,3)
        BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) = BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) + MATMUL(matrix,vector)
        matrix=dFdU_minus(:,:,p,q,k  ); vector=dUdUvol_minus(:,:,k  ,3)
        BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) = BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) + MATMUL(matrix,vector)
      END IF
#else
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + dFdU_plus( :,:,p,q,k-1)
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + dFdU_minus(:,:,p,q,k-1)
      BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + dFdU_minus(:,:,p,q,k  )
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) + dFdU_plus( :,:,p,q,k  )
#endif
    END DO !k
  END DO !p
END DO !q
#endif

END SUBROUTINE FVVolIntJac
#endif /* FV_ENABLED */


SUBROUTINE FinalizeJac_Ex()
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_Jac_ex_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(LL_minus)
SDEALLOCATE(LL_plus)
SDEALLOCATE(l_mp)
SDEALLOCATE(LL_mp)
SDEALLOCATE(nVec)
SDEALLOCATE(Surf)
#if PARABOLIC
SDEALLOCATE(R_Minus)
SDEALLOCATE(R_Plus)
SDEALLOCATE(JacLiftingFlux)
#endif /*PARABOLIC*/
#if FV_ENABLED && FV_RECONSTRUCT
SDEALLOCATE(UPrim_extended)
SDEALLOCATE(FV_sdx_XI_extended)
SDEALLOCATE(FV_sdx_ETA_extended)
#if PP_dim == 3
SDEALLOCATE(FV_sdx_ZETA_extended)
#endif
#endif
END SUBROUTINE FinalizeJac_Ex

END MODULE MOD_Jac_Ex
