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
USE MOD_Jac_br2            ,ONLY: BuildnVecTangSurf,Build_BR2_SurfTerms
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
  ALLOCATE(JacLiftingFlux(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_NZ,6))
  CALL Build_BR2_SurfTerms()
END IF
#endif /*PARABOLIC*/

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT BLOCK JACOBIAN DONE!'
END SUBROUTINE InitJac_Ex

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
USE MOD_Jac_br2                   ,ONLY:FillJacLiftingFlux
#endif
#if FV_ENABLED
USE MOD_FV_Vars                   ,ONLY:FV_Elems
#endif
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
  CALL FillJacLiftingFlux(t,iElem)
END IF
#endif /*PARABOLIC*/
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
    IF(FVEM.EQ.0)THEN
      CALL DGVolIntGradJac(BJ,iElem)               !d(F^v)/dQ
#if FV_ENABLED
    ELSE
      CALL FVVolIntGradJac(BJ,iElem)               !d(F^v)/dQ
#endif
    END IF
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
#if EQNSYSNR==2 && PARABOLIC
USE MOD_Precond_Vars  ,ONLy: EulerPrecond
USE MOD_Lifting_Vars  ,ONLY: gradUx,gradUy,gradUz
USE MOD_Jacobian      ,ONLY: EvalDiffFluxJacobian 
USE MOD_DG_Vars       ,ONLY: nDOFElem
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif /*EDDYVISCOSITY*/
#endif /*EQNSYSNR && PARABOLIC*/
USE MOD_Mesh_Vars     ,ONLY: Metrics_fTilde,Metrics_gTilde
#if PP_dim==3
USE MOD_Mesh_Vars     ,ONLY: Metrics_hTilde
#endif
USE MOD_Implicit_Vars ,ONLY: nDOFVarElem
#ifdef SPLIT_DG
USE MOD_DG_Vars       ,ONLY: DVolSurf
USE MOD_Jac_Split     ,ONLY: Jac_Split
#if EQNSYSNR==2 && PARABOLIC
USE MOD_DG_Vars       ,ONLY: D_hat
#endif
#else
USE MOD_DG_Vars       ,ONLY: D_hat,U,UPrim
#endif
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
INTEGER                                                 :: mm,nn,oo
INTEGER                                                 :: s,r1,r2,ll,vn1,vn2
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJacTilde,gJacTilde
#if PP_dim==3
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: hJacTilde
#endif
#ifdef SPLIT_DG
INTEGER                                                 :: tt
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJacVisc,gJacVisc
#if PP_dim==3
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: hJacVisc
#endif
#else
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: fJac,gJac,hJac
#endif
#if PP_dim==3
INTEGER                                                 :: r3
#endif
#if EQNSYSNR==2 && PARABOLIC
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: fJac_loc,gJac_loc,hJac_loc
#endif /*EQNSYSNR && PARABOLIC*/
!===================================================================================================================================
vn1=PP_nVar*(PP_N+1)
vn2=vn1*(PP_N+1)
#ifndef SPLIT_DG
! dF^hyp/DU
CALL EvalAdvFluxJacobian(U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem),fJac,gJac,hJac)
#endif

#if EQNSYSNR==2 && PARABOLIC
!No CALL of EvalDiffFlux Jacobian for EQNSYSNR==1, since the diffusive flux is not depending on U
IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
  ! dF^visc/DU
  CALL EvalDiffFluxJacobian(nDOFElem,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem) &
                            ,gradUx(:,:,:,:,iElem) &
                            ,gradUy(:,:,:,:,iElem) &
                            ,gradUz(:,:,:,:,iElem) &
                            ,fJac_loc,gJac_loc,hJac_loc &
#if EDDYVISCOSITY
                            ,muSGS(:,:,:,:,iElem)  &
#endif
                            )
#ifndef SPLIT_DG
  fJac=fJac+fJac_loc                                         
  gJac=gJac+gJac_loc                                         
#if PP_dim==3
  hJac=hJac+hJac_loc                                         
#endif/*PP_dim*/
#endif /*SPLIT_DG*/
END IF !EulerPrecond
#endif /*EQNSYSNR && PARABOLIC*/

s=0
DO oo=0,PP_NZ
  DO nn=0,PP_N
    DO mm=0,PP_N
#ifdef SPLIT_DG

#if EQNSYSNR==2 && PARABOLIC
      IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
        ! weak formulation for parabolic terms
        fJacVisc(:,:) = (fJac_loc(:,:,mm,nn,oo)*Metrics_fTilde(1,mm,nn,oo,iElem,0)  &
                        +gJac_loc(:,:,mm,nn,oo)*Metrics_fTilde(2,mm,nn,oo,iElem,0)  &
#if PP_dim==3
                        +hJac_loc(:,:,mm,nn,oo)*Metrics_fTilde(3,mm,nn,oo,iElem,0)  &
#endif
                         )

        gJacVisc(:,:) = (fJac_loc(:,:,mm,nn,oo)*Metrics_gTilde(1,mm,nn,oo,iElem,0)  &
                        +gJac_loc(:,:,mm,nn,oo)*Metrics_gTilde(2,mm,nn,oo,iElem,0)  &
#if PP_dim==3
                        +hJac_loc(:,:,mm,nn,oo)*Metrics_gTilde(3,mm,nn,oo,iElem,0)  &
#endif
                         )

#if PP_dim==3
        hJacVisc(:,:) =  fJac_loc(:,:,mm,nn,oo)*Metrics_hTilde(1,mm,nn,oo,iElem,0)   &
                        +gJac_loc(:,:,mm,nn,oo)*Metrics_hTilde(2,mm,nn,oo,iElem,0)  &
                        +hJac_loc(:,:,mm,nn,oo)*Metrics_hTilde(3,mm,nn,oo,iElem,0) 
#endif
      END IF
#endif /*EQNSYSNR && PARABOLIC*/


      ! strong split form (modified surface terms!)
      r1=           vn1*nn+vn2*oo
      r2=mm*PP_nVar       +vn2*oo
#if PP_dim==3
      r3=mm*PP_nVar+vn1*nn
#endif
      ! Case i!=m (i is the loop variable, similar to weak form)
      DO ll=0,PP_N
        
#if EQNSYSNR==2 && PARABOLIC
        IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
          BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,mm)*fJacVisc(:,:) 
          BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,nn)*gJacVisc(:,:) 
#if PP_dim==3
          BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,oo)*hJacVisc(:,:) 
#endif
        END IF
#endif /*EQNSYSNR && PARABOLIC*/

        CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,ll,nn,oo,iElem),UPrim(:,ll,nn,oo,iElem), &
                       Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,ll,nn,oo,iElem,0),fJacTilde(:,:))
        BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + DVolSurf(mm,ll)*fJacTilde(:,:) 

        CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,ll,oo,iElem),UPrim(:,mm,ll,oo,iElem), &
                       Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,ll,oo,iElem,0),gJacTilde(:,:))
        BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + DVolSurf(nn,ll)*gJacTilde(:,:) 

#if PP_dim==3
        CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,nn,ll,iElem),UPrim(:,mm,nn,ll,iElem), &
                       Metrics_hTilde(:,mm,nn,oo,iElem,0),Metrics_hTilde(:,mm,nn,ll,iElem,0),hJacTilde(:,:))
        BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) + DVolSurf(oo,ll)*hJacTilde(:,:) 
#endif

        r1=r1+PP_nVar
        r2=r2+vn1
#if PP_dim==3
        r3=r3+vn2
#endif
      END DO !ll
      ! Case i=m (influence on main diagonal of block jacobian), additional entries from split formulation
      DO tt=0,PP_N
        CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,tt,nn,oo,iElem),UPrim(:,tt,nn,oo,iElem), &
                       Metrics_fTilde(:,mm,nn,oo,iElem,0),Metrics_fTilde(:,tt,nn,oo,iElem,0),fJacTilde(:,:))
        BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + DVolSurf(tt,mm)*fJacTilde(:,:)

        CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,tt,oo,iElem),UPrim(:,mm,tt,oo,iElem), &
                       Metrics_gTilde(:,mm,nn,oo,iElem,0),Metrics_gTilde(:,mm,tt,oo,iElem,0),gJacTilde(:,:))
        BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + DVolSurf(tt,nn)*gJacTilde(:,:)

#if PP_dim==3
        CALL Jac_Split(U(:,mm,nn,oo,iElem),UPrim(:,mm,nn,oo,iElem),U(:,mm,nn,tt,iElem),UPrim(:,mm,nn,tt,iElem), &
                       Metrics_hTilde(:,mm,nn,oo,iElem,0),Metrics_hTilde(:,mm,nn,tt,iElem,0),hJacTilde(:,:))
        BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) + DVolSurf(tt,oo)*hJacTilde(:,:)
#endif
      END DO
#else /*SPLIT_DG*/

      ! weak formulation
      fJacTilde(:,:) = (fJac(:,:,mm,nn,oo)*Metrics_fTilde(1,mm,nn,oo,iElem,0)  &
                       +gJac(:,:,mm,nn,oo)*Metrics_fTilde(2,mm,nn,oo,iElem,0)  &
#if PP_dim==3
                       +hJac(:,:,mm,nn,oo)*Metrics_fTilde(3,mm,nn,oo,iElem,0)  &
#endif
                       )

      gJacTilde(:,:) = (fJac(:,:,mm,nn,oo)*Metrics_gTilde(1,mm,nn,oo,iElem,0)  &
                       +gJac(:,:,mm,nn,oo)*Metrics_gTilde(2,mm,nn,oo,iElem,0)  &
#if PP_dim==3
                       +hJac(:,:,mm,nn,oo)*Metrics_gTilde(3,mm,nn,oo,iElem,0)  &
#endif
                       )

#if PP_dim==3
      hJacTilde(:,:) = fJac(:,:,mm,nn,oo)*Metrics_hTilde(1,mm,nn,oo,iElem,0)   &
                       +gJac(:,:,mm,nn,oo)*Metrics_hTilde(2,mm,nn,oo,iElem,0)  &
                       +hJac(:,:,mm,nn,oo)*Metrics_hTilde(3,mm,nn,oo,iElem,0) 
#endif
      r1=           vn1*nn+vn2*oo
      r2=mm*PP_nVar       +vn2*oo
#if PP_dim==3
      r3=mm*PP_nVar+vn1*nn
#endif
      DO ll=0,PP_N
        BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,mm)*fJacTilde(:,:) 
        BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,nn)*gJacTilde(:,:) 
#if PP_dim==3
        BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) + D_hat(ll,oo)*hJacTilde(:,:) 
#endif
        r1=r1+PP_nVar
        r2=r2+vn1
#if PP_dim==3
        r3=r3+vn2
#endif
      END DO !ll
#endif /*SPLIT_DG*/
      s=s+PP_nVar
    END DO !mm
  END DO !nn
END DO !oo 

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
!> volume integral: the total derivative of the viscous flux with resprect to U:
!>                  dF^v/DU_cons = dF^v/dQ_prim* DQ_prim/DU_prim* DU_prim/DU_cons + dF^v/DU_cons
!>                                       |              |                |              |
!>                              FluxGradJacobian     Lifting      dPrimTempdCons  (already done in DGVolIntJac) 
!===================================================================================================================================
SUBROUTINE  DGVolIntGradJac(BJ,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Jac_br2       ,ONLY: JacLifting_VolInt
USE MOD_Precond_Vars  ,ONLY: NoFillIn
USE MOD_DG_Vars       ,ONLY: D_hat,U,UPrim,nDOFElem
USE MOD_Mesh_Vars     ,ONLY: Metrics_fTilde,Metrics_gTilde
#if PP_dim==3
USE MOD_Mesh_Vars     ,ONLY: Metrics_hTilde
#endif
USE MOD_Implicit_Vars ,ONLY: nDOFVarElem
USE MOD_GradJacobian  ,ONLY: EvalFluxGradJacobian
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars ,ONLY: muSGS
#endif /*EDDYVISCOSITY*/
USE MOD_Jacobian      ,ONLY: dPrimTempdCons
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
INTEGER                                                                    :: i,j,k,l,mm,nn,oo
INTEGER                                                                    :: r,s,vn1,vn2
REAL,DIMENSION(PP_nVarPrim,PP_nVar)                                        :: PrimConsJac
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)              :: fJacQx,gJacQx,hJacQx
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)              :: fJacQy,gJacQy,hJacQy
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)              :: fJacQz,gJacQz,hJacQz
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: fJacTilde,gJacTilde
REAL,DIMENSION(PP_nVar    ,PP_nVar)                                        :: fJac,gJac
#if PP_dim==3
REAL,DIMENSION(PP_nVar    ,PP_nVar)                                        :: hJac
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: hJacTilde
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,0:PP_N,3)      :: JacLifting_Z
#endif
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim):: JacLifting_X,JacLifting_Y
!===================================================================================================================================
vn1=PP_nVar*(PP_N+1)
vn2=vn1*(PP_N+1)
! Calculation of the derivative of the diffusive flux with respect to gradU
! dF/dQ
CALL EvalFluxGradJacobian(nDOFElem,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem) &
                          ,fJacQx,fJacQy,fJacQz &
                          ,gJacQx,gJacQy,gJacQz &
                          ,hJacQx,hJacQy,hJacQz &
#if EDDYVISCOSITY
                          ,muSGS(:,:,:,:,iElem) &
#endif
                         )

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
                          +gJacTilde(:,:)*Metrics_fTilde(2,i,j,k,iElem,0) &
#if PP_dim==3
                          +hJacTilde(:,:)*Metrics_fTilde(3,i,j,k,iElem,0) &
#endif
                          )
      gJacQx(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_gTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_gTilde(2,i,j,k,iElem,0) &
#if PP_dim==3
                          +hJacTilde(:,:)*Metrics_gTilde(3,i,j,k,iElem,0) &
#endif
                          )
#if PP_dim==3
      hJacQx(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_hTilde(1,i,j,k,iElem,0) &
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
      hJacQy(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_hTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_hTilde(2,i,j,k,iElem,0) &
                          +hJacTilde(:,:)*Metrics_hTilde(3,i,j,k,iElem,0) &
                          ) 
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
      hJacQz(:,:,i,j,k) = (fJacTilde(:,:)*Metrics_hTilde(1,i,j,k,iElem,0) &
                          +gJacTilde(:,:)*Metrics_hTilde(2,i,j,k,iElem,0) &
                          +hJacTilde(:,:)*Metrics_hTilde(3,i,j,k,iElem,0) &
                          )
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
      CALL dPrimTempdCons(UPrim(:,mm,nn,oo,iElem),PrimConsJac(:,:))
      IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
        DO j=0,PP_N
          fJacTilde(:,:)= ( MATMUL(fJacQx(:,:,mm,j,oo) , JacLifting_X(:,:,mm,j,oo,nn,2) )  &
                          + MATMUL(fJacQy(:,:,mm,j,oo) , JacLifting_Y(:,:,mm,j,oo,nn,2) )  &
#if PP_dim==3
                          + MATMUL(fJacQz(:,:,mm,j,oo) , JacLifting_Z(:,:,mm,j,oo,nn,2) )  &
#endif
                          )
          fJac(:,:) = MATMUL(fJacTilde(:,:),PrimConsJac(:,:))
          DO i=0,PP_N
            r=PP_nVar*i+vn1*j+vn2*oo
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(i,mm)*fJac(:,:)
          END DO !i
        END DO !j
        DO i=0,PP_N
          gJacTilde(:,:)= ( MATMUL(gJacQx(:,:,i,nn,oo) , JacLifting_X(:,:,i,nn,oo,mm,1) )  &
                          + MATMUL(gJacQy(:,:,i,nn,oo) , JacLifting_Y(:,:,i,nn,oo,mm,1) )  &
#if PP_dim==3
                          + MATMUL(gJacQz(:,:,i,nn,oo) , JacLifting_Z(:,:,i,nn,oo,mm,1) )  &
#endif
                          )
          gJac(:,:) = MATMUL(gJacTilde(:,:),PrimConsJac(:,:))
          DO j=0,PP_N
            r=PP_nVar*i+vn1*j+vn2*oo
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(j,nn)*gJac(:,:)
          END DO !j
        END DO !i
#if PP_dim==3
        DO k=0,PP_N
          fJacTilde(:,:)= ( MATMUL(fJacQx(:,:,mm,nn,k) , JacLifting_X(:,:,mm,nn,k,oo,3) )  &
                          + MATMUL(fJacQy(:,:,mm,nn,k) , JacLifting_Y(:,:,mm,nn,k,oo,3) )  &
                          + MATMUL(fJacQz(:,:,mm,nn,k) , JacLifting_Z(:,:,mm,nn,k,oo,3) ))
          fJac(:,:) = MATMUL(fJacTilde(:,:),PrimConsJac(:,:))
          DO i=0,PP_N
            r=PP_nVar*i+vn1*nn+vn2*k
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(i,mm)*fJac(:,:)
          END DO !i
        END DO !k 
        DO i=0,PP_N
          hJacTilde(:,:)= ( MATMUL(hJacQx(:,:,i,nn,oo) , JacLifting_X(:,:,i,nn,oo,mm,1) )  &
                          + MATMUL(hJacQy(:,:,i,nn,oo) , JacLifting_Y(:,:,i,nn,oo,mm,1) )  &
                          + MATMUL(hJacQz(:,:,i,nn,oo) , JacLifting_Z(:,:,i,nn,oo,mm,1) )) 
          hJac(:,:) = MATMUL(hJacTilde(:,:),PrimConsJac(:,:))
          DO k=0,PP_N
            r=PP_nVar*i+vn1*nn+vn2*k
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(k,oo)*hJac(:,:)
          END DO !k
        END DO !i 
        DO k=0,PP_NZ
          gJacTilde(:,:)= ( MATMUL(gJacQx(:,:,mm,nn,k) , JacLifting_X(:,:,mm,nn,k,oo,3) )  &
                           +MATMUL(gJacQy(:,:,mm,nn,k) , JacLifting_Y(:,:,mm,nn,k,oo,3) )  &
                           +MATMUL(gJacQz(:,:,mm,nn,k) , JacLifting_Z(:,:,mm,nn,k,oo,3) ))
          gJac(:,:) = MATMUL(gJacTilde(:,:),PrimConsJac(:,:))
          DO j=0,PP_N
            r=PP_nVar*mm+vn1*j+vn2*k
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(j,nn)*gJac(:,:)
          END DO !j
        END DO !k 
        DO j=0,PP_N
          hJacTilde(:,:)= ( MATMUL(hJacQx(:,:,mm,j,oo) , JacLifting_X(:,:,mm,j,oo,nn,2) )  &
                          + MATMUL(hJacQy(:,:,mm,j,oo) , JacLifting_Y(:,:,mm,j,oo,nn,2) )  &
                          + MATMUL(hJacQz(:,:,mm,j,oo) , JacLifting_Z(:,:,mm,j,oo,nn,2) ))
          hJac(:,:) = MATMUL(hJacTilde(:,:),PrimConsJac(:,:))
          DO k=0,PP_NZ
            r=PP_nVar*mm+vn1*j+vn2*k
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(k,oo)*hJac(:,:)
          END DO !k
        END DO !j 
#endif
      END IF !NoFillIn
      DO l=0,PP_N
        fJacTilde(:,:)= ( MATMUL(fJacQx(:,:,l,nn,oo) , JacLifting_X(:,:,l,nn,oo,mm,1) )  &
                        + MATMUL(fJacQy(:,:,l,nn,oo) , JacLifting_Y(:,:,l,nn,oo,mm,1) )  &
#if PP_dim==3
                        + MATMUL(fJacQz(:,:,l,nn,oo) , JacLifting_Z(:,:,l,nn,oo,mm,1) )  &
#endif
                        )
        fJac(:,:) = MATMUL(fJacTilde(:,:),PrimConsJac(:,:))
        DO i=0,PP_N
          r=PP_nVar*i+vn1*nn+vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(i,l)*fJac(:,:)
        END DO !i
        gJacTilde(:,:)= ( MATMUL(gJacQx(:,:,mm,l,oo) , JacLifting_X(:,:,mm,l,oo,nn,2) )  &
                        + MATMUL(gJacQy(:,:,mm,l,oo) , JacLifting_Y(:,:,mm,l,oo,nn,2) )  &
#if PP_dim==3
                        + MATMUL(gJacQz(:,:,mm,l,oo) , JacLifting_Z(:,:,mm,l,oo,nn,2) )  &
#endif
                       )
        gJac(:,:) = MATMUL(gJacTilde(:,:),PrimConsJac(:,:))
        DO j=0,PP_N
          r=PP_nVar*mm+vn1*j+vn2*oo
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + D_hat(j,l)*gJac(:,:)
        END DO !j
#if PP_dim==3
        hJacTilde(:,:)= ( MATMUL(hJacQx(:,:,mm,nn,l) , JacLifting_X(:,:,mm,nn,l,oo,3) )  &
                        + MATMUL(hJacQy(:,:,mm,nn,l) , JacLifting_Y(:,:,mm,nn,l,oo,3) )  &
                        + MATMUL(hJacQz(:,:,mm,nn,l) , JacLifting_Z(:,:,mm,nn,l,oo,3) ))
        hJac(:,:) = MATMUL(hJacTilde(:,:),PrimConsJac(:,:))
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

END SUBROUTINE DGVolIntGradJac
#endif

#if FV_ENABLED
!===================================================================================================================================
!> Volume integral Jacobian of the convective flux for FV elements.
!> Each finite volume subcell requires the solution of a Riemann problem. Hence, the derivation is done using a finite difference
!> for the derivation of the Riemann problem (same is always done for the surface fluxes). It also takes into account the additional
!> dependencies caused by the 2nd order reconstruction procedure. Here, a primitive reconstruction is assumed resulting in a back
!> and forth transformation from conservative to primitive. The additional dependencies are stored in dUdUvol_minus/plus and then
!> multiplied with the finite difference of the Riemann flux. Finally the derivatives are assembled at the right position of
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
#if EQNSYSNR==2 && PARABOLIC
USE MOD_Precond_Vars        ,ONLy: EulerPrecond
USE MOD_Lifting_Vars        ,ONLY: gradUx,gradUy,gradUz
USE MOD_Jacobian            ,ONLY: EvalDiffFluxJacobian 
USE MOD_DG_Vars             ,ONLY: U,nDOFElem
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars       ,ONLY: muSGS
#endif /*EDDYVISCOSITY*/
#endif /*EQNSYSNR && PARABOLIC*/
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
#if EQNSYSNR==2 && PARABOLIC
REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)    :: fJac_visc,gJac_visc,hJac_visc
REAL,DIMENSION(PP_nVar,PP_nVar)                          :: Jac_Visc_plus,Jac_Visc_minus
#endif /*EQNSYSNR && PARABOLIC*/
!===================================================================================================================================
vn1 = PP_nVar*(PP_N+1)
vn2 = vn1*(PP_N+1)
dfDU_plus  = 0
dFdU_minus = 0

#if EQNSYSNR==2 && PARABOLIC
!No CALL of EvalDiffFlux Jacobian for EQNSYSNR==1, since the diffusive flux is not depending on U
IF(EulerPrecond.EQV..FALSE.) THEN !Euler Precond = False
  ! dF^visc/DU
  CALL EvalDiffFluxJacobian(nDOFElem,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem) &
                            ,gradUx(:,:,:,:,iElem) &
                            ,gradUy(:,:,:,:,iElem) &
                            ,gradUz(:,:,:,:,iElem) &
                            ,fJac_visc,gJac_visc,hJac_visc &
#if EDDYVISCOSITY
                            ,muSGS(:,:,:,:,iElem)  &
#endif
                            )
END IF !EulerPrecond
#endif /*EQNSYSNR && PARABOLIC*/

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
#if EQNSYSNR==2 && PARABOLIC
      IF(EulerPrecond.EQV..FALSE.) THEN
        Jac_Visc_plus  = 0.5*(FV_NormVecXi(1,p,q,i,iElem)*fJac_visc(:,:,i-1,p,q) + &
                              FV_NormVecXi(2,p,q,i,iElem)*gJac_visc(:,:,i-1,p,q) + &
                              FV_NormVecXi(3,p,q,i,iElem)*hJac_visc(:,:,i-1,p,q))
        Jac_Visc_minus = 0.5*(FV_NormVecXi(1,p,q,i,iElem)*fJac_visc(:,:,i  ,p,q) + &
                              FV_NormVecXi(2,p,q,i,iElem)*gJac_visc(:,:,i  ,p,q) + &
                              FV_NormVecXi(3,p,q,i,iElem)*hJac_visc(:,:,i  ,p,q))

        BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_SurfElemXi_sw(p,q,i,iElem) * Jac_Visc_plus
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_SurfElemXi_sw(p,q,i,iElem) * Jac_Visc_minus
        BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) - FV_SurfElemXi_sw(p,q,i,iElem) * Jac_Visc_minus
        BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) - FV_SurfElemXi_sw(p,q,i,iElem) * Jac_Visc_plus
      END IF
#endif /*EQNSYSNR && PARABOLIC*/
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
#if EQNSYSNR==2 && PARABOLIC
      IF(EulerPrecond.EQV..FALSE.) THEN
        Jac_Visc_plus  = 0.5*(FV_NormVecEta(1,p,q,j,iElem)*fJac_visc(:,:,p,j-1,q) + &
                              FV_NormVecEta(2,p,q,j,iElem)*gJac_visc(:,:,p,j-1,q) + &
                              FV_NormVecEta(3,p,q,j,iElem)*hJac_visc(:,:,p,j-1,q))
        Jac_Visc_minus = 0.5*(FV_NormVecEta(1,p,q,j,iElem)*fJac_visc(:,:,p,j  ,q) + &
                              FV_NormVecEta(2,p,q,j,iElem)*gJac_visc(:,:,p,j  ,q) + &
                              FV_NormVecEta(3,p,q,j,iElem)*hJac_visc(:,:,p,j  ,q))

        BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_SurfElemEta_sw(p,q,j,iElem) * Jac_Visc_plus
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_SurfElemEta_sw(p,q,j,iElem) * Jac_Visc_minus
        BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) - FV_SurfElemEta_sw(p,q,j,iElem) * Jac_Visc_minus
        BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) - FV_SurfElemEta_sw(p,q,j,iElem) * Jac_Visc_plus
      END IF
#endif /*EQNSYSNR && PARABOLIC*/
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
#if EQNSYSNR==2 && PARABOLIC
      IF(EulerPrecond.EQV..FALSE.) THEN
        Jac_Visc_plus  = 0.5*(FV_NormVecZeta(1,p,q,k,iElem)*fJac_visc(:,:,p,q,k-1) + &
                              FV_NormVecZeta(2,p,q,k,iElem)*gJac_visc(:,:,p,q,k-1) + &
                              FV_NormVecZeta(3,p,q,k,iElem)*hJac_visc(:,:,p,q,k-1))
        Jac_Visc_minus = 0.5*(FV_NormVecZeta(1,p,q,k,iElem)*fJac_visc(:,:,p,q,k  ) + &
                              FV_NormVecZeta(2,p,q,k,iElem)*gJac_visc(:,:,p,q,k  ) + &
                              FV_NormVecZeta(3,p,q,k,iElem)*hJac_visc(:,:,p,q,k  ))

        BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_SurfElemZeta_sw(p,q,k,iElem) * Jac_Visc_plus
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_SurfElemZeta_sw(p,q,k,iElem) * Jac_Visc_minus
        BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) - FV_SurfElemZeta_sw(p,q,k,iElem) * Jac_Visc_minus
        BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) - FV_SurfElemZeta_sw(p,q,k,iElem) * Jac_Visc_plus
      END IF
#endif /*EQNSYSNR && PARABOLIC*/
    END DO !k
  END DO !p
END DO !q
#endif

END SUBROUTINE FVVolIntJac

#if PARABOLIC
!===================================================================================================================================
!> volume integral: the total derivative of the viscous flux with resprect to U:
!>                  dF^v/DU_cons = dF^v/dQ_prim* DQ_prim/DU_prim* DU_prim/DU_cons + dF^v/DU_cons
!>                                       |              |                |              |
!>                           FluxGradJacobian  FV-Reconstruction    dPrimTempdCons  (already done in FVVolIntJac) 
!===================================================================================================================================
SUBROUTINE  FVVolIntGradJac(BJ,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Jac_br2       ,ONLY: JacLifting_VolInt
USE MOD_Precond_Vars  ,ONLY: NoFillIn
USE MOD_DG_Vars       ,ONLY: U,UPrim,nDOFElem
USE MOD_Implicit_Vars ,ONLY: nDOFVarElem
USE MOD_GradJacobian  ,ONLY: EvalFluxGradJacobian
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars ,ONLY: muSGS
#endif /*EDDYVISCOSITY*/
USE MOD_Jacobian      ,ONLY: dPrimTempdCons
USE MOD_FV_Vars       ,ONLY: FV_NormVecXi,FV_NormVecEta,FV_SurfElemXi_sw,FV_SurfElemEta_sw
#if PP_dim==3
USE MOD_FV_Vars       ,ONLY: FV_NormVecZeta,FV_SurfElemZeta_sw
#endif
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
INTEGER                                                                    :: p,q,i,j,k
INTEGER                                                                    :: vn1,vn2,r,s,m,l,rr,ss,rrr,sss
REAL,DIMENSION(PP_nVarPrim,PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ)              :: PrimConsJac
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)              :: fJacQx,gJacQx,hJacQx
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)              :: fJacQy,gJacQy,hJacQy
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)              :: fJacQz,gJacQz,hJacQz
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: JacQx_plus,JacQx_minus,JacQy_plus,JacQy_minus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: fJacTilde_plus,gJacTilde_plus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: fJacTilde_minus,gJacTilde_minus
#if PP_dim==3
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: JacQz_plus,JacQz_minus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: hJacTilde_minus,hJacTilde_plus
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,0:PP_N,3)      :: JacLifting_Z
#endif
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim):: JacLifting_X,JacLifting_Y
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                                    :: JacTilde
REAL,DIMENSION(PP_nVar    ,PP_nVar)                                        :: Jac
!===================================================================================================================================
vn1=PP_nVar*(PP_N+1)
vn2=vn1*(PP_N+1)
! Calculation of derivative of gradUPrim with respect to UPrim_vol
CALL JacLifting_VolInt(1,iElem,JacLifting_X) !d(Q^1)/dU(1:3) (3 sums)
CALL JacLifting_VolInt(2,iElem,JacLifting_Y) !d(Q^2)/dU(1:3)
#if PP_dim==3
CALL JacLifting_VolInt(3,iElem,JacLifting_Z) !d(Q^3)/dU(1:3)
#endif

! Calculation of the derivative of the diffusive flux with respect to gradUPrim
! dF/dQ
CALL EvalFluxGradJacobian(nDOFElem,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem) &
                          ,fJacQx,fJacQy,fJacQz &
                          ,gJacQx,gJacQy,gJacQz &
                          ,hJacQx,hJacQy,hJacQz &
#if EDDYVISCOSITY
                          ,muSGS(:,:,:,:,iElem) &
#endif
                         )

! Calculation of the derivative of UPrim_vol with respect to UCons_vol
DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  CALL dPrimTempdCons(UPrim(:,i,j,k,iElem),PrimConsJac(:,:,i,j,k))
END DO; END DO; END DO! i,j,k=0,PP_N

! === Xi-Direction ================================================================================================================
DO p=0,PP_N
  DO q=0,PP_NZ
    DO i=1,PP_N
      fJacTilde_plus =fJacQx(:,:,i-1,p,q)
      gJacTilde_plus =gJacQx(:,:,i-1,p,q)
      fJacTilde_minus=fJacQx(:,:,i  ,p,q)
      gJacTilde_minus=gJacQx(:,:,i  ,p,q)
#if PP_dim==3
      hJacTilde_plus =hJacQx(:,:,i-1,p,q)
      hJacTilde_minus=hJacQx(:,:,i  ,p,q)
#endif
      JacQx_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecXi(1,p,q,i,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecXi(2,p,q,i,iElem) &
#if PP_dim==3
                         +hJacTilde_plus( :,:)*FV_NormVecXi(3,p,q,i,iElem) &
#endif
                         )*FV_SurfElemXi_sw(p,q,i,iElem)
      JacQx_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecXi(1,p,q,i,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecXi(2,p,q,i,iElem) &
#if PP_dim==3
                         +hJacTilde_minus(:,:)*FV_NormVecXi(3,p,q,i,iElem) &
#endif
                         )*FV_SurfElemXi_sw(p,q,i,iElem)
      fJacTilde_plus =fJacQy(:,:,i-1,p,q)
      gJacTilde_plus =gJacQy(:,:,i-1,p,q)
      fJacTilde_minus=fJacQy(:,:,i  ,p,q)
      gJacTilde_minus=gJacQy(:,:,i  ,p,q)
#if PP_dim==3
      hJacTilde_plus =hJacQy(:,:,i-1,p,q)
      hJacTilde_minus=hJacQy(:,:,i  ,p,q)
#endif
      JacQy_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecXi(1,p,q,i,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecXi(2,p,q,i,iElem) &
#if PP_dim==3
                         +hJacTilde_plus( :,:)*FV_NormVecXi(3,p,q,i,iElem) &
#endif
                         )*FV_SurfElemXi_sw(p,q,i,iElem)
      JacQy_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecXi(1,p,q,i,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecXi(2,p,q,i,iElem) &
#if PP_dim==3
                         +hJacTilde_minus(:,:)*FV_NormVecXi(3,p,q,i,iElem) &
#endif
                         )*FV_SurfElemXi_sw(p,q,i,iElem)
#if PP_dim==3
      fJacTilde_plus =fJacQz(:,:,i-1,p,q)
      gJacTilde_plus =gJacQz(:,:,i-1,p,q)
      fJacTilde_minus=fJacQz(:,:,i  ,p,q)
      gJacTilde_minus=gJacQz(:,:,i  ,p,q)
      hJacTilde_plus =hJacQz(:,:,i-1,p,q)
      hJacTilde_minus=hJacQz(:,:,i  ,p,q)
      JacQz_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecXi(1,p,q,i,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecXi(2,p,q,i,iElem) &
                         +hJacTilde_plus( :,:)*FV_NormVecXi(3,p,q,i,iElem) &
                         )*FV_SurfElemXi_sw(p,q,i,iElem)
      JacQz_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecXi(1,p,q,i,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecXi(2,p,q,i,iElem) &
                         +hJacTilde_minus(:,:)*FV_NormVecXi(3,p,q,i,iElem) &
                         )*FV_SurfElemXi_sw(p,q,i,iElem)
#endif
      ! assemble preconditioner
      r = vn1*p + vn2*q + PP_nVar*(i-1)
      s = vn1*p + vn2*q + PP_nVar*i
      m = vn1*p + vn2*q + PP_nVar*(i-2)
      l = vn1*p + vn2*q + PP_nVar*(i+1)
      IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
        CALL Assemble_FVVolIntGradJac(i,r,s,m,l,PrimConsJac(:,:,:,p,q),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                      JacLifting_X(:,:,:,p,q,:,2),JacLifting_Y(:,:,:,p,q,:,2),                        &
#if PP_dim==3
                                      fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,:,p,q,:,2),                           &
#endif
                                      BJ,.FALSE.)
        rr  = vn1*(p-1) + vn2*q + PP_nVar*(i-1)
        rrr = vn1*(p+1) + vn2*q + PP_nVar*(i-1)
        ss  = vn1*(p-1) + vn2*q + PP_nVar*i
        sss = vn1*(p+1) + vn2*q + PP_nVar*i
        IF(i.GE.1)THEN
          IF(p.GE.1)THEN ! fill dependency on lower left dof
            ! (i-1) Flux w.r.t (i-1,p-1) State
            JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,p,q,i-1,2) )  &
                           + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,p,q,i-1,2) )  &
#if PP_dim==3
                           + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,p,q,i-1,2) )  &
#endif
                           )
            Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-1,p-1,q))
            BJ(r+1:r+PP_nVar,rr+1:rr+PP_nVar) = BJ(r+1:r+PP_nVar,rr+1:rr+PP_nVar) + Jac
            ! (i) Flux w.r.t (i-1,p-1) State
            BJ(s+1:s+PP_nVar,rr+1:rr+PP_nVar) = BJ(s+1:s+PP_nVar,rr+1:rr+PP_nVar) - Jac
          END IF
          IF(p.LT.PP_N)THEN ! fill dependency on upper left dof
            ! (i-1) Flux w.r.t (i-1,p+1) State
            JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,p,q,i+1,2) )  &
                           + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,p,q,i+1,2) )  &
#if PP_dim==3
                           + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,p,q,i+1,2) )  &
#endif
                           )
            Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-1,p+1,q))
            BJ(r+1:r+PP_nVar,rrr+1:rrr+PP_nVar) = BJ(r+1:r+PP_nVar,rrr+1:rrr+PP_nVar) + Jac
            ! (i) Flux w.r.t (i-1,p+1) State
            BJ(s+1:s+PP_nVar,rrr+1:rrr+PP_nVar) = BJ(s+1:s+PP_nVar,rrr+1:rrr+PP_nVar) - Jac
          END IF
        END IF
        IF(p.GE.1)THEN ! fill dependency on lower right dof
          ! (i-1) Flux w.r.t (i,p-1) State
          JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,p,q,i-1,2) )  &
                         + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,p,q,i-1,2) )  &
#if PP_dim==3
                         + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,p,q,i-1,2) )  &
#endif
                         )
          Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i,p-1,q))
          BJ(r+1:r+PP_nVar,ss+1:ss+PP_nVar) = BJ(r+1:r+PP_nVar,ss+1:ss+PP_nVar) + Jac
          ! (i) Flux w.r.t (i,p-1) State
          BJ(s+1:s+PP_nVar,ss+1:ss+PP_nVar) = BJ(s+1:s+PP_nVar,ss+1:ss+PP_nVar) - Jac
        END IF
        IF(p.LT.PP_N)THEN ! fill dependency on upper right dof
          ! (i-1) Flux w.r.t (i,p+1) State
          JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,p,q,i+1,2) )  &
                         + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,p,q,i+1,2) )  &
#if PP_dim==3
                         + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,p,q,i+1,2) )  &
#endif
                         )
          Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i,p+1,q))
          BJ(r+1:r+PP_nVar,sss+1:sss+PP_nVar) = BJ(r+1:r+PP_nVar,sss+1:sss+PP_nVar) + Jac
          ! (i) Flux w.r.t (i,p+1) State
          BJ(s+1:s+PP_nVar,sss+1:sss+PP_nVar) = BJ(s+1:s+PP_nVar,sss+1:sss+PP_nVar) - Jac
        END IF
#if PP_dim==3
        !rrr = vn1*p + vn2*(q+1) + PP_nVar*(i-1)
        !sss = vn1*p + vn2*(q+1) + PP_nVar*i
        !rr  = vn1*p + vn2*(q-1) + PP_nVar*(i-1)
        !ss  = vn1*p + vn2*(q-1) + PP_nVar*i
        CALL Assemble_FVVolIntGradJac(i,r,s,m,l,PrimConsJac(:,:,:,p,q),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                      JacLifting_X(:,:,:,p,q,:,3),JacLifting_Y(:,:,:,p,q,:,3),                        &
                                      fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,:,p,q,:,3),                           &
                                      BJ,.FALSE.)
#endif
      END IF
      CALL Assemble_FVVolIntGradJac(i,r,s,m,l,PrimConsJac(:,:,:,p,q),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                    JacLifting_X(:,:,:,p,q,:,1),JacLifting_Y(:,:,:,p,q,:,1),                        &
#if PP_dim==3
                                    fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,:,p,q,:,1),                           &
#endif
                                    BJ,.TRUE.)
    END DO !i
  END DO !q
END DO !p
! === Eta-Direction ===============================================================================================================
DO p=0,PP_N
  DO q=0,PP_NZ
    DO j=1,PP_N
      fJacTilde_plus =fJacQx(:,:,p,j-1,q)
      gJacTilde_plus =gJacQx(:,:,p,j-1,q)
      fJacTilde_minus=fJacQx(:,:,p,j  ,q)
      gJacTilde_minus=gJacQx(:,:,p,j  ,q)
#if PP_dim==3
      hJacTilde_plus =hJacQx(:,:,p,j-1,q)
      hJacTilde_minus=hJacQx(:,:,p,j  ,q)
#endif
      JacQx_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecEta(1,p,q,j,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecEta(2,p,q,j,iElem) &
#if PP_dim==3
                         +hJacTilde_plus( :,:)*FV_NormVecEta(3,p,q,j,iElem) &
#endif
                         )*FV_SurfElemEta_sw(p,q,j,iElem)
      JacQx_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecEta(1,p,q,j,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecEta(2,p,q,j,iElem) &
#if PP_dim==3
                         +hJacTilde_minus(:,:)*FV_NormVecEta(3,p,q,j,iElem) &
#endif
                         )*FV_SurfElemEta_sw(p,q,j,iElem)
      fJacTilde_plus =fJacQy(:,:,p,j-1,q)
      gJacTilde_plus =gJacQy(:,:,p,j-1,q)
      fJacTilde_minus=fJacQy(:,:,p,j  ,q)
      gJacTilde_minus=gJacQy(:,:,p,j  ,q)
#if PP_dim==3                     
      hJacTilde_plus =hJacQy(:,:,p,j-1,q)
      hJacTilde_minus=hJacQy(:,:,p,j  ,q)
#endif
      JacQy_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecEta(1,p,q,j,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecEta(2,p,q,j,iElem) &
#if PP_dim==3
                         +hJacTilde_plus( :,:)*FV_NormVecEta(3,p,q,j,iElem) &
#endif
                         )*FV_SurfElemEta_sw(p,q,j,iElem)
      JacQy_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecEta(1,p,q,j,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecEta(2,p,q,j,iElem) &
#if PP_dim==3
                         +hJacTilde_minus(:,:)*FV_NormVecEta(3,p,q,j,iElem) &
#endif
                         )*FV_SurfElemEta_sw(p,q,j,iElem)
#if PP_dim==3
      fJacTilde_plus =fJacQz(:,:,p,j-1,q)
      gJacTilde_plus =gJacQz(:,:,p,j-1,q)
      fJacTilde_minus=fJacQz(:,:,p,j  ,q)
      gJacTilde_minus=gJacQz(:,:,p,j  ,q)
      hJacTilde_plus =hJacQz(:,:,p,j-1,q)
      hJacTilde_minus=hJacQz(:,:,p,j  ,q)
      JacQz_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecEta(1,p,q,j,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecEta(2,p,q,j,iElem) &
                         +hJacTilde_plus( :,:)*FV_NormVecEta(3,p,q,j,iElem) &
                         )*FV_SurfElemEta_sw(p,q,j,iElem)
      JacQz_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecEta(1,p,q,j,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecEta(2,p,q,j,iElem) &
                         +hJacTilde_minus(:,:)*FV_NormVecEta(3,p,q,j,iElem) &
                         )*FV_SurfElemEta_sw(p,q,j,iElem)
#endif
      ! assemble preconditioner
      r = vn1*(j-1) + vn2*q + PP_nVar*p
      s = vn1*j     + vn2*q + PP_nVar*p
      m = vn1*(j-2) + vn2*q + PP_nVar*p
      l = vn1*(j+1) + vn2*q + PP_nVar*p
      IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
        CALL Assemble_FVVolIntGradJac(j,r,s,m,l,PrimConsJac(:,:,p,:,q),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                      JacLifting_X(:,:,p,:,q,:,1),JacLifting_Y(:,:,p,:,q,:,1),                        &
#if PP_dim==3
                                      fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,p,:,q,:,1),                           &
#endif
                                      BJ,.FALSE.)
#if PP_dim==3
        CALL Assemble_FVVolIntGradJac(j,r,s,m,l,PrimConsJac(:,:,p,:,q),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                      JacLifting_X(:,:,p,:,q,:,3),JacLifting_Y(:,:,p,:,q,:,3),                        &
                                      fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,p,:,q,:,3),                           &
                                      BJ,.FALSE.)
#endif
      END IF
      CALL Assemble_FVVolIntGradJac(j,r,s,m,l,PrimConsJac(:,:,p,:,q),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                    JacLifting_X(:,:,p,:,q,:,2),JacLifting_Y(:,:,p,:,q,:,2),                        &
#if PP_dim==3
                                    fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,p,:,q,:,2),                           &
#endif
                                    BJ,.TRUE.)
    END DO !j
  END DO !q
END DO !p
#if PP_dim==3
! === Zeta-Direction ==============================================================================================================
DO p=0,PP_N
  DO q=0,PP_N
    DO k=1,PP_N
      fJacTilde_plus =fJacQx(:,:,p,q,k-1)
      gJacTilde_plus =gJacQx(:,:,p,q,k-1)
      fJacTilde_minus=fJacQx(:,:,p,q,k  )
      gJacTilde_minus=gJacQx(:,:,p,q,k  )
      hJacTilde_plus =hJacQx(:,:,p,q,k-1)
      hJacTilde_minus=hJacQx(:,:,p,q,k  )
      JacQx_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecZeta(1,p,q,k,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecZeta(2,p,q,k,iElem) &
                         +hJacTilde_plus( :,:)*FV_NormVecZeta(3,p,q,k,iElem) &
                         )*FV_SurfElemZeta_sw(p,q,k,iElem)
      JacQx_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecZeta(1,p,q,k,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecZeta(2,p,q,k,iElem) &
                         +hJacTilde_minus(:,:)*FV_NormVecZeta(3,p,q,k,iElem) &
                         )*FV_SurfElemZeta_sw(p,q,k,iElem)
      fJacTilde_plus =fJacQy(:,:,p,q,k-1)
      gJacTilde_plus =gJacQy(:,:,p,q,k-1)
      fJacTilde_minus=fJacQy(:,:,p,q,k  )
      gJacTilde_minus=gJacQy(:,:,p,q,k  )
      hJacTilde_plus =hJacQy(:,:,p,q,k-1)
      hJacTilde_minus=hJacQy(:,:,p,q,k  )
      JacQy_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecZeta(1,p,q,k,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecZeta(2,p,q,k,iElem) &
                         +hJacTilde_plus( :,:)*FV_NormVecZeta(3,p,q,k,iElem) &
                         )*FV_SurfElemZeta_sw(p,q,k,iElem)
      JacQy_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecZeta(1,p,q,k,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecZeta(2,p,q,k,iElem) &
                         +hJacTilde_minus(:,:)*FV_NormVecZeta(3,p,q,k,iElem) &
                         )*FV_SurfElemZeta_sw(p,q,k,iElem)
      fJacTilde_plus =fJacQz(:,:,p,q,k-1)
      gJacTilde_plus =gJacQz(:,:,p,q,k-1)
      fJacTilde_minus=fJacQz(:,:,p,q,k  )
      gJacTilde_minus=gJacQz(:,:,p,q,k  )
      hJacTilde_plus =hJacQz(:,:,p,q,k-1)
      hJacTilde_minus=hJacQz(:,:,p,q,k  )
      JacQz_plus  =  0.5*(fJacTilde_plus( :,:)*FV_NormVecZeta(1,p,q,k,iElem) &
                         +gJacTilde_plus( :,:)*FV_NormVecZeta(2,p,q,k,iElem) &
                         +hJacTilde_plus( :,:)*FV_NormVecZeta(3,p,q,k,iElem) &
                         )*FV_SurfElemZeta_sw(p,q,k,iElem)
      JacQz_minus = 0.5*(fJacTilde_minus(:,:)*FV_NormVecZeta(1,p,q,k,iElem) &
                         +gJacTilde_minus(:,:)*FV_NormVecZeta(2,p,q,k,iElem) &
                         +hJacTilde_minus(:,:)*FV_NormVecZeta(3,p,q,k,iElem) &
                         )*FV_SurfElemZeta_sw(p,q,k,iElem)
      ! assemble preconditioner
      r = vn1*q + vn2*(k-1) + PP_nVar*p
      s = vn1*q + vn2*k     + PP_nVar*p
      m = vn1*q + vn2*(k-2) + PP_nVar*p
      l = vn1*q + vn2*(k+1) + PP_nVar*p
      IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
        CALL Assemble_FVVolIntGradJac(k,r,s,m,l,PrimConsJac(:,:,p,q,:),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                      JacLifting_X(:,:,:,:,p,q,1),JacLifting_Y(:,:,:,:,p,q,1),                        &
                                      fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,:,:,p,q,1),                           &
                                      BJ,.FALSE.)
        CALL Assemble_FVVolIntGradJac(k,r,s,m,l,PrimConsJac(:,:,p,q,:),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                      JacLifting_X(:,:,:,:,p,q,2),JacLifting_Y(:,:,:,:,p,q,2),                        &
                                      fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,:,:,p,q,2),                           &
                                      BJ,.FALSE.)
      END IF                          
      CALL Assemble_FVVolIntGradJac(k,r,s,m,l,PrimConsJac(:,:,p,q,:),JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus, &
                                    JacLifting_X(:,:,:,:,p,q,3),JacLifting_Y(:,:,:,:,p,q,3),                        &
                                    fJacQz_plus,fJacQy_minus,JacLifting_Z(:,:,:,:,p,q,3),                           &
                                    BJ,.TRUE.)
    END DO !k
  END DO !q
END DO !p
#endif
END SUBROUTINE FVVolIntGradJac

!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE Assemble_FVVolIntGradJac(i,r,s,m,l,PrimConsJac,JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus,JacLifting_X,JacLifting_Y, &
#if PP_dim==3
                                    JacQz_plus,JacQz_minus,JacLifting_Z, &
#endif
                                    BJ,isAligned)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Implicit_Vars ,ONLY: nDOFVarElem
USE MOD_Precond_Vars  ,ONLY: NoFillIn
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                               :: i,r,s,m,l
REAL,DIMENSION(PP_nVarPrim,PP_nVar    ,0:PP_N)       ,INTENT(IN) :: PrimConsJac
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)              ,INTENT(IN) :: JacQx_plus,JacQy_plus,JacQx_minus,JacQy_minus
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N),INTENT(IN) :: JacLifting_X,JacLifting_Y
#if PP_dim==3
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)              ,INTENT(IN) :: JacQz_plus,JacQz_minus
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N),INTENT(IN) :: JacLifting_X,JacLifting_Y
#endif
LOGICAL,INTENT(IN)                                               :: isAligned !< .TRUE. if e.g. flux over XI-sides w.r.t. gradients
                                                                              !> in x-direction are considered. This leads to a
                                                                              !> five point stencil in one dimension (4 points for
                                                                              !> the considered flux over a specific sub-cell side)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: BJ(1:nDOFVarElem,1:nDOFVarElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim)                          :: JacTilde
REAL,DIMENSION(PP_nVar    ,PP_nVar)                              :: Jac
!===================================================================================================================================
!----------- derivatives w.r.t. (i-2) -------------!
IF((i.GT.1).AND.isAligned)THEN
  ! (i-1) Flux w.r.t (i-2) State
  JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,i-2) )  &
                 + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,i-2) )  &
#if PP_dim==3
                 + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,i-2) )  &
#endif
                 )
  Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-2))
  BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) = BJ(r+1:r+PP_nVar,m+1:m+PP_nVar) + Jac
  IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
    ! (i) Flux w.r.t. (i-2) State
    JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,i-2) )  &
                   + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,i-2) )  &
#if PP_dim==3
                   + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,i-2) )  &
#endif
                   )
    Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-2))
    BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) = BJ(s+1:s+PP_nVar,m+1:m+PP_nVar) - Jac
  END IF
END IF
!----------- derivatives w.r.t. (i-1) -------------!
! (i-1) Flux w.r.t (i-1) State
JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,i-1) )  &
               + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,i-1) )  &
#if PP_dim==3
               + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,i-1) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-1))
BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + Jac
JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,i-1) )  &
               + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,i-1) )  &
#if PP_dim==3
               + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,i-1) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-1))
BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + Jac
! (i) Flux w.r.t (i-1) State
JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,i-1) )  &
               + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,i-1) )  &
#if PP_dim==3
               + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,i-1) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-1))
BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) - Jac
JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,i-1) )  &
               + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,i-1) )  &
#if PP_dim==3
               + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,i-1) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i-1))
BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) - Jac
!----------- derivatives w.r.t. (i) -------------!
! (i-1) Flux w.r.t (i) State
JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,i) )  &
               + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,i) )  &
#if PP_dim==3
               + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,i) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i))
BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + Jac
JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,i) )  &
               + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,i) )  &
#if PP_dim==3
               + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,i) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i))
BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + Jac
! (i) Flux w.r.t (i) State
JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,i) )  &
               + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,i) )  &
#if PP_dim==3
               + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,i) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i))
BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) - Jac
JacTilde(:,:)= ( MATMUL(JacQx_plus , JacLifting_X(:,:,i-1,i) )  &
               + MATMUL(JacQy_plus , JacLifting_Y(:,:,i-1,i) )  &
#if PP_dim==3
               + MATMUL(JacQz_plus , JacLifting_Z(:,:,i-1,i) )  &
#endif
               )
Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i))
BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) - Jac
!----------- derivatives w.r.t. (i+1,p,q) -------------!
IF((i.LT.PP_N).AND.isAligned)THEN
  IF(NoFillIn.EQV..FALSE.) THEN !NoFillIn has the same sparsity as the EulerPrecond
    ! (i-1) Flux w.r.t (i+1) State
    JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,i+1) )  &
                   + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,i+1) )  &
#if PP_dim==3
                   + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,i+1) )  &
#endif
                   )
    Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i+1))
    BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) = BJ(r+1:r+PP_nVar,l+1:l+PP_nVar) + Jac
  END IF
  ! (i) Flux w.r.t (i+1) State
  JacTilde(:,:)= ( MATMUL(JacQx_minus, JacLifting_X(:,:,i  ,i+1) )  &
                 + MATMUL(JacQy_minus, JacLifting_Y(:,:,i  ,i+1) )  &
#if PP_dim==3
                 + MATMUL(JacQz_minus, JacLifting_Z(:,:,i  ,i+1) )  &
#endif
                 )
  Jac(:,:) = MATMUL(JacTilde(:,:),PrimConsJac(:,:,i+1))
  BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) = BJ(s+1:s+PP_nVar,l+1:l+PP_nVar) - Jac
END IF

END SUBROUTINE Assemble_FVVolIntGradJac
#endif /*PARABOLIC*/
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
#if PARABOLIC
SDEALLOCATE(nVec)
SDEALLOCATE(Surf)
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
