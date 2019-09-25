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
#include "eos.h"

!===================================================================================================================================
!> This module contains routines required for the jacobian of the br2 lifting scheme.
!===================================================================================================================================
MODULE MOD_Jac_br2
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE FillJacLiftingFlux
  MODULE PROCEDURE FillJacLiftingFlux
END INTERFACE

INTERFACE JacLifting_VolInt
  MODULE PROCEDURE JacLifting_VolInt
END INTERFACE

INTERFACE dQOuter
  MODULE PROCEDURE dQOuter
END INTERFACE

INTERFACE dQInner
  MODULE PROCEDURE dQInner
END INTERFACE

PUBLIC::FillJacLiftingFlux,JacLifting_VolInt,BuildnVecTangSurf,Build_BR2_SurfTerms,dQOuter,dQInner
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Derivative of the numerical flux h* of the gradient system with respect to U_jk
!> h*= 0.5*(U_xi + U_neighbour_(-xi)) (mean value of the face values)
!> => h*-U_xi=-0.5*U_xi + 0.5*U_neighbour_(-xi)
!> => d(h*-U_xi)_jk(prim)/dU_jk(prim) = -0.5
!===================================================================================================================================
SUBROUTINE FillJacLiftingFlux(t,iElem)
! MODULES
USE MOD_Globals
USE MOD_Jac_ex_Vars               ,ONLY:JacLiftingFlux,Surf
USE MOD_GetBoundaryFlux_fd        ,ONLY:Lifting_GetBoundaryFlux_FD
USE MOD_Mesh_Vars                 ,ONLY:nBCSides,ElemToSide,S2V2
USE MOD_Mesh_Vars                 ,ONLY:NormVec,TangVec1,TangVec2,SurfElem
USE MOD_DG_Vars                   ,ONLY:UPrim_master
USE MOD_Mesh_Vars                 ,ONLY:Face_xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL   ,INTENT(IN) :: t
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
    CALL Lifting_GetBoundaryFlux_FD(SideID,t,JacLiftingFlux(:,:,:,:,iLocSide),UPrim_master, &
                                    SurfElem,Face_xGP,NormVec,TangVec1,TangVec2,S2V2(:,:,:,0,iLocSide)) !flip=0 for BCSide
  ELSE
    JacLiftingFlux(:,:,:,:,iLocSide)=0.
    DO iVar=1,PP_nVarPrim
      JacLiftingFlux(iVar,iVar,:,:,iLocSide)=-0.5*Surf(:,:,iLocSide,iElem)
    END DO !iVar
  END IF !SideID
END DO!iLocSide
END SUBROUTINE FillJacLiftingFlux

!===================================================================================================================================
!> Computes the Volume gradient Jacobian of the BR2 scheme dQprim/dUprim (Q= Grad U)
!> Normal vectors are supposed to point outwards!
!===================================================================================================================================
SUBROUTINE JacLifting_VolInt(dir,iElem,JacLifting)
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
REAL,INTENT(OUT)                             :: JacLifting(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: i,j,k,ll
INTEGER                                      :: iVar
REAL                                         :: delta(1:PP_nVarPrim,1:PP_nVarPrim) 
!===================================================================================================================================
delta=0.
DO iVar=1,PP_nVarPrim
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

!===================================================================================================================================
!> used for BR2: normal vectors, outward pointing and sorted in ijk element fashion!!! 
!> The usual NormVec is only outward poiting for the master Sides.
!> nVec are also outward pointing for the slave Sides
!===================================================================================================================================
SUBROUTINE BuildnVecTangSurf()
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

!===================================================================================================================================
!> used for BR2: normal vectors, outward pointing and sorted in ijk element fashion!!! 
!===================================================================================================================================
SUBROUTINE Build_BR2_SurfTerms()
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
REAL,INTENT(OUT)                                    :: dQ_dUVolInner(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_NZ,6,0:PP_N)
#else
REAL,INTENT(OUT)                                    :: dQ_dUVolInner(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_NZ,2:5,0:PP_N)
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
REAL                                                :: r(0:PP_N,0:PP_NZ)
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

END MODULE MOD_Jac_br2
#endif /*PARABOLIC*/
