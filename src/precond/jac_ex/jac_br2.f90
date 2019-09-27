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
INTEGER            :: iVar,iLocSide,SideID,flip
REAL               :: mp
!===================================================================================================================================
#if PP_dim==3
DO iLocSide=1,6
#else
DO iLocSide=2,5
#endif 
  SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
  flip  =ElemToSide(E2S_FLIP   ,ilocSide,iElem)
  IF (SideID.LE.nBCSides) THEN !BCSides
    CALL Lifting_GetBoundaryFlux_FD(SideID,t,JacLiftingFlux(:,:,:,:,iLocSide),UPrim_master, &
                                    SurfElem,Face_xGP,NormVec,TangVec1,TangVec2,S2V2(:,:,:,0,iLocSide)) !flip=0 for BCSide
  ELSE
    JacLiftingFlux(:,:,:,:,iLocSide)=0.
    IF(flip.EQ.0)THEN
      mp = -1.
    ELSE
      mp = 1.
    END IF
    DO iVar=1,PP_nVarPrim
      JacLiftingFlux(iVar,iVar,:,:,iLocSide) = mp*0.5*Surf(:,:,iLocSide,iElem)
    END DO !iVar
  END IF !SideID
END DO!iLocSide
END SUBROUTINE FillJacLiftingFlux

!===================================================================================================================================
!> Computes the Volume gradient Jacobian of the BR2 scheme dQprim/dUprim (Q= Grad U)
!> Normal vectors are supposed to point outwards!
!===================================================================================================================================
SUBROUTINE JacLifting_VolInt(dir,iElem,UPrim,JacLifting)
! MODULES
USE MOD_Jac_Ex_Vars        ,ONLY: LL_minus,LL_plus
USE MOD_Jac_Ex_Vars        ,ONLY: JacLiftingFlux 
USE MOD_DG_Vars            ,ONLY: D,UPrim_master,UPrim_slave
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,sJ   ! metrics
#if PP_dim==3
USE MOD_Mesh_Vars          ,ONLY: Metrics_hTilde
#endif
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: ElemToSide,S2V2,Normvec
USE MOD_Jacobian           ,ONLY: dConsdPrimTemp,dPrimTempdCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                           :: iElem
INTEGER,INTENT(IN)                           :: dir
REAL,INTENT(IN)                              :: UPrim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                             :: JacLifting(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: i,j,k,ll
INTEGER                                      :: iVar
REAL,DIMENSION(1:PP_nVarPrim,1:PP_nVarPrim)  :: delta,SurfVol_PrimJac_plus,SurfVol_PrimJac_minus
REAL,DIMENSION(1:PP_nVar    ,1:PP_nVarPrim)  :: ConsPrimJac
REAL,DIMENSION(1:PP_nVarPrim,1:PP_nVar    )  :: PrimConsJac
INTEGER :: pq_p(2),pq_m(2),SideID_p,SideID_m,flip
REAL    :: mp_plus,mp_minus
!===================================================================================================================================
! fill volume jacobian of lifting flux (flux is primitive state vector -> identiy matrix)
delta=0.
DO iVar=1,PP_nVarPrim
  delta(iVar,iVar)=1.
END DO

JacLifting=0.
DO ll=0,PP_N
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        ! ++++++++++++++++++++++++++++ Variant without special built SurfElem and NormVec ++++++++++++++++++++++++++++++++++++++++
        ! XI sides
        CALL dConsdPrimTemp(UPrim(:,ll,j,k),ConsPrimJac)
        SideID_m=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
        flip=ElemToSide(    E2S_FLIP   ,XI_MINUS,iElem)
        pq_m=S2V2(:,j,k,flip           ,XI_MINUS)
        IF(flip.EQ.0)THEN
          mp_minus = -1.
          CALL dPrimTempdCons(UPrim_master(:,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
        ELSE
          mp_minus = 1.
          CALL dPrimTempdCons(UPrim_slave(:,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
        END IF
        SurfVol_PrimJac_minus = MATMUL(PrimConsJac,ConsPrimJac)

        SideID_p =ElemToSide(E2S_SIDE_ID,XI_PLUS ,iElem)
        flip=ElemToSide(     E2S_FLIP   ,XI_PLUS ,iElem)
        pq_p=S2V2(:,j,k,flip            ,XI_PLUS)
        IF(flip.EQ.0)THEN
          mp_plus = -1.
          CALL dPrimTempdCons(UPrim_master(:,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
        ELSE
          mp_plus = 1.
          CALL dPrimTempdCons(UPrim_slave(:,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
        END IF
        SurfVol_PrimJac_plus = MATMUL(PrimConsJac,ConsPrimJac)

        JacLifting(:,:,i,j,k,ll,1) = JacLifting(:,:,i,j,k,ll,1) + sJ(i,j,k,iElem,0)*                       &
                                     ( D(i,ll)*Metrics_fTilde(dir,ll,j,k,iElem,0)*delta(:,:)               &
                                      +LL_plus( i,ll)*NormVec(dir,pq_p(1),pq_p(2),0,SideID_p)*             &
                                       MATMUL(JacLiftingFlux(:,:,j,k,XI_PLUS),SurfVol_PrimJac_plus(:,:))   &
                                      +LL_minus(i,ll)*NormVec(dir,pq_m(1),pq_m(2),0,SideID_m)*             &
                                       MATMUL(JacLiftingFlux(:,:,j,k,XI_MINUS),SurfVol_PrimJac_minus(:,:)))
        ! ETA sides
        CALL dConsdPrimTemp(UPrim(:,i,ll,k),ConsPrimJac)
        SideID_m=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
        flip=ElemToSide(    E2S_FLIP   ,ETA_MINUS,iElem)
        pq_m=S2V2(:,i,k,flip           ,ETA_MINUS)
        IF(flip.EQ.0)THEN
          mp_minus = -1.
          CALL dPrimTempdCons(UPrim_master(:,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
        ELSE
          mp_minus = 1.
          CALL dPrimTempdCons(UPrim_slave(:,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
        END IF
        SurfVol_PrimJac_minus = MATMUL(PrimConsJac,ConsPrimJac)

        SideID_p =ElemToSide(E2S_SIDE_ID,ETA_PLUS ,iElem)
        flip=ElemToSide(     E2S_FLIP   ,ETA_PLUS ,iElem)
        pq_p=S2V2(:,i,k,flip            ,ETA_PLUS)
        IF(flip.EQ.0)THEN
          mp_plus = -1.
          CALL dPrimTempdCons(UPrim_master(:,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
        ELSE
          mp_plus = 1.
          CALL dPrimTempdCons(UPrim_slave(:,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
        END IF
        SurfVol_PrimJac_plus = MATMUL(PrimConsJac,ConsPrimJac)

        JacLifting(:,:,i,j,k,ll,2) = JacLifting(:,:,i,j,k,ll,2) + sJ(i,j,k,iElem,0)*                        &
                                     ( D(j,ll)*Metrics_gTilde(dir,i,ll,k,iElem,0)*delta(:,:)                &
                                      +LL_plus( j,ll)*NormVec(dir,pq_p(1),pq_p(2),0,SideID_p)*              &
                                       MATMUL(JacLiftingFlux(:,:,i,k,ETA_PLUS),SurfVol_PrimJac_plus(:,:))   &
                                      +LL_minus(j,ll)*NormVec(dir,pq_m(1),pq_m(2),0,SideID_m)*              &
                                       MATMUL(JacLiftingFlux(:,:,i,k,ETA_MINUS),SurfVol_PrimJac_minus(:,:)))
#if PP_dim==3
        ! ZETA sides
        CALL dConsdPrimTemp(UPrim(:,i,j,ll),ConsPrimJac)
        SideID_m=ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
        flip=ElemToSide(    E2S_FLIP   ,ZETA_MINUS,iElem)
        pq_m=S2V2(:,i,j,flip           ,ZETA_MINUS)
        IF(flip.EQ.0)THEN
          mp_minus = -1.
          CALL dPrimTempdCons(UPrim_master(:,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
        ELSE
          mp_minus = 1.
          CALL dPrimTempdCons(UPrim_slave(:,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
        END IF
        SurfVol_PrimJac_minus = MATMUL(PrimConsJac,ConsPrimJac)

        SideID_p =ElemToSide(E2S_SIDE_ID,ZETA_PLUS ,iElem)
        flip=ElemToSide(     E2S_FLIP   ,ZETA_PLUS ,iElem)
        pq_p=S2V2(:,i,j,flip            ,ZETA_PLUS)
        IF(flip.EQ.0)THEN
          mp_plus = -1.
          CALL dPrimTempdCons(UPrim_master(:,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
        ELSE
          mp_plus = 1.
          CALL dPrimTempdCons(UPrim_slave(:,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
        END IF
        SurfVol_PrimJac_plus = MATMUL(PrimConsJac,ConsPrimJac)

        JacLifting(:,:,i,j,k,ll,3) = JacLifting(:,:,i,j,k,ll,3) + sJ(i,j,k,iElem,0)*                         &
                                     ( D(k,ll)*Metrics_hTilde(dir,i,j,ll,iElem,0)*delta(:,:)                 &
                                      +LL_plus( k,ll)*NormVec(dir,pq_p(1),pq_p(2),0,SideID_p)*               &
                                       MATMUL(JacLiftingFlux(:,:,i,j,ZETA_PLUS),SurfVol_PrimJac_plus(:,:))   &
                                      +LL_minus(k,ll)*NormVec(dir,pq_m(1),pq_m(2),0,SideID_m)*               &
                                       MATMUL(JacLiftingFlux(:,:,i,j,ZETA_MINUS),SurfVol_PrimJac_minus(:,:)))
#endif

        ! ++++++++++++++++++++++++++++ Variant using special built SurfElem and NormVec ++++++++++++++++++++++++++++++++++++++++

        !JacLifting(:,:,i,j,k,ll,1) = JacLifting(:,:,i,j,k,ll,1) + sJ(i,j,k,iElem,0)*                                    &
                                     !( D(i,ll)*Metrics_fTilde(dir,ll,j,k,iElem,0)*delta(:,:)                            &
                                      !+LL_plus( i,ll)*nVec(dir,j,k,XI_PLUS ,iElem)*JacLiftingFlux(:,:,j,k,XI_PLUS)      &
                                      !+LL_minus(i,ll)*nVec(dir,j,k,XI_MINUS,iElem)*JacLiftingFlux(:,:,j,k,XI_MINUS))

        !JacLifting(:,:,i,j,k,ll,2) = JacLifting(:,:,i,j,k,ll,2) + sJ(i,j,k,iElem,0)*                                    &
                                     !( D(j,ll)*Metrics_gTilde(dir,i,ll,k,iElem,0)*delta(:,:)                            &
                                      !+LL_plus( j,ll)*nVec(dir,i,k,ETA_PLUS ,iElem)*JacLiftingFlux(:,:,i,k,ETA_PLUS)    &
                                      !+LL_minus(j,ll)*nVec(dir,i,k,ETA_MINUS,iElem)*JacLiftingFlux(:,:,i,k,ETA_MINUS))
!#if PP_dim==3
        !JacLifting(:,:,i,j,k,ll,3) = JacLifting(:,:,i,j,k,ll,3) + sJ(i,j,k,iElem,0)*                                    &
                                     !( D(k,ll)*Metrics_hTilde(dir,i,j,ll,iElem,0)*delta(:,:)                            &
                                      !+LL_plus( k,ll)*nVec(dir,i,j,ZETA_PLUS ,iElem)*JacLiftingFlux(:,:,i,j,ZETA_PLUS)  &
                                      !+LL_minus(k,ll)*nVec(dir,i,j,ZETA_MINUS,iElem)*JacLiftingFlux(:,:,i,j,ZETA_MINUS))
!#endif
      END DO !i
    END DO !j
  END DO !k
END DO !ll
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
INTEGER :: iElem,p,q,l,iLocSide,SideID,Flip
REAL    :: RFace(0:PP_N,0:PP_NZ)
#if USE_MPI
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
USE MOD_DG_Vars                   ,ONLY: D,UPrim,UPrim_master,UPrim_slave
USE MOD_Lifting_Vars              ,ONLY: etaBR2
USE MOD_Jacobian                  ,ONLY: dConsdPrimTemp,dPrimTempdCons
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
REAL                                                :: r (0:PP_N,0:PP_NZ),mp
REAL                                                :: UPrim_face(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                                                :: dQ_dUVolInner_loc(0:PP_N,0:PP_NZ,0:PP_N)
REAL                                                :: ConsPrimJac(    1:PP_nVar    ,1:PP_nVarPrim)
REAL                                                :: PrimConsJac(    1:PP_nVarPrim,1:PP_nVar    )
REAL                                                :: SurfVol_PrimJac(1:PP_nVarPrim,1:PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N)
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
        UPrim_face(:,jk(1),jk(2)) = UPrim_master(:,p,q,SideID)
      END DO !p
    END DO !q
    mp = 1.
  ELSE !slave
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Plus(dir,p,q,SideID)
        UPrim_face(:,jk(1),jk(2)) = UPrim_slave(:,p,q,SideID)
      END DO !p
    END DO !q
    mp = -1.
  END IF !Flip=0

  dQ_dUVolInner_loc=0.
  SELECT CASE(iLocSide)

  CASE(XI_MINUS,XI_PLUS)
    DO mm=0,PP_N
      DO k=0,PP_NZ
        DO j=0,PP_N
          dQ_dUVolInner_loc(j,k,mm)   = dQ_dUVolInner_loc(j,k,mm) + etaBR2 * r(j,k) * l_mp(mm,iLocSide)
          CALL dConsdPrimTemp(UPrim(:,mm,j,k,iElem),ConsPrimJac)
          CALL dPrimTempdCons(UPrim_face(:,j,k),PrimConsJac)
          SurfVol_PrimJac(:,:,j,k,mm) = MATMUL(PrimConsJac,ConsPrimJac)
        END DO !j
      END DO !k
    END DO !mm
  CASE(ETA_MINUS,ETA_PLUS)
    DO nn=0,PP_N
      DO k=0,PP_NZ
        DO i=0,PP_N
          dQ_dUVolInner_loc(i,k,nn) = dQ_dUVolInner_loc(i,k,nn) + etaBR2 * r(i,k) * l_mp(nn,iLocSide)
          CALL dConsdPrimTemp(UPrim(:,i,nn,k,iElem),ConsPrimJac)
          CALL dPrimTempdCons(UPrim_face(:,i,k),PrimConsJac)
          SurfVol_PrimJac(:,:,i,k,nn) = MATMUL(PrimConsJac,ConsPrimJac)
        END DO !i
      END DO !k
    END DO !nn
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
    DO oo=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          dQ_dUVolInner_loc(i,j,oo) = dQ_dUVolInner_loc(i,j,oo) + etaBR2 * r(i,j) * l_mp(oo,iLocSide)
          CALL dConsdPrimTemp(UPrim(:,i,j,oo,iElem),ConsPrimJac)
          CALL dPrimTempdCons(UPrim_face(:,i,j),PrimConsJac)
          SurfVol_PrimJac(:,:,i,j,oo) = MATMUL(PrimConsJac,ConsPrimJac)
        END DO !i
      END DO !j
    END DO !oo
#endif
  END SELECT
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO mm=0,PP_N
        dQ_dUVolInner(:,:,j,k,iLocSide,mm)=mp*MATMUL(SurfVol_PrimJac(:,:,j,k,mm),JacLiftingFlux(:,:,j,k,iLocSide))* &
                                           dq_dUVolinner_loc(j,k,mm)
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
USE MOD_DG_Vars                   ,ONLY: UPrim,UPrim_master,UPrim_slave
USE MOD_Jacobian                  ,ONLY: dConsdPrimTemp,dPrimTempdCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: dir
INTEGER,INTENT(IN)                                  :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#if PP_dim == 3
REAL,INTENT(OUT)                                    :: dQ_dUVolOuter(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_NZ,6,0:PP_N)
#else
REAL,INTENT(OUT)                                    :: dQ_dUVolOuter(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: iLocSide,p,q,i,j,k,mm,nn
#if PP_dim == 3
INTEGER                                             :: oo
#endif
INTEGER                                             :: SideID,Flip,jk(2)
REAL                                                :: r(0:PP_N,0:PP_NZ)
REAL                                                :: UPrim_face(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                                                :: ConsPrimJac(    1:PP_nVar    ,1:PP_nVarPrim)
REAL                                                :: PrimConsJac(    1:PP_nVarPrim,1:PP_nVar    )
REAL                                                :: SurfVol_PrimJac(1:PP_nVarPrim,1:PP_nVarPrim)
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
        UPrim_face(:,jk(1),jk(2)) = UPrim_master(:,p,q,SideID)
      END DO !p
    END DO !q
  ELSE
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Minus(dir,p,q,SideID)
        UPrim_face(:,jk(1),jk(2)) = UPrim_slave(:,p,q,SideID)
      END DO !p
    END DO !q
  END IF !Flip=0

  SELECT CASE(iLocSide)
  CASE(XI_MINUS,XI_PLUS)
    DO mm=0,PP_N
      DO k=0,PP_NZ
        DO j=0,PP_N
          CALL dConsdPrimTemp(UPrim(:,mm,j,k,iElem),ConsPrimJac)
          CALL dPrimTempdCons(UPrim_face(:,j,k),PrimConsJac)
          SurfVol_PrimJac = MATMUL(PrimConsJac,ConsPrimJac)
          dQ_dUVolOuter(:,:,j,k,iLocSide,mm) = dQ_dUVolOuter(:,:,j,k,iLocSide,mm) + 0.5*etaBR2 * r(j,k) * &
                                               l_mp(mm,iLocSide)*Surf(j,k,iLocSide,iElem)*SurfVol_PrimJac
        END DO !j
      END DO !k
    END DO !mm
  CASE(ETA_MINUS,ETA_PLUS)
    DO nn=0,PP_N
      DO k=0,PP_NZ
        DO i=0,PP_N
          CALL dConsdPrimTemp(UPrim(:,i,nn,k,iElem),ConsPrimJac)
          CALL dPrimTempdCons(UPrim_face(:,i,k),PrimConsJac)
          SurfVol_PrimJac = MATMUL(PrimConsJac,ConsPrimJac)
          dQ_dUVolOuter(:,:,i,k,iLocSide,nn) = dQ_dUVolOuter(:,:,i,k,iLocSide,nn) + 0.5*etaBR2 * r(i,k) * &
                                               l_mp(nn,iLocSide)*Surf(i,k,iLocSide,iElem)*SurfVol_PrimJac

        END DO !i
      END DO !k
    END DO !nn
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
    DO oo=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          CALL dConsdPrimTemp(UPrim(:,i,j,oo,iElem),ConsPrimJac)
          CALL dPrimTempdCons(UPrim_face(:,i,j),PrimConsJac)
          SurfVol_PrimJac = MATMUL(PrimConsJac,ConsPrimJac)
          dQ_dUVolOuter(:,:,i,j,iLocSide,oo) = dQ_dUVolOuter(:,:,i,j,iLocSide,oo) + 0.5*etaBR2 * r(i,j) * &
                                               l_mp(oo,iLocSide)*Surf(i,j,iLocSide,iElem)*SurfVol_PrimJac

        END DO !i
      END DO !j
    END DO !oo
#endif
  END SELECT
END DO !iLocSide

END SUBROUTINE dQOuter
END MODULE MOD_Jac_br2
#endif /*PARABOLIC*/
