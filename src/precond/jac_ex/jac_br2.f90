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
!> This module contains routines required for the Jacobian of the br2 lifting scheme assuming non-conservative lifting.
!> Attention:
!> If br1 lifting scheme is chosen, etabr2 is set to 1. In br1 lifting the surface gradient at a specific side is depending on
!> UPrim in the volume and UPrim_face on all sides. For br2 it is depending only on UPrim in the volume and UPrim_face of the 
!> current face itself. For the preconditioner this additional dependency is neglected: No additional fillin
!> (in comparison to br2) is done. For Gauss-Lobatto nodes this approach gives correct resutls except for for the dofs at the
!> cornes as here the dependency of the dof on multiple surfaces is done seperately.
!===================================================================================================================================
MODULE MOD_Jac_br2
#if PARABOLIC
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

INTERFACE JacLifting_VolGrad
  MODULE PROCEDURE JacLifting_VolGrad
END INTERFACE

INTERFACE dQOuter
  MODULE PROCEDURE dQOuter
END INTERFACE

INTERFACE dQInner
  MODULE PROCEDURE dQInner
END INTERFACE

PUBLIC::FillJacLiftingFlux,JacLifting_VolGrad,Build_BR2_SurfTerms,dQOuter,dQInner
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
USE MOD_PreProc
USE MOD_Jac_ex_Vars               ,ONLY:JacLiftingFlux
USE MOD_GetBoundaryFlux_fd        ,ONLY:Lifting_GetBoundaryFlux_FD
USE MOD_Mesh_Vars                 ,ONLY:nBCSides,ElemToSide,S2V2
USE MOD_Mesh_Vars                 ,ONLY:NormVec,TangVec1,TangVec2,SurfElem,firstInnerSide,MortarType,MortarInfo,FS2M
USE MOD_DG_Vars                   ,ONLY:UPrim_master
USE MOD_Mesh_Vars                 ,ONLY:Face_xGP
USE MOD_Jac_Ex_MortarLifting      ,ONLY:Jacobian_MortarLifting
USE MOD_Jacobian                  ,ONLY: dConsdPrimTemp,dPrimTempdCons
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
INTEGER :: iVar,iLocSide,SideID,flip,p,q,pq(2)
REAL    :: mp
INTEGER :: nMortars,iMortarSide,iMortar
INTEGER :: sideMap(1:2,0:PP_N,0:PP_NZ)
REAL    :: DLiftingFluxMortar_DUInner(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,4)
REAL    :: JacLiftingFlux_tmp(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ)
REAL    :: PrimConsJac(PP_nVarPrim,PP_nVar)
REAL    :: ConsPrimJac(PP_nVar,PP_nVarPrim)
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
    IF (SideID.LT.firstInnerSide) THEN
      ! This is a (big) mortar side
      nMortars=MERGE(4,2,MortarType(1,SideID).EQ.1)
      ! Index in mortar side list
      iMortarSide=MortarType(2,SideID)
      DLiftingFluxMortar_DUInner = 0.
    ELSE
      nMortars = 1
      sideMap = S2V2(:,:,:,flip,iLocSide)
    END IF
    DO iMortar = 1, nMortars
      JacLiftingFlux(:,:,:,:,iLocSide)=0.
      IF (nMortars.GT.1) THEN
        SideID= MortarInfo(MI_SIDEID,iMortar,iMortarSide)
        sideMap = FS2M(:,:,:,MortarInfo(MI_FLIP,iMortar,iMortarSide))
      END IF
      IF(flip.EQ.0)THEN
        mp = -1.
      ELSE
        mp = 1.
      END IF
      DO iVar=1,PP_nVarPrim
        DO q=0,PP_NZ; DO p=0,PP_N
          pq(:) = sideMap(:,p,q)
          JacLiftingFlux(iVar,iVar,pq(1),pq(2),iLocSide) = mp*0.5*SurfElem(p,q,0,SideID)
        END DO; END DO ! p,q=0,PP_N
      END DO !iVar
      IF (nMortars.GT.1) THEN
        DO q=0,PP_NZ; DO p=0,PP_N
          pq(:) = sideMap(:,p,q)
          CALL dPrimTempdCons(UPrim_master(:,pq(1),pq(2),SideID),PrimConsJac)
          DLiftingFluxMortar_DUInner(:,:,p,q,iMortar) = MATMUL(JacLiftingFlux(:,:,p,q,iLocSide),PrimConsJac)
        END DO; END DO ! p,q=0,PP_N
      END IF
    END DO ! iMortar = 1, nMortars
    IF (nMortars.GT.1) THEN
      CALL Jacobian_MortarLifting(MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                                  DLiftingFluxMortar_DUInner,JacLiftingFlux_tmp(:,:,:,:))
      DO q=0,PP_NZ; DO p=0,PP_N
        pq(:)  = S2V2(:,p,q,flip,iLocSide)
        CALL dConsdPrimTemp(UPrim_master(:,pq(1),pq(2),ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),ConsPrimJac)
        JacLiftingFlux(:,:,p,q,iLocSide) = MATMUL(JacLiftingFlux_tmp(:,:,p,q),ConsPrimJac)
      END DO; END DO ! p,q=0,PP_N
    END IF
  END IF !SideID
END DO!iLocSide
END SUBROUTINE FillJacLiftingFlux

!===================================================================================================================================
!> Computes the Volume gradient Jacobian of the BR2 scheme dQprim/dUprim (Q= Grad U). This consists of the volume and surface
!> integral of the lifting procedure. This get's complicated since we interpolate the volume CONSERVATIVE variables onto the
!> surface, then do a ConsToPrim and then compute the surface flux in the lifting equation. Thus, this dependency needs to be taken
!> into account.
!> dQprim^vol(Uprim^vol,Uprim^surf)/dUprim^vol = dQprim^vol/dUprim^vol +  dQprim^vol/dUprim^surf * dUprim^surf/dUprim^vol
!>                                                   LiftingVolInt            LiftingSurfint
!> Taking the prolongation into account:
!> dUprim^surf/dUprim^vol = dUprim^surf/dUcons^surf * dUcons^surf/dUcons^vol * dUcons^vol/dUprim^vol
!>                            ConsToPrimSurf             ProlongToFace             PrimToConsVol
!===================================================================================================================================
SUBROUTINE JacLifting_VolGrad(dir,iElem,JacLifting)
! MODULES
USE MOD_PreProc
USE MOD_Jac_Ex_Vars        ,ONLY: LL_minus,LL_plus
USE MOD_Jac_Ex_Vars        ,ONLY: JacLiftingFlux 
USE MOD_DG_Vars            ,ONLY: D
#if (PP_NodeType==1)
USE MOD_DG_Vars            ,ONLY: UPrim
USE MOD_DG_Vars            ,ONLY: UPrim_master,UPrim_slave
#endif
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,sJ
#if PP_dim==3
USE MOD_Mesh_Vars          ,ONLY: Metrics_hTilde
#endif
USE MOD_Mesh_Vars          ,ONLY: ElemToSide,S2V2,Normvec
USE MOD_Jacobian           ,ONLY: dConsdPrimTemp,dPrimTempdCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                           :: iElem   !< current element index
INTEGER,INTENT(IN)                           :: dir     !< current physical direction
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                             :: JacLifting(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
                                                !< Jacobian of volume gradients in direction dir w.r.t. primitive volume solution
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: i,j,k,ll
INTEGER                                      :: iVar
REAL,DIMENSION(1:PP_nVarPrim,1:PP_nVarPrim)  :: delta,SurfVol_PrimJac_plus,SurfVol_PrimJac_minus
#if (PP_NodeType==1)
REAL,DIMENSION(1:PP_nVar    ,1:PP_nVarPrim)  :: ConsPrimJac
REAL,DIMENSION(1:PP_nVarPrim,1:PP_nVar    )  :: PrimConsJac
#endif
INTEGER                                      :: pq_p(2),pq_m(2),SideID_p,SideID_m,flip
INTEGER                                      :: sideIDs(2,3),sideIDPlus,sideIDMinus,refDir
!===================================================================================================================================
! fill volume Jacobian of lifting flux (flux is primitive state vector -> identity matrix)
delta=0.
DO iVar=1,PP_nVarPrim
  delta(iVar,iVar)=1.
END DO

! Array that tells us which sides to work with for each direction
sideIDs = RESHAPE((/XI_MINUS,XI_PLUS,ETA_MINUS,ETA_PLUS,ZETA_MINUS,ZETA_PLUS/),(/2,3/))

JacLifting=0.
! ll index: Loop variable along the (two or) three one-dimensional lines (XI, ETA, ZETA). We derive w.r.t. the DOFs along those
! lines.
DO ll=0,PP_N
  ! i,j,k index: The DOF of the volume gradient we are deriving
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        ! Loop over XI,ETA,ZETA direction
        DO refDir=1,PP_dim
          sideIDMinus = sideIDs(1,refDir)
          sideIDPlus  = sideIDs(2,refDir)
#if (PP_NodeType==1)
          ! dUcons^vol/dUprim^vol
          SELECT CASE(refDir)
          CASE(1)
            CALL dConsdPrimTemp(UPrim(:,ll,j,k,iElem),ConsPrimJac)
          CASE(2)
            CALL dConsdPrimTemp(UPrim(:,i,ll,k,iElem),ConsPrimJac)
          CASE(3)
            CALL dConsdPrimTemp(UPrim(:,i,j,ll,iElem),ConsPrimJac)
          END SELECT
#endif
          !!!!!!! Minux side !!!!!!
          SideID_m=ElemToSide(E2S_SIDE_ID,sideIDMinus,iElem)
          flip=ElemToSide(    E2S_FLIP   ,sideIDMinus,iElem)
          SELECT CASE(refDir)
          CASE(1)
            pq_m=S2V2(:,j,k,flip         ,sideIDMinus)
          CASE(2)
            pq_m=S2V2(:,i,k,flip         ,sideIDMinus)
          CASE(3)
            pq_m=S2V2(:,i,j,flip         ,sideIDMinus)
          END SELECT
#if (PP_NodeType==1)
          ! dUprim^surf/dUcons^surf
          IF(flip.EQ.0)THEN
            CALL dPrimTempdCons(UPrim_master(:,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
          ELSE
            CALL dPrimTempdCons(UPrim_slave( :,pq_m(1),pq_m(2),SideID_m),PrimConsJac)
          END IF
          ! Since the ProlongToFace is independant of iVar, we can combine the two EOS calls into one
          SurfVol_PrimJac_minus = MATMUL(PrimConsJac,ConsPrimJac)
#elif (PP_NodeType==2)
          ! For GL, since the surface points are directly volume points, the dependency reduces to the identity matrix
          SurfVol_PrimJac_minus = delta
#endif
          !!!!!!! Plus side !!!!!!
          SideID_p =ElemToSide(E2S_SIDE_ID,sideIDPlus ,iElem)
          flip=ElemToSide(     E2S_FLIP   ,sideIDPlus ,iElem)
          SELECT CASE(refDir)
          CASE(1)
            pq_p=S2V2(:,j,k,flip          ,sideIDPlus)
          CASE(2)
            pq_p=S2V2(:,i,k,flip          ,sideIDPlus)
          CASE(3)
            pq_p=S2V2(:,i,j,flip          ,sideIDPlus)
          END SELECT
#if (PP_NodeType==1)
          ! dUprim^surf/dUcons^surf
          IF(flip.EQ.0)THEN
            CALL dPrimTempdCons(UPrim_master(:,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
          ELSE
            CALL dPrimTempdCons(UPrim_slave( :,pq_p(1),pq_p(2),SideID_p),PrimConsJac)
          END IF
          ! Since the ProlongToFace is independant of iVar, we can combine the two EOS calls into one
          SurfVol_PrimJac_plus = MATMUL(PrimConsJac,ConsPrimJac)
#elif (PP_NodeType==2)
          ! For GL, since the surface points are directly volume points, the dependency reduces to the identity matrix
          SurfVol_PrimJac_plus = delta
#endif
          ! Assemble lifting Jacobian with surface and volume contribution
          SELECT CASE(refDir)
          CASE(1)
            JacLifting(:,:,i,j,k,ll,1) = sJ(i,j,k,iElem,0)*(D(i,ll)*Metrics_fTilde(dir,ll,j,k,iElem,0)*delta(:,:)& !LiftingVolInt
                                          +LL_plus( i,ll)*NormVec(dir,pq_p(1),pq_p(2),0,SideID_p)*               & !LiftingSurfInt
                                           MATMUL(JacLiftingFlux(:,:,j,k,XI_PLUS),SurfVol_PrimJac_plus(:,:))     & !Prolongation
                                          +LL_minus(i,ll)*NormVec(dir,pq_m(1),pq_m(2),0,SideID_m)*               & !LiftingSurfInt
                                           MATMUL(JacLiftingFlux(:,:,j,k,XI_MINUS),SurfVol_PrimJac_minus(:,:)))    !Prolongation
          CASE(2)
            JacLifting(:,:,i,j,k,ll,2) = sJ(i,j,k,iElem,0)*(D(j,ll)*Metrics_gTilde(dir,i,ll,k,iElem,0)*delta(:,:)&
                                          +LL_plus( j,ll)*NormVec(dir,pq_p(1),pq_p(2),0,SideID_p)*               &
                                           MATMUL(JacLiftingFlux(:,:,i,k,ETA_PLUS),SurfVol_PrimJac_plus(:,:))    &
                                          +LL_minus(j,ll)*NormVec(dir,pq_m(1),pq_m(2),0,SideID_m)*               &
                                           MATMUL(JacLiftingFlux(:,:,i,k,ETA_MINUS),SurfVol_PrimJac_minus(:,:)))
#if PP_dim==3
          CASE(3)
            JacLifting(:,:,i,j,k,ll,3) = sJ(i,j,k,iElem,0)*(D(k,ll)*Metrics_hTilde(dir,i,j,ll,iElem,0)*delta(:,:)&
                                          +LL_plus( k,ll)*NormVec(dir,pq_p(1),pq_p(2),0,SideID_p)*               &
                                           MATMUL(JacLiftingFlux(:,:,i,j,ZETA_PLUS),SurfVol_PrimJac_plus(:,:))   &
                                          +LL_minus(k,ll)*NormVec(dir,pq_m(1),pq_m(2),0,SideID_m)*               &
                                           MATMUL(JacLiftingFlux(:,:,i,j,ZETA_MINUS),SurfVol_PrimJac_minus(:,:)))
#endif
          END SELECT
        END DO ! refDir=1,PP_dim
      END DO !i
    END DO !j
  END DO !k
END DO !ll
END SUBROUTINE JacLifting_VolGrad


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
CALL StartReceiveMPIData(R_Plus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RPlus(:,RECV),SendID=2)   ! Recv MINE/Geo: slave->master
CALL StartSendMPIData(   R_Plus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RPlus(:,SEND),SendID=2)   ! SEND YOUR/Geo: slave->master
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_RPlus) 

CALL StartReceiveMPIData(R_Minus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RMinus(:,RECV),SendID=1) ! Recv YOUR/Geo: master->slave
CALL StartSendMPIData(R_Minus,3*(PP_N+1)**(PP_dim-1),1,nSides,MPIRequest_RMinus(:,SEND),SendID=1)    ! SEND MINE/Geo: master->slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_RMinus) 
#endif

END SUBROUTINE Build_BR2_SurfTerms

!===================================================================================================================================
!> Contains the required derivatives of the BR2 scheme for my own element of both the volume and the surface integral part of the
!> lifting.
!> ONLY THE DERIVATIVE OF Q_INNER !!!!!
!> Computation is done for one element!
!===================================================================================================================================
SUBROUTINE dQInner(dir,iElem,dQInner_dUVol,dQVol_dU)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                 ,ONLY: ElemToSide,S2V2
USE MOD_Jac_Ex_Vars               ,ONLY: l_mp
USE MOD_Jac_Ex_Vars               ,ONLY: R_Minus,R_Plus
USE MOD_Jac_Ex_Vars               ,ONLY: JacLiftingFlux 
USE MOD_Mesh_Vars                 ,ONLY: Metrics_fTilde,Metrics_gTilde,sJ
#if PP_dim==3
USE MOD_Mesh_Vars                 ,ONLY: Metrics_hTilde
#endif
USE MOD_DG_Vars                   ,ONLY: D
USE MOD_DG_Vars                   ,ONLY: UPrim_master,UPrim_slave
#if PP_Lifting==2
USE MOD_Lifting_Vars              ,ONLY: etaBR2
#endif
USE MOD_Jacobian                  ,ONLY: dConsdPrimTemp,dPrimTempdCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: dir   !< considered direction (1,2,3)
INTEGER,INTENT(IN) :: iElem !< considered element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#if PP_dim == 3
REAL,INTENT(OUT) :: dQInner_dUVol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N) !< Jacobian of surface integral of the lifting w.r.t.
                                                                               !< conservative solution, only from my element!
#else
REAL,INTENT(OUT) :: dQInner_dUVol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif
REAL,INTENT(OUT) :: dQVol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)              !< Jacobian of volume integral of the lifting w.r.t.
                                                                               !< conservative solution
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iLocSide,p,q,i,j,k,mm,nn,ll
#if PP_dim == 3
INTEGER           :: oo
#endif
INTEGER           :: SideID,Flip,jk(2)
REAL              :: r(0:PP_N,0:PP_NZ),mp
REAL              :: dQInner_dUVol_loc(0:PP_N,0:PP_NZ,0:PP_N)
REAL              :: UPrim_face(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL              :: PrimConsJac(PP_nVarPrim,1:PP_nVar,0:PP_N,0:PP_NZ)
#if PP_Lifting==1
REAL,PARAMETER    :: etaBR2=1.
#endif
!===================================================================================================================================
dQInner_dUVol=0.

! Derivative of the volume integral of the lifting procedure w.r.t. the lifting variables - since this consists of only the local
! gradients for BR2, this reduced to the D matrix including the metrics.
DO ll=0,PP_N
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        dQVol_dU(i,j,k,ll,1) =  sJ(i,j,k,iElem,0) * D(i,ll)*Metrics_fTilde(dir,ll,j,k,iElem,0)
        dQVol_dU(i,j,k,ll,2) =  sJ(i,j,k,iElem,0) * D(j,ll)*Metrics_gTilde(dir,i,ll,k,iElem,0)
#if PP_dim==3
        dQVol_dU(i,j,k,ll,3) =  sJ(i,j,k,iElem,0) * D(k,ll)*Metrics_hTilde(dir,i,j,ll,iElem,0)
#endif
      END DO !i
    END DO !j
  END DO !k
END DO !ll

! Computation of the derivative of the surface integral part
#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  ! Flip either slave or master solution into my own volume system, mp considers that master and slave fluxes have opposing signs
  ! r considers the metric terms on the surface and the integration => this gives us the actual surface integral.
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

  dQInner_dUVol_loc=0.
  SELECT CASE(iLocSide)

  ! Now compute the derivatives for the sides. We additionaly take the derivative of the prolongation into account by using the
  ! array l_mp, and the dependency of the primitive surface variables from the conservative ones
  CASE(XI_MINUS,XI_PLUS)
    DO mm=0,PP_N
      DO k=0,PP_NZ
        DO j=0,PP_N
          dQInner_dUVol_loc(j,k,mm)   = dQInner_dUVol_loc(j,k,mm) + etaBR2 * r(j,k) * l_mp(mm,iLocSide)
          CALL dPrimTempdCons(UPrim_face(:,j,k),PrimConsJac(:,:,j,k))
        END DO !j
      END DO !k
    END DO !mm
  CASE(ETA_MINUS,ETA_PLUS)
    DO nn=0,PP_N
      DO k=0,PP_NZ
        DO i=0,PP_N
          dQInner_dUVol_loc(i,k,nn) = dQInner_dUVol_loc(i,k,nn) + etaBR2 * r(i,k) * l_mp(nn,iLocSide)
          CALL dPrimTempdCons(UPrim_face(:,i,k),PrimConsJac(:,:,i,k))
        END DO !i
      END DO !k
    END DO !nn
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
    DO oo=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          dQInner_dUVol_loc(i,j,oo) = dQInner_dUVol_loc(i,j,oo) + etaBR2 * r(i,j) * l_mp(oo,iLocSide)
          CALL dPrimTempdCons(UPrim_face(:,i,j),PrimConsJac(:,:,i,j))
        END DO !i
      END DO !j
    END DO !oo
#endif
  END SELECT
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO mm=0,PP_N
        dQInner_dUVol(:,:,j,k,iLocSide,mm)=mp*MATMUL(JacLiftingFlux(:,:,j,k,iLocSide),PrimConsJac(:,:,j,k))* &
                                           dQInner_dUVol_loc(j,k,mm)
      END DO
    END DO !p
  END DO !q
END DO !iLocSide

END SUBROUTINE dQInner

!===================================================================================================================================
!> Contains the required derivatives of the BR2 scheme for the neighbouring element: only the surface integral part of the lifting.
!> ONLY THE DERIVATIVE OF Q_OUTER !!!!!
!> computation is done for one element!
!===================================================================================================================================
SUBROUTINE dQOuter(dir,iElem,dQOuter_dUVol)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                 ,ONLY: ElemToSide,S2V2,nBCSides,SurfElem,firstInnerSide,MortarType,MortarInfo,FS2M
USE MOD_Jac_Ex_Vars               ,ONLY: l_mp
USE MOD_Jac_Ex_Vars               ,ONLY: R_Minus,R_Plus
#if PP_Lifting==2
USE MOD_Lifting_Vars              ,ONLY: etaBR2
#endif
USE MOD_DG_Vars                   ,ONLY: UPrim_master,UPrim_slave
USE MOD_Jacobian                  ,ONLY: dConsdPrimTemp,dPrimTempdCons
#if FV_ENABLED
USE MOD_FV_Vars                   ,ONLY: FV_Elems_master,FV_Elems_slave
#endif
USE MOD_Jac_Ex_MortarScalar       ,ONLY: Jacobian_MortarScalar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: dir   !< considered direction (1,2,3)
INTEGER,INTENT(IN) :: iElem !< considered element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#if PP_dim == 3
REAL,INTENT(OUT) :: dQOuter_dUVol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)!< Jacobian of surface integral of the lifting w.r.t.
                                                                              !< conservative solution, only from neigbour element!
#else
REAL,INTENT(OUT) :: dQOuter_dUVol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iLocSide,p,q,i,j,k,mm,nn
#if PP_dim == 3
INTEGER           :: oo
#endif
INTEGER           :: SideID,Flip,jk(2)
REAL              :: r(0:PP_N,0:PP_NZ)
REAL              :: UPrim_face(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL              :: PrimConsJac(PP_nVarPrim,PP_nVar)
#if PP_Lifting==1
REAL,PARAMETER    :: etaBR2=1.
#endif
#if FV_ENABLED
INTEGER           :: FVSide
#endif
REAL              :: DQOuter_DUSurf(1,1,0:PP_N,0:PP_NZ)
REAL              :: DQOuterMortar_DUSurf(1,1,0:PP_N,0:PP_NZ,4)
INTEGER           :: nMortars,iMortarSide,iMortar,pq(2)
INTEGER           :: sideMap(1:2,0:PP_N,0:PP_NZ)
!===================================================================================================================================
dQOuter_dUVol=0.

! Computation of the derivative of the surface integral part
#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  IF(SideID.LE.nBCSides) CYCLE  !for boundary conditions, dQ_dUVol=0.
  IF (SideID.LT.firstInnerSide) THEN
    ! This is a (big) mortar side
    nMortars=MERGE(4,2,MortarType(1,SideID).EQ.1)
    ! Index in mortar side list
    iMortarSide=MortarType(2,SideID)
    DQOuterMortar_DUSurf = 0.
  ELSE
    nMortars = 1
  END IF
  ! Flip either slave or master solution into my own volume system
  ! r considers the metric terms on the surface and the integration => this gives us the actual surface integral.
  Flip  =ElemToSide(E2S_FLIP,iLocSide,iElem)
  IF(Flip.EQ.0)THEN
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Plus(dir,p,q,SideID)
        UPrim_face(:,jk(1),jk(2)) = UPrim_master(:,p,q,SideID)
      END DO !p
    END DO !q
#if FV_ENABLED
    FVSide = FV_Elems_slave(SideID)
#endif
  ELSE
    DO q=0,PP_NZ
      DO p=0,PP_N
        jk(:)=S2V2(:,p,q,Flip,iLocSide)
        r(jk(1),jk(2))=R_Minus(dir,p,q,SideID)
        UPrim_face(:,jk(1),jk(2)) = UPrim_slave(:,p,q,SideID)
      END DO !p
    END DO !q
#if FV_ENABLED
    FVSide = FV_Elems_master(SideID)
#endif
  END IF !Flip=0

  ! First, compute the derivatives of the gradients w.r.t. the surface DOFs (each gradient only to the solution at the same point)
  DO iMortar=1,nMortars
    IF (nMortars.GT.1) THEN
      ! We overwrite the SideID with the ID of the SMALL mortar side!
      SideID = MortarInfo(MI_SIDEID,iMortar,iMortarSide)
      ! For mortar sides, we transform the side-based Jacobians of the small sides into the side system of the big side. The
      ! transformation into the volume system is done later.
      sideMap = FS2M(:,:,:,MortarInfo(MI_FLIP,iMortar,iMortarSide))
    ELSE
      ! For non-mortar sides, directly transform into volume system
      sideMap = S2V2(:,:,:,flip,iLocSide)
    END IF
    DO q=0,PP_NZ
      DO p=0,PP_N
        pq = sideMap(:,p,q)
        dQOuter_dUSurf(:,:,p,q) = 0.5*etaBR2*SurfElem(pq(1),pq(2),0,SideID)
      END DO !p
    END DO !q
    IF (nMortars.GT.1) THEN
      ! Store the Jacobians for each small side, later combine them into the Jacobian for the big side
      DQOuterMortar_DUSurf(:,:,:,:,iMortar) = dQOuter_dUSurf
    END IF
  END DO ! iMortar
  ! If this is a big mortar side, we combine the Jacobians from the small sides to the Jacobians on the big side and flip into the
  ! volume system
  IF (nMortars.GT.1) THEN
    CALL Jacobian_MortarScalar(MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                               DQOuterMortar_DUSurf,dQOuter_dUSurf)
    !dQOuter_dUSurf = 0.

  ELSEIF (MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)).LT.0) THEN 
    ! This is a small mortar. Since the neighbouring element thus is a big mortar, the surface integral in the BR2 scheme has a
    ! different derivative than a conforming element. We can't compute that influence at the moment, thus set it to zero.
    ! We would need to know which small mortar we are (1,2,3,4 and in what mortar type) to select the right dependencies, but this
    ! information is not stored anywhere.
    dQOuter_dUSurf = 0.
  END IF


  ! Now compute the derivatives for the sides w.r.t. volume DOFs. The influence of the surface integral for all volume DOFs is
  ! stored in the array l_mp, the influence of the metric terms in r. Also take the dependency of the primitive surface variables
  ! from the conservative ones into account.
  SELECT CASE(iLocSide)
  CASE(XI_MINUS,XI_PLUS)
#if FV_ENABLED
    IF(FVSide.EQ.0)THEN
#endif
      DO mm=0,PP_N
        DO k=0,PP_NZ
          DO j=0,PP_N
            CALL dPrimTempdCons(UPrim_face(:,j,k),PrimConsJac)
            dQOuter_dUVol(:,:,j,k,iLocSide,mm) = r(j,k)*l_mp(mm,iLocSide)*PrimConsJac*DQOuter_dUSurf(1,1,j,k)
          END DO !j
        END DO !k
      END DO !mm
#if FV_ENABLED
    END IF
#endif
  CASE(ETA_MINUS,ETA_PLUS)
#if FV_ENABLED
    IF(FVSide.EQ.0)THEN
#endif
      DO nn=0,PP_N
        DO k=0,PP_NZ
          DO i=0,PP_N
            CALL dPrimTempdCons(UPrim_face(:,i,k),PrimConsJac)
            dQOuter_dUVol(:,:,i,k,iLocSide,nn) = r(i,k)*l_mp(nn,iLocSide)*PrimConsJac*DQOuter_dUSurf(1,1,i,k)
          END DO !i
        END DO !k
      END DO !nn
#if FV_ENABLED
    END IF
#endif
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
#if FV_ENABLED
    IF(FVSide.EQ.0)THEN
#endif
      DO oo=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            CALL dPrimTempdCons(UPrim_face(:,i,j),PrimConsJac)
            dQOuter_dUVol(:,:,i,j,iLocSide,oo) = r(i,j)*l_mp(oo,iLocSide)*PrimConsJac*DQOuter_dUSurf(1,1,i,j)
          END DO !i
        END DO !j
      END DO !oo
#if FV_ENABLED
    END IF
#endif
#endif
  END SELECT
END DO !iLocSide

END SUBROUTINE dQOuter
#endif /*PARABOLIC*/
END MODULE MOD_Jac_br2
