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
!> Contains the computation of the local Jacobian of the surface integral
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
INTERFACE JacSurfInt
  MODULE PROCEDURE JacSurfInt
END INTERFACE

PUBLIC::JacSurfInt
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Contains the computation of the local Jacobian of the surface integral.
!> Computation is done for one element!
!>
!> Dependency of the surface integral:
!>---------------------------------------------------------------------------------------------------- 
!> DF(F^v,F^h)/DU |_surfInt = DF^h/DU |_surfInt + DF^v/DU |_surfInt
!> 
!> DF^h/DU |_surfInt = DF^h/DF^h_surf * DF^h_surf/DU
!> DF^h/DU |_surfInt = DF^h/DF^h_surf * DF^h_surf/DU_surf * DU_surf/DU
!>                                       => Df_DUinner
!> 
!> DF^v/DU |_surfInt = DF^v/DF^v_surf * DF^v_surf/DU
!> DF^v/DU |_surfInt = DF^v/DF^v_surf * DF^v(U,Q(U))_surf/DU
!> DF^v/DU |_surfInt = DF^v/DF^v_surf * (DF^v(U,Q(U))_surf/DU + DF^v(U,Q(U))_surf/DQ * DQ/DU)
!> DF^v/DU |_surfInt = DF^v/DF^v_surf * (DF^v(U,Q(U))_surf/DU_surf * DU_surf/DU + DF^v(U,Q(U))_surf/DQ_surf * DQ_surf/DU)
!> 
!> Derivatives of diffusive numerical flux function w.r.t. solution and gradients:
!> DF^v(U,Q(U))_surf/DU_surf = D(0.5*(F^v_L+ F^v_R)/DU_surf = 0.5*D((F^v_L)/DU_surf
!>                                                                 => Df_DUInner
!> DF^v(U,Q(U))_surf/DQ_surf = D(0.5*(F^v_L+ F^v_R)/DQ_surf = 0.5*(D(F^v_L)/DQ_surf_L + D(F^v_R)/DQ_surf_R)
!>                                                                   =>Df_dQxInner        =>Df_DQxOuter
!> 
!> Derivative of lifting procedure or FV gradient reconstruction:
!> DQ_surf_L/DU = DQ_surf_L|_VolInt/DU + DQ_surf_L|_SurfInt/DU
!> DQ_surf_L/DU = DQ_surf_L|_VolInt/DU_prim * DU_prim/DU + DQ_surf_L|_SurfInt/DU
!> DQ_surf_L/DU = DQ_surf_L|_VolInt/DU_prim * DU_prim/DU + DQ_surf_L|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU
!>                       => dQxvol_dU        =>PrimConsJac       => dQxInner              => PrimConsJac     => Prolongation
!> 
!> DQ_surf_R/DU = DQ_surf_R|_VolInt/DU + DQ_surf_R|_SurfInt/DU
!> 
!> DQ_surf_R/DU = DQ_surf_R|_VolInt/DU + DQ_surf_R|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU
!>                          =0                => dQxOuter              =>PrimConsJac       => Prolongation
!>---------------------------------------------------------------------------------------------------- 
!>
!> The surface integral depends both on the hyperbolic and the viscous fluxes.
!> For the hyperbolic part we need:
!>  * Derivative of the Riemann solver (computed using FD approach) w.r.t. the surface solution, computed in normal system and
!>    pre-multiplied with SurfElem
!>  * Application of surface integral over those fluxes, and dependency of surface solution w.r.t. volume DOFs => both considered in
!>    LL_plus and LL_minus arrays
!> For the diffusive part we need:
!>  * Analytical derivation of surface flux w.r.t. surface solution
!>  * Analytical derivation of surface flux w.r.t. surface gradients
!>  * Transform those into the normal system and pre-multiply with SurfElem
!>  * Add dependency on solution to hyperbolic part
!>  * Now for the dependeny on the surface gradients. Since the surface flux function of the lifting in my neighbouring cells 
!>    depends on the solution in my cell, the gradients from my neighbouring cells also depend on my solution! Since the viscous
!>    flux then depends on those gradients, we get an additional dependency from the outer gradients.
!>    * dependency of the lifting volume integral w.r.t. volume DOFs
!>    * dependency of the inner surface gradients w.r.t. primitive volume DOFs
!>    * dependency of the outer surface gradients w.r.t. primitive volume DOFs
!>---------------------------------------------------------------------------------------------------- 
!> After calculating the required derivatives of the flux (Df_DUinner,Df_dQxInner,Df_DQxOuter) the assembling is called either for
!> the DG case or the FV case.
!===================================================================================================================================
SUBROUTINE JacSurfInt(t,BJ,iElem)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars                   ,ONLY: U_master,U_slave,UPrim_master,UPrim_slave
USE MOD_Mesh_Vars                 ,ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars                 ,ONLY: nBCSides,ElemToSide,S2V2,firstInnerSide,MortarType,MortarInfo,FS2M
USE MOD_Mesh_Vars                 ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_Implicit_Vars             ,ONLY: nDOFVarElem
USE MOD_Riemann_FD                ,ONLY: Riemann_FD
USE MOD_GetBoundaryFlux_FD        ,ONLY: GetBoundaryFlux_FD
#if PARABOLIC
USE MOD_Jac_Ex_br2                ,ONLY: dQOuter,dQInner
USE MOD_Jacobian                  ,ONLY: dPrimTempdCons
USE MOD_DG_Vars                   ,ONLY: nDOFFace
USE MOD_Precond_Vars              ,ONLY: HyperbolicPrecond
USE MOD_Lifting_Vars              ,ONLY: gradUx_master,gradUx_slave
USE MOD_Lifting_Vars              ,ONLY: gradUy_master,gradUy_slave
USE MOD_Lifting_Vars              ,ONLY: gradUz_master,gradUz_slave
USE MOD_Jacobian                  ,ONLY: EvalDiffFluxJacobian,EvalFluxGradJacobian
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars             ,ONLY: muSGS_master,muSGS_slave
#endif /*EDDYVISCOSITY*/
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV_Vars                   ,ONLY: FV_Elems_Sum
USE MOD_FV_Vars                   ,ONLY: FV_Elems_master,FV_Elems_slave,FV_Elems
#if FV_RECONSTRUCT
USE MOD_GetBoundaryFlux           ,ONLY: GetBoundaryState
#endif
#endif
USE MOD_Jac_Ex_MortarU            ,ONLY: Jacobian_MortarU
#if PARABOLIC
USE MOD_Jac_Ex_MortarGrad         ,ONLY: Jacobian_MortarGrad
#endif
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
INTEGER                                                   :: iLocSide,SideID,Flip,FVElem,FVSide,FVSum
REAL                                                      :: signum
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ)                    :: USideL,USideR
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_NZ)                :: UPrimSideL,UPrimSideR
#if PP_dim == 3
REAL                                                      :: Df_DUInner(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2,6)
#else
REAL                                                      :: Df_DUInner(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2,2:5)
#endif
#if PARABOLIC
INTEGER                                                   :: jk(2),p,q,i
REAL,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_NZ)             :: gradUxSideL,gradUySideL,gradUzSideL,gradUxSideR,gradUySideR,gradUzSideR
REAL,DIMENSION(PP_nVar,PP_nVarPrim,0:PP_N,0:PP_NZ)        :: fJacQx,fJacQz,fJacQy,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz
REAL,DIMENSION(PP_nVar,PP_nVar    ,0:PP_N,0:PP_NZ,2)      :: fJac,gJac,hJac
#if EDDYVISCOSITY
REAL,DIMENSION(1,0:PP_N,0:PP_NZ)                          :: muSGSL,muSGSR
#endif
#if PP_dim==3
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_N,1:6) :: Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter, &
                                                             Df_dQzInner,Df_dQzOuter
#else
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_NZ,2:5):: Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter
#endif /*PP_dim*/
#endif /*PARABOLIC*/
! Mortars
INTEGER                                                   :: nMortars,iMortarSide,iMortar,flip_small
INTEGER                                                   :: sideMap(1:2,0:PP_N,0:PP_NZ)
REAL                                                      :: DfMortar_DUInner(PP_nVar,PP_nVar,0:PP_N,0:PP_NZ,2,4)
#if PARABOLIC
REAL,DIMENSION(PP_nVar,PP_nVarPrim,0:PP_N,0:PP_NZ,4)      :: DfMortar_DQxInner,DfMortar_DQyInner,DfMortar_DQxOuter,DfMortar_DQyOuter
#if PP_dim==3
REAL,DIMENSION(PP_nVar,PP_nVarPrim,0:PP_N,0:PP_NZ,4)      :: DfMortar_DQzInner,DfMortar_DQzOuter
#endif /*PP_dim*/
#endif /*PARABOLIC*/
!===================================================================================================================================
#if FV_ENABLED
FVElem = FV_Elems(iElem)
#else
FVElem = 0
#endif

#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
  Flip=ElemToSide(E2S_FLIP,iLocSide,iElem)
#if FV_ENABLED
  FVSide = MAX(FV_Elems_master(SideID),FV_Elems_slave(SideID))
  FVSum  = FV_Elems_Sum(SideID)
#else
  FVSide = 0
  FVSum  = 0
#endif

  IF (SideID.GT.nBCSides) THEN !InnerSide
    IF ((SideID.LT.firstInnerSide).OR.((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN
      ! This is a (big) mortar side
      nMortars=MERGE(4,2,MortarType(1,SideID).EQ.1)
      ! Index in mortar side list
      iMortarSide=MortarType(2,SideID)
      DfMortar_DUInner = 0.
    ELSE
      nMortars = 1
      sideMap = S2V2(:,:,:,flip,iLocSide)
      flip_small = flip
    END IF
    ! For mortars: Loop over small mortar sides
    DO iMortar=1,nMortars
      IF (nMortars.GT.1) THEN
        ! We overwrite the SideID with the ID of the SMALL mortar side!
        SideID = MortarInfo(MI_SIDEID,iMortar,iMortarSide)
        ! For mortar sides, we transform the side-based Jacobians of the small sides into the side system of the big side. The
        ! transformation into the volume system is done later, while we do it in Riemann_FD for non-mortar sides.
        flip_small = MortarInfo(MI_FLIP,iMortar,iMortarSide)
        sideMap = FS2M(:,:,:,flip_small)
      END IF

      ! We always want to take the derivative w.r.t. the DOF on our own side. Depending on whether we are master or slave, this means
      ! we need to input different arrays into our routines that calculate the derivative always w.r.t. U_L!
      IF (flip_small.EQ.0) THEN
        UPrimSideL = UPrim_master(:,:,:,SideID);USideL = U_master(:,:,:,SideID)
        UPrimSideR = UPrim_slave( :,:,:,SideID);USideR = U_slave( :,:,:,SideID)
#if PARABOLIC
        gradUxSideL = gradUx_master(:,:,:,SideID);gradUxSideR = gradUx_slave(:,:,:,SideID)
        gradUySideL = gradUy_master(:,:,:,SideID);gradUySideR = gradUy_slave(:,:,:,SideID)
        gradUzSideL = gradUz_master(:,:,:,SideID);gradUzSideR = gradUz_slave(:,:,:,SideID)
#if EDDYVISCOSITY
        muSGSL = muSGS_master(:,:,:,SideID)
        muSGSR = muSGS_slave(:,:,:,SideID)
#endif
#endif
        signum = 1.
      ELSE
        UPrimSideL = UPrim_slave( :,:,:,SideID);USideL = U_slave( :,:,:,SideID)
        UPrimSideR = UPrim_master(:,:,:,SideID);USideR = U_master(:,:,:,SideID)
#if PARABOLIC
        gradUxSideL = gradUx_slave(:,:,:,SideID);gradUxSideR = gradUx_master(:,:,:,SideID)
        gradUySideL = gradUy_slave(:,:,:,SideID);gradUySideR = gradUy_master(:,:,:,SideID)
        gradUzSideL = gradUz_slave(:,:,:,SideID);gradUzSideR = gradUz_master(:,:,:,SideID)
#if EDDYVISCOSITY
        muSGSL = muSGS_slave(:,:,:,SideID)
        muSGSR = muSGS_master(:,:,:,SideID)
#endif
#endif
        signum = -1.
      END IF
      ! d(f*_ad)_jk/dU_master_jk, with SurfaceIntegral already considered!
      CALL Riemann_FD(Df_DUInner(:,:,:,:,:,iLocSide),USideL,USideR,UPrimSideL,UPrimSideR,                                                &
                                 signum*NormVec(:,:,:,FVSide,SideID),signum*tangVec1(:,:,:,FVSide,SideID),tangVec2(:,:,:,FVSide,SideID), &
                                 SurfElem(:,:,FVSide,SideID),sideMap,FVSum,FVElem,FVSide)
#if PARABOLIC
      IF (.NOT.(HyperbolicPrecond)) THEN 
        ! Call the analytical Jacobian (viscous flux w.r.t. conservative variables)
        ! d(f_diff)_jk/dU_master_jk
        CALL EvalDiffFluxJacobian(nDOFFace,USideL,UPrimSideL,gradUxSideL,gradUySideL,gradUzSideL, &
                                  fJac(:,:,:,:,1),gJac(:,:,:,:,1),hJac(:,:,:,:,1)                 &
#if EDDYVISCOSITY
                                 ,muSGSL                                                          &
#endif
                                  )
#if FV_ENABLED
        ! Call the analytical Jacobian (viscous flux w.r.t. conservative variables), but for the neighbouring DOF
        CALL EvalDiffFluxJacobian(nDOFFace,USideR,UPrimSideR,gradUxSideR,gradUySideR,gradUzSideR, &
                                  fJac(:,:,:,:,2),gJac(:,:,:,:,2),hJac(:,:,:,:,2)                 &
#if EDDYVISCOSITY
                                 ,muSGSR                                                          &
#endif
      
                                 )
#endif
        ! Derivative of the diffusive flux with respect to gradU_L
        CALL EvalFluxGradJacobian(nDOFFace,USideL,UPrimSideL, &
                                  fJacQx,fJacQy,fJacQz,       &
                                  gJacQx,gJacQy,gJacQz,       &
                                  hJacQx,hJacQy,hJacQz        &
#if EDDYVISCOSITY
                                 ,muSGSL                      &
#endif
                                  )
        DO q=0,PP_NZ
          DO p=0,PP_N
            jk(:)=sideMap(:,p,q)
            ! BR1/2 use central fluxes for the viscous flux, so f*_diff = 0.5*(f_diff^L+f_diff^R).
            ! Transform the diffusive flux Jacobians into the normal system for the suface integral, and pre-multiply with 
            ! the SurfElem here (already done in Riemann_FD for the hyperbolic parts).
            DO i=1,FVElem+1
              ! Direct dependency of the viscous flux from the solution variables on the side.
              Df_DUInner(:,:,jk(1),jk(2),i,iLocSide)= Df_DUInner(:,:,jk(1),jk(2),i,iLocSide) +                  &
                                                      0.5*( fJac(:,:,p,q,i)*signum*NormVec(1,p,q,FVSide,SideID) &
                                                           +gJac(:,:,p,q,i)*signum*NormVec(2,p,q,FVSide,SideID) &
#if PP_dim==3
                                                           +hJac(:,:,p,q,i)*signum*NormVec(3,p,q,FVSide,SideID) &
#endif
                                                           )*SurfElem(p,q,FVSide,SideID)
            END DO
            ! Df_dQxInner = d(f,g,h)_diff/dQ_x * NormVec
            ! Dependencies of the surface flux w.r.t. the inner surface gradients in each direction
            Df_dQxInner(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQx(:,:,p,q)*signum*NormVec(1,p,q,FVSide,SideID) &
                                                         +gJacQx(:,:,p,q)*signum*NormVec(2,p,q,FVSide,SideID) &
#if PP_dim==3
                                                         +hJacQx(:,:,p,q)*signum*NormVec(3,p,q,FVSide,SideID) &
#endif
                                                        )*SurfElem(p,q,FVSide,SideID)
            Df_dQyInner(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQy(:,:,p,q)*signum*NormVec(1,p,q,FVSide,SideID) &
                                                         +gJacQy(:,:,p,q)*signum*NormVec(2,p,q,FVSide,SideID) &
#if PP_dim==3
                                                         +hJacQy(:,:,p,q)*signum*NormVec(3,p,q,FVSide,SideID) &
#endif
                                                        )*SurfElem(p,q,FVSide,SideID)
#if PP_dim==3
            Df_dQzInner(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQz(:,:,p,q)*signum*NormVec(1,p,q,FVSide,SideID) &
                                                         +gJacQz(:,:,p,q)*signum*NormVec(2,p,q,FVSide,SideID) &
                                                         +hJacQz(:,:,p,q)*signum*NormVec(3,p,q,FVSide,SideID) )*SurfElem(p,q,FVSide,SideID)
#endif
          END DO !p
        END DO !q
        ! Evaluate Jacobians of diffusive fluxes w.r.t. the outer gradients
        CALL  EvalFluxGradJacobian(nDOFFace,USideR,UPrimSideR, &
                                   fJacQx,fJacQy,fJacQz,       &
                                   gJacQx,gJacQy,gJacQz,       &
                                   hJacQx,hJacQy,hJacQz        &
#if EDDYVISCOSITY
                                  ,muSGSR                      &
#endif
                                   )
        IF((FVSum.EQ.1).OR.(FVSum.EQ.2))THEN
          ! For mixed interfaces, we ignore the contribution from the outer gradients
          Df_DQxOuter(:,:,:,:,iLocSide)=0.
          Df_DQyOuter(:,:,:,:,iLocSide)=0.
#if PP_dim==3
          Df_DQzOuter(:,:,:,:,iLocSide)=0.
#endif /*PP_dim*/
        ELSE
          DO q=0,PP_NZ
            DO p=0,PP_N
              ! Again, transform into normal system and pre-multiply with the SurfElem
              jk(:)=sideMap(:,p,q)
              ! Df_dQxOuter = d(f,g,h)_diff/dQ_x * NormVec
              ! Dependencies of the surface flux w.r.t. the outer surface gradients in each direction
              Df_dQxOuter(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQx(:,:,p,q)*signum*NormVec(1,p,q,FVSide,SideID) &
                                                           +gJacQx(:,:,p,q)*signum*NormVec(2,p,q,FVSide,SideID) &
#if PP_dim==3
                                                           +hJacQx(:,:,p,q)*signum*NormVec(3,p,q,FVSide,SideID) &
#endif
                                                          )*SurfElem(p,q,FVSide,SideID)
              Df_dQyOuter(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQy(:,:,p,q)*signum*NormVec(1,p,q,FVSide,SideID) &
                                                           +gJacQy(:,:,p,q)*signum*NormVec(2,p,q,FVSide,SideID) &
#if PP_dim==3
                                                           +hJacQy(:,:,p,q)*signum*NormVec(3,p,q,FVSide,SideID) &
#endif
                                                          )*SurfElem(p,q,FVSide,SideID)
#if PP_dim==3
              Df_dQzOuter(:,:,jk(1),jk(2),iLocSide)=  0.5*( fJacQz(:,:,p,q)*signum*NormVec(1,p,q,FVSide,SideID) &
                                                           +gJacQz(:,:,p,q)*signum*NormVec(2,p,q,FVSide,SideID) &
                                                           +hJacQz(:,:,p,q)*signum*NormVec(3,p,q,FVSide,SideID) )*SurfElem(p,q,FVSide,SideID)
#endif
            END DO !p
          END DO !q
        END IF !FVSum
      END IF ! HyperbolicPrecond==False
#endif /*PARABOLIC*/
      IF (nMortars.GT.1) THEN
        ! Store the Jacobians for each small side, later combine them into the Jacobian for the big side
        DfMortar_DUInner(:,:,:,:,1,iMortar)=Df_DUInner(:,:,:,:,1,iLocSide)
#if FV_ENABLED && FV_RECONSTRUCT
        DfMortar_DUInner(:,:,:,:,2,iMortar)=Df_DUInner(:,:,:,:,2,iLocSide)
#endif
#if PARABOLIC
        DfMortar_dQxInner(:,:,:,:,iMortar)=Df_DQxInner(:,:,:,:,iLocSide)
        DfMortar_dQyInner(:,:,:,:,iMortar)=Df_DQyInner(:,:,:,:,iLocSide)
        DfMortar_dQxOuter(:,:,:,:,iMortar)=Df_DQxOuter(:,:,:,:,iLocSide)
        DfMortar_dQyOuter(:,:,:,:,iMortar)=Df_DQyOuter(:,:,:,:,iLocSide)
#if PP_dim==3
        DfMortar_dQzInner(:,:,:,:,iMortar)=Df_DQzInner(:,:,:,:,iLocSide)
        DfMortar_dQzOuter(:,:,:,:,iMortar)=Df_DQzOuter(:,:,:,:,iLocSide)
#endif
#endif /*PARABOLIC*/
      END IF
    END DO ! nMortars

    IF (nMortars.GT.1) THEN
      ! Combine the Jacobians on the small sides into the Jacobian on the big side, also flip into volume system
      CALL Jacobian_MortarU(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                            DfMortar_DUInner(:,:,:,:,1,:),Df_DUInner(:,:,:,:,1,iLocSide))
#if FV_ENABLED && FV_RECONSTRUCT
      CALL Jacobian_MortarU(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                            DfMortar_DUInner(:,:,:,:,2,:),Df_DUInner(:,:,:,:,2,iLocSide))
#endif
#if PARABOLIC
      CALL Jacobian_MortarGrad(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                               DfMortar_dQxInner,Df_DQxInner(:,:,:,:,iLocSide))
      CALL Jacobian_MortarGrad(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                               DfMortar_dQyInner,Df_DQyInner(:,:,:,:,iLocSide))
      CALL Jacobian_MortarGrad(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                               DfMortar_dQxOuter,Df_DQxOuter(:,:,:,:,iLocSide))
      CALL Jacobian_MortarGrad(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                               DfMortar_dQyOuter,Df_DQyOuter(:,:,:,:,iLocSide))
#if PP_dim==3
      CALL Jacobian_MortarGrad(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                               DfMortar_dQzInner,Df_DQzInner(:,:,:,:,iLocSide))
      CALL Jacobian_MortarGrad(FVSide,MortarType(1,ElemToSide(E2S_SIDE_ID,ilocSide,iElem)),S2V2(:,:,:,flip,iLocSide), &
                               DfMortar_dQzOuter,Df_DQzOuter(:,:,:,:,iLocSide))
#endif
#endif /*PARABOLIC*/
    END IF ! mortar
  ELSE !Boundary
#if FV_ENABLED
    FVSide = FV_Elems_master(SideID)
#endif
#if PARABOLIC
    Df_DQxOuter(:,:,:,:,iLocSide)=0.
    Df_DQyOuter(:,:,:,:,iLocSide)=0.
#if PP_dim==3
    Df_DQzOuter(:,:,:,:,iLocSide)=0.
#endif /*PP_dim*/
#endif /*parabolic*/
    !df*_Boundary_jk/dU_jk*surfElem
    !df*_Boundary_jk/dQ_jk*surfElem
    CALL GetBoundaryFlux_FD(SideID,t,Df_DUinner(:,:,:,:,:,iLocSide),U_master(:,:,:,SideID),UPrim_master(:,:,:,SideID), &
#if PARABOLIC
                            gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),       &
                            Df_DQxInner(:,:,:,:,iLocSide),                                                             &
                            Df_DQyInner(:,:,:,:,iLocSide),                                                             &
#if PP_dim==3
                            Df_DQzInner(:,:,:,:,iLocSide),                                                             &
#endif
#endif /*PARABOLIC*/
                            SurfElem(:,:,FVSide,SideID),Face_xGP(:,:,:,FVSide,SideID),                                 &
                            NormVec(:,:,:,FVSide,SideID),TangVec1(:,:,:,FVSide,SideID),TangVec2(:,:,:,FVSide,SideID),  &
                            S2V2(:,:,:,Flip,iLocSide))
#if FV_ENABLED && FV_RECONSTRUCT
    ! Set UPrim_slave at boundaries
    CALL GetBoundaryState(SideID,t,PP_N,UPrim_slave(:,:,:,SideID),UPrim_master(:,:,:,SideID),                        &
                          NormVec(:,:,:,FVSide,SideID),TangVec1(:,:,:,FVSide,SideID),TangVec2(:,:,:,FVSide,SideID),  &
                          Face_xGP(:,:,:,FVSide,SideID))
#endif
  END IF!SideID
END DO !iLocSide

!Assembling of the preconditioner
!BJ=d(f*_adv+f*_diff)_jk/dU_mno
#if FV_ENABLED
IF(FVElem.EQ.0)THEN ! DG Element
#endif
  CALL Assemble_JacSurfInt_DG(iElem,Df_DUinner,                                &
#if PARABOLIC
                              Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter, &
#if PP_dim==3
                              Df_dQzInner,Df_dQzOuter,                         &
#endif
#endif
                              BJ)
#if FV_ENABLED
ELSE
  CALL Assemble_JacSurfInt_FV(iElem,Df_DUinner,                                &
#if PARABOLIC
                              Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter, &
#if PP_dim==3
                              Df_dQzInner,Df_dQzOuter,                         &
#endif
#endif
                              BJ)
END IF
#endif

END SUBROUTINE JacSurfInt

!===================================================================================================================================
!> This routine takes the derivatives of the surface flux w.r.t. the solution and the gradients, calculates the derivatives of the
!> lifting procedure and assembles the block-Jacobi preconditioner with the influence of the surface integral.
!> For more details about the derivation of the different contributions see comment above subroutine JacSurfInt.
!===================================================================================================================================
SUBROUTINE Assemble_JacSurfInt_DG(iElem,Df_DUInner,                                &
#if PARABOLIC
                                  Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter, &
#if PP_dim==3
                                  Df_dQzInner,Df_dQzOuter,                         &
#endif
#endif
                                  BJ)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Jac_Ex_Vars               ,ONLY: LL_minus, LL_plus
USE MOD_Implicit_Vars             ,ONLY: nDOFVarElem
#if PARABOLIC
USE MOD_Jac_Ex_br2                ,ONLY: dQOuter,dQInner
USE MOD_Jacobian                  ,ONLY: dPrimTempdCons
USE MOD_Precond_Vars              ,ONLY: NoFillIn
USE MOD_DG_Vars                   ,ONLY: L_Hatminus,L_Hatplus,UPrim
USE MOD_Precond_Vars              ,ONLY: HyperbolicPrecond
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                                      :: iElem      !< index of current element
#if PP_dim == 3
REAL,DIMENSION(1:PP_nVar,1:PP_nVar    ,0:PP_N,0:PP_NZ,2,1:6),INTENT(IN) :: Df_DUInner !< Jacobi of f_surf w.r.t. U
#if PARABOLIC
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_N   ,1:6),INTENT(IN) :: Df_dQxInner,Df_dQyInner,Df_dQzInner !< Jacobi of f_surf
                                                                                      !> w.r.t. inner gradients
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_N   ,1:6),INTENT(IN) :: Df_dQxOuter,Df_dQyOuter,Df_dQzOuter !< Jacobi of f_surf
                                                                                      !> w.r.t. gradients of neighbour
#endif
#else
REAL,DIMENSION(1:PP_nVar,1:PP_nVar    ,0:PP_N,0:PP_NZ,2,2:5),INTENT(IN) :: Df_DUInner
#if PARABOLIC
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_NZ  ,2:5),INTENT(IN) :: Df_dQxInner,Df_dQyInner
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_NZ  ,2:5),INTENT(IN) :: Df_dQxOuter,Df_dQyOuter
#endif
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(nDOFVarElem,nDOFVarElem),INTENT(INOUT)                   :: BJ         !< block-Jacobian of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: mm,nn,oo,r,s,vn1,vn2,ll
#if PARABOLIC
INTEGER         :: i,j,l
REAL            :: Df_dQ_minus(      1:PP_nVar,1:PP_nVar)
REAL            :: Df_dQ_plus(       1:PP_nVar,1:PP_nVar)
REAL            :: Df_dQ_plus_Tilde( 1:PP_nVar,1:PP_nVarPrim)
REAL            :: Df_dQ_minus_Tilde(1:PP_nVar,1:PP_nVarPrim)
REAL            :: dQxvol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
REAL            :: dQyvol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
REAL            :: PrimConsJac(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#if PP_dim==3
INTEGER         :: k
REAL            :: dQxInner_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL            :: dQyInner_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL            :: dQzInner_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL            :: dQxOuter_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL            :: dQyOuter_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL            :: dQzOuter_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,6,0:PP_N)
REAL            :: dQzvol_dU(0:PP_N,0:PP_N,0:PP_N,0:PP_N,3)
#else
REAL            :: dQxInner_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
REAL            :: dQyInner_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
REAL            :: dQxOuter_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
REAL            :: dQyOuter_dUvol(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif /*PP_dim*/
#endif /*PARABOLIC*/
!===================================================================================================================================
! Helper variables to quickly build the one-dimensional mapping: ind = iVar+PP_nVar*i+vn1*j+vn2*k for DOF(iVar,i,j,k)
vn1 = PP_nVar * (PP_N + 1)
vn2 = vn1 * (PP_N +1)
!Assembling of the preconditioner
!BJ=d(f*_adv+f*_diff)_jk/dU_mno


! The Jacobians of the surface fluxes have already been multiplied by the surface element. Df_DUinner considers both the
! dependency of the diffusive and the hyperbolic flux w.r.t. the inner solution.
! The matrix LL_plus/minus now takes into account how the surface solution is depending on the considered volume DOF and the
! derivative of the surface integral itself (consists of prolongation L and integration L_hat).

! Loop over all volume DOFs
DO oo = 0,PP_NZ
  DO nn = 0,PP_N
    DO mm = 0,PP_N
      ! One-dimensional index for the DOF (mm,nn,oo)
      s = vn2 * oo + vn1 * nn + PP_nVar * mm
      ! Loop over XI,ETA and ZETA lines
      DO ll = 0,PP_N
        ! Dependency of flux on XI line with index ll on volume DOF mm,nn,oo
        r = PP_nVar*ll + vn1*nn + vn2*oo
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                          + LL_Minus(ll,mm)*Df_DUinner(:,:,nn,oo,1,XI_MINUS  )  &
                                          + LL_Plus (ll,mm)*Df_DUinner(:,:,nn,oo,1,XI_PLUS   )
        ! Dependency of flux on ETA line with index ll on volume DOF mm,nn,oo
        r = PP_nVar*mm + vn1*ll + vn2*oo
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                          + LL_Minus(ll,nn)*Df_DUinner(:,:,mm,oo,1,ETA_MINUS )  &
                                          + LL_Plus (ll,nn)*Df_DUinner(:,:,mm,oo,1,ETA_PLUS  )
#if PP_dim==3
        ! Dependency of flux on ZETA line with index ll on volume DOF mm,nn,oo
        r = PP_nVar*mm + vn1*nn + vn2*ll
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                          + LL_Minus(ll,oo)*Df_DUinner(:,:,mm,nn,1,ZETA_MINUS)  &
                                          + LL_Plus (ll,oo)*Df_DUinner(:,:,mm,nn,1,ZETA_PLUS )
#endif
      END DO ! ll
    END DO ! nn
  END DO ! mm 
END DO ! oo

#if PARABOLIC
IF (.NOT.HyperbolicPrecond) THEN
  ! Compute the following derivatives (depending on the lifting scheme): 
  !  * dependency of inner surface gradients w.r.t. primitive volume DOFs dQInner_dUvol
  !  * dependency of the lifting volume integral w.r.t. the volume DOFs dQVol_dUvol
  !  * dependency of the outer surface gradients w.r.t. the volume DOFs dQOuter_dUvol
  CALL dQInner(1,iElem,dQxInner_dUvol,dQxVol_dU)
  CALL dQInner(2,iElem,dQyInner_dUvol,dQyVol_dU)
  CALL dQOuter(1,iElem,dQxOuter_dUvol)
  CALL dQOuter(2,iElem,dQyOuter_dUvol)
#if PP_dim==3
  CALL dQInner(3,iElem,dQzInner_dUvol,dQzVol_dU)
  CALL dQOuter(3,iElem,dQzOuter_dUvol)
#endif

  ! What is left now:
  !   DF^v/DF^v_surf * DF^v(U,Q(U))_surf/DQ_surf * DQ_surf/DU
  ! = DF^v/DF^v_surf * 0.5*(D(F^v_L)/DQ_surf_L * DQ_surf_L/DU + D(F^v_R)/DQ_surf_R * DQ_surf_R/DU)
  ! We split DQ_surf_L/DU:
  ! DQ_surf_L/DU = DQ_surf_L|_VolInt/DU_prim * DU_prim/DU + DQ_surf_L|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU
  ! We first compute the term:
  ! 0.5*(D(F^v_L)/DQ_surf_L * DQ_surf_L|_VolInt/DU_prim * DU_prim/DU
  !  => Df_dQxInner              => dQxVol_dU             => PrimConsJac
  ! How is the viscous flux in the volume depending on the volume DOFs via the volume integral of the lifting?
  DO oo = 0,PP_NZ
    DO nn = 0,PP_N
      DO mm = 0,PP_N
        ! XI-direction--------------------------------------------------------------------------------------------------------------
        CALL dPrimTempdCons(UPrim(:,mm,nn,oo,iElem),PrimConsJac(:,:,mm,nn,oo))
        s = vn2 * oo + vn1 * nn + PP_nVar * mm
        ! Dependeny of all viscous XI-fluxes along the XI line w.r.t. my DOFs
        DO l=0,PP_N
          Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,nn,oo,XI_MINUS  ) * dQxVol_dU(l,nn,oo,mm,1) &
                                  +Df_dQyInner(:,:,nn,oo,XI_MINUS  ) * dQyVol_dU(l,nn,oo,mm,1) & 
#if PP_dim==3
                                  +Df_dQzInner(:,:,nn,oo,XI_MINUS  ) * dQzVol_dU(l,nn,oo,mm,1) &
#endif
                                  )
          df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          Df_dQ_plus_Tilde(:,:) = (Df_dQxInner(:,:,nn,oo,XI_PLUS   ) * dQxvol_dU(l,nn,oo,mm,1) &
                                  +Df_dQyInner(:,:,nn,oo,XI_PLUS   ) * dQyvol_dU(l,nn,oo,mm,1) & 
#if PP_dim==3
                                  +Df_dQzInner(:,:,nn,oo,XI_PLUS   ) * dQzvol_dU(l,nn,oo,mm,1) &
#endif
                                  )        
          df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          ! Surface integral and prolongation
          DO i = 0,PP_N
            r = PP_nVar * i + vn1 * nn+ vn2 * oo 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                            + LL_Minus(i,l)* df_dQ_minus(:,:) &
                                            + LL_Plus(i,l) * df_dQ_plus(:,:)
          END DO !i
        END DO ! l
        IF(.NOT.NoFillIn) THEN
          ! Dependeny of all viscous XI-fluxes along the ETA line w.r.t. my DOFs
          DO j = 0,PP_N
            Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,j,oo,XI_MINUS  ) * dQxVol_dU(mm,j,oo,nn,2) &
                                    +Df_dQyInner(:,:,j,oo,XI_MINUS  ) * dQyVol_dU(mm,j,oo,nn,2) & 
#if PP_dim==3
                                    +Df_dQzInner(:,:,j,oo,XI_MINUS  ) * dQzVol_dU(mm,j,oo,nn,2) &
#endif
                                    )
            df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            df_dQ_plus_Tilde(:,:) = (Df_dQxInner(:,:,j,oo,XI_PLUS   ) * dQxvol_dU(mm,j,oo,nn,2) &
                                    +Df_dQyInner(:,:,j,oo,XI_PLUS   ) * dQyvol_dU(mm,j,oo,nn,2) &
#if PP_dim==3
                                    +Df_dQzInner(:,:,j,oo,XI_PLUS   ) * dQzvol_dU(mm,j,oo,nn,2) &
#endif
                                    )
            df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            ! Surface integral and prolongation
            DO i = 0,PP_N
              r = PP_nVar * i + vn1 * j+ vn2 * oo 
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                              + LL_Minus(i,mm)*df_dQ_minus(:,:) & 
                                              + LL_Plus(i,mm) * df_dQ_plus(:,:) 
            END DO ! i
          END DO ! j
#if PP_dim==3
          ! Dependeny of all viscous XI-fluxes along the ZETA line w.r.t. my DOFs
          DO k = 0,PP_NZ
            Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,nn,k,XI_MINUS  ) * dQxVol_dU(mm,nn,k,oo,3) &
                                    +Df_dQyInner(:,:,nn,k,XI_MINUS  ) * dQyVol_dU(mm,nn,k,oo,3) &
                                    +Df_dQzInner(:,:,nn,k,XI_MINUS  ) * dQzVol_dU(mm,nn,k,oo,3) &
                                    )
            df_dQ_minus(:,:)= MATMUL(Df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            df_dQ_plus_Tilde(:,:) = (Df_dQxInner(:,:,nn,k,XI_PLUS   ) * dQxvol_dU(mm,nn,k,oo,3) &
                                    +Df_dQyInner(:,:,nn,k,XI_PLUS   ) * dQyvol_dU(mm,nn,k,oo,3) &
                                    +Df_dQzInner(:,:,nn,k,XI_PLUS   ) * dQzvol_dU(mm,nn,k,oo,3) &
                                    )
            df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            ! Surface integral and prolongation
            DO i = 0,PP_N
              r = PP_nVar * i + vn1 * nn+ vn2 * k 
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                              + LL_Minus(i,mm)*df_dQ_minus(:,:) &
                                              + LL_Plus(i,mm) * df_dQ_plus(:,:)
            END DO ! i
          END DO ! k
#endif
          ! ETA-direction-------------------------------------------------------------------------------------------------------------
          ! Dependeny of all viscous ETA-fluxes along the XI line w.r.t. my DOFs
          DO i = 0,PP_N
            Df_dQ_minus_Tilde(:,:)= ( Df_dQxInner(:,:,i,oo,ETA_MINUS ) * dQxvol_dU(i,nn,oo,mm,1) &
                                     +Df_dQyInner(:,:,i,oo,ETA_MINUS ) * dQyvol_dU(i,nn,oo,mm,1) &
#if PP_dim==3
                                     +Df_dQzInner(:,:,i,oo,ETA_MINUS ) * dQzvol_dU(i,nn,oo,mm,1) &
#endif
                                    )
            df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            df_dQ_plus_Tilde(:,:) = ( Df_dQxInner(:,:,i,oo,ETA_PLUS  ) * dQxvol_dU(i,nn,oo,mm,1) &
                                     +Df_dQyInner(:,:,i,oo,ETA_PLUS  ) * dQyvol_dU(i,nn,oo,mm,1) &
#if PP_dim==3
                                     +Df_dQzInner(:,:,i,oo,ETA_PLUS  ) * dQzvol_dU(i,nn,oo,mm,1) &
#endif
                                    )
            df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            DO j = 0,PP_N
              r = PP_nVar * i + vn1 * j+ vn2 * oo 
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                              + LL_Minus(j,nn)*df_dQ_minus(:,:) & 
                                              + LL_Plus(j,nn) * df_dQ_plus(:,:)
            END DO ! j
          END DO ! i
        END IF !NoFillIn
        ! Dependeny of all viscous ETA-fluxes along the ETA line w.r.t. my DOFs
        DO l = 0,PP_N
          df_dQ_minus_Tilde(:,:)= ( Df_dQxInner(:,:,mm,oo,ETA_MINUS ) * dQxvol_dU(mm,l,oo,nn,2) &
                                   +Df_dQyInner(:,:,mm,oo,ETA_MINUS ) * dQyvol_dU(mm,l,oo,nn,2) &
#if PP_dim==3
                                   +Df_dQzInner(:,:,mm,oo,ETA_MINUS ) * dQzvol_dU(mm,l,oo,nn,2) &
#endif
                                  )
          df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus_Tilde(:,:)= ( Df_dQxInner(:,:,mm,oo,ETA_PLUS  ) * dQxvol_dU(mm,l,oo,nn,2) &
                                  +Df_dQyInner(:,:,mm,oo,ETA_PLUS  ) * dQyvol_dU(mm,l,oo,nn,2) & 
#if PP_dim==3
                                  +Df_dQzInner(:,:,mm,oo,ETA_PLUS  ) * dQzvol_dU(mm,l,oo,nn,2) &
#endif
                                 )
          df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO j=0,PP_N
            r = PP_nVar * mm + vn1 * j+ vn2 * oo
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                            + LL_Minus(j,l) *df_dQ_minus(:,:) &
                                            + LL_Plus(j,l)  * df_dQ_plus(:,:) 
          END DO !j
        END DO ! l
#if PP_dim==3
        IF(.NOT.NoFillIn) THEN
          ! Dependeny of all viscous ETA-fluxes along the ZETA line w.r.t. my DOFs
          DO k = 0,PP_NZ
            df_dQ_minus_Tilde(:,:)=  (Df_dQxInner(:,:,mm,k,ETA_MINUS ) * dQxvol_dU(mm,nn,k,oo,3) &
                                     +Df_dQyInner(:,:,mm,k,ETA_MINUS ) * dQyvol_dU(mm,nn,k,oo,3) & 
                                     +Df_dQzInner(:,:,mm,k,ETA_MINUS ) * dQzvol_dU(mm,nn,k,oo,3) &
                                     )
            df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            df_dQ_plus_Tilde(:,:)= ( Df_dQxInner(:,:,mm,k,ETA_PLUS  ) * dQxvol_dU(mm,nn,k,oo,3) &
                                    +Df_dQyInner(:,:,mm,k,ETA_PLUS  ) * dQyvol_dU(mm,nn,k,oo,3) &
                                    +Df_dQzInner(:,:,mm,k,ETA_PLUS  ) * dQzvol_dU(mm,nn,k,oo,3) &
                                   )
            df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            DO j = 0,PP_N
              r = PP_nVar * mm + vn1 * j+ vn2 * k 
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                              + LL_Minus(j,nn)*df_dQ_minus(:,:) &
                                              + LL_Plus(j,nn) * df_dQ_plus(:,:) 
            END DO ! j
          END DO ! k
          ! ZETA-direction------------------------------------------------------------------------------------------------------------
          ! Dependeny of all viscous ZETA-fluxes along the XI line w.r.t. my DOFs
          DO i = 0,PP_N
            df_dQ_minus_Tilde(:,:)=  (Df_dQxInner(:,:,i,nn,ZETA_MINUS) * dQxvol_dU(i,nn,oo,mm,1) &
                                     +Df_dQyInner(:,:,i,nn,ZETA_MINUS) * dQyvol_dU(i,nn,oo,mm,1) &
                                     +Df_dQzInner(:,:,i,nn,ZETA_MINUS) * dQzvol_dU(i,nn,oo,mm,1) &
                                     )
            df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            
            df_dQ_plus_Tilde(:,:) =  (Df_dQxInner(:,:,i,nn,ZETA_PLUS ) * dQxvol_dU(i,nn,oo,mm,1) &
                                     +Df_dQyInner(:,:,i,nn,ZETA_PLUS ) * dQyvol_dU(i,nn,oo,mm,1) &
                                     +Df_dQzInner(:,:,i,nn,ZETA_PLUS ) * dQzvol_dU(i,nn,oo,mm,1) &
                                     )
            df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            DO k = 0,PP_NZ
              r = PP_nVar * i + vn1 * nn+ vn2 * k 
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                              + LL_Minus(k,oo)*df_dQ_minus(:,:) & 
                                              + LL_Plus(k,oo) * df_dQ_plus(:,:) 
            END DO ! k
          END DO ! i
          ! Dependeny of all viscous ZETA-fluxes along the ETA line w.r.t. my DOFs
          DO j = 0,PP_N
            df_dQ_minus_Tilde(:,:)=  (Df_dQxInner(:,:,mm,j,ZETA_MINUS) * dQxvol_dU(mm,j,oo,nn,2) &
                                     +Df_dQyInner(:,:,mm,j,ZETA_MINUS) * dQyvol_dU(mm,j,oo,nn,2) &
                                     +Df_dQzInner(:,:,mm,j,ZETA_MINUS) * dQzvol_dU(mm,j,oo,nn,2) &
                                     )
            df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            df_dQ_plus_Tilde(:,:) =  (Df_dQxInner(:,:,mm,j,ZETA_PLUS ) * dQxvol_dU(mm,j,oo,nn,2) &
                                     +Df_dQyInner(:,:,mm,j,ZETA_PLUS ) * dQyvol_dU(mm,j,oo,nn,2) &
                                     +Df_dQzInner(:,:,mm,j,ZETA_PLUS ) * dQzvol_dU(mm,j,oo,nn,2) &
                                     )
            df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
            DO k = 0,PP_NZ
              r = PP_nVar * mm + vn1 * j+ vn2 * k 
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                              + LL_Minus(k,oo)*df_dQ_minus(:,:) &
                                              + LL_Plus(k,oo) * df_dQ_plus(:,:)
            END DO ! k
          END DO ! j
        END IF !NoFillIn
        ! Dependeny of all viscous ZETA-fluxes along the ZETA line w.r.t. my DOFs
        DO l = 0,PP_N
          df_dQ_minus_Tilde(:,:)=  (Df_dQxInner(:,:,mm,nn,ZETA_MINUS) * dQxvol_dU(mm,nn,l,oo,3) &
                                   +Df_dQyInner(:,:,mm,nn,ZETA_MINUS) * dQyvol_dU(mm,nn,l,oo,3) &
                                   +Df_dQzInner(:,:,mm,nn,ZETA_MINUS) * dQzvol_dU(mm,nn,l,oo,3) &
                                   )
          df_dQ_minus(:,:)= MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          df_dQ_plus_Tilde(:,:)=  (Df_dQxInner(:,:,mm,nn,ZETA_PLUS ) * dQxvol_dU(mm,nn,l,oo,3) &
                                  +Df_dQyInner(:,:,mm,nn,ZETA_PLUS ) * dQyvol_dU(mm,nn,l,oo,3) &
                                  +Df_dQzInner(:,:,mm,nn,ZETA_PLUS ) * dQzvol_dU(mm,nn,l,oo,3) &
                                  )
          df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,mm,nn,oo))
          DO k=0,PP_NZ
            r = PP_nVar * mm + vn1 * nn+ vn2 * k 
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)  &
                                            + LL_Minus(k,l) * df_dQ_minus(:,:) & 
                                            + LL_Plus(k,l)  * df_dQ_plus(:,:) 
          END DO !k
        END DO ! l
#endif
      END DO ! nn
    END DO ! mm 
  END DO ! oo

! What is left now:
! DF^v/DF^v_surf * 0.5*(D(F^v_L)/DQ_surf_L * DQ_surf_L|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU 
!                     + D(F^v_R)/DQ_surf_R * DQ_surf_R/DU)
! Since: 
! DQ_surf_R/DU = DQ_surf_R|_VolInt/DU + DQ_surf_R|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU
!                         = 0
!   DF^v/DF^v_surf * 0.5*(D(F^v_L)/DQ_surf_L * DQ_surf_L|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU 
!                       + D(F^v_R)/DQ_surf_R * DQ_surf_R/DU)
! = DF^v/DF^v_surf * 0.5*(D(F^v_L)/DQ_surf_L * DQ_surf_L|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU 
!                           => Df_dQxInner          => dQxInner_dUVol
!                       + D(F^v_R)/DQ_surf_R * DQ_surf_R|_SurfInt/DU_Lprim * DU_Lprim/DU_Lcons * DU_Lcons/DU)
!                           => Df_dQxOuter          => dQxOuter_dUVol
! The DQxInner/Outer_DUvol already consider the prolongation DU_Lcons/DU and the DU_Lprim/DU_Lcons on the surface!
! How is the viscous flux in the volume depending on the volume DOFs via the surface integral of the lifting for both 
! my element and my neighbouring element?
  DO oo = 0,PP_NZ
    DO nn = 0,PP_N
      DO mm = 0,PP_N
        s = vn2 * oo + vn1 * nn + PP_nVar * mm
        df_dQ_minus(:,:)= (MATMUL(Df_dQxinner(:,:,nn,oo,XI_MINUS  ) ,dQxinner_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
                          +MATMUL(Df_dQyinner(:,:,nn,oo,XI_MINUS  ) ,dQyinner_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
#if PP_dim==3
                          +MATMUL(Df_dQzinner(:,:,nn,oo,XI_MINUS  ) ,dQzinner_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
#endif
                          +MATMUL(Df_dQxOuter(:,:,nn,oo,XI_MINUS  ) ,dQxOuter_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
                          +MATMUL(Df_dQyOuter(:,:,nn,oo,XI_MINUS  ) ,dQyOuter_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) & 
#if PP_dim==3
                          +MATMUL(Df_dQzOuter(:,:,nn,oo,XI_MINUS  ) ,dQzOuter_dUvol(:,:,nn,oo,XI_MINUS  ,mm) ) &
#endif
                          )
        df_dQ_plus(:,:) = (MATMUL(Df_dQxinner(:,:,nn,oo,XI_PLUS   ) , dQxinner_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
                          +MATMUL(Df_dQyinner(:,:,nn,oo,XI_PLUS   ) , dQyinner_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
#if PP_dim==3
                          +MATMUL(Df_dQzinner(:,:,nn,oo,XI_PLUS   ) , dQzinner_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
#endif
                          +MATMUL(Df_dQxOuter(:,:,nn,oo,XI_PLUS   ) , dQxOuter_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
                          +MATMUL(Df_dQyOuter(:,:,nn,oo,XI_PLUS   ) , dQyOuter_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
#if PP_dim==3
                          +MATMUL(Df_dQzOuter(:,:,nn,oo,XI_PLUS   ) , dQzOuter_dUvol(:,:,nn,oo,XI_PLUS   ,mm) ) &
#endif
                          )
        DO l = 0,PP_N
          r = PP_nVar * l + vn1 * nn+ vn2 * oo 
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)  &
                                          + L_HatMinus(l) * df_dQ_minus(:,:) &
                                          + L_HatPlus(l)  * df_dQ_plus(:,:) 
        END DO ! l
        df_dQ_minus(:,:)= (MATMUL(Df_dQxinner(:,:,mm,oo,ETA_MINUS ) , dQxinner_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) &
                          +MATMUL(Df_dQyinner(:,:,mm,oo,ETA_MINUS ) , dQyinner_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) & 
#if PP_dim==3
                          +MATMUL(Df_dQzinner(:,:,mm,oo,ETA_MINUS ) , dQzinner_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) &
#endif
                          +MATMUL(Df_dQxOuter(:,:,mm,oo,ETA_MINUS ) , dQxOuter_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) &
                          +MATMUL(Df_dQyOuter(:,:,mm,oo,ETA_MINUS ) , dQyOuter_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) &
#if PP_dim==3
                          +MATMUL(Df_dQzOuter(:,:,mm,oo,ETA_MINUS ) , dQzOuter_dUvol(:,:,mm,oo,ETA_MINUS ,nn) ) &
#endif
                          )
        df_dQ_plus(:,:) = (MATMUL(Df_dQxinner(:,:,mm,oo,ETA_PLUS  ) , dQxinner_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) &
                          +MATMUL(Df_dQyinner(:,:,mm,oo,ETA_PLUS  ) , dQyinner_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) & 
#if PP_dim==3
                          +MATMUL(Df_dQzinner(:,:,mm,oo,ETA_PLUS  ) , dQzinner_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) &
#endif
                          +MATMUL(Df_dQxOuter(:,:,mm,oo,ETA_PLUS  ) , dQxOuter_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) &
                          +MATMUL(Df_dQyOuter(:,:,mm,oo,ETA_PLUS  ) , dQyOuter_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) &
#if PP_dim==3                                                                           
                          +MATMUL(Df_dQzOuter(:,:,mm,oo,ETA_PLUS  ) , dQzOuter_dUvol(:,:,mm,oo,ETA_PLUS  ,nn) ) &
#endif
                          )
        DO l = 0,PP_N
          r = PP_nVar * mm + vn1 * l+ vn2 * oo 
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)  &
                                          + L_HatMinus(l) * df_dQ_minus(:,:) &
                                          + L_HatPlus(l)  * df_dQ_plus(:,:)
        END DO ! l
#if PP_dim==3
        df_dQ_minus(:,:)= (MATMUL(Df_dQxinner(:,:,mm,nn,ZETA_MINUS) , dQxinner_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                          +MATMUL(Df_dQyinner(:,:,mm,nn,ZETA_MINUS) , dQyinner_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                          +MATMUL(Df_dQzinner(:,:,mm,nn,ZETA_MINUS) , dQzinner_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                          +MATMUL(Df_dQxOuter(:,:,mm,nn,ZETA_MINUS) , dQxOuter_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                          +MATMUL(Df_dQyOuter(:,:,mm,nn,ZETA_MINUS) , dQyOuter_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                          +MATMUL(Df_dQzOuter(:,:,mm,nn,ZETA_MINUS) , dQzOuter_dUvol(:,:,mm,nn,ZETA_MINUS,oo) ) &
                          )
        df_dQ_plus(:,:)= (MATMUL(Df_dQxinner(:,:,mm,nn,ZETA_PLUS ) , dQxinner_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                         +MATMUL(Df_dQyinner(:,:,mm,nn,ZETA_PLUS ) , dQyinner_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                         +MATMUL(Df_dQzinner(:,:,mm,nn,ZETA_PLUS ) , dQzinner_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                         +MATMUL(Df_dQxOuter(:,:,mm,nn,ZETA_PLUS ) , dQxOuter_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                         +MATMUL(Df_dQyOuter(:,:,mm,nn,ZETA_PLUS ) , dQyOuter_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                         +MATMUL(Df_dQzOuter(:,:,mm,nn,ZETA_PLUS ) , dQzOuter_dUvol(:,:,mm,nn,ZETA_PLUS ,oo) ) &
                         )
        DO l = 0,PP_N
          r = PP_nVar * mm + vn1 * nn+ vn2 * l 
          BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)  &
                                          + L_HatMinus(l) * df_dQ_minus(:,:) &
                                          + L_HatPlus(l)  * df_dQ_plus(:,:) 
        END DO ! l
#endif
      END DO ! nn
    END DO ! mm 
  END DO ! oo
END IF !HyperbolicPrecond
#endif /*PARABOLIC*/
END SUBROUTINE Assemble_JacSurfInt_DG

#if FV_ENABLED
!===================================================================================================================================
!> This routine takes the derivatives of the surface flux w.r.t. the solution and the gradients, calculates the derivatives of the
!> FV reconstruction procedure and assembles the block-Jacobi preconditioner with the influence of the surface integral.
!> For more details about the derivation of the different contributions see comment above subroutine JacSurfInt.
!===================================================================================================================================
SUBROUTINE Assemble_JacSurfInt_FV(iElem,Df_DUInner,                                &
#if PARABOLIC
                                  Df_dQxInner,Df_dQyInner,Df_dQxOuter,Df_dQyOuter, &
#if PP_dim==3
                                  Df_dQzInner,Df_dQzOuter,                         &
#endif
#endif
                                  BJ)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars         ,ONLY: nDOFVarElem
USE MOD_FV_Vars               ,ONLY: FV_w_inv
#if PARABOLIC
USE MOD_Jac_Ex_Reconstruction ,ONLY: JacFVGradients_Vol,JacFVGradients_nb
USE MOD_Jacobian              ,ONLY: dPrimTempdCons
USE MOD_Precond_Vars          ,ONLY: NoFillIn
USE MOD_DG_Vars               ,ONLY: UPrim
USE MOD_Precond_Vars          ,ONLY: HyperbolicPrecond
#endif
#if FV_RECONSTRUCT
USE MOD_DG_Vars               ,ONLY: UPrim_master,UPrim_slave
USE MOD_Mesh_Vars             ,ONLY: ElemToSide,S2V2,firstInnerSide,nBCSides,firstMortarMPISide,lastMortarMPISide
USE MOD_Jac_Ex_Vars           ,ONLY: UPrim_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended
USE MOD_Jac_Ex_Reconstruction ,ONLY: FV_Reconstruction_Derivative_Surf
USE MOD_FV_Vars               ,ONLY: FV_dx_master,FV_dx_slave
#if PP_dim == 3
USE MOD_Jac_Ex_Vars           ,ONLY: FV_sdx_ZETA_extended
#endif
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                                      :: iElem      !< index of current element
#if PP_dim == 3
REAL,DIMENSION(1:PP_nVar,1:PP_nVar    ,0:PP_N,0:PP_NZ,2,1:6),INTENT(IN) :: Df_DUInner !< Jacobi of f_surf w.r.t. U
#if PARABOLIC
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_N   ,1:6),INTENT(IN) :: Df_dQxInner,Df_dQyInner,Df_dQzInner !< Jacobi of f_surf
                                                                                      !> w.r.t. inner gradients
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_N   ,1:6),INTENT(IN) :: Df_dQxOuter,Df_dQyOuter,Df_dQzOuter !< Jacobi of f_surf
                                                                                      !> w.r.t. gradients of neighbour
#endif
#else
REAL,DIMENSION(1:PP_nVar,1:PP_nVar    ,0:PP_N,0:PP_NZ,2,2:5),INTENT(IN) :: Df_DUInner
#if PARABOLIC
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_NZ  ,2:5),INTENT(IN) :: Df_dQxInner,Df_dQyInner
REAL,DIMENSION(1:PP_nVar,1:PP_nVarPrim,0:PP_N,0:PP_NZ  ,2:5),INTENT(IN) :: Df_dQxOuter,Df_dQyOuter
#endif
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(nDOFVarElem,nDOFVarElem),INTENT(INOUT)                   :: BJ         !< block-Jacobian of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: mm,nn,oo,r,s,vn1,vn2
INTEGER                     :: SideID,Flip
#if PARABOLIC
INTEGER                     :: p,q
REAL                        :: Df_dQ_minus(      1:PP_nVar,1:PP_nVar)
REAL                        :: Df_dQ_plus(       1:PP_nVar,1:PP_nVar)
REAL                        :: Df_dQ_plus_Tilde( 1:PP_nVar,1:PP_nVarPrim)
REAL                        :: Df_dQ_minus_Tilde(1:PP_nVar,1:PP_nVarPrim)
REAL                        :: dQxvol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
REAL                        :: dQyvol_dU(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
REAL                        :: PrimConsJac(PP_nVarPrim,PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#if PP_dim==3
REAL                        :: dQxOuter_dUvol(0:PP_N,0:PP_NZ,6,0:PP_N)
REAL                        :: dQyOuter_dUvol(0:PP_N,0:PP_NZ,6,0:PP_N)
REAL                        :: dQzOuter_dUvol(0:PP_N,0:PP_NZ,6,0:PP_N)
REAL                        :: dQzvol_dU(0:PP_N,0:PP_N,0:PP_N,0:PP_N,3)
#else
REAL                        :: dQxOuter_dUvol(0:PP_N,0:PP_NZ,2:5,0:PP_N)
REAL                        :: dQyOuter_dUvol(0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif /*PP_dim*/
#endif /*PARABOLIC*/
#if FV_RECONSTRUCT
REAL,DIMENSION(PP_nVarPrim) :: UPrim_plus,UPrim_minus,UPrim_plus_nb,UPrim_minus_nb
REAL                        :: dUdUvol_plus( PP_nVar,PP_nVar,PP_N:PP_N+1,1:2)
REAL                        :: dUdUvol_minus(PP_nVar,PP_nVar,-1:0,2:3)
INTEGER                     :: pq(2)
REAL                        :: FV_dx_L,FV_dx_R,FV_dx_L_nb,FV_dx_R_nb
LOGICAL                     :: Mortar_minus,Mortar_plus
#if PARABOLIC
INTEGER                     :: ss(2)
#endif
#endif
!===================================================================================================================================
! Helper variables to quickly build the one-dimensional mapping: ind = iVar+PP_nVar*i+vn1*j+vn2*k for DOF(iVar,i,j,k)
vn1 = PP_nVar * (PP_N + 1)
vn2 = vn1 * (PP_N +1)
!Assembling of the preconditioner
!BJ=d(f*_adv+f*_diff)_jk/dU_mno

! The Jacobians of the surface fluxes have already been multiplied by the surface element. Df_DUinner considers both the
! dependency of the diffusive and the hyperbolic flux w.r.t. the inner solution.
! The first part calculates DF/DF_surf * DF_surf/DU_surf * DU_surf/DU for the XI/ETA/ZETA direction
!                           =>FV_w_inv   => Df_DUinner     =>dUdUvol_minus/plus (derivative of reconstruction)
! XI-direction-------------------------------------------------------------------------------------------------------------------
DO oo = 0,PP_NZ
  DO nn = 0,PP_N
#if FV_RECONSTRUCT
    SideID=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
    Flip=ElemToSide(  E2S_FLIP   ,XI_MINUS,iElem)
    pq=S2V2(:,nn,oo,flip         ,XI_MINUS)
    IF (((SideID.GT.nBCSides).AND.(SideID.LT.firstInnerSide)).OR. &
        ((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN ! big mortar side
      Mortar_minus = .TRUE.
    ELSE
      Mortar_minus = .FALSE.
      IF(Flip.EQ.0)THEN
        UPrim_minus    = UPrim_master(:,pq(1),pq(2),SideID)
        UPrim_minus_nb = UPrim_slave( :,pq(1),pq(2),SideID)
        FV_dx_L        = FV_dx_master(1,pq(1),pq(2),SideID)
        IF(SideID.GE.firstInnerSide)THEN
          FV_dx_R_nb     = FV_dx_slave( 1,pq(1),pq(2),SideID)
        ELSE ! BC
          FV_dx_R_nb     = FV_dx_master(1,pq(1),pq(2),SideID)
        END IF
      ELSE
        UPrim_minus    = UPrim_slave( :,pq(1),pq(2),SideID)
        UPrim_minus_nb = UPrim_master(:,pq(1),pq(2),SideID)
        FV_dx_L        = FV_dx_slave( 1,pq(1),pq(2),SideID)
        FV_dx_R_nb     = FV_dx_master(1,pq(1),pq(2),SideID)
      END IF
    END IF
    SideID=ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem)
    Flip=ElemToSide(  E2S_FLIP   ,XI_PLUS,iElem)
    pq=S2V2(:,nn,oo,flip         ,XI_PLUS)
    IF (((SideID.GT.nBCSides).AND.(SideID.LT.firstInnerSide)).OR. &
        ((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN ! big mortar side
      Mortar_plus = .TRUE.
    ELSE
      Mortar_plus = .FALSE.
      IF(Flip.EQ.0)THEN
        UPrim_plus    = UPrim_master(:,pq(1),pq(2),SideID)
        UPrim_plus_nb = UPrim_slave( :,pq(1),pq(2),SideID)
        FV_dx_R       = FV_dx_master(1,pq(1),pq(2),SideID)
        IF(SideID.GE.firstInnerSide)THEN
          FV_dx_L_nb    = FV_dx_slave( 1,pq(1),pq(2),SideID)
        ELSE ! BC
          FV_dx_L_nb    = FV_dx_master(1,pq(1),pq(2),SideID)
        END IF
      ELSE
        UPrim_plus    = UPrim_slave( :,pq(1),pq(2),SideID)
        UPrim_plus_nb = UPrim_master(:,pq(1),pq(2),SideID)
        FV_dx_R       = FV_dx_slave( 1,pq(1),pq(2),SideID)
        FV_dx_L_nb    = FV_dx_master(1,pq(1),pq(2),SideID)
      END IF
    END IF
    CALL FV_Reconstruction_Derivative_Surf(FV_sdx_XI_extended(nn,oo,:,iElem),FV_dx_L,FV_dx_R,              & 
                                           FV_dx_L_nb,FV_dx_R_nb,UPrim_plus,UPrim_minus,                   &
                                           UPrim_plus_nb,UPrim_minus_nb,UPrim_extended(:,:,nn,oo,iElem),   &
                                           Mortar_minus,Mortar_plus,dUdUvol_plus(:,:,:,:),dUdUvol_minus(:,:,:,:))
    !-------------------Derivatives at MINUS side with respect to volume dofs----------------------------------------------------
    s = vn2*oo + vn1*nn
    ! direct dependency of prolongated value (of current element, minus) to volume dofs next to interface
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,nn,oo,1,XI_MINUS),dUdUvol_minus(:,:,0,2)) 
    ! dependency of prolongated value (of neighbouring element, plus) to volume dofs next to interface
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,nn,oo,2,XI_MINUS),dUdUvol_minus(:,:,-1,3)) 
    ! dependency of prolongated value (of current element, minus) to volume dofs neighbouring the dofs next to interface
    r = vn2*oo + vn1*nn + PP_nVar
    BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,nn,oo,1,XI_MINUS),dUdUvol_minus(:,:,0,3)) 
    !-------------------Derivatives at PLUS side with respect to volume dofs-----------------------------------------------------
    s = vn2*oo + vn1*nn + PP_nVar*PP_N
    ! direct dependency of prolongated value (of current element, plus) to volume dofs next to interface
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,nn,oo,1,XI_PLUS),dUdUvol_plus(:,:,PP_N,2))
    ! dependency of prolongated value (of neighbouring element, minus) to volume dofs next to interface
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,nn,oo,2,XI_PLUS),dUdUvol_plus(:,:,PP_N+1,1))
    ! dependency of prolongated value (of current element, plus) to volume dofs neighbouring the dofs next to interface
    r = vn2*oo + vn1*nn + PP_nVar*(PP_N-1)
    BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,nn,oo,1,XI_PLUS),dUdUvol_plus(:,:,PP_N,1))
#else
    !--------------------only direct dependencies of prolongated (copied) values to volume dofs next to interface----------------
    s = vn2*oo + vn1*nn
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*Df_DUinner(:,:,nn,oo,1,XI_MINUS) 
    s = vn2*oo + vn1*nn + PP_nVar*PP_N
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*Df_DUinner(:,:,nn,oo,1,XI_PLUS)
#endif
  END DO !nn
END DO !oo
! ETA-direction-------------------------------------------------------------------------------------------------------------------
DO oo = 0,PP_NZ
  DO mm = 0,PP_N
#if FV_RECONSTRUCT
    SideID=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
    Flip=ElemToSide(  E2S_FLIP   ,ETA_MINUS,iElem)
    pq=S2V2(:,mm,oo,flip         ,ETA_MINUS)
    IF (((SideID.GT.nBCSides).AND.(SideID.LT.firstInnerSide)).OR. &
        ((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN ! big mortar side
      Mortar_minus = .TRUE.
    ELSE
      Mortar_minus = .FALSE.
      IF(Flip.EQ.0)THEN
        UPrim_minus    = UPrim_master(:,pq(1),pq(2),SideID)
        UPrim_minus_nb = UPrim_slave( :,pq(1),pq(2),SideID)
        FV_dx_L        = FV_dx_master(1,pq(1),pq(2),SideID)
        IF(SideID.GT.nBCSides)THEN
          FV_dx_R_nb     = FV_dx_slave( 1,pq(1),pq(2),SideID)
        ELSE !BC
          FV_dx_R_nb     = FV_dx_master(1,pq(1),pq(2),SideID)
        END IF
      ELSE
        UPrim_minus    = UPrim_slave( :,pq(1),pq(2),SideID)
        UPrim_minus_nb = UPrim_master(:,pq(1),pq(2),SideID)
        FV_dx_L        = FV_dx_slave( 1,pq(1),pq(2),SideID)
        FV_dx_R_nb     = FV_dx_master(1,pq(1),pq(2),SideID)
      END IF
    END IF
    SideID=ElemToSide(E2S_SIDE_ID,ETA_PLUS,iElem)
    Flip=ElemToSide(  E2S_FLIP   ,ETA_PLUS,iElem)
    pq=S2V2(:,mm,oo,flip         ,ETA_PLUS)
    IF (((SideID.GT.nBCSides).AND.(SideID.LT.firstInnerSide)).OR. &
        ((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN ! big mortar side
      Mortar_plus = .TRUE.
    ELSE
      Mortar_plus = .FALSE.
      IF(Flip.EQ.0)THEN
        UPrim_plus    = UPrim_master(:,pq(1),pq(2),SideID)
        UPrim_plus_nb = UPrim_slave( :,pq(1),pq(2),SideID)
        FV_dx_R       = FV_dx_master(1,pq(1),pq(2),SideID)
        IF(SideID.GE.firstInnerSide)THEN
          FV_dx_L_nb    = FV_dx_slave( 1,pq(1),pq(2),SideID)
        ELSE !BC
          FV_dx_L_nb    = FV_dx_master(1,pq(1),pq(2),SideID)
        END IF
      ELSE
        UPrim_plus    = UPrim_slave( :,pq(1),pq(2),SideID)
        UPrim_plus_nb = UPrim_master(:,pq(1),pq(2),SideID)
        FV_dx_R       = FV_dx_slave( 1,pq(1),pq(2),SideID)
        FV_dx_L_nb    = FV_dx_master(1,pq(1),pq(2),SideID)
      END IF
    END IF
    CALL FV_Reconstruction_Derivative_Surf(FV_sdx_ETA_extended(mm,oo,:,iElem),FV_dx_L,FV_dx_R,             & 
                                           FV_dx_L_nb,FV_dx_R_nb,UPrim_plus,UPrim_minus,                   &
                                           UPrim_plus_nb,UPrim_minus_nb,UPrim_extended(:,mm,:,oo,iElem),   &
                                           Mortar_minus,Mortar_plus,dUdUvol_plus(:,:,:,:),dUdUvol_minus(:,:,:,:))
    s = vn2*oo + PP_nVar*mm
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,oo,1,ETA_MINUS),dUdUvol_minus(:,:,0,2))
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,oo,2,ETA_MINUS),dUdUvol_minus(:,:,-1,3))
    r = vn2*oo + PP_nVar*mm + vn1*1
    BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,oo,1,ETA_MINUS),dUdUvol_minus(:,:,0,3))
    s = vn2*oo + PP_nVar*mm + vn1*PP_N 
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,oo,1,ETA_PLUS),dUdUvol_plus(:,:,PP_N,2))
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,oo,2,ETA_PLUS),dUdUvol_plus(:,:,PP_N+1,1))
    r = vn2*oo + PP_nVar*mm + vn1*(PP_N-1) 
    BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,oo,1,ETA_PLUS),dUdUvol_plus(:,:,PP_N,1))
#else
    s = vn2*oo + PP_nVar*mm
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*Df_DUinner(:,:,mm,oo,1,ETA_MINUS) 
    s = vn2*oo + PP_nVar*mm + vn1*PP_N 
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*Df_DUinner(:,:,mm,oo,1,ETA_PLUS)
#endif
  END DO
END DO
#if PP_dim==3
! ZETA-direction------------------------------------------------------------------------------------------------------------------
DO nn = 0,PP_N
  DO mm = 0,PP_N
#if FV_RECONSTRUCT
    SideID=ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
    Flip=ElemToSide(  E2S_FLIP   ,ZETA_MINUS,iElem)
    pq=S2V2(:,mm,nn,flip         ,ZETA_MINUS)
    IF (((SideID.GT.nBCSides).AND.(SideID.LT.firstInnerSide)).OR. &
        ((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN ! big mortar side
      Mortar_minus = .TRUE.
    ELSE
      Mortar_minus = .FALSE.
      IF(Flip.EQ.0)THEN
        UPrim_minus    = UPrim_master(:,pq(1),pq(2),SideID)
        UPrim_minus_nb = UPrim_slave( :,pq(1),pq(2),SideID)
        FV_dx_L        = FV_dx_master(1,pq(1),pq(2),SideID)
        IF(SideID.GE.firstInnerSide)THEN
          FV_dx_R_nb     = FV_dx_slave( 1,pq(1),pq(2),SideID)
        ELSE
          FV_dx_R_nb     = FV_dx_master(1,pq(1),pq(2),SideID)
        END IF
      ELSE
        UPrim_minus    = UPrim_slave( :,pq(1),pq(2),SideID)
        UPrim_minus_nb = UPrim_master(:,pq(1),pq(2),SideID)
        FV_dx_L        = FV_dx_slave( 1,pq(1),pq(2),SideID)
        FV_dx_R_nb     = FV_dx_master(1,pq(1),pq(2),SideID)
      END IF
    END IF
    SideID=ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
    Flip=ElemToSide(  E2S_FLIP   ,ZETA_PLUS,iElem)
    pq=S2V2(:,mm,nn,flip         ,ZETA_PLUS)
    IF (((SideID.GT.nBCSides).AND.(SideID.LT.firstInnerSide)).OR. &
        ((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN ! big mortar side
      Mortar_plus = .TRUE.
    ELSE
      Mortar_plus = .FALSE.
      IF(Flip.EQ.0)THEN
        UPrim_plus    = UPrim_master(:,pq(1),pq(2),SideID)
        UPrim_plus_nb = UPrim_slave( :,pq(1),pq(2),SideID)
        FV_dx_R       = FV_dx_master(1,pq(1),pq(2),SideID)
        IF(SideID.GE.firstInnerSide)THEN
          FV_dx_L_nb    = FV_dx_slave( 1,pq(1),pq(2),SideID)
        ELSE
          FV_dx_L_nb    = FV_dx_master(1,pq(1),pq(2),SideID)
        END IF
      ELSE
        UPrim_plus    = UPrim_slave( :,pq(1),pq(2),SideID)
        UPrim_plus_nb = UPrim_master(:,pq(1),pq(2),SideID)
        FV_dx_R       = FV_dx_slave( 1,pq(1),pq(2),SideID)
        FV_dx_L_nb    = FV_dx_master(1,pq(1),pq(2),SideID)
      END IF
    END IF
    CALL FV_Reconstruction_Derivative_Surf(FV_sdx_ZETA_extended(mm,nn,:,iElem),FV_dx_L,FV_dx_R,            & 
                                           FV_dx_L_nb,FV_dx_R_nb,UPrim_plus,UPrim_minus,                   &
                                           UPrim_plus_nb,UPrim_minus_nb,UPrim_extended(:,mm,nn,:,iElem),   &
                                           Mortar_minus,Mortar_plus,dUdUvol_plus(:,:,:,:),dUdUvol_minus(:,:,:,:))
    s = vn1*nn + PP_nVar*mm
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,nn,1,ZETA_MINUS),dUdUvol_minus(:,:,0,2))
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,nn,2,ZETA_MINUS),dUdUvol_minus(:,:,-1,3))
    r = vn1*nn + PP_nVar*mm + vn2*1
    BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,nn,1,ZETA_MINUS),dUdUvol_minus(:,:,0,3))
    s = vn1*nn + PP_nVar*mm + vn2*PP_NZ
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,nn,1,ZETA_PLUS),dUdUvol_plus(:,:,PP_N,2))
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,nn,2,ZETA_PLUS),dUdUvol_plus(:,:,PP_N+1,1))
    r = vn1*nn + PP_nVar*mm + vn2*(PP_NZ-1)
    BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) = BJ(s+1:s+PP_nVar,r+1:r+PP_nVar) &
                                      + FV_w_inv*MATMUL(Df_DUinner(:,:,mm,nn,1,ZETA_PLUS),dUdUvol_plus(:,:,PP_N,1))
#else
    s = vn1*nn + PP_nVar*mm
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*Df_DUinner(:,:,mm,nn,1,ZETA_MINUS) 
    s = vn1*nn + PP_nVar*mm + vn2*PP_NZ
    BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) = BJ(s+1:s+PP_nVar,s+1:s+PP_nVar) &
                                      + FV_w_inv*Df_DUinner(:,:,mm,nn,1,ZETA_PLUS)
#endif
  END DO
END DO
#endif

#if PARABOLIC
IF (.NOT.HyperbolicPrecond) THEN
  ! Compute the following derivatives (depending on the FV gradient reconstruction): 
  !  * dependency of the volume gradients w.r.t. the volume DOFs dQVol_dUvol
  !  * dependency of the outer surface gradients w.r.t. the volume DOFs dQOuter_dUvol
  ! Note that dQInner_dUvol is not required as surface gradients are the same as the volume gradients of the outer layer. Hence,
  ! DQ_surf_L|_SurfInt/DU_Lprim is zero and only DQ_surf_L|_VolInt/DU_prim is present.
  CALL JacFVGradients_Vol(1,iElem,dQxVol_dU) !d(Q^1)/dU(1:3) (3 sums)
  CALL JacFVGradients_Vol(2,iElem,dQyVol_dU) !d(Q^1)/dU(1:3) (3 sums)
  CALL JacFVGradients_nb( 1,iElem,dQxOuter_dUvol)
  CALL JacFVGradients_nb( 2,iElem,dQyOuter_dUvol)
#if PP_dim==3
  CALL JacFVGradients_Vol(3,iElem,dQzVol_dU) !d(Q^1)/dU(1:3) (3 sums)
  CALL JacFVGradients_nb( 3,iElem,dQzOuter_dUvol)
#endif

  ! What is left now:
  !   DF^v/DF^v_surf * DF^v(U,Q(U))_surf/DQ_surf * DQ_surf/DU
  ! = DF^v/DF^v_surf * 0.5*(D(F^v_L)/DQ_surf_L * DQ_surf_L/DU + D(F^v_R)/DQ_surf_R * DQ_surf_R/DU)
  ! We first compute the dependecies of F^v_L on U:
  ! 0.5*(D(F^v_L)/DQ_surf_L * DQ_surf_L/DU_prim * DU_prim/DU
  !  => Df_dQxInner           => dQxVol_dU      => PrimConsJac
  ! How is the viscous flux in the volume depending on the volume DOFs via the left part of the flux
  !  * SplitXI_ETA
  ! Than we compute the dependencies of F^v_R on U:
  ! 0.5* D(F^v_R)/DQ_surf_R * DQ_surf_R/DU_prim * DU_prim/DU)
  !  => Df_dQxOuter           => dQxOuter_dUvol  => PrimConsJac
  ! How is the viscous flux in the volume depending on the volume DOFs via the right part of the flux
  ! Here, we do both parts together:
  !   DF^v/DF^v_surf * DF^v(U,Q(U))_surf/DQ_surf * DQ_surf/DU
  ! = DF^v/DF^v_surf * (Df_dQxInner*dQxVol_dU + Df_dQxOuter*dQxOuter_dUvol) * PrimConsJac
  ! === XI-Direction =============================================================================================================
  DO p=0,PP_N
    DO q=0,PP_NZ
      ! XI_MINUS
      r = vn1*p + vn2*q ! index of flux
      CALL dPrimTempdCons(UPrim(:,0,p,q,iElem),PrimConsJac(:,:,0,p,q))
      Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,XI_MINUS) * dQxVol_dU(0,p,q,0,1) &
                              +Df_dQxOuter(:,:,p,q,XI_MINUS) * dQxOuter_dUvol(p,q,XI_MINUS,0) &
                              +Df_dQyInner(:,:,p,q,XI_MINUS) * dQyVol_dU(0,p,q,0,1) &
                              +Df_dQyOuter(:,:,p,q,XI_MINUS) * dQyOuter_dUvol(p,q,XI_MINUS,0) &
#if PP_dim==3
                              +Df_dQzInner(:,:,p,q,XI_MINUS) * dQzVol_dU(0,p,q,0,1) &
                              +Df_dQzOuter(:,:,p,q,XI_MINUS) * dQzOuter_dUvol(p,q,XI_MINUS,0) &
#endif
                              )
      df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,0,p,q))
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_w_inv*df_dQ_minus

      s = vn1*p + vn2*q + PP_nVar*1 ! index of point right of flux
      CALL dPrimTempdCons(UPrim(:,1,p,q,iElem),PrimConsJac(:,:,1,p,q))
      Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,XI_MINUS) * dQxVol_dU(0,p,q,1,1) &
                              +Df_dQyInner(:,:,p,q,XI_MINUS) * dQyVol_dU(0,p,q,1,1) & 
#if PP_dim==3
                              +Df_dQzInner(:,:,p,q,XI_MINUS) * dQzVol_dU(0,p,q,1,1) &
#endif
                              )
      df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,1,p,q))
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_w_inv*df_dQ_minus

      IF(NoFillIn.EQV..FALSE.)THEN 
        ! dependency on gradient in eta-direction
        ss(1) = vn1*(p-1) + vn2*q ! index of point lower of flux (ETA)
        ss(2) = vn1*(p+1) + vn2*q ! index of point upper of flux (ETA)
        CALL Assemble_FVSurfIntGradJac(p,r,ss,UPrim(:,0,:,q,iElem),Df_dQxInner(:,:,p,q,XI_MINUS),Df_dQyInner(:,:,p,q,XI_MINUS), &
                                       dQxVol_dU(0,p,q,:,2),dQyVol_dU(0,p,q,:,2),                                               &
#if PP_dim==3
                                       Df_dQzInner(:,:,p,q,XI_MINUS),dQzVol_dU(0,p,q,:,2),                                      &
#endif
                                       BJ)
#if PP_dim==3
        ! dependency on gradient in zeta-direction
        ss(1) = vn1*p + vn2*(q-1) ! index of point lower of flux (ZETA)
        ss(2) = vn1*p + vn2*(q+1) ! index of point upper of flux (ZETA)
        CALL Assemble_FVSurfIntGradJac(q,r,ss,UPrim(:,0,p,:,iElem),Df_dQxInner(:,:,p,q,XI_MINUS),Df_dQyInner(:,:,p,q,XI_MINUS), &
                                       dQxVol_dU(0,p,q,:,3),dQyVol_dU(0,p,q,:,3),                                               &
                                       Df_dQzInner(:,:,p,q,XI_MINUS),dQzVol_dU(0,p,q,:,3),                                      &
                                       BJ)
#endif
      END IF

      ! XI_PLUS
      r = vn1*p + vn2*q +PP_nVar*PP_N ! index of flux
      CALL dPrimTempdCons(UPrim(:,PP_N,p,q,iElem),PrimConsJac(:,:,PP_N,p,q))
      Df_dQ_plus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,XI_PLUS) * dQxVol_dU(PP_N,p,q,PP_N,1) &
                             -Df_dQxOuter(:,:,p,q,XI_PLUS) * dQxOuter_dUvol(p,q,XI_PLUS,PP_N) &
                             +Df_dQyInner(:,:,p,q,XI_PLUS) * dQyVol_dU(PP_N,p,q,PP_N,1) &
                             -Df_dQyOuter(:,:,p,q,XI_PLUS) * dQyOuter_dUvol(p,q,XI_PLUS,PP_N) &
#if PP_dim==3
                             +Df_dQzInner(:,:,p,q,XI_PLUS) * dQzVol_dU(PP_N,p,q,PP_N,1) &
                             -Df_dQzOuter(:,:,p,q,XI_PLUS) * dQzOuter_dUvol(p,q,XI_PLUS,PP_N) &
#endif
                             )
      df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,PP_N,p,q))
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_w_inv*df_dQ_plus

      s = vn1*p + vn2*q + PP_nVar*(PP_N-1) ! index of point left of flux
      CALL dPrimTempdCons(UPrim(:,PP_N-1,p,q,iElem),PrimConsJac(:,:,PP_N-1,p,q))
      Df_dQ_plus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,XI_PLUS) * dQxVol_dU(PP_N,p,q,PP_N-1,1) &
                             +Df_dQyInner(:,:,p,q,XI_PLUS) * dQyVol_dU(PP_N,p,q,PP_N-1,1) & 
#if PP_dim==3
                             +Df_dQzInner(:,:,p,q,XI_PLUS) * dQzVol_dU(PP_N,p,q,PP_N-1,1) &
#endif
                             )
      df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,PP_N-1,p,q))
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_w_inv*df_dQ_plus

      IF(NoFillIn.EQV..FALSE.)THEN
        ! dependency on gradient in eta-direction
        ss(1) = vn1*(p-1) + vn2*q + PP_nVar*PP_N ! index of point lower of flux (ETA)
        ss(2) = vn1*(p+1) + vn2*q + PP_nVar*PP_N ! index of point upper of flux (ETA)
        CALL Assemble_FVSurfIntGradJac(p,r,ss,UPrim(:,PP_N,:,q,iElem),Df_dQxInner(:,:,p,q,XI_PLUS),Df_dQyInner(:,:,p,q,XI_PLUS), &
                                       dQxVol_dU(PP_N,p,q,:,2),dQyVol_dU(PP_N,p,q,:,2),                                          &
#if PP_dim==3
                                       Df_dQzInner(:,:,p,q,XI_PLUS),dQzVol_dU(PP_N,p,q,:,2),                                     &
#endif
                                       BJ)
#if PP_dim==3
        ! dependency on gradient in zeta-direction
        ss(1) = vn1*p + vn2*(q-1) + PP_nVar*PP_N ! index of point lower of flux (ZETA)
        ss(2) = vn1*p + vn2*(q+1) + PP_nVar*PP_N ! index of point upper of flux (ZETA)
        CALL Assemble_FVSurfIntGradJac(q,r,ss,UPrim(:,PP_N,p,:,iElem),Df_dQxInner(:,:,p,q,XI_PLUS),Df_dQyInner(:,:,p,q,XI_PLUS), &
                                       dQxVol_dU(PP_N,p,q,:,3),dQyVol_dU(PP_N,p,q,:,3),                                          &
                                       Df_dQzInner(:,:,p,q,XI_PLUS),dQzVol_dU(PP_N,p,q,:,3),                                     &
                                       BJ)
#endif
      END IF
    END DO !q
  END DO !p
  ! === ETA-Direction ============================================================================================================
  DO p=0,PP_N
    DO q=0,PP_NZ
      ! ETA_MINUS
      r = PP_nVar * p + vn2 * q  ! index of flux
      CALL dPrimTempdCons(UPrim(:,p,0,q,iElem),PrimConsJac(:,:,p,0,q))
      Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ETA_MINUS) * dQxVol_dU(p,0,q,0,2) &
                              +Df_dQxOuter(:,:,p,q,ETA_MINUS) * dQxOuter_dUvol(p,q,ETA_MINUS,0) &
                              +Df_dQyInner(:,:,p,q,ETA_MINUS) * dQyVol_dU(p,0,q,0,2) &
                              +Df_dQyOuter(:,:,p,q,ETA_MINUS) * dQyOuter_dUvol(p,q,ETA_MINUS,0) &
#if PP_dim==3                                                                                                       
                              +Df_dQzInner(:,:,p,q,ETA_MINUS) * dQzVol_dU(p,0,q,0,2) &
                              +Df_dQzOuter(:,:,p,q,ETA_MINUS) * dQzOuter_dUvol(p,q,ETA_MINUS,0) &
#endif
                              )
      df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,p,0,q))
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_w_inv*df_dQ_minus

      s = vn1*1 + vn2*q + PP_nVar*p ! index of point right of flux
      CALL dPrimTempdCons(UPrim(:,p,1,q,iElem),PrimConsJac(:,:,p,1,q))
      Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ETA_MINUS) * dQxVol_dU(p,0,q,1,2) &
                              +Df_dQyInner(:,:,p,q,ETA_MINUS) * dQyVol_dU(p,0,q,1,2) & 
#if PP_dim==3                                                                     
                              +Df_dQzInner(:,:,p,q,ETA_MINUS) * dQzVol_dU(p,0,q,1,2) &
#endif
                              )
      df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,p,1,q))
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_w_inv*df_dQ_minus

      IF(NoFillIn.EQV..FALSE.)THEN
        ! dependency on gradient in xi-direction
        ss(1) = PP_nVar * (p-1) + vn2 * q ! index of point lower of flux (XI)
        ss(2) = PP_nVar * (p+1) + vn2 * q ! index of point upper of flux (XI)
        CALL Assemble_FVSurfIntGradJac(p,r,ss,UPrim(:,:,0,q,iElem),Df_dQxInner(:,:,p,q,ETA_MINUS),Df_dQyInner(:,:,p,q,ETA_MINUS), &
                                       dQxVol_dU(p,0,q,:,1),dQyVol_dU(p,0,q,:,1),                                                 &
#if PP_dim==3
                                       Df_dQzInner(:,:,p,q,ETA_MINUS),dQzVol_dU(p,0,q,:,1),                                       &
#endif
                                       BJ)
#if PP_dim==3
        ! dependency on gradient in zeta-direction
        ss(1) = PP_nVar * p + vn2 * (q-1) ! index of point lower of flux (ZETA)
        ss(2) = PP_nVar * p + vn2 * (q+1) ! index of point upper of flux (ZETA)
        CALL Assemble_FVSurfIntGradJac(q,r,ss,UPrim(:,p,0,:,iElem),Df_dQxInner(:,:,p,q,ETA_MINUS),Df_dQyInner(:,:,p,q,ETA_MINUS), &
                                       dQxVol_dU(p,0,q,:,3),dQyVol_dU(p,0,q,:,3),                                                 &
                                       Df_dQzInner(:,:,p,q,ETA_MINUS),dQzVol_dU(p,0,q,:,3),                                       &
                                       BJ)
#endif
      END IF

      ! ETA_PLUS
      r = PP_nVar * p + vn1 * PP_N + vn2 * q ! index of flux
      CALL dPrimTempdCons(UPrim(:,p,PP_N,q,iElem),PrimConsJac(:,:,p,PP_N,q))
      Df_dQ_plus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ETA_PLUS) * dQxVol_dU(p,PP_N,q,PP_N,2) &
                             -Df_dQxOuter(:,:,p,q,ETA_PLUS) * dQxOuter_dUvol(p,q,ETA_PLUS,PP_N) &
                             +Df_dQyInner(:,:,p,q,ETA_PLUS) * dQyVol_dU(p,PP_N,q,PP_N,2) &
                             -Df_dQyOuter(:,:,p,q,ETA_PLUS) * dQyOuter_dUvol(p,q,ETA_PLUS,PP_N) &
#if PP_dim==3
                             +Df_dQzInner(:,:,p,q,ETA_PLUS) * dQzVol_dU(p,PP_N,q,PP_N,2) &
                             -Df_dQzOuter(:,:,p,q,ETA_PLUS) * dQzOuter_dUvol(p,q,ETA_PLUS,PP_N) &
#endif
                             )
      df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,p,PP_N,q))
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_w_inv*df_dQ_plus

      s = vn1*(PP_N-1) + vn2*q + PP_nVar*p ! index of point left of flux
      CALL dPrimTempdCons(UPrim(:,p,PP_N-1,q,iElem),PrimConsJac(:,:,p,PP_N-1,q))
      Df_dQ_plus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ETA_PLUS) * dQxVol_dU(p,PP_N,q,PP_N-1,2) &
                             +Df_dQyInner(:,:,p,q,ETA_PLUS) * dQyVol_dU(p,PP_N,q,PP_N-1,2) & 
#if PP_dim==3                                                                      
                             +Df_dQzInner(:,:,p,q,ETA_PLUS) * dQzVol_dU(p,PP_N,q,PP_N-1,2) &
#endif
                             )
      df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,p,PP_N-1,q))
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_w_inv*df_dQ_plus

      IF(NoFillIn.EQV..FALSE.)THEN
        ! dependency on gradient in xi-direction
        ss(1) = PP_nVar * (p-1) + vn2 * q + vn1 * PP_N ! index of point lower of flux (XI)
        ss(2) = PP_nVar * (p+1) + vn2 * q + vn1 * PP_N ! index of point upper of flux (XI)
        CALL Assemble_FVSurfIntGradJac(p,r,ss,UPrim(:,:,PP_N,q,iElem),Df_dQxInner(:,:,p,q,ETA_PLUS),Df_dQyInner(:,:,p,q,ETA_PLUS), &
                                       dQxVol_dU(p,PP_N,q,:,1),dQyVol_dU(p,PP_N,q,:,1),                                            &
#if PP_dim==3
                                       Df_dQzInner(:,:,p,q,ETA_PLUS),dQzVol_dU(p,PP_N,q,:,1),                                      &
#endif
                                       BJ)
#if PP_dim==3
        ! dependency on gradient in zeta-direction
        ss(1) = PP_nVar * p + vn2 * (q-1) + vn1 * PP_N ! index of point lower of flux (ZETA)
        ss(2) = PP_nVar * p + vn2 * (q+1) + vn1 * PP_N ! index of point upper of flux (ZETA)
        CALL Assemble_FVSurfIntGradJac(q,r,ss,UPrim(:,p,PP_N,:,iElem),Df_dQxInner(:,:,p,q,ETA_PLUS),Df_dQyInner(:,:,p,q,ETA_PLUS), &
                                       dQxVol_dU(p,PP_N,q,:,3),dQyVol_dU(p,PP_N,q,:,3),                                            &
                                       Df_dQzInner(:,:,p,q,ETA_PLUS),dQzVol_dU(p,PP_N,q,:,3),                                      &
                                       BJ)
#endif
      END IF
    END DO !q
  END DO !p
#if PP_dim==3
  ! === ZETA-Direction ===========================================================================================================
  DO p=0,PP_N
    DO q=0,PP_N
      ! ZETA_MINUS
      r = PP_nVar * p + vn1 * q  ! index of flux
      CALL dPrimTempdCons(UPrim(:,p,q,0,iElem),PrimConsJac(:,:,p,q,0))
      Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ZETA_MINUS) * dQxVol_dU(p,q,0,0,3) &
                              +Df_dQxOuter(:,:,p,q,ZETA_MINUS) * dQxOuter_dUvol(p,q,ZETA_MINUS,0) &
                              +Df_dQyInner(:,:,p,q,ZETA_MINUS) * dQyVol_dU(p,q,0,0,3) &
                              +Df_dQyOuter(:,:,p,q,ZETA_MINUS) * dQyOuter_dUvol(p,q,ZETA_MINUS,0) &
                              +Df_dQzInner(:,:,p,q,ZETA_MINUS) * dQzVol_dU(p,q,0,0,3) &
                              +Df_dQzOuter(:,:,p,q,ZETA_MINUS) * dQzOuter_dUvol(p,q,ZETA_MINUS,0) &
                              )
      df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,p,q,0))
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_w_inv*df_dQ_minus

      s = vn2*1 + vn1*q + PP_nVar*p ! index of point right of flux
      CALL dPrimTempdCons(UPrim(:,p,q,1,iElem),PrimConsJac(:,:,p,q,1))
      Df_dQ_minus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ZETA_MINUS) * dQxVol_dU(p,q,0,1,3) &
                              +Df_dQyInner(:,:,p,q,ZETA_MINUS) * dQyVol_dU(p,q,0,1,3) & 
                              +Df_dQzInner(:,:,p,q,ZETA_MINUS) * dQzVol_dU(p,q,0,1,3) &
                              )
      df_dQ_minus(:,:) = MATMUL(df_dQ_minus_Tilde(:,:),PrimConsJac(:,:,p,q,1))
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_w_inv*df_dQ_minus

      IF(NoFillIn.EQV..FALSE.)THEN
        ! dependency on gradient in xi-direction
        ss(1) = PP_nVar * (p-1) + vn1 * q ! index of point lower of flux (XI)
        ss(2) = PP_nVar * (p+1) + vn1 * q ! index of point upper of flux (XI)
        CALL Assemble_FVSurfIntGradJac(p,r,ss,UPrim(:,:,q,0,iElem),Df_dQxInner(:,:,p,q,ZETA_MINUS),Df_dQyInner(:,:,p,q,ZETA_MINUS), &
                                       dQxVol_dU(p,q,0,:,1),dQyVol_dU(p,q,0,:,1),                                                   &
                                       Df_dQzInner(:,:,p,q,ZETA_MINUS),dQzVol_dU(p,q,0,:,1),                                        &
                                       BJ)
        ! dependency on gradient in eta-direction
        ss(1) = PP_nVar * p + vn1 * (q-1) ! index of point lower of flux (ETA)
        ss(2) = PP_nVar * p + vn1 * (q+1) ! index of point upper of flux (ETA)
        CALL Assemble_FVSurfIntGradJac(q,r,ss,UPrim(:,p,:,0,iElem),Df_dQxInner(:,:,p,q,ZETA_MINUS),Df_dQyInner(:,:,p,q,ZETA_MINUS), &
                                       dQxVol_dU(p,q,0,:,2),dQyVol_dU(p,q,0,:,2),                                                   &
                                       Df_dQzInner(:,:,p,q,ZETA_MINUS),dQzVol_dU(p,q,0,:,2),                                        &
                                       BJ)
      END IF

      ! ZETA_PLUS
      r = PP_nVar * p + vn2 * PP_N + vn1 * q ! index of flux
      CALL dPrimTempdCons(UPrim(:,p,q,PP_N,iElem),PrimConsJac(:,:,p,q,PP_N))
      Df_dQ_plus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ZETA_PLUS) * dQxVol_dU(p,q,PP_N,PP_N,3) &
                             -Df_dQxOuter(:,:,p,q,ZETA_PLUS) * dQxOuter_dUvol(p,q,ZETA_PLUS,PP_N) &
                             +Df_dQyInner(:,:,p,q,ZETA_PLUS) * dQyVol_dU(p,q,PP_N,PP_N,3) &
                             -Df_dQyOuter(:,:,p,q,ZETA_PLUS) * dQyOuter_dUvol(p,q,ZETA_PLUS,PP_N) &
                             +Df_dQzInner(:,:,p,q,ZETA_PLUS) * dQzVol_dU(p,q,PP_N,PP_N,3) &
                             -Df_dQzOuter(:,:,p,q,ZETA_PLUS) * dQzOuter_dUvol(p,q,ZETA_PLUS,PP_N) &
                             )
      df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,p,q,PP_N))
      BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) = BJ(r+1:r+PP_nVar,r+1:r+PP_nVar) + FV_w_inv*df_dQ_plus

      s = vn2*(PP_N-1) + vn1*q + PP_nVar*p ! index of point left of flux
      CALL dPrimTempdCons(UPrim(:,p,q,PP_N-1,iElem),PrimConsJac(:,:,p,q,PP_N-1))
      Df_dQ_plus_Tilde(:,:)= (Df_dQxInner(:,:,p,q,ZETA_PLUS) * dQxVol_dU(p,q,PP_N,PP_N-1,3) &
                             +Df_dQyInner(:,:,p,q,ZETA_PLUS) * dQyVol_dU(p,q,PP_N,PP_N-1,3) & 
                             +Df_dQzInner(:,:,p,q,ZETA_PLUS) * dQzVol_dU(p,q,PP_N,PP_N-1,3) &
                             )
      df_dQ_plus(:,:) = MATMUL(df_dQ_plus_Tilde(:,:),PrimConsJac(:,:,p,q,PP_N-1))
      BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + FV_w_inv*df_dQ_plus

      IF(NoFillIn.EQV..FALSE.)THEN
        ! dependency on gradient in xi-direction
        ss(1) = PP_nVar * (p-1) + vn1 * q + vn2 * PP_N ! index of point lower of flux (XI)
        ss(2) = PP_nVar * (p+1) + vn1 * q + vn2 * PP_N ! index of point upper of flux (XI)
        CALL Assemble_FVSurfIntGradJac(p,r,ss,UPrim(:,:,q,PP_N,iElem),Df_dQxInner(:,:,p,q,ZETA_PLUS),Df_dQyInner(:,:,p,q,ZETA_PLUS), &
                                       dQxVol_dU(p,q,PP_N,:,1),dQyVol_dU(p,q,PP_N,:,1),                                              &
                                       Df_dQzInner(:,:,p,q,ZETA_PLUS),dQzVol_dU(p,q,PP_N,:,1),                                       &
                                       BJ)
        ! dependency on gradient in eta-direction
        ss(1) = PP_nVar * p + vn1 * (q-1) + vn2 * PP_N ! index of point lower of flux (ETA)
        ss(2) = PP_nVar * p + vn1 * (q+1) + vn2 * PP_N ! index of point upper of flux (ETA)
        CALL Assemble_FVSurfIntGradJac(q,r,ss,UPrim(:,p,:,PP_N,iElem),Df_dQxInner(:,:,p,q,ZETA_PLUS),Df_dQyInner(:,:,p,q,ZETA_PLUS), &
                                       dQxVol_dU(p,q,PP_N,:,2),dQyVol_dU(p,q,PP_N,:,2),                                              &
                                       Df_dQzInner(:,:,p,q,ZETA_PLUS),dQzVol_dU(p,q,PP_N,:,2),                                       &
                                       BJ)
      END IF
    END DO !q
  END DO !p
#endif
END IF !HyperbolicPrecond
#endif /*PARABOLIC*/
END SUBROUTINE Assemble_JacSurfInt_FV

#if PARABOLIC
!===================================================================================================================================
!> This routine assembles the dependency of the viscous flux on the volume solution via the surface integral, and the gradients.
!> It is called only for the cross dependencies: e.g. XI-Flux dependency via ETA-gradient
!>   DF^v/DF^v_surf * DF^v(U,Q(U))_surf/DQ_surf * DQ_surf/DU
!> = DF^v/DF^v_surf * (Df_dQxInner*dQxVol_dU + Df_dQxOuter*dQxOuter_dUvol) * PrimConsJac
!>     => FV_w_inv      =>df_dQ     =>JacQ                    => 0         =>PrimConsJac
!> dQxOuter_dUvol is zero as e.g. ETA gradient on neighbour element does not depend on current element for XI-sides
!> _________________________
!>             |
!>             |
!>       x     |    x
!>       :     |    :
!> ______:_____|____:_______
!>  dep. :     |    :dependency
!>       :    F_XI  :
!>  Q_nb_ETA   =>   Q_ETA
!>       :     |    :
!>  dep. :     |    :dependency
!> ______:_____|____:_______
!>       :     |    :
!>       :     |    :
!>       x     |    x
!>             |
!> ____________|____________
!===================================================================================================================================
SUBROUTINE Assemble_FVSurfIntGradJac(i,r,ss,UPrim,df_dQx,df_dQy,JacQx,JacQy, &
#if PP_dim==3
                                    df_dQz,JacQz, &
#endif
                                    BJ)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Implicit_Vars ,ONLY: nDOFVarElem
USE MOD_Jacobian      ,ONLY: dPrimTempdCons
USE MOD_FV_Vars       ,ONLY: FV_w_inv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: i              !< index in flux direction
INTEGER,INTENT(IN)                                 :: r              !< one dimensional index of considered flux
INTEGER,INTENT(IN)                                 :: ss(2)          !< one dimensional index of considered DOFs: 1: below, 2:above
REAL,DIMENSION(PP_nVarPrim,0:PP_N)     ,INTENT(IN) :: UPrim          !< primitive solution along i index
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim),INTENT(IN) :: df_dQx,df_dQy  !< flux Jacobian w.r.t gradients in x-,y-direction
REAL,DIMENSION(0:PP_N)                             :: JacQx,JacQy    !< Jacobian of gradients in x-,y-direction w.r.t. solution
                                                                     !> along i index
#if PP_dim==3
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim),INTENT(IN) :: df_dQz         !< flux Jacobian w.r.t gradients in z-direction
REAL,DIMENSION(0:PP_N)                             :: JacQz          !< Jacobian of gradients in z-dir w.r.t. solution along i index
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:nDOFVarElem,1:nDOFVarElem),INTENT(INOUT) :: BJ      !< block-Jacobian of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                 :: ii,j
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim) :: JacTilde
REAL,DIMENSION(PP_nVar    ,PP_nVar)     :: Jac
REAL,DIMENSION(PP_nVarPrim,PP_nVar)     :: PrimConsJac
!===================================================================================================================================
j = 1
DO ii=i-1,i+1,2
  IF(ii.LT.0)THEN
    j = j+1
    CYCLE
  END IF
  IF(ii.GT.PP_N) CYCLE
  CALL dPrimTempdCons(UPrim(:,ii),PrimConsJac)
  JacTilde(:,:) = (Df_dQx * JacQx(ii) &
                  +Df_dQy * JacQy(ii) & 
#if PP_dim==3
                  +Df_dQz * JacQz(ii) &
#endif
                  )
  Jac = MATMUL(JacTilde(:,:),PrimConsJac)
  BJ(r+1:r+PP_nVar,ss(j)+1:ss(j)+PP_nVar) = BJ(r+1:r+PP_nVar,ss(j)+1:ss(j)+PP_nVar) + FV_w_inv*Jac
  j = j+1
END DO
END SUBROUTINE Assemble_FVSurfIntGradJac
#endif
#endif

END MODULE MOD_JacSurfInt
