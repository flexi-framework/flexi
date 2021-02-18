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

!==================================================================================================================================
!> Contains routines that fill the inter element and boundary surface fluxes for the BR1 scheme.
!==================================================================================================================================
MODULE MOD_Lifting_FillFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE Lifting_FillFlux
  MODULE PROCEDURE Lifting_FillFlux
END INTERFACE

INTERFACE Lifting_FillFlux_BC
  MODULE PROCEDURE Lifting_FillFlux_BC
END INTERFACE

INTERFACE Lifting_FillFlux_NormVec
  MODULE PROCEDURE Lifting_FillFlux_NormVec
END INTERFACE

PUBLIC::Lifting_FillFlux,Lifting_FillFlux_BC,Lifting_FillFlux_NormVec
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> \brief Computes the untransformed BR1 surface fluxes used to compute the gradient in direction "dir",
!> surfelem contribution is considered as well.
!>
!> Fills the untransformed interior surface fluxes for the BR1 scheme. The numerical flux in the
!> BR1 lifting is simply taken as the arithmetic mean of the solution.
!> The physical flux will be multiplied by the surface element contribution and in a later routine by the normal vector to
!> transform into reference space,
!> see e.g. "Explicit discontinuous Galerkin methods for unsteady problems" (Hindenlang et al. 2012) for details.
!> For the strong form, in the surface integral the inner solution is substracted form the numerical flux. Since the numerical
!> flux is \f$ \frac{1}{2} (U^+ + U^-) \f$ and the inner solution is simply \f$ U^- \f$ for the master side
!> , the surface flux will become \f$ \frac{1}{2} (U^+ + U^-) - U^- = \frac{1}{2} (U^+ - U^-) \f$ for the strong form.
!>
!> The flux is filled for the master side, the contribution for the slave side (which is different because the inner solution
!> is equal to \f$ U^+ \f$) is taken into account in the SurfInt routine.
!==================================================================================================================================
PPURE SUBROUTINE Lifting_FillFlux(ULiftface_master,ULiftface_slave,Flux,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars,    ONLY: doWeakLifting
USE MOD_Mesh_Vars,       ONLY: SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,   lastInnerSide
USE MOD_Mesh_Vars,       ONLY: firstMPISide_MINE,lastMPISide_MINE
#if FV_ENABLED
USE MOD_FV_Vars         ,ONLY: FV_Elems_Sum,FV_sVdm
USE MOD_ChangeBasisByDim,ONLY: ChangeBasisSurf
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides                                             !< = .TRUE. only MINE MPISides are filled,
                                                                             !< =.FALSE. InnerSides
REAL,INTENT(INOUT) :: ULiftface_master(PP_nVarLifting,0:PP_N,0:PP_NZ,nSides) !< Solution on master sides
REAL,INTENT(INOUT) :: ULiftface_slave( PP_nVarLifting,0:PP_N,0:PP_NZ,nSides) !< Solution on slave sides
REAL,INTENT(INOUT) :: Flux(            PP_nVarLifting,0:PP_N,0:PP_NZ,nSides) !< Untransformed lifting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID,sig
#if FV_ENABLED
REAL               :: UPrim_glob(1:PP_nVarLifting,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide
   lastSideID =  lastInnerSide
END IF

! for strong form subtract solution from the inside
sig = MERGE(1, -1, doWeakLifting)

DO SideID=firstSideID,lastSideID
#if FV_ENABLED
  SELECT CASE(FV_Elems_Sum(SideID))
  CASE(0) ! both DG
#endif
    Flux(:,:,:,SideID) = sig*ULiftface_master(:,:,:,SideID) + ULiftface_slave(:,:,:,SideID)
#if FV_ENABLED
  CASE(1) ! master=FV, slave=DG
    CALL ChangeBasisSurf(PP_nVarLifting,PP_N,PP_N,FV_sVdm,ULiftface_master(:,:,:,SideID),ULift_glob)
    Flux(:,:,:,SideID) = sig*UPrim_glob + ULiftface_slave(:,:,:,SideID)
  CASE(2) ! master=DG, slave=FV
    CALL ChangeBasisSurf(PP_nVarLifting,PP_N,PP_N,FV_sVdm,ULiftface_slave(:,:,:,SideID),ULift_glob)
    Flux(:,:,:,SideID) = sig*ULiftface_master(:,:,:,SideID) + UPrim_glob
  CASE(3) ! both FV
    CYCLE
  END SELECT
#endif

  ! BR1 uses arithmetic mean value of states for the Riemann flux
  ! Transformation with NormVec is applied in separate routine
  DO q=0,PP_NZ; DO p=0,PP_N
    Flux(:,p,q,SideID) = 0.5*Flux(:,p,q,SideID) * SurfElem(p,q,0,SideID)
  END DO; END DO
END DO ! SideID

END SUBROUTINE Lifting_FillFlux



!==================================================================================================================================
!> \brief Computes the BR1 surface fluxes for boundary conditions, surfelem contribution is considered as well.
!>
!> Fills the Untransformed boundary surface fluxes for the BR1 scheme. The numerical flux in the
!> BR1 lifting is simply taken as the arithmetic mean of the solution.
!> This routine calls the equation system dependant routine Lifting_GetBoundaryFlux which will fill the flux depending on the
!> boundary condition that has to be applied. The Lifting_GetBoundaryFlux routine will also differentiate between weak and
!> strong form and already multiply the flux by the surface element.
!>
!> The multiplication by the normal vector is done in a separate routine.
!>
!> The flux is filled for the master side, the contribution for the slave side (which is different because the inner solution
!> is equal to \f$ U^+ \f$) is taken into account in the SurfInt routine.
!==================================================================================================================================
SUBROUTINE Lifting_FillFlux_BC(t,UPrim_master,Flux)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,TangVec1,TangVec2,Face_xGP,SurfElem,nSides,nBCSides
USE MOD_GetBoundaryFlux, ONLY: Lifting_GetBoundaryFlux
#if FV_ENABLED
USE MOD_FV_Vars         ,ONLY: FV_Elems_master
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: t                                                  !< Current solution time
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim   ,0:PP_N,0:PP_NZ,nSides) !< Primitive solution from the inside
REAL,INTENT(INOUT) :: Flux(        PP_nVarLifting,0:PP_N,0:PP_NZ,nSides) !< Untransformed lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID
!==================================================================================================================================
! compute untransformed flux for boundary sides, normvec is applied later
DO SideID=1,nBCSides
  IF (FV_Elems_master(SideID).GT.0) CYCLE
  CALL Lifting_GetBoundaryFlux(SideID,t,UPrim_master(:,:,:,SideID),Flux(:,:,:,SideID),&
       NormVec(:,:,:,0,SideID),TangVec1(:,:,:,0,SideID),TangVec2(:,:,:,0,SideID),Face_xGP(:,:,:,0,SideID),SurfElem(:,:,0,SideID))
END DO
END SUBROUTINE Lifting_FillFlux_BC


!==================================================================================================================================
!> \brief Multiplies the untransformed flux (already weighted by the surface element) by the normal vector to complete the
!>        transformation into reference space.
!>
!> The untransformed fluxes is multiplied for each gradient direction by the corresponding normal vector,
!> which is the second step of the transformation. We end up with the lifting fluxes for X/Y/Z direction.
!==================================================================================================================================
PPURE SUBROUTINE Lifting_FillFlux_NormVec(Flux,FluxX,FluxY,FluxZ,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,nSides,MortarType
USE MOD_Mesh_Vars,       ONLY: firstMPISide_YOUR,lastMPISide_MINE,lastMPISide_YOUR
#if FV_ENABLED
USE MOD_FV_Vars,         ONLY: FV_Elems_master,FV_Elems_Sum
USE MOD_Mesh_Vars,       ONLY: nBCSides
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides                                 !< = .TRUE. MPI_YOUR
                                                                 !< =.FALSE. BC+Inner+MPI_Mine
REAL,INTENT(IN)    :: Flux( PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides) !< Untransformed Lifting boundary flux
REAL,INTENT(INOUT) :: FluxX(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides) !< Lifting boundary flux in x direction
REAL,INTENT(INOUT) :: FluxY(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides) !< Lifting boundary flux in y direction
REAL,INTENT(INOUT) :: FluxZ(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides) !< Lifting boundary flux in z direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: firstSideID,lastSideID,SideID,p,q
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! apply normal vector to YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID =  lastMPISide_YOUR
ELSE
  ! apply normal vector to all MINE fluxes
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
#if FV_ENABLED
  IF (FV_Elems_master(SideID).GT.0.AND.SideID.LE.nBCSides) CYCLE ! FV BC already filled
  IF (FV_Elems_Sum(SideID).EQ.3)                           CYCLE ! FV/FV already filled
#endif
  IF(MortarType(1,SideID).GT.0)                            CYCLE ! no big mortars
  DO q=0,PP_NZ; DO p=0,PP_N
    FluxX(:,p,q,SideID)=Flux(:,p,q,SideID)*NormVec(1,p,q,0,SideID)
    FluxY(:,p,q,SideID)=Flux(:,p,q,SideID)*NormVec(2,p,q,0,SideID)
#if (PP_dim==3)
    FluxZ(:,p,q,SideID)=Flux(:,p,q,SideID)*NormVec(3,p,q,0,SideID)
#else
    FluxZ(:,p,q,SideID)=0.
#endif
  END DO; END DO
END DO

END SUBROUTINE Lifting_FillFlux_NormVec

END MODULE MOD_Lifting_FillFlux
#endif /*PARABOLIC*/
