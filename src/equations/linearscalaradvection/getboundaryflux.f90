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

!==================================================================================================================================
!> Routines to provide boundary conditions for the domain. Fills the boundary part of the fluxes list.
!==================================================================================================================================
MODULE MOD_GetBoundaryFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

INTERFACE GetBoundaryState
  MODULE PROCEDURE GetBoundaryState
END INTERFACE

#if FV_ENABLED
#if FV_RECONSTRUCT
INTERFACE GetBoundaryFVgradient
  MODULE PROCEDURE GetBoundaryFVgradient
END INTERFACE
#endif
#endif

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

#if PARABOLIC
INTERFACE Lifting_GetBoundaryFlux
  MODULE PROCEDURE Lifting_GetBoundaryFlux
END INTERFACE
PUBLIC :: Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/

PUBLIC :: InitBC
PUBLIC :: GetBoundaryFlux
PUBLIC :: GetBoundaryState
#if FV_ENABLED
#if FV_RECONSTRUCT
PUBLIC :: GetBoundaryFVgradient
#endif
#endif
PUBLIC :: FinalizeBC
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Initialize boundary conditions. Read parameters and sort boundary conditions by types.
!> Call boundary condition specific init routines.
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: BCData,nBCByType,BCSideID
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,nBCs
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState
INTEGER :: MaxBCState,MaxBCStateGlobal
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
  CALL CollectiveStop(__STAMP__,&
    "InitBC not ready to be called or already called.")
END IF
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
END DO
MaxBCStateGLobal=MaxBCState
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

! Allocate buffer array to store temp data for all BC sides
ALLOCATE(BCData(PP_nVar,0:PP_N,0:PP_N,nBCSides))
BCData=0.

! Initialize State File Boundary condition
DO i=1,nBCs
  locType =BoundaryType(i,BC_TYPE)
END DO

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i) nBCByType(i)=nBCByType(i)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i)THEN
      nBCByType(i)=nBCByType(i)+1
      BCSideID(i,nBCByType(i))=iSide
    END IF
  END DO
END DO

END SUBROUTINE InitBC

!==================================================================================================================================
!> Computes the boundary state for the different boundary conditions.
!==================================================================================================================================
SUBROUTINE GetBoundaryState(SideID,t,Nloc,UPrim_boundary,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BoundaryType,BC
USE MOD_ExactFunc    ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: SideID                                          !< ID of current side
REAL,INTENT(IN)         :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)      :: Nloc    !< polynomial degree
REAL,INTENT(IN)         :: UPrim_master(  PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< inner surface solution
REAL,INTENT(IN)         :: NormVec(                 3,0:Nloc,0:ZDIM(Nloc)) !< normal surface vectors
REAL,INTENT(IN)         :: TangVec1(                3,0:Nloc,0:ZDIM(Nloc)) !< tangent surface vectors 1
REAL,INTENT(IN)         :: TangVec2(                3,0:Nloc,0:ZDIM(Nloc)) !< tangent surface vectors 2
REAL,INTENT(IN)         :: Face_xGP(                3,0:Nloc,0:ZDIM(Nloc)) !< positions of surface flux points
REAL,INTENT(OUT)        :: UPrim_boundary(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< resulting boundary state

! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: p,q
INTEGER                 :: BCType
!===================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
SELECT CASE(BCType)
CASE(2) !Exact function or refstate
  DO q=0,ZDIM(Nloc); DO p=0,Nloc
    CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q),UPrim_boundary(:,p,q))
  END DO; END DO
CASE(1) !Periodic already filled!
  CALL Abort(__STAMP__, &
      "GetBoundaryState called for periodic side!")
CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in linearscalaradvection/getboundaryflux.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryState

!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(SideID,t,Nloc,Flux,UPrim_master,          &
#if PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BC,BoundaryType
USE MOD_Riemann      ,ONLY: GetFlux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                   :: SideID
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: Nloc    !< polynomial degree
REAL,INTENT(IN)                      :: UPrim_master( PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< inner surface solution
#if PARABOLIC
                                                                           !> inner surface solution gradients in x/y/z-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx_master,gradUy_master,gradUz_master
#endif /*PARABOLIC*/
                                                                           !> normal and tangential vectors on surfaces
REAL,DIMENSION(      3,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: NormVec,TangVec1,TangVec2
REAL,DIMENSION(      3,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: Face_xGP   !< positions of surface flux points
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: Flux       !< resulting boundary fluxes
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: BCType,BCState
REAL                                 :: UPrim_boundary(PP_nVarPrim,0:Nloc,0:Nloc)
!==================================================================================================================================
BCType =BoundaryType(BC(SideID),BC_TYPE)
BCState=BoundaryType(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!
CASE(2) !Exact function or refstate
  CALL GetBoundaryState(SideID,t,Nloc,UPrim_boundary,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
  CALL GetFlux(Nloc,Flux,UPrim_master,UPrim_boundary,    &
#if PARABOLIC
               gradUx_master,gradUy_master,gradUz_master,&
               gradUx_master,gradUy_master,gradUz_master,&
#endif /*PARABOLIC*/
               NormVec,TangVec1,TangVec2,doBC=.TRUE.)
CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in linearscalaradvection/getboundaryflux.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryFlux

#if FV_ENABLED
#if FV_RECONSTRUCT
!==================================================================================================================================
!> Computes the gradient at a boundary for FV subcells.
!==================================================================================================================================
SUBROUTINE GetBoundaryFVgradient(SideID,t,gradU,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP,sdx_Face)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: Abort
USE MOD_Mesh_Vars     ,ONLY: BoundaryType,BC
USE MOD_TestCase      ,ONLY: GetBoundaryFVgradientTestcase
USE MOD_Exactfunc     ,ONLY: ExactFunc
USE MOD_Equation_Vars ,ONLY: IniExactFunc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(IN)   :: t
REAL,INTENT(IN)   :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL,INTENT(OUT)  :: gradU       (PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: NormVec (              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: TangVec1(              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: TangVec2(              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: Face_xGP(              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: sdx_Face(                0:PP_N,0:PP_NZ,3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: p,q
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVarPrim)
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

IF (BCType.LT.0) THEN ! testcase boundary condition
  CALL GetBoundaryFVgradientTestcase(SideID,t,gradU,UPrim_master)
ELSE
  SELECT CASE(BCType)
  CASE(2) ! exact BC = Dirichlet BC !!
    ! Determine the exact BC state
    DO q=0,PP_NZ; DO p=0,PP_N
      CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q),UPrim_boundary)
      gradU(:,p,q) = (UPrim_master(:,p,q) - UPrim_boundary) * sdx_Face(p,q,3)
    END DO ; END DO

  CASE(1) !Periodic already filled!
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in linearscalaradvection/getboundaryfvgradient.f90!')
  END SELECT ! BCType
END IF
END SUBROUTINE GetBoundaryFVgradient
#endif
#endif


#if PARABOLIC
!==================================================================================================================================
!> Computes the boundary fluxes for the lifting procedure for a given mesh face (defined by SideID).
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(SideID,t,UPrim_master,Flux,NormVec,TangVec1,TangVec2,Face_xGP,SurfElem)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BoundaryType,BC
USE MOD_Lifting_Vars ,ONLY: doWeakLifting
USE MOD_TestCase     ,ONLY: Lifting_GetBoundaryFluxTestcase
USE MOD_Exactfunc    ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(IN)   :: t                                       !< current time (provided by time integration scheme)
REAL,INTENT(IN)   :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)  :: Flux(     PP_nVarLifting,0:PP_N,0:PP_NZ) !< lifting boundary flux
REAL,INTENT(IN)   :: NormVec (              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: TangVec1(              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: TangVec2(              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: Face_xGP(              3,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)   :: SurfElem(                0:PP_N,0:PP_NZ)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q
INTEGER           :: BCType,BCState
REAL              :: UPrim_boundary(PP_nVarPrim,0:PP_N,0:PP_NZ)
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

IF (BCType.LT.0) THEN ! testcase boundary conditions
  CALL Lifting_GetBoundaryFluxTestcase(SideID,t,UPrim_master,Flux)
ELSE
  SELECT CASE(BCType)
  CASE(2)
    DO q=0,PP_NZ; DO p=0,PP_N
      CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q),UPrim_boundary(:,p,q))
    END DO; END DO
    Flux=0.5*(UPrim_master+UPrim_boundary)
  CASE(1) !Periodic already filled!
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in navierstokes/getboundaryflux.f90!')
  END SELECT

  !in case lifting is done in strong form
  IF(.NOT.doWeakLifting) Flux=Flux-UPrim_master

  DO q=0,PP_NZ; DO p=0,PP_N
    Flux(:,p,q)=Flux(:,p,q)*SurfElem(p,q)
  END DO; END DO
END IF

END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/



!==================================================================================================================================
!> Finalize boundary conditions
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,BCSideID,nBCByType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC

END MODULE MOD_GetBoundaryFlux
