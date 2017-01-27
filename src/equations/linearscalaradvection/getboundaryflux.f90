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
! Private Part --------------------------------------------------------------------------------------------------------------------
INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

#if FV_ENABLED
INTERFACE GetBoundaryFVgradient
  MODULE PROCEDURE GetBoundaryFVgradient
END INTERFACE
#endif

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

#if PARABOLIC
INTERFACE Lifting_GetBoundaryFlux
  MODULE PROCEDURE Lifting_GetBoundaryFlux
END INTERFACE

! Public Part ---------------------------------------------------------------------------------------------------------------------

PUBLIC :: Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/
PUBLIC :: InitBC
PUBLIC :: GetBoundaryFlux
#if FV_ENABLED
PUBLIC :: GetBoundaryFVgradient
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
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,iError)
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
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(SideID,t,Nloc,Flux,UPrim_master,                   &
#if PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BC,nBCs,BoundaryType,nSides
USE MOD_Exactfunc    ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
USE MOD_Riemann      ,ONLY: GetFlux
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems_master
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                   :: SideID  
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: Nloc    !< polynomial degree
REAL,INTENT(IN)                      :: UPrim_master( PP_nVarPrim,0:Nloc,0:Nloc,1:nSides) !< inner surface solution
#if PARABOLIC
                                                                           !> inner surface solution gradients in x/y/z-direction
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,1:nSides),INTENT(IN)  :: gradUx_master,gradUy_master,gradUz_master
#endif /*PARABOLIC*/
                                                                           !> normal and tangential vectors on surfaces
REAL,DIMENSION(      3,0:Nloc,0:Nloc,0:FV_ENABLED,1:nSides),INTENT(IN)  :: NormVec,TangVec1,TangVec2
REAL,DIMENSION(      3,0:Nloc,0:Nloc,0:FV_ENABLED,1:nSides),INTENT(IN)  :: Face_xGP   !< positions of surface flux points
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc,1:nSides),INTENT(OUT) :: Flux       !< resulting boundary fluxes
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q
INTEGER                              :: BCType,BCState,nBCLoc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc):: U_Face_loc
INTEGER                              :: FVEM=0 ! FV_Elems_master
!==================================================================================================================================
DO iBC=1,nBCs
  !IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  !nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2) !Exact function or refstate
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    !DO iSide=1,nBCLoc
#if FV_ENABLED
      FVEM = FV_Elems_master(SideID)
#endif
      !IF(BCState.EQ.0)THEN
        DO q=0,Nloc; DO p=0,Nloc
          CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q,FVEM,SideID),U_Face_loc(:,p,q))
        END DO; END DO
      !ELSE
        !DO q=0,Nloc; DO p=0,Nloc
          !U_Face_loc(:,p,q) = RefStateCons(BCState,:)
        !END DO; END DO
      !END IF
      CALL GetFlux(Nloc,Flux(:,:,:,SideID),UPrim_master(1:PP_nVar,:,:,SideID),U_Face_loc, &
#if PARABOLIC
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
                 gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID),gradUz_master(:,:,:,SideID),&
#endif /*PARABOLIC*/
                 NormVec(:,:,:,FVEM,SideID),TangVec1(:,:,:,FVEM,SideID),TangVec2(:,:,:,FVEM,SideID),doBC=.TRUE.)
    !END DO

  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in linearscalaradvection/getboundaryflux.f90!')
  END SELECT ! BCType
END DO

END SUBROUTINE GetBoundaryFlux

#if FV_ENABLED
SUBROUTINE GetBoundaryFVgradient(t,gradU,U_inner)
!==================================================================================================================================
!==================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: Abort
USE MOD_PreProc
USE MOD_Exactfunc     ,ONLY: ExactFunc
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Equation_Vars ,ONLY: nBCByType,BCSideID
USE MOD_Mesh_Vars     ,ONLY: nSides,nBCs,BoundaryType
USE MOD_Mesh_Vars     ,ONLY: Face_xGP
USE MOD_FV_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: t
REAL,INTENT(IN)     :: U_inner(PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)    :: gradU(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: p,q
REAL    :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
REAL    :: UPrim_inner   (PP_nVarPrim,0:PP_N,0:PP_N)
REAL    :: UPrim_Face_loc(PP_nVarPrim,0:PP_N,0:PP_N)
INTEGER        :: iBC,iSide,SideID
INTEGER        :: BCType,BCState,nBCLoc
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
 
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2) ! exact BC = Dirichlet BC !!
    ! Determine the exact BC state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q,1,SideID),U_Face_loc(:,p,q))
      END DO ; END DO
      gradU(:,:,:,SideID) = U_inner(:,:,:,SideID) - U_Face_loc
    END DO
    
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in linearscalaradvection/getboundaryfvgradient.f90!',999,999.)
  END SELECT ! BCType
END DO
END SUBROUTINE GetBoundaryFVgradient
#endif




#if PARABOLIC
!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!>              SUBROUTINE CalcSurfInt
!> Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(t,U_master,Flux)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars    ,ONLY: Face_xGP,SurfElem
USE MOD_Exactfunc    ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
USE MOD_Lifting_Vars ,ONLY: doWeakLifting
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems_master
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                      :: t                                    !< current time (provided by time integration scheme)
REAL,INTENT(IN)                      :: U_master(PP_nVar,0:PP_N,0:PP_N,nBCSides)
REAL,INTENT(OUT)                     :: Flux(PP_nVar,0:PP_N,0:PP_N,nBCSides) !< lifting boundary flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N):: U_Face_loc
INTEGER                              :: FVEM=0 ! FV_Elems_master
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2)
    !IF(BCState.EQ.0)THEN
      ! BCState specifies refstate to be used, if 0 then use iniexactfunc
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
#if FV_ENABLED
        FVEM = FV_Elems_master(SideID)
#endif
        DO q=0,PP_N; DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q,FVEM,SideID),U_Face_loc(:,p,q))
        END DO; END DO
        Flux(:,:,:,SideID)=0.5*(U_master(:,:,:,SideID)+U_Face_loc)
      END DO
    !ELSE
      !DO iSide=1,nBCLoc
        !SideID=BCSideID(iBC,iSide)
        !DO q=0,PP_N; DO p=0,PP_N
          !Flux(:,p,q,SideID)=0.5*(U_master(:,p,q,SideID)+RefStateCons(BCState,:))
        !END DO; END DO
      !END DO
    !END IF
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in navierstokes/getboundaryflux.f90!')
  END SELECT
END DO ! iBC

IF(.NOT.doWeakLifting)THEN
  !in case lifting is done in strong form
  Flux(:,:,:,1:nBCSides)=Flux(:,:,:,1:nBCSides)-U_master(    :,:,:,1:nBCSides)
END IF

DO iSide=1,nBCSides
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,iSide)=Flux(:,p,q,iSide)*SurfElem(p,q,0,iSide)
  END DO; END DO
END DO ! iSide
END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/



!==================================================================================================================================
!> Initialize boundary conditions
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,nBCByType,BCSideID
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
