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
!> Contains routines to predict the solution at the new stage
!===================================================================================================================================
MODULE MOD_Predictor
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE Predictor
  MODULE PROCEDURE Predictor
END INTERFACE

PUBLIC:: InitPredictor,FillInitPredictor,Predictor,StorePredictor,FinalizePredictor
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Allocate global variables required for predictor
!===================================================================================================================================
SUBROUTINE InitPredictor(TimeDiscMethod)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools,          ONLY: GETINT
USE MOD_Implicit_Vars,        ONLY: Upast,PredictorType,PredictorOrder,t_old,U_predictor,Ut_old,Un_old
USE MOD_Mesh_Vars,            ONLY: nElems
USE MOD_TimeDisc_Vars,        ONLY: nRKStages,TimeDiscType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: TimeDiscMethod !< name of time discretization to be used
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF(TimeDiscType.EQ.'ESDIRK')THEN
  PredictorType=GETINT('PredictorType','0')

  SELECT CASE(PredictorType)
  CASE(0) ! nothing to do
  CASE(1) ! use right hand side as predictor
    ALLOCATE(U_predictor(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE(2) ! Lagrange polynomial predictor
    PredictorOrder=GETINT('PredictorOrder','1')
    ALLOCATE(t_old(0:PredictorOrder,2:nRKStages))
    ALLOCATE(Upast(      1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,0:PredictorOrder,2:nRKStages))
    ALLOCATE(U_predictor(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE(3) ! predictor using dense output polynomial for extrapolation
    SELECT CASE (TRIM(TimeDiscMethod))
    CASE('esdirk3-4')
      PredictorOrder = 2
    CASE('esdirk4-6')
      PredictorOrder = 3
    CASE DEFAULT
      CALL abort(__STAMP__,'No dense output extrapolation available for chosen time discretization!',PredictorType,999.)
    END SELECT
    ALLOCATE(Ut_old(     1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,1:nRKStages))
    ALLOCATE(Un_old(     1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
    ALLOCATE(U_predictor(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE DEFAULT
    CALL abort(__STAMP__,'PredictorType not implemented!',PredictorType,999.)
  END SELECT
END IF
END SUBROUTINE InitPredictor

!===================================================================================================================================
!> Initially fills stage values so that the predictor works also at first timestep
!===================================================================================================================================
SUBROUTINE FillInitPredictor(t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars,  ONLY: t_old,Upast,PredictorType,PredictorOrder,PredictorCounter,Un_old,Ut_old
USE MOD_DG_Vars,        ONLY: U
USE MOD_TimeDisc_Vars,  ONLY: nRKStages,RKb_denseout
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: i,iStage
!===================================================================================================================================
SELECT CASE(PredictorType)
CASE(2)
  DO i=0,PredictorOrder
    DO iStage=2,nRKStages
      t_old(i,iStage) = t
      Upast(:,:,:,:,:,i,iStage) = U
    END DO
  END DO
  PredictorCounter = -1
CASE(3)
  IF(ALLOCATED(RKb_denseout))THEN
    Ut_old = 0.
    Un_old = U
  ELSE
    CALL abort(__STAMP__, &
        'Predictor Type 3 (extrapolation via dense output formula) not available for chosen Runge-Kutta method or predictor order')
  END IF
END SELECT

END SUBROUTINE FillInitPredictor

!===================================================================================================================================
!> predicts the new Stage-value to decrease computational time
!===================================================================================================================================
SUBROUTINE Predictor(tStage,iStage)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,          ONLY: U
USE MOD_Mesh_Vars,        ONLY: nElems
USE MOD_Implicit_Vars,    ONLY: Upast,t_old,PredictorType,PredictorOrder,PredictorCounter,U_predictor,Un_old,Ut_old,LinSolverRHS
USE MOD_TimeDisc_Vars,    ONLY: dt,dt_old,RKb_denseout,RKc_implicit,nRKStages
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: tStage
INTEGER,INTENT(IN)           :: iStage
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: i,j,k,ii,jj,iElem
REAL                         :: Lagrange(0:PredictorOrder)
!===================================================================================================================================
SELECT CASE(PredictorType)
CASE(0)
  ! trival guess
  ! nothing to do, because old stage value is used as a prediction
CASE(1) ! use RHS as Predictor
  U_predictor=LinSolverRHS
CASE(2)
  IF(PredictorCounter.GT.PredictorOrder)THEN ! do not build predictor if not a sufficient amount of states are present
    ! extrapolate chimera boundary conditions from old values (Lagrange polynomial)
    ! build Lagrange polynomial with degree PredictorOrder
    Lagrange(:)=1.
    DO j=0,PredictorOrder
      DO i=0,PredictorOrder
        IF(i.NE.j)THEN 
          Lagrange(j) = Lagrange(j)*((tStage-t_old(i,iStage))/(t_old(j,iStage)-t_old(i,iStage)))
        END IF
      END DO ! i
    END DO ! j
    ! build predictor
    DO iElem=1,nElems
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        U_predictor(:,i,j,k,iElem) = 0.
        DO jj=0,PredictorOrder
          U_predictor(:,i,j,k,iElem) = U_predictor(:,i,j,k,iElem) + Upast(:,i,j,k,iElem,jj,iStage)*Lagrange(jj)
        END DO
      END DO; END DO; END DO! i,j,k=0,PP_N
    END DO ! iElem
  ELSE
    U_predictor = U
  END IF ! PredictorCounter
CASE(3) ! predictor using dense output polynomial for extrapolation
  ! Kennedey, Carpenter: Additive Runge-Kutta Schemes for Convection-Diffusion-Reaction Equtions
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      U_predictor(:,i,j,k,iElem) = Un_old(:,i,j,k,iElem)
      DO ii=1,nRKStages
        DO jj=1,PredictorOrder
          U_predictor(:,i,j,k,iElem) = U_predictor(:,i,j,k,iElem) + dt_old * &
                                       Ut_old(:,i,j,k,iElem,ii)*RKb_denseout(jj,ii)*(1.+RKc_implicit(iStage)*dt/dt_old)**jj
        END DO
      END DO
    END DO; END DO; END DO! i,j,k=0,PP_N
  END DO ! iElem
END SELECT

END SUBROUTINE Predictor

!===================================================================================================================================
!> Stores information of old stages required for extrapolation
!===================================================================================================================================
SUBROUTINE StorePredictor(Ut_implicit,Un,tStage,iStage)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,        ONLY: nElems
USE MOD_DG_Vars,          ONLY: U
USE MOD_Implicit_Vars,    ONLY: PredictorType,t_old,Upast,PredictorOrder,Ut_old,Un_old,PredictorCounter
USE MOD_TimeDisc_Vars,    ONLY: nRKStages
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Ut_implicit(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,1:nRKStages) 
REAL,INTENT(IN)              :: Un(         1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) 
REAL,INTENT(IN)              :: tStage
INTEGER,INTENT(IN)           :: iStage
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: j
!===================================================================================================================================
! 1 - old stage, 0 - actual stage
SELECT CASE(PredictorType)
CASE(2)
  IF(nRKStages.EQ.iStage) PredictorCounter = PredictorCounter + 1
  ! shift old values one time level
  DO j=PredictorOrder,1,-1
    t_old(j,iStage) = t_old(j-1,iStage)
    Upast(:,:,:,:,:,j,iStage) = Upast(:,:,:,:,:,j-1,iStage)
  END DO
  !fill new predictor time level
  t_old(0,iStage) = tStage
  Upast(:,:,:,:,:,0,iStage)=U
CASE(3)
  IF(iStage.EQ.nRKStages)THEN
    Ut_old = Ut_implicit
    Un_old = Un
  END IF
END SELECT

END SUBROUTINE StorePredictor

!===================================================================================================================================
!> Deallocate global variables required for predictor
!===================================================================================================================================
SUBROUTINE FinalizePredictor()
! MODULES
USE MOD_Implicit_vars,    ONLY:Upast,t_old,U_predictor,Ut_old,Un_old
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(Upast)
SDEALLOCATE(t_old)
SDEALLOCATE(U_predictor)
SDEALLOCATE(Ut_old)
SDEALLOCATE(Un_old)
END SUBROUTINE FinalizePredictor

END MODULE MOD_Predictor
