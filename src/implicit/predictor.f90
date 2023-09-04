!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
!> Contains routines to predict the solution at the new stage.
!> Four different approaches are implemented:
!>  * Simply use the current solution (no predictor)
!>  * Use the right hand side as a predictor
!>  * Use a polynomial extrapolation of varying order. The extrapolation is build from values of the same stage from old time steps
!>  * Use the dense output polynomial for extrapolation, only available for certain time integration schemes
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

INTERFACE InitPredictor
  MODULE PROCEDURE InitPredictor
END INTERFACE

INTERFACE FillInitPredictor
  MODULE PROCEDURE FillInitPredictor
END INTERFACE

INTERFACE Predictor
  MODULE PROCEDURE Predictor
END INTERFACE

INTERFACE PredictorStoreValues
  MODULE PROCEDURE PredictorStoreValues
END INTERFACE

INTERFACE FinalizePredictor
  MODULE PROCEDURE FinalizePredictor
END INTERFACE

PUBLIC:: InitPredictor,FillInitPredictor,Predictor,PredictorStoreValues,FinalizePredictor
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Allocate global variables and read in parameters required for predictor
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
  PredictorType=GETINT('PredictorType')

  SELECT CASE(PredictorType)
  CASE(0) ! nothing to do
  CASE(1) ! use right hand side as predictor
    ALLOCATE(U_predictor(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE(2) ! Lagrange polynomial predictor
    PredictorOrder=GETINT('PredictorOrder')
    ! For the extrapolation, we need values from old time steps. The extrapolaton will be done for each stage separately, so we need
    ! to store the solution and time instance at each stage. How many time steps we need to store depends on the polynomial order.
    ALLOCATE(t_old(0:PredictorOrder,2:nRKStages))
    ALLOCATE(Upast(      1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,0:PredictorOrder,2:nRKStages))
    ALLOCATE(U_predictor(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE(3) ! predictor using dense output polynomial for extrapolation
    ! The dense output can be used to calculate the solution at an arbitrary time instance inside of a time step. We use this to get
    ! an extrapolated value. Needs specific coefficients for each scheme, only some are available.
    SELECT CASE (TRIM(TimeDiscMethod))
    CASE('esdirk3-4')
      PredictorOrder = 2
    CASE('esdirk4-6')
      PredictorOrder = 3
    CASE DEFAULT
      CALL CollectiveStop(__STAMP__,'No dense output extrapolation available for chosen time discretization!',PredictorType,999.)
    END SELECT
    ! For the extrapolation we need the old solution Un and the old Uts from all stages of the last time step.
    ALLOCATE(Ut_old(     1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,1:nRKStages))
    ALLOCATE(Un_old(     1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
    ALLOCATE(U_predictor(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'PredictorType not implemented!',PredictorType,999.)
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
USE MOD_Implicit_Vars,  ONLY: t_old,Upast,PredictorType,PredictorOrder,Un_old,Ut_old
USE MOD_DG_Vars,        ONLY: U
USE MOD_TimeDisc_Vars,  ONLY: nRKStages
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
CASE(3)
  Ut_old = 0.
  Un_old = U
END SELECT

END SUBROUTINE FillInitPredictor

!===================================================================================================================================
!> Predicts the new stage-value to decrease computational time
!===================================================================================================================================
SUBROUTINE Predictor(tStage,iStage)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,          ONLY: U
USE MOD_Mesh_Vars,        ONLY: nElems
USE MOD_Implicit_Vars,    ONLY: Upast,t_old,PredictorType,PredictorOrder,U_predictor,Un_old,Ut_old,LinSolverRHS
USE MOD_TimeDisc_Vars,    ONLY: dt,dt_old,RKb_denseout,RKc_implicit,nRKStages,iter
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
  IF((iter-1).GT.PredictorOrder)THEN ! do not build predictor if not a sufficient amount of states are present
    ! build Lagrange polynomial with degree PredictorOrder
    ! the nodes are equal to the old stage times
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
    ! If not enough values are stored yet, simply use the current solution
    U_predictor = U
  END IF ! Enough stored values
CASE(3) ! predictor using dense output polynomial for extrapolation
  ! Kennedey, Carpenter: Additive Runge-Kutta Schemes for Convection-Diffusion-Reaction Equations
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      U_predictor(:,i,j,k,iElem) = Un_old(:,i,j,k,iElem)
      DO ii=1,nRKStages
        DO jj=1,PredictorOrder
          U_predictor(:,i,j,k,iElem) = U_predictor(:,i,j,k,iElem) + dt_old * &
                                       Ut_old(:,i,j,k,iElem,ii)*RKb_denseout(jj,ii)*(1.+RKc_implicit(iStage)*dt/dt_old)**jj
        END DO ! jj=1,PredictorOrder
      END DO ! ii=1,nRKStages
    END DO; END DO; END DO! i,j,k=0,PP_N
  END DO ! iElem
END SELECT

END SUBROUTINE Predictor

!===================================================================================================================================
!> Stores information of old stages required for extrapolation, get's called after a stage is finished.
!===================================================================================================================================
SUBROUTINE PredictorStoreValues(Ut,Un,tStage,iStage)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,        ONLY: nElems
USE MOD_DG_Vars,          ONLY: U
USE MOD_Implicit_Vars,    ONLY: PredictorType,t_old,Upast,PredictorOrder,Ut_old,Un_old
USE MOD_TimeDisc_Vars,    ONLY: nRKStages
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,1:nRKStages)  !< Time derivative, not all stages may be
                                                                                          !< filled!
REAL,INTENT(IN)              :: Un(         1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)     !< Current solution at begin of time step
REAL,INTENT(IN)              :: tStage                                                    !< Stage time
INTEGER,INTENT(IN)           :: iStage                                                    !< Current RK stage
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: j
!===================================================================================================================================
SELECT CASE(PredictorType)
CASE(2)
  ! shift old values one time level
  DO j=PredictorOrder,1,-1
    t_old(j,iStage) = t_old(j-1,iStage)
    Upast(:,:,:,:,:,j,iStage) = Upast(:,:,:,:,:,j-1,iStage)
  END DO
  ! fill new predictor time level
  t_old(0,iStage) = tStage
  Upast(:,:,:,:,:,0,iStage)=U
CASE(3)
  IF(iStage.EQ.nRKStages)THEN
    Ut_old = Ut
    Un_old = Un
  END IF
END SELECT

END SUBROUTINE PredictorStoreValues

!===================================================================================================================================
!> Deallocate global variables required for predictor
!===================================================================================================================================
SUBROUTINE FinalizePredictor()
! MODULES
USE MOD_Implicit_Vars
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
