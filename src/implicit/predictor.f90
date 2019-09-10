#include "flexi.h"

MODULE MOD_Predictor
!===================================================================================================================================
! Contains routines to predict the solution at the new stage
!===================================================================================================================================
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

SUBROUTINE InitPredictor()
!===================================================================================================================================
! Allocate global variable 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools,          ONLY:GETINT
USE MOD_Implicit_Vars,        ONLY:U_stage_old,PredictorType,PredictorOrder,t_old
USE MOD_Mesh_Vars,            ONLY:nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
PredictorType=GETINT('PredictorType','0')

SELECT CASE(PredictorType)
CASE(0)
CASE(1) ! linear predictor
  ALLOCATE(t_old(0:1))
  ALLOCATE(U_stage_old(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:1))
CASE(2) ! Lagrange polynomial predictor
  PredictorOrder=GETINT('PredictorOrder','0')
  ALLOCATE(t_old(0:PredictorOrder))
  ALLOCATE(U_stage_old(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:PredictorOrder))
CASE DEFAULT
  CALL abort(&
__STAMP__ &
,'PredictorType not implemented!',PredictorType,999.)
END SELECT

END SUBROUTINE InitPredictor

SUBROUTINE FillInitPredictor(t)
!===================================================================================================================================
! predicts the new Stage-value to decrease computational time
!===================================================================================================================================
! MODULES
USE MOD_Implicit_Vars  ,ONLY: t_old,U_stage_old,PredictorType,PredictorOrder
USE MOD_DG_Vars        ,ONLY: U

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: i
!===================================================================================================================================
SELECT CASE(PredictorType)
CASE(1)
  t_old(:)=t
  U_stage_old(:,:,:,:,:,1) = U(:,:,:,:,:)
  U_stage_old(:,:,:,:,:,0) = U(:,:,:,:,:)
CASE(2)
  DO i=0,PredictorOrder
    t_old(i)=t
    U_stage_old(:,:,:,:,:,i) = U
  END DO
END SELECT


END SUBROUTINE FillInitPredictor

SUBROUTINE Predictor
!===================================================================================================================================
! predicts the new Stage-value to decrease computational time
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,          ONLY: U
USE MOD_Implicit_Vars    ,ONLY: U_stage_old,t_old,PredictorType,PredictorOrder,PredictorCounter,timeStage
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                         :: sdt_old
INTEGER                      :: i,j
REAL                         :: Lagrange(0:PredictorOrder)
!===================================================================================================================================
SELECT CASE(PredictorType)
CASE(0)
  ! trival guess
  ! nothing to do, because old stage value is used as a prediction
CASE(1)
  ! extrapolate linearly U at the old stages
  U_stage_old(:,:,:,:,:,0)=U
  IF(ABS(t_old(0)-t_old(1)).GT.1.E-15)THEN
    sdt_old = 1./(t_old(0)-t_old(1))
    U = U + (timeStage-t_old(0))*sdt_old * (U - U_stage_old(:,:,:,:,:,1))
  END IF
CASE(2)
  IF(PredictorCounter.GT.PredictorOrder)THEN ! do not build predictor if not a sufficient amount of states are present
    ! extrapolate U stages from old values (Lagrange polynomial)
    ! build Lagrange polynomial with degree PredictorOrder
    U_stage_old(:,:,:,:,:,0)=U
    Lagrange(:)=1.
    DO j=0,PredictorOrder
      DO i=0,PredictorOrder
        IF(i.NE.j)THEN 
          Lagrange(j) = Lagrange(j)*((timeStage-t_old(i))/(t_old(j)-t_old(i)))
        END IF
      END DO
    END DO
    ! build predictor
    U = 0.
    DO j=0,PredictorOrder
      U = U + U_stage_old(:,:,:,:,:,j)*Lagrange(j)
    END DO
  END IF
  PredictorCounter = PredictorCounter + 1
END SELECT

END SUBROUTINE Predictor


SUBROUTINE StorePredictor
!===================================================================================================================================
! predicts the new Stage-value to decrease computational time
!===================================================================================================================================
! MODULES
USE MOD_Implicit_Vars    ,ONLY: PredictorType,t_old,U_stage_old,PredictorOrder,timeStage

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: j
!===================================================================================================================================
! 1 - old stage, 0 - actual stage
SELECT CASE(PredictorType)
CASE(1)
  t_old(1) = t_old(0)
  t_old(0) = timeStage
  U_stage_old(:,:,:,:,:,1) = U_stage_old(:,:,:,:,:,0)
CASE(2)
  ! shift old values one time level
  DO j=PredictorOrder,1,-1
    t_old(j) = t_old(j-1)
    U_stage_old(:,:,:,:,:,j) = U_stage_old(:,:,:,:,:,j-1)
  END DO
  !fill new predictor time level
  t_old(0) = timeStage
END SELECT

END SUBROUTINE StorePredictor


SUBROUTINE FinalizePredictor()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Implicit_vars,    ONLY:U_stage_old,t_old
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(U_stage_old)
SDEALLOCATE(t_old)
END SUBROUTINE FinalizePredictor

END MODULE MOD_Predictor
