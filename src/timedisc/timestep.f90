!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
!> Module for the GTS Temporal discretization
!==================================================================================================================================
MODULE MOD_TimeStep
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE SetTimeStep
  MODULE PROCEDURE SetTimeStep
END INTERFACE

INTERFACE TimeStepByLSERKW2
  MODULE PROCEDURE TimeStepByLSERKW2
END INTERFACE

INTERFACE TimeStepByLSERKK3
  MODULE PROCEDURE TimeStepByLSERKK3
END INTERFACE

INTERFACE TimeStepByESDIRK
  MODULE PROCEDURE TimeStepByESDIRK
END INTERFACE

! > Dummy interface for time step function pointer
ABSTRACT INTERFACE
  SUBROUTINE TimeIntegrator(t)
    REAL,INTENT(INOUT) :: t
  END SUBROUTINE
END INTERFACE

PROCEDURE(TimeIntegrator),POINTER :: TimeStep !< pointer to timestepping routine, depends on td

PUBLIC :: TimeStep
PUBLIC :: SetTimeStep
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Set the timestep pointer to the specified type of timestep routine
!===================================================================================================================================
SUBROUTINE SetTimeStep(TimeDiscType)
! MODULES
USE MOD_Globals   ,ONLY: Abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: TimeDiscType   !< type of timestep required, i.e. number of registers
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (ASSOCIATED(TimeStep)) CALL Abort(__STAMP__, &
   'TimeStep pointer is already associated! Dereference the pointer before associating it with a new time step routine!')

SELECT CASE(TimeDiscType)
  CASE('LSERKW2')
    TimeStep=>TimeStepByLSERKW2
  CASE('LSERKK3')
    TimeStep=>TimeStepByLSERKK3
  CASE('ESDIRK')
    TimeStep=>TimeStepByESDIRK
  CASE DEFAULT
    CALL Abort(__STAMP__, 'Unknown timestep routine!')
END SELECT

END SUBROUTINE SetTimeStep


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE TimeStepByLSERKW2(t)
! MODULES
USE MOD_PreProc
USE MOD_Vector
USE MOD_DG            ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars       ,ONLY: U,Ut,nTotalU
USE MOD_PruettDamping ,ONLY: TempFilterTimeDeriv
USE MOD_TimeDisc_Vars ,ONLY: dt,Ut_tmp,RKA,RKb,RKc,nRKStages,CurrentStage
#if FV_ENABLED
USE MOD_Indicator     ,ONLY: CalcIndicator
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
USE MOD_FV_Switching  ,ONLY: FV_Switch
USE MOD_FV_Vars       ,ONLY: FV_toDGinRK
#endif /*FV_ENABLED==1*/
#if PP_LIMITER
USE MOD_PPLimiter     ,ONLY: PPLimiter
USE MOD_Filter_Vars   ,ONLY: DoPPLimiter
#endif /*PP_LIMITER*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: t                                     !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: b_dt(1:nRKStages)
REAL     :: tStage
INTEGER  :: iStage
!===================================================================================================================================

! Premultiply with dt
b_dt = RKb*dt

DO iStage = 1,nRKStages
  ! NOTE: perform timestep in rk
  CurrentStage = iStage
  IF (CurrentStage.EQ.1) THEN; tStage = t
  ELSE                       ; tStage = t+RKc(CurrentStage)*dt
  END IF

  ! Update gradients
  CALL DGTimeDerivative_weakForm(tStage)

  IF (iStage.EQ.1) THEN
    CALL VCopy(nTotalU,Ut_tmp,Ut)                        !Ut_tmp = Ut
  ELSE
    CALL VAXPBY(nTotalU,Ut_tmp,Ut,ConstOut=-RKA(iStage)) !Ut_tmp = Ut - Ut_tmp*RKA (iStage)
  END IF
  CALL VAXPBY(nTotalU,U,Ut_tmp,   ConstIn =b_dt(iStage)) !U       = U + Ut_tmp*b_dt(iStage)

#if FV_ENABLED
  ! Time needs to be evaluated at the next step because time integration was already performed
  ASSOCIATE(tFV => MERGE(t+dt,t,iStage.EQ.nRKStages))
  CALL CalcIndicator(U,tFV)
  END ASSOCIATE
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
  ! NOTE: Apply switch and update FV_Elems
  CALL FV_Switch(U,Ut_tmp,AllowToDG=(iStage.EQ.nRKStages .OR. FV_toDGinRK))
#endif /*FV_ENABLED==1*/
#if PP_LIMITER
  IF(DoPPLimiter) CALL PPLimiter()
#endif /*PP_LIMITER*/
END DO

END SUBROUTINE TimeStepByLSERKW2


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration:  3 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE TimeStepByLSERKK3(t)
! MODULES
USE MOD_PreProc
USE MOD_Vector
USE MOD_DG           ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars      ,ONLY: U,Ut,nTotalU
USE MOD_TimeDisc_Vars,ONLY: dt,UPrev,S2,RKdelta,RKg1,RKg2,RKg3,RKb,RKc,nRKStages,CurrentStage
#if FV_ENABLED
USE MOD_Indicator     ,ONLY: CalcIndicator
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
USE MOD_FV_Switching  ,ONLY: FV_Switch
USE MOD_FV_Vars       ,ONLY: FV_toDGinRK
#endif /*FV_ENABLED==1*/
#if PP_LIMITER
USE MOD_PPLimiter    ,ONLY: PPLimiter
USE MOD_Filter_Vars  ,ONLY: DoPPLimiter
#endif /*PP_LIMITER*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: t                                     !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: b_dt(1:nRKStages)
REAL     :: tStage
INTEGER  :: iStage
!===================================================================================================================================

! Nomenclature:
! S1 == U, S2 == S2, S3 == UPrev

! Premultiply with dt
b_dt = RKb*dt

DO iStage = 1,nRKStages
  ! NOTE: perform timestep in rk
  CurrentStage = iStage
  IF (CurrentStage.EQ.1) THEN; tStage = t
  ELSE                       ; tStage = t+RKc(CurrentStage)*dt
  END IF

  ! Update gradients
  CALL DGTimeDerivative_weakForm(tStage)

  IF (iStage.EQ.1) THEN
    CALL VCopy(nTotalU,UPrev,U)                                          !UPrev = U
    CALL VCopy(nTotalU,S2,U)                                             !S2    = U
  ELSE
    CALL VAXPBY(nTotalU,S2,U,ConstIn=RKdelta(iStage))                    !S2    = S2 + U*RKdelta(iStage)
    CALL VAXPBY(nTotalU,U,S2,ConstOut=RKg1(iStage),ConstIn=RKg2(iStage)) !U     = RKg1(iStage)*U + RKg2(iStage)*S2
    CALL VAXPBY(nTotalU,U,UPrev,ConstIn=RKg3(iStage))                    !U     = U + RKg3(ek)*UPrev
  END IF
  CALL VAXPBY(nTotalU,U,Ut,ConstIn=b_dt(iStage))                         !U     = U + Ut*b_dt(iStage)

#if FV_ENABLED
  ! Time needs to be evaluated at the next step because time integration was already performed
  ASSOCIATE(tFV => MERGE(t+dt,t,iStage.EQ.nRKStages))
  CALL CalcIndicator(U,tFV)
  ! NOTE: Apply switch and update FV_Elems
  END ASSOCIATE
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
  CALL FV_Switch(U,Uprev,S2,AllowToDG=(iStage.EQ.nRKStages .OR. FV_toDGinRK))
#endif /*FV_ENABLED==1*/
#if PP_LIMITER
  IF(DoPPLimiter) CALL PPLimiter()
#endif /*PP_LIMITER*/
END DO

END SUBROUTINE TimeStepByLSERKK3

!===================================================================================================================================
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> ESDIRK time integrator with RKA/b/c from Butcher tableau:
!>
!> Un=U
!> Calculation of the explicit terms for every stage i=1,...,s of the ESDIRK
!> RHS_(i-1) = Un + dt * sum(j=1,i-1) a_ij  Ut(t^n + c_j delta t^n, U_j)
!> Call Newton for searching the roots of the function
!> F(U) = U - RHS(i-1) - a_ii * delta t * Ut(t^n + c_i * delta t^n, U) = 0
!===================================================================================================================================
SUBROUTINE TimeStepByESDIRK(t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG                ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars           ,ONLY: U,Ut
USE MOD_Implicit          ,ONLY: Newton
USE MOD_Implicit_Vars     ,ONLY: LinSolverRHS,adaptepsNewton,epsNewton,nDOFVarProc,nGMRESIterdt,NewtonConverged,nInnerGMRES
USE MOD_Mathtools         ,ONLY: GlobalVectorDotProduct
USE MOD_Mesh_Vars         ,ONLY: nElems
#if USE_PRECOND
USE MOD_Precond           ,ONLY: BuildPrecond
USE MOD_Precond_Vars      ,ONLY: PrecondIter
#endif /*USE_PRECOND*/
USE MOD_Predictor         ,ONLY: Predictor,PredictorStoreValues
USE MOD_TimeDisc_Vars     ,ONLY: dt,nRKStages,RKA_implicit,RKc_implicit,iter,CFLScale,CFLScale_Readin
USE MOD_TimeDisc_Vars     ,ONLY: RKb_implicit,RKb_embedded,safety,ESDIRK_gamma
#if PARABOLIC
USE MOD_TimeDisc_Vars     ,ONLY: DFLScale,DFLScale_Readin
#endif
#if FV_ENABLED
USE MOD_Indicator         ,ONLY: CalcIndicator
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
USE MOD_FV_Switching      ,ONLY: FV_Switch
#endif /*FV_ENABLED==1*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: t   !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: Ut_implicit(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,1:nRKStages) ! temporal variable for Ut_implicit
REAL    :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
INTEGER :: iStage,iCounter
REAL    :: tStage
REAL    :: delta_embedded(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)          ! difference between solution obtained with
                                                                             ! full order scheme and embedded scheme
!===================================================================================================================================
!CALL DGTimeDerivative_weakForm(t)! has to be called before preconditioner to fill U_master/slave ! already called in timedisc
#if USE_PRECOND
IF ((iter==0).OR.(MOD(iter,PrecondIter)==0)) CALL BuildPrecond(t,ESDIRK_gamma,dt)
#endif /*USE_PRECOND*/
tStage                   = t
Un                       = U
Ut_implicit(:,:,:,:,:,1) = Ut
DO iStage=2,nRKStages
  IF (NewtonConverged) THEN
    ! Time of current stage
    tStage = tStage + RKc_implicit(iStage)*dt
    ! Compute RHS for linear solver
    LinSolverRHS=Un
    DO iCounter=1,iStage-1
      LinSolverRHS = LinSolverRHS + dt*(RKA_implicit(iStage,iCounter)*Ut_implicit(:,:,:,:,:,iCounter))
    END DO
    ! Get predictor of u^s+1
    CALL Predictor(tStage,iStage)
    ! Solve to new stage
    CALL Newton(tStage,RKA_implicit(iStage,iStage))
    ! Store old values for use in next stages
    !CALL DGTimeDerivative_weakForm(tStage) ! already set in last Newton iteration
    Ut_implicit(:,:,:,:,:,iStage)=Ut
    ! Store predictor
    CALL PredictorStoreValues(Ut_implicit,Un,tStage,iStage)
  END IF
END DO

! Adaptive Newton tolerance, see: Kennedy,Carpenter: Additive Runge-Kutta Schemes for Convection-Diffusion-Reaction Equations
IF (adaptepsNewton.AND.NewtonConverged) THEN
  delta_embedded = 0.
  DO iStage=1,nRKStages
    delta_embedded = delta_embedded + (RKb_implicit(iStage)-RKb_embedded(iStage)) * Ut_implicit(:,:,:,:,:,iStage)
  END DO
  CALL GlobalVectorDotProduct(delta_embedded,delta_embedded,nDOFVarProc,epsNewton)
  epsNewton = (MIN(dt*SQRT(epsNewton)/safety,1E-3))
#if DEBUG
  SWRITE(UNIT_stdOut,'(A,ES16.7)') 'epsNewton = ',epsNewton
#endif
END IF

IF (NewtonConverged) THEN
  ! increase timestep size until target CFLScale is reached
  CFLScale = MIN(CFLScale_Readin,1.05*CFLScale)
#if PARABOLIC
  DFLScale = MIN(DFLScale_Readin,1.05*DFLScale)
#endif
#if FV_ENABLED
  ! Time needs to be evaluated at the next step
  CALL CalcIndicator(U,t+dt)
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
  ! NOTE: Apply switch and update FV_Elems
  CALL FV_Switch(U,AllowToDG=.TRUE.)
#endif /*FV_ENABLED==1*/
ELSE
  ! repeat current timestep with decreased timestep size
  U = Un
  t = t-dt
  CFLScale = 0.5*CFLScale
#if PARABOLIC
  DFLScale = 0.5*DFLScale
#endif
  NewtonConverged = .TRUE.
  IF (CFLScale(0).LT.0.01*CFLScale_Readin(0)) THEN
    CALL Abort(__STAMP__, &
    'Newton not converged with GMRES Iterations of last Newton step and CFL reduction',nInnerGMRES,CFLScale(0)/CFLScale_Readin(0))
  END IF
  SWRITE(UNIT_stdOut,'(A)') 'Attention: Timestep failed, repeating with dt/2!'
END IF

nGMRESIterdt = 0
END SUBROUTINE TimeStepByESDIRK

END MODULE MOD_TimeStep
