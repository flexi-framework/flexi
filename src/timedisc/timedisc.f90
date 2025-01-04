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
MODULE MOD_TimeDisc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: TimeDisc
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> GTS Temporal discretization
!==================================================================================================================================
SUBROUTINE TimeDisc()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze             ,ONLY: Analyze
USE MOD_Analyze_Vars        ,ONLY: analyze_dt,tWriteData,WriteData_dt
USE MOD_AnalyzeEquation_Vars,ONLY: doCalcTimeAverage
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
USE MOD_BaseFlow            ,ONLY: UpdateBaseFlow
USE MOD_BaseFlow_Vars       ,ONLY: doBaseFlow
USE MOD_DG                  ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars             ,ONLY: U
USE MOD_Equation_Vars       ,ONLY: StrVarNames
USE MOD_HDF5_Output         ,ONLY: WriteState
USE MOD_IO_HDF5             ,ONLY:
USE MOD_Mesh_Vars           ,ONLY: MeshFile,nGlobalElems
USE MOD_Output              ,ONLY: Visualize,PrintStatusLine
USE MOD_Overintegration     ,ONLY: Overintegration
USE MOD_Overintegration_Vars,ONLY: OverintegrationType
USE MOD_Predictor           ,ONLY: FillInitPredictor
USE MOD_RecordPoints        ,ONLY: RecordPoints
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc
USE MOD_Restart_Vars        ,ONLY: DoRestart,RestartTime,FlushInitialState
USE MOD_TestCase            ,ONLY: AnalyzeTestCase,CalcForcing
USE MOD_TimeDisc_Functions  ,ONLY: InitTimeStep,UpdateTimeStep,AnalyzeTimeStep
USE MOD_TimeStep            ,ONLY: TimeStep
USE MOD_TestCase_Vars       ,ONLY: doTCSource
USE MOD_TimeDisc_Vars       ,ONLY: iter,iter_analyze,maxIter
USE MOD_TimeDisc_Vars       ,ONLY: t,tStart,tEnd,dt,tAnalyze
USE MOD_TimeDisc_Vars       ,ONLY: TimeDiscType
USE MOD_TimeDisc_Vars       ,ONLY: doAnalyze,doFinalize,writeCounter,nCalcTimestep
USE MOD_TimeDisc_Vars       ,ONLY: time_start
USE MOD_TimeAverage         ,ONLY: CalcTimeAverage
#if FV_ENABLED
USE MOD_Indicator           ,ONLY: CalcIndicator
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
USE MOD_FV_Switching        ,ONLY: FV_FillIni,FV_Switch,FV_Info
#elif FV_ENABLED == 2 || FV_ENABLED == 3
USE MOD_FV_Blending         ,ONLY: FV_Info
#endif /*FV_ENABLED == 1*/
#if PP_LIMITER
USE MOD_PPLimiter           ,ONLY: PPLimiter,PPLimiter_Info
USE MOD_Filter_Vars         ,ONLY: DoPPLimiter
#endif /*PP_LIMITER*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: TimeArray(8)              !< Array for system time
!==================================================================================================================================

SWRITE(UNIT_stdOut,'(132("-"))')

! Get system time
CALL DATE_AND_TIME(values=TimeArray)
SWRITE(UNIT_stdOut,'(A,I2.2,A1,I2.2,A1,I4.4,A1,I2.2,A1,I2.2,A1,I2.2)') &
  ' Sys date  :    ',timeArray(3),'.',timeArray(2),'.',timeArray(1),' ',timeArray(5),':',timeArray(6),':',timeArray(7)

! write number of grid cells and dofs only once per computation
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#GridCells : ',REAL(nGlobalElems)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs      : ',REAL(nGlobalElems)*REAL((PP_N+1)**PP_dim)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#Procs     : ',REAL(nProcessors)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs/Proc : ',REAL(nGlobalElems*(PP_N+1)**PP_dim/nProcessors)

! Fill correct initial time
t = MERGE(RestartTime,tStart,DoRestart)

! NOTE: Set initial variables
! t is set in init solution
tWriteData        = MIN(t+WriteData_dt,tEnd)
tAnalyze          = MIN(t+analyze_dt  ,tEnd)
iter              = 0
iter_analyze      = 0
writeCounter      = 0
nCalcTimestep     = 0
doAnalyze         = .FALSE.
doFinalize        = .FALSE.

! do initial filtering of solution (has to be done after DG and FV are filled -> after restart)
! --- Perform some preparational steps ---  overintegrate solution first time
SELECT CASE(OverintegrationType)
  CASE (1)
    CALL Overintegration(U)
  CASE (2)
    CALL ApplyJacobianCons(U,toPhysical=.FALSE.,FVE=0)
    CALL Overintegration(U)
END SELECT

! Indicator (and FV_FillIni) need to be called for the DG operator below to return a valid solution
#if FV_ENABLED
CALL CalcIndicator(U,t)
#endif
#if FV_ENABLED == 1
IF(.NOT.DoRestart)  CALL FV_FillIni()
#endif
! initial update of baseflow
IF (doBaseFlow) CALL UpdateBaseFlow(0.)

! Do first RK stage of first timestep to fill gradients
CALL DGTimeDerivative_weakForm(t)

! Compute indicator again (and switch) since first indicator call above did not include gradients (needed by e.g. Ducros indicator)
#if FV_ENABLED
CALL CalcIndicator(U,t)
#endif
#if FV_ENABLED == 1
CALL FV_Switch(U,AllowToDG=.FALSE.)
#endif /* FV_ENABLED == 1 */
#if PP_LIMITER
IF(DoPPLimiter) CALL PPLimiter()
#endif /*PP_LIMITER*/

IF(.NOT.DoRestart) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WRITING INITIAL SOLUTION:'
ELSEIF (FlushInitialState) THEN
  SWRITE(UNIT_stdOut,'(A)') ' REWRITING SOLUTION:'
END IF

! Call Testcase analyze routines for the initial solution
CALL AnalyzeTestCase(t,.FALSE.)

! Write the state at time=0, i.e. the initial condition
CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,FutureTime=tWriteData,isErrorFile=.FALSE.)
CALL Visualize(t,U)

! compute initial timestep
CALL InitTimeStep()

! Run initial analyze
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' Errors of initial solution:'
! print initial analyze
CALL Analyze(t,iter)

! compute initial timestep
CALL InitTimeStep()

! Fill recordpoints buffer (initialization/restart)
IF(RP_onProc) CALL RecordPoints(PP_nVar,StrVarNames,iter,t,.TRUE.)

#if FV_ENABLED
CALL FV_Info(INT(1, KIND=DP))
#endif /*FV_ENABLED*/
#if PP_LIMITER
CALL PPLimiter_Info(INT(1, KIND=DP))
#endif /*PP_LIMITER*/

SWRITE(UNIT_stdOut,'(A)') ' CALCULATION RUNNING...'

IF(TimeDiscType.EQ.'ESDIRK') CALL FillInitPredictor(t)
CALL CPU_TIME(time_start)

! Run computation
DO
  IF (doBaseFlow)       CALL UpdateBaseFlow(dt)
  ! Update time step
  CALL UpdateTimeStep()

  IF(doCalcTimeAverage) CALL CalcTimeAverage(.FALSE.,dt,t)
  IF(doTCSource)        CALL CalcForcing(t,dt)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Perform Timestep using a global time stepping routine, attention: only RK3 has time dependent BC
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL TimeStep(t)
  iter         = iter         + 1
  iter_analyze = iter_analyze + 1
  t            = t            + dt

  ! Perform analysis at the end of the RK loop
  CALL AnalyzeTimeStep()

  CALL PrintStatusLine(t,dt,tStart,tEnd,iter,maxIter)

  IF(doFinalize) EXIT
END DO

END SUBROUTINE TimeDisc

END MODULE MOD_TimeDisc
