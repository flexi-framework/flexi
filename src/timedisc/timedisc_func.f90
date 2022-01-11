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
!> Module for the GTS Temporal discretization
!==================================================================================================================================
MODULE MOD_TimeDisc_Functions
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersTimeDisc
  MODULE PROCEDURE DefineParametersTimeDisc
END INTERFACE

INTERFACE InitTimeDisc
  MODULE PROCEDURE InitTimeDisc
END INTERFACE

! INTERFACE InitTimeStep
!   MODULE PROCEDURE InitTimeStep
! END INTERFACE

! INTERFACE UpdateTimeStep
!   MODULE PROCEDURE UpdateTimeStep
! END INTERFACE

! INTERFACE AnalyzeTimeStep
!   MODULE PROCEDURE AnalyzeTimeStep
! END INTERFACE

INTERFACE FinalizeTimeDisc
  MODULE PROCEDURE FinalizeTimeDisc
END INTERFACE

PUBLIC :: DefineParametersTimeDisc
PUBLIC :: InitTimeDisc
PUBLIC :: InitTimeStep
PUBLIC :: UpdateTimeStep
PUBLIC :: AnalyzeTimeStep
PUBLIC :: FinalizeTimeDisc
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersTimeDisc()
! MODULES
USE MOD_ReadInTools         ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("TimeDisc")
CALL prms%CreateStringOption('TimeDiscMethod', "Specifies the type of time-discretization to be used, e.g. the name of&
                                               & a specific Runge-Kutta scheme. Possible values:\n"//&
                                               "  * standardrk3-3\n  * carpenterrk4-5\n  * niegemannrk4-14\n"//&
                                               "  * toulorgerk4-8c\n  * toulorgerk3-7c\n  * toulorgerk4-8f\n"//&
                                               "  * ketchesonrk4-20\n  * ketchesonrk4-18\n  * eulerimplicit\n"//&
                                               "  * cranknicolson2-2\n  * esdirk2-3\n  * esdirk3-4\n"//&
                                               "  * esdirk4-6" , value='CarpenterRK4-5')
CALL prms%CreateRealOption(  'TEnd',           "End time of the simulation (mandatory).")
CALL prms%CreateRealOption(  'CFLScale',       "Scaling factor for the theoretical CFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateRealOption(  'DFLScale',       "Scaling factor for the theoretical DFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateIntOption(   'maxIter',        "Stop simulation when specified number of timesteps has been performed.", value='-1')
CALL prms%CreateIntOption(   'NCalcTimeStepMax',"Compute dt at least after every Nth timestep.", value='1')
END SUBROUTINE DefineParametersTimeDisc


!==================================================================================================================================
!> Get information for end time and max time steps from ini file
!==================================================================================================================================
SUBROUTINE InitTimeDisc()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars         ,ONLY: NFilter,FilterType
USE MOD_IO_HDF5             ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Overintegration_Vars,ONLY: NUnder
USE MOD_Predictor           ,ONLY: InitPredictor
USE MOD_ReadInTools         ,ONLY: GETREAL,GETINT,GETSTR
USE MOD_StringTools         ,ONLY: LowCase,StripSpaces
USE MOD_TimeDisc_Vars       ,ONLY: b_dt,CFLScale,DFLScale,dtElem,dt,tend
USE MOD_TimeDisc_Vars       ,ONLY: Ut_tmp,UPrev,S2
USE MOD_TimeDisc_Vars       ,ONLY: maxIter,nCalcTimeStepMax
USE MOD_TimeDisc_Vars       ,ONLY: SetTimeDiscCoefs,TimeStep,TimeDiscName,TimeDiscType,TimeDiscInitIsDone,nRKStages
USE MOD_TimeStep            ,ONLY: TimeStepByLSERKW2,TimeStepByLSERKK3,TimeStepByESDIRK
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255):: TimeDiscMethod
INTEGER           :: NEff
!==================================================================================================================================
IF(TimeDiscInitIsDone)THEN
   SWRITE(UNIT_stdOut,'(A)') "InitTimeDisc already called."
   RETURN
END IF

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TIMEDISC...'

TimeDiscMethod   = GETSTR('TimeDiscMethod','Carpenter RK4-5')
maxIter          = GETINT('maxIter','-1')
nCalcTimeStepMax = GETINT('nCalcTimeStepMax','1')

CALL StripSpaces(TimeDiscMethod)
CALL LowCase(TimeDiscMethod)
SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: '//TRIM(TimeDiscName)

CALL SetTimeDiscCoefs(TimeDiscMethod)
SELECT CASE(TimeDiscType)
  CASE('LSERKW2')
    TimeStep=>TimeStepByLSERKW2
    ALLOCATE(Ut_tmp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE('LSERKK3')
    TimeStep=>TimeStepByLSERKK3
    ALLOCATE(S2   (1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) &
            ,UPrev(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE('ESDIRK')
    ! Implicit time integration
    TimeStep=>TimeStepByESDIRK
    ! Predictor for Newton
    CALL InitPredictor(TimeDiscMethod)
END SELECT

ALLOCATE(b_dt(1:nRKStages))

! Read the end time TEnd from ini file
TEnd     = GETREAL('TEnd')
! Read the normalized CFL number
CFLScale = GETREAL('CFLScale')
#if PARABOLIC
! Read the normalized DFL number
DFLScale = GETREAL('DFLScale')
#endif /*PARABOLIC*/
NEff     = MIN(PP_N,NFilter,NUnder)
IF(FilterType.GT.2) NEff = PP_N  !LAF,HESTHAVEN no timestep effect
CALL FillCFL_DFL(NEff,PP_N)
! Set timestep to a large number
ALLOCATE(dtElem(nElems))
dt       = HUGE(1.)
dtElem   = 0.

CALL AddToElemData(ElementOut,'dt',dtElem)

TimeDiscInitIsDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT TIMEDISC DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitTimeDisc


!===================================================================================================================================
!> Initial time step calculation for new timedisc-loop
!===================================================================================================================================
SUBROUTINE InitTimeStep()
! MODULES
USE MOD_Globals
USE MOD_CalcTimeStep        ,ONLY: CalcTimeStep
USE MOD_TimeDisc_Vars       ,ONLY: t,dt,dt_minOld,b_dt,RKb
USE MOD_TimeDisc_Vars       ,ONLY: ViscousTimeStep,CalcTimeStart,nCalcTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: errType
!===================================================================================================================================
dt            = CalcTimeStep(errType)
nCalcTimestep = 0
dt_minOld     = -999.
IF (errType.NE.0) CALL abort(__STAMP__,&
#if EQNSYSNR == 3
  'Error: (1) density, (2) convective / (3) viscous timestep / muTilde (4) is NaN. Type/time:',errType,t)
#else
  'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)
#endif

! Premultiply with dt
b_dt = RKb*dt

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,ES16.7)') ' Initial Timestep  : ', dt
IF(ViscousTimeStep.AND.MPIRoot) WRITE(UNIT_stdOut,'(A)') ' Viscous timestep dominates! '

CalcTimeStart = FLEXITIME()

END SUBROUTINE InitTimeStep


!===================================================================================================================================
!> Update time step at the beginning of each timedisc loop
!===================================================================================================================================
SUBROUTINE UpdateTimeStep()
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars        ,ONLY: tWriteData
USE MOD_CalcTimeStep        ,ONLY: CalcTimeStep
USE MOD_HDF5_Output         ,ONLY: WriteState
USE MOD_Mesh_Vars           ,ONLY: MeshFile
USE MOD_TimeDisc_Vars       ,ONLY: t,tAnalyze,tEnd,dt,dt_min,dt_minOld,b_dt,RKb
USE MOD_TimeDisc_Vars       ,ONLY: iter,maxIter,nCalcTimeStep,nCalcTimeStepMax
USE MOD_TimeDisc_Vars       ,ONLY: doAnalyze,doFinalize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: errType
!===================================================================================================================================

! Return if no timestep update requested in this iteration
IF (nCalcTimestep.GE.1) THEN
  nCalcTimestep = nCalcTimestep-1
  RETURN
END IF

dt_min(DT_MIN)     = CALCTIMESTEP(errType)
dt_min(DT_ANALYZE) = tAnalyze-t             ! Time to next analysis, put in extra variable so number does not change due to numerical errors
dt_min(DT_END)     = tEnd    -t             ! Do the same for end time
dt                 = MINVAL(dt_min)

IF (dt.EQ.dt_min(DT_ANALYZE)) doAnalyze  = .TRUE.
IF (dt.EQ.dt_min(DT_END    )) doFinalize = .TRUE.
dt                 = MINVAL(dt_min,MASK=dt_min.GT.0)

nCalcTimestep = MIN(FLOOR(ABS(LOG10(ABS(dt_minOld/dt-1.)**2.*100.+EPSILON(0.)))),nCalcTimeStepMax) - 1
dt_minOld     = dt
IF (errType.NE.0) THEN
  CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,FutureTime=tWriteData,isErrorFile=.TRUE.)
  CALL abort(__STAMP__,&
#if EQNSYSNR == 3
  'Error: (1) density, (2) convective / (3) viscous timestep / muTilde (4) is NaN. Type/time:',errType,t)
#else
  'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)
#endif
END IF

! Increase time step if the NEXT time step would be smaller than dt/100
IF(dt_min(DT_ANALYZE)-dt.LT.dt/100.0 .AND. dt_min(DT_ANALYZE).GT.0) THEN; dt = dt_min(DT_ANALYZE); doAnalyze  = .TRUE.; END IF
! Increase time step if the LAST time step would be smaller than dt/100
IF(    dt_min(DT_END)-dt.LT.dt/100.0 .AND. dt_min(DT_END    ).GT.0) THEN; dt = dt_min(DT_END)    ; doFinalize = .TRUE.; END IF

! Premultiply with dt
b_dt = RKb*dt

IF (iter.EQ.maxIter) THEN
  tEnd=t; tAnalyze=t; tWriteData=t
  doAnalyze=.TRUE.; doFinalize=.TRUE.
END IF

END SUBROUTINE UpdateTimeStep


!===================================================================================================================================
!> Update time step at the beginning of each timedisc loop
!===================================================================================================================================
SUBROUTINE AnalyzeTimeStep(writeCounter)
! MODULES
USE MOD_Globals
USE MOD_Analyze             ,ONLY: Analyze
USE MOD_Analyze_Vars        ,ONLY: analyze_dt,WriteData_dt,tWriteData,nWriteData
USE MOD_AnalyzeEquation_Vars,ONLY: doCalcTimeAverage
USE MOD_DG                  ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars             ,ONLY: U
USE MOD_Equation_Vars       ,ONLY: StrVarNames
USE MOD_HDF5_Output         ,ONLY: WriteState,WriteBaseFlow
USE MOD_Mesh_Vars           ,ONLY: MeshFile
USE MOD_Output              ,ONLY: Visualize,PrintAnalyze,PrintStatusLine
USE MOD_PruettDamping       ,ONLY: TempFilterTimeDeriv
USE MOD_RecordPoints        ,ONLY: RecordPoints,WriteRP
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc
USE MOD_Sponge_Vars         ,ONLY: CalcPruettDamping
USE MOD_TestCase            ,ONLY: AnalyzeTestCase
USE MOD_TestCase_Vars       ,ONLY: nAnalyzeTestCase
USE MOD_TimeAverage         ,ONLY: CalcTimeAverage
USE MOD_TimeDisc_Vars       ,ONLY: t,dt,tAnalyze,tEnd,CalcTimeStart
USE MOD_TimeDisc_Vars       ,ONLY: Ut_tmp,iter,iter_analyze,nCalcTimestep
USE MOD_TimeDisc_Vars       ,ONLY: doAnalyze,doFinalize
#if FV_ENABLED
USE MOD_FV                  ,ONLY: FV_Info,FV_Switch
USE MOD_Indicator           ,ONLY: CalcIndicator
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT) :: writeCounter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if FV_ENABLED
CALL CalcIndicator(U,t)
CALL FV_Switch(U,Ut_temp)
#endif
! Call DG operator to fill face data, fluxes, gradients for analyze
CALL DGTimeDerivative_weakForm(t)

! Call your analysis routine for your testcase here.
IF((MOD(iter,INT(nAnalyzeTestCase,KIND=8)).EQ.0).OR.doAnalyze) CALL AnalyzeTestCase(t)
! Evaluate recordpoints
IF(RP_onProc) CALL RecordPoints(PP_nVar,StrVarNames,iter,t,doAnalyze)
! Update Pruett filter base flow
IF(CalcPruettDamping) CALL TempFilterTimeDeriv(U,dt)

! Analyze and output now
IF(doAnalyze)THEN
  CALL PrintAnalyze(dt)
#if FV_ENABLED
  ! Summation has one more iter step
  CALL FV_Info(iter_analyze+1)
#endif
#if PP_LIMITER
    CALL PPLimiter_Info(iter_analyze+1)
#endif

  ! Visualize data and write solution
  writeCounter = writeCounter+1
  IF((writeCounter.EQ.nWriteData).OR.doFinalize)THEN
    ! Write various derived data
    IF(doCalcTimeAverage) CALL CalcTimeAverage(.TRUE.,dt,t)
    IF(RP_onProc)         CALL WriteRP(PP_nVar,StrVarNames,t,.TRUE.)
    IF(CalcPruettDamping) CALL WriteBaseFlow(TRIM(MeshFile),t)
    ! Write state file
    ! NOTE: this should be last in the series, so we know all previous data
    ! has been written correctly when the state file is present
    tWriteData = MIN(tAnalyze+WriteData_dt,tEnd)
    CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,FutureTime=tWriteData,isErrorFile=.FALSE.)
    ! Visualize data
    CALL Visualize(t,U)
    writeCounter = 0
  END IF

  ! do analysis
  CALL Analyze(t,iter)
  iter_analyze  = 0
  CalcTimeStart = FLEXITIME()
  tAnalyze      = MIN(tAnalyze+analyze_dt,tEnd)

  ! Disable analyze for next time step
  doAnalyze     = .FALSE.
END IF

END SUBROUTINE AnalyzeTimeStep



!===================================================================================================================================
!> Scaling of the CFL number, from paper GASSNER, KOPRIVA, "A comparision of the Gauss and Gauss-Lobatto
!> Discontinuous Galerkin Spectral Element Method for Wave Propagation Problems" .
!> For N=1-10 input CFLscale can now be (nearly) 1. and will be scaled adequately depending on
!> polynomial degree N, NodeType and TimeDisc method.
!===================================================================================================================================
SUBROUTINE FillCFL_DFL(Nin_CFL,Nin_DFL)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY:CFLScale,CFLScale_Readin,CFLScaleAlpha
#if PARABOLIC
USE MOD_TimeDisc_Vars       ,ONLY:DFLScale,DFLScale_Readin,DFLScaleAlpha,RelativeDFL
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_TimeDisc_Vars       ,ONLY:CFLScaleFV
#endif /*FV*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nin_CFL !< input polynomial degree for advection terms
INTEGER,INTENT(IN) :: Nin_DFL !< input polynomial degree for viscous terms
                              !< for overintegration via Filtering, the gradients and thus the visscous flux remains of order N+1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: alpha
#if !(PARABOLIC)
INTEGER            :: dummy
#endif
!===================================================================================================================================
! CFL in DG depends on the polynomial degree
! Runge-Kutta methods
alpha       = CFLScaleAlpha(MIN(15,Nin_CFL))
CFLScale(0) = CFLScale(0)*alpha
#if FV_ENABLED
CFLScale(1) = CFLScale(1)*CFLScaleFV/(PP_N+1.) ! equidistant distribution
#endif
IF((Nin_CFL.GT.15).OR.(CFLScale(0).GT.alpha))THEN
  SWRITE(UNIT_stdOut,'(132("!"))')
  SWRITE(UNIT_stdOut,'(A)')'Warning: The chosen CFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_stdOut,'(132("!"))')
END IF
!scale with 2N+1
CFLScale(0) = CFLScale(0)/(2.*Nin_CFL+1.)
SWRITE(UNIT_stdOut,'(A,2ES16.7)') '   CFL (DG/FV):',CFLScale
CFLScale_Readin = CFLScale

#if PARABOLIC
!########################### DFL ########################################
! DFL in DG depends on the polynomial degree
! since DFl is only on real axis, stability numbers are defined for RK3 and then scaled for RK4

alpha       = DFLScaleAlpha(MIN(10,Nin_DFL))*RelativeDFL
DFLScale(0) = DFLScale(0)*alpha
#if FV_ENABLED
DFLScale(1) = DFLScale(1)*DFLScaleAlpha(1)*RelativeDFL/(PP_N+1.)**2
#endif
IF((Nin_DFL.GT.10).OR.(DFLScale(0).GT.alpha))THEN
  SWRITE(UNIT_stdOut,'(132("!"))')
  SWRITE(UNIT_stdOut,'(A)')'Warning: The chosen DFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_stdOut,'(132("!"))')
END IF
DFLScale(0) = DFLScale(0)/(2.*Nin_DFL+1.)**2
SWRITE(UNIT_stdOut,'(A,2ES16.7)') '   DFL (DG/FV):',DFLScale
DFLScale_Readin = DFLScale
#else
dummy = Nin_DFL ! prevent compile warning
#endif /*PARABOLIC*/
END SUBROUTINE fillCFL_DFL


!==================================================================================================================================
!> Finalizes variables necessary for timedisc subroutines
!==================================================================================================================================
SUBROUTINE FinalizeTimeDisc()
! MODULES
USE MOD_TimeDisc_Vars
IMPLICIT NONE
!==================================================================================================================================
TimeDiscInitIsDone = .FALSE.
SDEALLOCATE(dtElem)
SDEALLOCATE(Ut_tmp)
SDEALLOCATE(S2)
SDEALLOCATE(UPrev)
SDEALLOCATE(RKA)
SDEALLOCATE(RKb)
SDEALLOCATE(RKc)
SDEALLOCATE(RKg1)
SDEALLOCATE(RKg2)
SDEALLOCATE(RKg3)
SDEALLOCATE(RKdelta)
NULLIFY(TimeStep)
END SUBROUTINE FinalizeTimeDisc

END MODULE MOD_TimeDisc_Functions
