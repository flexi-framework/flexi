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
!> Basic routines performing an analysis of the solution valid for all testcases
!==================================================================================================================================
MODULE MOD_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: DefineParametersAnalyze
PUBLIC:: InitAnalyze
PUBLIC:: Analyze
PUBLIC:: FinalizeAnalyze
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyze()
! MODULES
USE MOD_ReadInTools,        ONLY: prms
USE MOD_AnalyzeEquation,    ONLY: DefineParametersAnalyzeEquation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Analyze")
CALL prms%CreateLogicalOption('CalcErrorNorms',  "Set true to compute L2 and LInf error norms at analyze step.",&
                                                 '.TRUE.')
CALL prms%CreateLogicalOption('AnalyzeToFile',   "Set true to output result of error norms to a file (CalcErrorNorms=T)",&
                                                 '.FALSE.')
CALL prms%CreateRealOption(   'analyze_dt',      "Specifies time intervall at which analysis routines are called.",&
                                                 '0.')
CALL prms%CreateIntOption(    'nWriteData',      "Intervall as multiple of analyze_dt at which HDF5 files "//&
                                                 "(e.g. State,TimeAvg,Fluc) are written.",&
                                                 '1')
CALL prms%CreateIntOption(    'NAnalyze',        "Polynomial degree at which analysis is performed (e.g. for L2 errors). "//&
                                                 "Default: 2*N.")
CALL prms%CreateIntOption(    'AnalyzeExactFunc',"Define exact function used for analyze (e.g. for computing L2 errors). "//&
                                                 "Default: Same as IniExactFunc")
CALL prms%CreateIntOption(    'AnalyzeRefState' ,"Define state used for analyze (e.g. for computing L2 errors). "//&
                                                 "Default: Same as IniRefState")
CALL prms%CreateLogicalOption('doMeasureFlops',  "Set true to measure flop count, if compiled with PAPI.",&
                                                 '.TRUE.')
CALL prms%CreateRealOption(   'PIDkill',         'Kill FLEXI if PID gets above this value (optional)',&
                                                 '-1.0')
CALL prms%CreateIntOption(    'NCalcPID'         ,'Compute PID after every Nth timestep.',&
                                                  '1')
CALL DefineParametersAnalyzeEquation()

END SUBROUTINE DefineParametersAnalyze


!==================================================================================================================================
!> Initializes variables necessary for analyze subroutines
!> - provides basic quantities like global domain volume, surface area of boundary conditions
!>   or precomputed surface and volume integration weights
!> - initializes other specific analysis and benchmarking routines
!==================================================================================================================================
SUBROUTINE InitAnalyze()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation,    ONLY: InitAnalyzeEquation
USE MOD_ReadInTools,        ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_StringTools,        ONLY: INTTOSTR
USE MOD_Interpolation_Vars, ONLY: xGP,wGP,wBary,InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: nBCs,SurfElem,nSides,AnalyzeSide,sJ,nElems
USE MOD_Equation_Vars,      ONLY: StrVarNames,IniExactFunc,IniRefState
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_Benchmarking,       ONLY: InitBenchmarking
USE MOD_Timedisc_Vars,      ONLY: TEnd
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                       :: lastLine(0:2*PP_nVar+5)
INTEGER                                    :: i,j,k,iSurf,iElem,iSide
LOGICAL                                    :: hasAnalyzeSides(nBCs)
! determine character length based on longest entry in 'StrVarNames'
! to this end, LEN_TRIM was previously called on entire string array, while standard only permits scalar strings as argument
! therefore, use implicit loop to call LEN_TRIM on each array entry, composing the integer array of string lengths on the fly
CHARACTER(MAXVAL( (/( LEN_TRIM(StrVarNames(i)), i=1,SIZE(StrVarNames(:)) )/) )+5)  :: VarNames(2*PP_nVar+5)
!==================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.AnalyzeInitIsDone) &
  CALL CollectiveStop(__STAMP__,'InitAnalyse not ready to be called or already called.')

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ANALYZE...'

! Get the various analysis/output variables
doCalcErrorNorms = GETLOGICAL('CalcErrorNorms')
doAnalyzeToFile  = GETLOGICAL('AnalyzeToFile')
AnalyzeExactFunc = GETINT('AnalyzeExactFunc',INTTOSTR(IniExactFunc))
AnalyzeRefState  = GETINT('AnalyzeRefState' ,INTTOSTR(IniRefState))

analyze_dt       = GETREAL('analyze_dt')
nWriteData       = GETINT('nWriteData')
NAnalyze         = GETINT('NAnalyze'   ,INTTOSTR(2*(PP_N+1)))
#if PP_dim == 3
NAnalyzeZ        = NAnalyze
#else
NAnalyzeZ        = 0
#endif
! If analyze_dt is set to 0 (default) or to a negative value, no analyze calls should be performed at all.
! To achieve this, analyze_dt is set to the final simulation time. This will prevent any calls of the analyze routine
! except at the beginning and the end of the simulation.
IF (analyze_dt.LE.0.) THEN
  analyze_dt = TEnd
  nWriteData = 1
END IF

WriteData_dt = analyze_dt*nWriteData

! precompute integration weights
ALLOCATE(wGPSurf(0:PP_N,0:PP_NZ),wGPVol(0:PP_N,0:PP_N,0:PP_NZ))
#if PP_dim == 3
DO j=0,PP_N; DO i=0,PP_N
  wGPSurf(i,j)  = wGP(i)*wGP(j)
END DO; END DO
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  wGPVol(i,j,k) = wGP(i)*wGP(j)*wGP(k)
END DO; END DO; END DO
#else
DO i=0,PP_N
  wGPSurf(i,0)  = wGP(i)
END DO
DO j=0,PP_N; DO i=0,PP_N
  wGPVol(i,j,0) = wGP(i)*wGP(j)
END DO; END DO
#endif

! precompute volume of the domain
ALLOCATE(ElemVol(nElems))
ElemVol=0.
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    ElemVol(iElem)=ElemVol(iElem)+wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO !i,j,k
END DO ! iElem
Vol=SUM(ElemVol)

! compute surface of each boundary
ALLOCATE(Surf(nBCs))
Surf=0.
hasAnalyzeSides=.FALSE.
DO iSide=1,nSides
  iSurf=AnalyzeSide(iSide)
  IF(iSurf.EQ.0) CYCLE
  hasAnalyzeSides(iSurf)=.TRUE.
  DO j=0,PP_NZ; DO i=0,PP_N
    Surf(iSurf)=Surf(iSurf)+wGPSurf(i,j)*SurfElem(i,j,0,iSide)
  END DO; END DO
END DO
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Vol ,1   ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Surf,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
! Communicate whether any processor has a surface at the respective boundary
CALL MPI_ALLREDUCE(MPI_IN_PLACE,hasAnalyzeSides,nBCs,MPI_LOGICAL,MPI_LOR,MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/
! Prevent division by 0 if a BC has no sides associated with it (e.g. periodic)
DO iSurf=1,nBCs
  IF (.NOT.hasAnalyzeSides(iSurf)) Surf(iSurf) = HUGE(1.)
END DO

! Initialize eval routines
CALL InitAnalyzeBasis(PP_N,NAnalyze,xGP,wBary)

IF(doAnalyzeToFile.AND.MPIRoot)THEN
  DO i=1,PP_nVar
    VarNames(i)                   = 'L2_'//TRIM(StrVarNames(i))
    VarNames(PP_nVar+1:PP_nVar+i) = 'LInf_'//TRIM(StrVarNames(i))
  END DO
  VarNames(2*PP_nVar+1:2*PP_nVar+5) = &
    [CHARACTER(9) :: 'timesteps','t_CPU','DOF','Ncells','nProcs'] ! gfortran hates mixed length arrays
  FileName_ErrNorm='out.'//TRIM(ProjectName)
  CALL InitOutputToFile(FileName_ErrNorm,'Analyze',2*PP_nVar+5,VarNames,lastLine)
  iterRestart    =MAX(lastLine(2*PP_nVar+1),0.)
  calcTimeRestart=MAX(lastLine(2*PP_nVar+2),0.)
END IF

CALL InitAnalyzeEquation()
CALL InitBenchmarking()

! Read in kill PID
PID_kill     = GETREAL('PIDkill')
IF (PID_kill.GT.0) THEN
  nCalcPIDMax  = GETINT('nCalcPID')
  PIDTimeStart = FLEXITIME()
END IF

AnalyzeInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A,ES18.9)')' Volume of computational domain : ',Vol
SWRITE(UNIT_stdOut,'(A)')       ' INIT ANALYZE DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitAnalyze


!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!> - Builds Vandermonde to interpolate the solution onto a Gauss-Lobatto mesh at a higher polynomial degree
!> - Precomputes volume interpolation weights
!==================================================================================================================================
SUBROUTINE InitAnalyzeBasis(N_in,Nloc,xGP,wBary)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars,       ONLY: wGPVolAnalyze,Vdm_GaussN_NAnalyze
USE MOD_Basis,              ONLY: InitializeVandermonde
USE MOD_Interpolation,      ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars, ONLY: NodeTypeGL
#if FV_ENABLED == 1
USE MOD_Analyze_Vars,       ONLY: FV_Vdm_NAnalyze
USE MOD_FV_Basis,           ONLY: FV_Build_X_w_BdryX
USE MOD_FV_Vars,            ONLY: FV_CellType
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: N_in                  !< input polynomial degree
INTEGER,INTENT(IN)               :: Nloc                  !< polynomial degree of analysis polynomial
REAL,INTENT(IN)                  :: xGP(0:N_in)           !< interpolation points
REAL,INTENT(IN)                  :: wBary(0:N_in)         !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: XiAnalyze(0:Nloc)
REAL                             :: wAnalyze (0:Nloc)
INTEGER                          :: i,j
#if PP_dim == 3
INTEGER                          :: k
#endif
#if FV_ENABLED == 1
REAL                             :: dummy1(0:Nloc),dummy2(0:Nloc)
REAL                             :: FV_BdryX(0:Nloc+1)
#endif
!==================================================================================================================================
ALLOCATE(wGPVolAnalyze(0:Nloc,0:Nloc,0:ZDIM(Nloc)),Vdm_GaussN_NAnalyze(0:Nloc,0:N_in))
CALL GetNodesAndWeights(Nloc,NodeTypeGL,XiAnalyze,wAnalyze)
CALL InitializeVandermonde(N_in,Nloc,wBary,xGP,XiAnalyze,Vdm_GaussN_NAnalyze)

#if PP_dim == 3
DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
  wGPVolAnalyze(i,j,k) = wAnalyze(i)*wAnalyze(j)*wAnalyze(k)
END DO; END DO; END DO
#else
DO j=0,Nloc; DO i=0,Nloc
  wGPVolAnalyze(i,j,0) = wAnalyze(i)*wAnalyze(j)
END DO; END DO
#endif

#if FV_ENABLED == 1
CALL FV_Build_X_w_BdryX(Nloc, dummy1, dummy2, FV_BdryX, FV_CellType)
ALLOCATE(FV_Vdm_NAnalyze(0:Nloc,0:N_in))
CALL InitializeVandermonde(N_in,Nloc,wBary,xGP,FV_BdryX,FV_Vdm_NAnalyze)
#endif

END SUBROUTINE InitAnalyzeBasis


!==================================================================================================================================
!> Controls analysis routines and is called at analyze time levels
!> - calls generic error norm computation
!> - calls equation system specific analysis
!==================================================================================================================================
SUBROUTINE Analyze(Time,iter)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,       ONLY: StartTime
USE MOD_PreProc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation,    ONLY: AnalyzeEquation
USE MOD_Benchmarking,       ONLY: Benchmarking
USE MOD_Mesh_Vars,          ONLY: nGlobalElems
USE MOD_Output,             ONLY: OutputToFile,PrintStatusLine
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_TimeDisc_Vars,      ONLY: dt,tStart,tEnd,maxIter
USE MOD_Restart_Vars,       ONLY: RestartTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                   !< current simulation time
INTEGER(KIND=DP),INTENT(IN)     :: iter                   !< current iteration
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)               :: formatStr
CHARACTER(LEN=255)              :: hilf
REAL                            :: CalcTime,RunTime
REAL                            :: L_Inf_Error(PP_nVar),L_2_Error(PP_nVar)
!==================================================================================================================================
! Graphical output
CalcTime = FLEXITIME()
RunTime  = CalcTime-StartTime
SWRITE(UNIT_stdOut,'(A14,ES16.7)')' Sim time   : ',Time

! Calculate error norms
IF(doCalcErrorNorms)THEN
  CALL CalcErrorNorms(Time,L_2_Error,L_Inf_Error)
  IF(MPIRoot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',PP_nVar,'ES18.9)'
    WRITE(UNIT_stdOut,formatStr)' L_2        : ',L_2_Error
    WRITE(UNIT_stdOut,formatStr)' L_inf      : ',L_Inf_Error
    IF(doAnalyzeToFile)THEN
      CALL OutputToFile(FileName_ErrNorm,(/Time/),(/2*PP_nVar+5,1/),  (/L_2_Error,L_Inf_Error, &
                        REAL(iter)+iterRestart,RunTime+CalcTimeRestart,              &
                        REAL(nGlobalElems*(PP_N+1)**3),REAL(nGlobalElems),REAL(nProcessors)/))
    END IF
  END IF !MPIRoot
END IF  ! ErrorNorms

CALL AnalyzeEquation(Time)
CALL AnalyzePerformance()
CALL Benchmarking()

IF(Time.GT.RestartTime) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  CALL PrintStatusLine(time,dt,tStart,tEnd,iter,maxIter,doETA=.TRUE.)
  WRITE(hilf,'(A,A,A)') 'RUNNING ',TRIM(ProjectName),'...'
  CALL DisplaySimulationTime(CalcTime, StartTime, hilf)
  ! SWRITE(UNIT_stdOut,'(132("-"))')
  ! SWRITE(UNIT_stdOut,*)
END IF
END SUBROUTINE Analyze


!==================================================================================================================================
!> Calculates current code performance
!==================================================================================================================================
SUBROUTINE AnalyzePerformance()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars,       ONLY: PID,PIDTimeStart,PIDTimeEnd,PID_kill,nCalcPID,nCalcPIDMax
USE MOD_Mesh_Vars,          ONLY: nGlobalElems
USE MOD_TimeDisc_Vars,      ONLY: nRKStages
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

IF (PID_kill.LE.0) RETURN

! Return if no timestep update requested in this iteration
IF (nCalcPID.GE.1) THEN
  nCalcPID = nCalcPID - 1
  RETURN
ELSE
  nCalcPID = nCalcPIDMax - 1
END IF

! Get calculation time per DOF
PIDTimeEnd = FLEXITIME()
PID        = (PIDTimeEnd-PIDTimeStart)*REAL(nProcessors)/(REAL(nGlobalElems)*REAL((PP_N+1)**PP_dim))/nRKStages

! Abort if PID is too high
IF (PID.GT.PID_kill) &
  CALL CollectiveStop(__STAMP__,'Aborting due to low performance, PID',RealInfo=PID)

PIDTimeStart = PIDTimeStart

END SUBROUTINE AnalyzePerformance


!==================================================================================================================================
!> Calculates L_infinfity and L_2 error norms of state variables using the analyze framework (GL points+weights)
!==================================================================================================================================
SUBROUTINE CalcErrorNorms(Time,L_2_Error,L_Inf_Error)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: Elem_xGP,sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Exactfunc,          ONLY: ExactFunc
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_Analyze_Vars,       ONLY: NAnalyze,NAnalyzeZ,Vdm_GaussN_NAnalyze
USE MOD_Analyze_Vars,       ONLY: wGPVolAnalyze,Vol,AnalyzeExactFunc,AnalyzeRefState
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_Vdm,FV_w
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                   !< current simulation time
REAL,INTENT(OUT)                :: L_2_Error(  PP_nVar)   !< L2 error of the solution
REAL,INTENT(OUT)                :: L_Inf_Error(PP_nVar)   !< LInf error of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem,k,l,m
REAL                            :: U_exact(PP_nVar)
REAL                            :: U_NAnalyze(   1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyzeZ)
REAL                            :: Coords_NAnalyze(      3,0:NAnalyze,0:NAnalyze,0:NAnalyzeZ)
REAL                            :: J_NAnalyze(           1,0:NAnalyze,0:NAnalyze,0:NAnalyzeZ)
REAL                            :: J_N(       1,0:PP_N,0:PP_N,0:PP_NZ)
REAL                            :: IntegrationWeight
#if FV_ENABLED
REAL                            :: U_DG(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                            :: U_FV(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================
! Calculate error norms
L_Inf_Error(:) = -HUGE(1.)
L_2_Error(  :) = 0.

! Interpolate values of Error-Grid from GP's
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) THEN ! FV Element
    DO m=0,PP_NZ; DO l=0,PP_N; DO k=0,PP_N
      CALL ExactFunc(AnalyzeExactFunc,time,Elem_xGP(1:3,k,l,m,iElem),U_DG(:,k,l,m),AnalyzeRefState)
    END DO; END DO; END DO ! k, l, m

    ! Project the DG solution to the FV subcells
    CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,U_DG(:,:,:,:),U_FV(:,:,:,:))

    ! Calculate the FV errors on the FV subcells
    DO m=0,PP_NZ; DO l=0,PP_N; DO k=0,PP_N
      L_Inf_Error = MAX(L_Inf_Error,ABS(U(:,k,l,m,iElem) - U_FV(:,k,l,m)))
      IntegrationWeight = FV_w(k)*FV_w(l)*FV_w(m)/sJ(k,l,m,iElem,1)
      ! To sum over the elements, we compute here the square of the L_2 error
      L_2_Error = L_2_Error+(U(:,k,l,m,iElem) - U_FV(:,k,l,m))*(U(:,k,l,m,iElem) - U_FV(:,k,l,m))*IntegrationWeight
    END DO; END DO; END DO ! k, l, m
  ELSE
#endif
   ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
   CALL ChangeBasisVolume(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Elem_xGP(1:3,:,:,:,iElem),Coords_NAnalyze(1:3,:,:,:))

   ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the Jacobian ;-)
   J_N(1,0:PP_N,0:PP_N,0:PP_NZ)=1./sJ(:,:,:,iElem,0)
   CALL ChangeBasisVolume(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N,J_NAnalyze)

   ! Interpolate the solution to the analyze grid
   CALL ChangeBasisVolume(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,U(1:PP_nVar,:,:,:,iElem),U_NAnalyze(1:PP_nVar,:,:,:))

  ! Calculate the DG errors on the analyze grid
   DO m=0,NAnalyzeZ; DO l=0,NAnalyze; DO k=0,NAnalyze
     CALL ExactFunc(AnalyzeExactFunc,time,Coords_NAnalyze(1:3,k,l,m),U_exact,AnalyzeRefState)
     L_Inf_Error = MAX(L_Inf_Error,abs(U_NAnalyze(:,k,l,m) - U_exact))
     IntegrationWeight = wGPVolAnalyze(k,l,m)*J_NAnalyze(1,k,l,m)
     ! To sum over the elements, we compute here the square of the L_2 error
     L_2_Error = L_2_Error+(U_NAnalyze(:,k,l,m) - U_exact)*(U_NAnalyze(:,k,l,m) - U_exact)*IntegrationWeight
    END DO; END DO; END DO ! k, l, m
#if FV_ENABLED
  END IF
#endif
END DO ! iElem=1,nElems

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,L_2_Error  ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,L_Inf_Error,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(L_2_Error  ,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(L_Inf_Error,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
END IF
#endif

! We normalize the L_2 Error with the volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/Vol)

END SUBROUTINE CalcErrorNorms


!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyze()
! MODULES
USE MOD_AnalyzeEquation,    ONLY: FinalizeAnalyzeEquation
USE MOD_Analyze_Vars,       ONLY: AnalyzeInitIsDone,wGPSurf,wGPVol,Surf,wGPVolAnalyze,Vdm_GaussN_NAnalyze,ElemVol
#if FV_ENABLED == 1
USE MOD_Analyze_Vars,       ONLY: FV_Vdm_NAnalyze
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL FinalizeAnalyzeEquation()
SDEALLOCATE(Vdm_GaussN_NAnalyze)
SDEALLOCATE(wGPVolAnalyze)
SDEALLOCATE(Surf)
SDEALLOCATE(wGPVol)
SDEALLOCATE(wGPSurf)
SDEALLOCATE(ElemVol)
#if FV_ENABLED == 1
SDEALLOCATE(FV_Vdm_NAnalyze)
#endif
AnalyzeInitIsDone = .FALSE.

END SUBROUTINE FinalizeAnalyze

END MODULE MOD_Analyze
