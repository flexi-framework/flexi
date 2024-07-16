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
!> Provides routines for visualization and ASCII output of time series. Initializes some variables used for HDF5 output.
!==================================================================================================================================
MODULE MOD_Output
! MODULES
USE MOD_ReadInTools
USE ISO_C_BINDING
IMPLICIT NONE

INTERFACE
  SUBROUTINE print_userblock(filename,inifilename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR
      CHARACTER(KIND=C_CHAR) :: filename(*)
      CHARACTER(KIND=C_CHAR) :: inifilename(*)
  END SUBROUTINE print_userblock
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

! Output format for state visualization
INTEGER,PARAMETER :: OUTPUTFORMAT_NONE         = 0
INTEGER,PARAMETER :: OUTPUTFORMAT_TECPLOT      = 1
INTEGER,PARAMETER :: OUTPUTFORMAT_TECPLOTASCII = 2
INTEGER,PARAMETER :: OUTPUTFORMAT_PARAVIEW     = 3

! Output format for ASCII data files
INTEGER,PARAMETER :: ASCIIOUTPUTFORMAT_CSV     = 0
INTEGER,PARAMETER :: ASCIIOUTPUTFORMAT_TECPLOT = 1

INTERFACE DefineParametersOutput
  MODULE PROCEDURE DefineParametersOutput
END INTERFACE

INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

INTERFACE PrintPercentage
  MODULE PROCEDURE PrintPercentage
END INTERFACE

INTERFACE PrintStatusLine
  MODULE PROCEDURE PrintStatusLine
END INTERFACE

INTERFACE Visualize
  MODULE PROCEDURE Visualize
END INTERFACE

INTERFACE InitOutputToFile
  MODULE PROCEDURE InitOutputToFile
END INTERFACE

INTERFACE OutputToFile
  MODULE PROCEDURE OutputToFile
END INTERFACE

INTERFACE FinalizeOutput
  MODULE PROCEDURE FinalizeOutput
END INTERFACE

PUBLIC:: InitOutput,PrintPercentage,PrintStatusLine,Visualize,InitOutputToFile,OutputToFile,FinalizeOutput
PUBLIC:: print_userblock
!==================================================================================================================================

PUBLIC::DefineParametersOutput
CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersOutput()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Output")
CALL prms%CreateIntOption(          'NVisu',       "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateIntOption(          'NOut',        "Polynomial degree at which solution is written. -1: NOut=N, >0: NOut", '-1')
CALL prms%CreateStringOption(       'ProjectName', "Name of the current simulation (mandatory).")
CALL prms%CreateLogicalOption(      'Logging',     "Write log files containing debug output.", '.FALSE.')
CALL prms%CreateLogicalOption(      'ErrorFiles',  "Write error files containing error output.", '.TRUE.')
CALL prms%CreateIntFromStringOption('OutputFormat',"File format for visualization: None, Tecplot, TecplotASCII, ParaView. "//&
                                                 " Note: Tecplot output is currently unavailable due to licensing issues.", 'None')
CALL addStrListEntry('OutputFormat','none',        OUTPUTFORMAT_NONE)
CALL addStrListEntry('OutputFormat','tecplot',     OUTPUTFORMAT_TECPLOT)
CALL addStrListEntry('OutputFormat','tecplotascii',OUTPUTFORMAT_TECPLOTASCII)
CALL addStrListEntry('OutputFormat','paraview',    OUTPUTFORMAT_PARAVIEW)
CALL prms%CreateIntFromStringOption('ASCIIOutputFormat',"File format for ASCII files, e.g. body forces: CSV, Tecplot."&
                                                       , 'CSV')
CALL addStrListEntry('ASCIIOutputFormat','csv',    ASCIIOUTPUTFORMAT_CSV)
CALL addStrListEntry('ASCIIOutputFormat','tecplot',ASCIIOUTPUTFORMAT_TECPLOT)
CALL prms%CreateLogicalOption(      'doPrintStatusLine','Print: percentage of time, ...', '.FALSE.')
CALL prms%CreateLogicalOption(      'WriteStateFiles','Write HDF5 state files. Disable this only for debugging issues. \n'// &
                                                      'NO SOLUTION WILL BE WRITTEN!', '.TRUE.')
CALL prms%CreateLogicalOption(      'WriteTimeAvgFiles','Write HDF5 time average files. Disable this only for debugging. \n'// &
                                                      'NO TIME AVERAGE FILES WILL BE WRITTEN!', '.TRUE.')
END SUBROUTINE DefineParametersOutput

!==================================================================================================================================
!> Initialize all output variables.
!==================================================================================================================================
SUBROUTINE InitOutput()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Output_Vars
USE MOD_ReadInTools       ,ONLY:GETSTR,GETLOGICAL,GETINT,GETINTFROMSTR
USE MOD_StringTools       ,ONLY:INTTOSTR
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone,NodeTypeVISU,NodeType
USE ISO_C_BINDING,         ONLY:C_NULL_CHAR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
CHARACTER(LEN=8)               :: StrDate
CHARACTER(LEN=10)              :: StrTime
CHARACTER(LEN=255)             :: LogFile
!==================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.OutputInitIsDone) &
  CALL CollectiveStop(__STAMP__,'InitOutput not ready to be called or already called.')

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

NVisu=GETINT('NVisu',INTTOSTR(PP_N))

! Gauss/Gl -> Visu : computation -> visualization
ALLOCATE(Vdm_GaussN_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVISU,Vdm_GaussN_NVisu)

! Output polynomial degree (to reduce storage,e.g.in case of overintegration)
NOut=GETINT('NOut') ! -1: PP_N(off), >0:custom
IF(NOut.EQ.-1)THEN
  NOut=PP_N
END IF
! Get Vandermonde for output changebasis
IF(NOut.NE.PP_N)THEN
  ALLOCATE(Vdm_N_NOut(0:NOut,0:PP_N))
  CALL GetVandermonde(PP_N,NodeType,NOut,NodeType,Vdm_N_NOut,modal=.TRUE.)
END IF

! Name for all output files
ProjectName = GETSTR('ProjectName')
Logging     = GETLOGICAL('Logging')
ErrorFiles  = GETLOGICAL('ErrorFiles')

doPrintStatusLine=GETLOGICAL("doPrintStatusLine")
WriteStateFiles=GETLOGICAL("WriteStateFiles")
IF (.NOT.WriteStateFiles) CALL PrintWarning("Write of state files disabled!")
WriteTimeAvgFiles=GETLOGICAL("WriteTimeAvgFiles")
IF (.NOT.WriteTimeAvgFiles) CALL PrintWarning("Write of time average files disabled!")

IF (MPIRoot) THEN
  ! prepare userblock file
  CALL print_userblock(TRIM(UserBlockTmpFile)//C_NULL_CHAR,TRIM(ParameterFile)//C_NULL_CHAR)
  INQUIRE(FILE=TRIM(UserBlockTmpFile),SIZE=userblock_total_len)
END IF

WRITE(ErrorFileName,'(A,A8,I6.6,A4)')TRIM(ProjectName),'_ERRORS_',myRank,'.out'

! Get output format for state visualization
OutputFormat = GETINTFROMSTR('OutputFormat')
! Get output format for ASCII data files
ASCIIOutputFormat = GETINTFROMSTR('ASCIIOutputFormat')

! Open file for logging
IF(Logging)THEN
  WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
  OPEN(UNIT     = UNIT_logOut,  &
       FILE     = LogFile,      &
       STATUS   = 'UNKNOWN',    &
       ACTION   = 'WRITE',      &
       POSITION = 'APPEND',     &
       IOSTAT   = OpenStat)
  CALL DATE_AND_TIME(StrDate,StrTime)
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,'(132("#"))')
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                      StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
END IF  ! Logging

OutputInitIsDone =.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitOutput


!==================================================================================================================================
!> Displays the current process of a given procedure
!==================================================================================================================================
SUBROUTINE PrintPercentage(NameOpt,percent)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars   ,ONLY: doPrintStatusLine
USE MOD_ReadInTools   ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: NameOpt        !< Option name
REAL,INTENT(IN)             :: percent        !< current progress percentage
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=20)    :: fmtName
!==================================================================================================================================

IF (.NOT.doPrintStatusLine) RETURN
IF (.NOT.MPIRoot)           RETURN

WRITE(fmtName,*) prms%maxNameLen
WRITE(UNIT_StdOut,'(A3,A'//TRIM(fmtName)//',A2)',ADVANCE='NO')' | ',TRIM(NameOpt),' | '
WRITE(UNIT_stdOut,'(A,A1,A)',ADVANCE='NO') REPEAT('=',MAX(CEILING(percent*(prms%maxValueLen+1)/100.),0)),'>',&
                                           REPEAT(' ',(prms%maxValueLen+1)-MAX(CEILING(percent*(prms%maxValueLen+1)/100.),0))
IF (percent.LT.100) THEN
  WRITE(UNIT_stdOut,'(A2,F6.2,A3,A1)',ADVANCE='NO' ) '|[',percent,'%]|',ACHAR(13) ! ACHAR(13) is carriage
ELSE
  WRITE(UNIT_stdOut,'(A2,F6.2,A3)'   ,ADVANCE='YES') '|[',percent,'%]|'
END IF

END SUBROUTINE PrintPercentage


!==================================================================================================================================
!> Displays the actual status of the simulation and counts the amount of FV elements
!==================================================================================================================================
SUBROUTINE PrintStatusLine(t,dt,tStart,tEnd,iter,maxIter,doETA)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: epsMach
USE MOD_PreProc
USE MOD_Output_Vars   ,ONLY: doPrintStatusLine
USE MOD_Restart_Vars  ,ONLY: DoRestart,RestartTime
USE MOD_TimeDisc_Vars ,ONLY: time_start
#if (FV_ENABLED == 1) || PP_LIMITER
USE MOD_Mesh_Vars     ,ONLY: nGlobalElems
#endif
#if FV_ENABLED == 1
USE MOD_Analyze_Vars  ,ONLY: totalFV_nElems
USE MOD_FV_Vars       ,ONLY: FV_Elems
#elif FV_ENABLED == 2 || FV_ENABLED == 3
USE MOD_Analyze_Vars  ,ONLY: FV_totalAlpha
USE MOD_FV_Vars       ,ONLY: FV_alpha
#endif /*FV_ENABLED*/
#if PP_LIMITER
USE MOD_Analyze_Vars, ONLY: totalPP_nElems
USE MOD_Filter_Vars , ONLY: PP_Elems
#endif /*PP_LIMITER*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)             :: t      !< current simulation time
REAL,INTENT(IN)             :: dt     !< current time step
REAL,INTENT(IN)             :: tStart !< start time of simulation
REAL,INTENT(IN)             :: tEnd   !< end time of simulation
INTEGER(KIND=8),INTENT(IN)  :: iter    !< current iteration
INTEGER(KIND=8),INTENT(IN)  :: maxIter !< end iteration of simulation
LOGICAL,INTENT(IN),OPTIONAL :: doETA !< flag to print ETA without carriage return
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL           :: doETA_loc
REAL              :: percent,percent_time,percent_iter,percent_ETA
REAL              :: time_remaining,mins,secs,hours,days
CHARACTER(LEN=3)  :: tmpString
#if !FV_ENABLED && PP_LIMITER
INTEGER,PARAMETER :: barWidth = 39
#elif (FV_ENABLED==1) &&  PP_LIMITER
INTEGER,PARAMETER :: barWidth = 27
#elif (FV_ENABLED==1) && !PP_LIMITER
INTEGER,PARAMETER :: barWidth = 39
#elif (FV_ENABLED==2 || FV_ENABLED == 3) &&  PP_LIMITER
INTEGER,PARAMETER :: barWidth = 20
#elif (FV_ENABLED==2 || FV_ENABLED == 3) && !PP_LIMITER
INTEGER,PARAMETER :: barWidth = 32
#else
INTEGER,PARAMETER :: barWidth = 51
#endif
#if FV_ENABLED == 1
INTEGER(KIND=8)   :: FVcounter
REAL              :: FV_percent
#elif FV_ENABLED == 2 || FV_ENABLED == 3
REAL              :: FV_alpha_range(2)
#endif /*FV_ENABLED*/
#if PP_LIMITER
INTEGER(KIND=8)   :: PPcounter
REAL              :: PP_percent
#endif /*PP_LIMITER*/
!==================================================================================================================================
#if FV_ENABLED == 1
FVcounter      = INT(SUM(FV_Elems),KIND=8)
totalFV_nElems = totalFV_nElems + FVcounter ! counter for output of FV amount during analyze
#elif FV_ENABLED == 2 || FV_ENABLED == 3
FV_alpha_range(1) = MINVAL(FV_alpha)
FV_alpha_range(2) = MAXVAL(FV_alpha)
FV_totalAlpha    = FV_totalAlpha + SUM(FV_alpha)
#endif /*FV_ENABLED*/
#if PP_LIMITER
PPcounter      = INT(SUM(PP_Elems),KIND=8)
totalPP_nElems = totalPP_nElems + PPcounter ! counter for output of FV amount during analyze
#endif

IF (PRESENT(doETA)) THEN
  doETA_loc = doETA
ELSE
  doETA_loc = .FALSE.
END IF

IF(.NOT.doPrintStatusLine .AND. .NOT.doETA_loc) RETURN

#if FV_ENABLED && USE_MPI
IF (MPIRoot) THEN
#if FV_ENABLED == 1
  CALL MPI_REDUCE(MPI_IN_PLACE    ,FVcounter       ,1,MPI_INTEGER8        ,MPI_SUM,0,MPI_COMM_FLEXI,iError)
#elif FV_ENABLED == 2 || FV_ENABLED == 3
  CALL MPI_REDUCE(MPI_IN_PLACE    ,FV_alpha_range(1),1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE    ,FV_alpha_range(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
#endif /*FV_ENABLED*/
ELSE
#if FV_ENABLED == 1
  CALL MPI_REDUCE(FVcounter       ,0               ,1,MPI_INTEGER8        ,MPI_SUM,0,MPI_COMM_FLEXI,iError)
#elif FV_ENABLED == 2 || FV_ENABLED == 3
  CALL MPI_REDUCE(FV_alpha_range(1),0               ,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(FV_alpha_range(2),0               ,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
#endif /*FV_ENABLED*/
END IF
#endif

#if PP_LIMITER && USE_MPI
IF (MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PPcounter,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(PPcounter,0           ,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif /*PP_LIMITER && USE_MPI*/

IF(MPIRoot)THEN
#ifdef INTEL
  OPEN(UNIT_stdOut,CARRIAGECONTROL='fortran')
#endif
  percent_time = (t-tStart) / (tEnd-tStart)
  percent_iter = REAL(iter) / REAL(maxIter)
  percent      = ABS(MAX(percent_time,percent_iter))
  ! ETA bar needs a tiny amount
  percent      = MAX(percent,TINY(1.))

  ! Calculate ETA with percent of current run
  ASSOCIATE(tBegin => MERGE(RestartTime,tStart,DoRestart))
  percent_ETA  = (t-tBegin) / (tEnd-tBegin)
  percent_ETA  = MAX(percent_ETA,percent_iter)
  END ASSOCIATE

  CALL CPU_TIME(time_remaining)
  time_remaining = time_remaining - time_start
  IF (percent_ETA.GT.0.0) time_remaining = time_remaining/percent_ETA - time_remaining
  percent = percent*100.
  secs = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  mins = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  hours = MOD(time_remaining,24.)
  days = time_remaining / 24

#if FV_ENABLED == 1
  FV_percent = REAL(FVcounter) / REAL(nGlobalElems) * 100.
  WRITE(UNIT_stdOut,'(F7.2,A5)',ADVANCE='NO') FV_percent, '% FV,'
#elif FV_ENABLED == 2 || FV_ENABLED == 3
  WRITE(UNIT_stdOut,'(A7,F5.3,A1,F5.3,A)',ADVANCE='NO') ' aFV = ',FV_alpha_range(1),'-',FV_alpha_range(2),','
#endif /*FV_ENABLED*/

#if PP_LIMITER
  PP_percent = REAL(PPcounter) / REAL(nGlobalElems) * 100.
  WRITE(UNIT_stdOut,'(F7.2,A5)',ADVANCE='NO') PP_percent, '% PP,'
#endif /*PP_LIMITER*/

  ! Status line or standard output
  IF (PRESENT(doETA)) THEN; tmpString = 'YES'
  ELSE                    ; tmpString = 'NO'
  ENDIF
END IF

IF (.NOT.MPIRoot) RETURN

IF (mins.LT.1 .AND. hours.EQ.0 .AND. days.EQ.0) THEN
    WRITE(UNIT_stdOut,'(A,E10.4,A,E10.4,A,A,A3,A,A1,A,A3,F6.2,A3,A1)',ADVANCE=tmpString)                               &
    '   Time = ', t,'  dt = ', dt, ' ', ' ETA [d:h:m]:<1 min remaining',' |',                                          &
    REPEAT('=',MAX(CEILING(percent*barWidth/100.)-1,0)),'>',REPEAT(' ',barWidth-MAX(CEILING(percent*barWidth/100.),0)),'| [',percent,'%] ',&
    ACHAR(13) ! ACHAR(13) is carriage return
  ELSE
    WRITE(UNIT_stdOut,'(A,E10.4,A,E10.4,A,A,I4,A1,I0.2,A1,I0.2,A1,I0.2,A7,A,A1,A,A3,F6.2,A3,A1)',ADVANCE=tmpString)    &
    '   Time = ', t,'  dt = ', dt, ' ', ' ETA [d:h:m]',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),' |',      &
    REPEAT('=',MAX(CEILING(percent*barWidth/100.)-1,0)),'>',REPEAT(' ',barWidth-MAX(CEILING(percent*barWidth/100.),0)),'| [',percent,'%] ',&
    ACHAR(13) ! ACHAR(13) is carriage return
END IF
#ifdef INTEL
 CLOSE(UNIT_stdOut)
#endif /*INTEL*/

END SUBROUTINE PrintStatusLine


!==================================================================================================================================
!> Displays the analysis of the simulation
!==================================================================================================================================
SUBROUTINE PrintAnalyze(dt)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY: SimulationEfficiency!,StartTime,WallTime
USE MOD_PreProc
USE MOD_Analyze_Vars        ,ONLY: PID
USE MOD_Implicit_Vars       ,ONLY: nGMRESIterGlobal,nNewtonIterGlobal
! USE MOD_Mesh_Vars           ,ONLY: nGlobalElems
! USE MOD_Restart_Vars        ,ONLY: RestartTime,RestartWallTime
! USE MOD_TimeDisc_Vars       ,ONLY: CalcTimeStart,CalcTimeEnd
USE MOD_TimeDisc_Vars       ,ONLY: TimeDiscType,ViscousTimeStep
USE MOD_TimeDisc_Vars       ,ONLY: iter!,t,iter_analyze,nRKStages
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: dt     !< current time step
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: timeArray(8)              !> Array for system time
! REAL    :: WallTimeEnd               !> wall time of simulation end
!==================================================================================================================================

! ! Get calculation time per DOF
! CalcTimeEnd          = FLEXITIME()
! WallTimeEnd          = CalcTimeEnd
! WallTime             = WallTimeEnd-StartTime
! SimulationEfficiency = (t-RestartTime)/((WallTimeEnd-RestartWallTime)*nProcessors/3600.) ! in [s] / [CPUh]
! PID                  = (CalcTimeEnd-CalcTimeStart)*REAL(nProcessors)/(REAL(nGlobalElems)*REAL((PP_N+1)**PP_dim)*REAL(iter_analyze))/nRKStages
! Get system time
CALL DATE_AND_TIME(values=timeArray)

IF (.NOT.MPIRoot) RETURN

WRITE(UNIT_stdOut,*)
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A,I2.2,A1,I2.2,A1,I4.4,A1,I2.2,A1,I2.2,A1,I2.2)') &
  ' Sys date   :    ',timeArray(3),'.',timeArray(2),'.',timeArray(1),' ',timeArray(5),':',timeArray(6),':',timeArray(7)
WRITE(UNIT_stdOut,'(A,ES12.5,A)')' CALCULATION TIME PER STAGE/DOF: [',PID,' sec ]'
! WRITE(UNIT_stdOut,'(A,ES12.5,A)')' EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [',SimulationEfficiency,' sec/h ]'
WRITE(UNIT_stdOut,'(A,ES12.5,A)')' EFFICIENCY: CALCULATION TIME [s]/[Core-h]: [',SimulationEfficiency,' sec/h ]'
WRITE(UNIT_stdOut,'(A,ES16.7)')  ' Timestep   : ',dt
IF(ViscousTimeStep) WRITE(UNIT_stdOut,'(A)')' Viscous timestep dominates! '
WRITE(UNIT_stdOut,'(A,ES16.7)')  '#Timesteps  : ',REAL(iter)
IF(TimeDiscType.EQ.'ESDIRK') THEN
  WRITE(UNIT_stdOut,'(A,ES16.7)')'#GMRES iter : ',REAL(nGMRESIterGlobal)
  WRITE(UNIT_stdOut,'(A,ES16.7)')'#Newton iter: ',REAL(nNewtonIterGlobal)
  nGMRESIterGlobal  = 0
  nNewtonIterGlobal = 0
END IF

END SUBROUTINE PrintAnalyze


!==================================================================================================================================
!> Supersample DG dataset at (equidistant) visualization points and output to file.
!> Currently only Paraview binary format is supported.
!> Tecplot support has been removed due to licensing issues (possible GPL incompatibility).
!==================================================================================================================================
SUBROUTINE Visualize(OutputTime,U)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars,    ONLY:StrVarNames
USE MOD_Output_Vars,      ONLY:ProjectName,OutputFormat
USE MOD_Mesh_Vars  ,      ONLY:Elem_xGP,nElems
USE MOD_Output_Vars,      ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_ChangeBasisByDim, ONLY:ChangeBasisVolume
USE MOD_VTK,              ONLY:WriteDataToVTK,WriteVTKMultiBlockDataSet
#if FV_ENABLED
USE MOD_FV_Vars,          ONLY:FV_Elems
#if FV_RECONSTRUCT
USE MOD_FV_Vars,          ONLY:FV_dx_XI_L,FV_dx_XI_R,FV_dx_ETA_L,FV_dx_ETA_R
USE MOD_FV_Vars,          ONLY:gradUxi,gradUeta
#if PP_dim == 3
USE MOD_FV_Vars,          ONLY:FV_dx_ZETA_L,FV_dx_ZETA_R
USE MOD_FV_Vars,          ONLY:gradUzeta
#endif
#endif
USE MOD_EOS,              ONLY:PrimToCons,ConsToPrim
USE MOD_FV_Basis
USE MOD_Indicator_Vars,   ONLY:IndValue
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)               :: OutputTime                                !< simulation time of output
REAL,INTENT(INOUT)            :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution vector to be visualized
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem,FV_iElem,DG_iElem,PP_nVar_loc,nFV_Elems,iVar
REAL,ALLOCATABLE,TARGET       :: Coords_NVisu(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET       :: U_NVisu(:,:,:,:,:)
REAL,POINTER                  :: Coords_NVisu_p(:,:,:,:,:)
REAL,POINTER                  :: U_NVisu_p(:,:,:,:,:)
CHARACTER(LEN=255)            :: FileString_DG
#if FV_ENABLED
CHARACTER(LEN=255)            :: FileString_FV
CHARACTER(LEN=255)            :: FileString_multiblock
INTEGER                       :: i,j,k,NVisu_FV,iii,jjj,kkk,ii,jj,kk
REAL                          :: UPrim(1:PP_nVarPrim)
REAL                          :: UPrim2(1:PP_nVarPrim)
REAL,ALLOCATABLE,TARGET       :: FV_Coords_NVisu(:,:,:,:,:)
REAL,POINTER                  :: FV_Coords_NVisu_p(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET       :: FV_U_NVisu(:,:,:,:,:)
REAL,POINTER                  :: FV_U_NVisu_p(:,:,:,:,:)
REAL,ALLOCATABLE              :: Vdm_GaussN_NVisu_FV(:,:)
#endif
CHARACTER(LEN=255),ALLOCATABLE:: StrVarNames_loc(:)
#if FV_ENABLED && FV_RECONSTRUCT
REAL                          :: dx,dy
#if PP_dim == 3
REAL                          :: dz
#endif
#endif
!==================================================================================================================================
IF(outputFormat.LE.0) RETURN
! Specify output names

nFV_Elems = 0
PP_nVar_loc=PP_nVar
#if FV_ENABLED
PP_nVar_loc=PP_nVar+2
DO iElem=1,nElems
  IF (FV_Elems(iElem).GT.0) nFV_Elems = nFV_Elems + 1
END DO
NVisu_FV = (PP_N+1)*2-1
ALLOCATE(FV_U_NVisu(PP_nVar_loc,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),1:nFV_Elems))
ALLOCATE(FV_Coords_NVisu(1:3,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),1:nFV_Elems))
ALLOCATE(Vdm_GaussN_NVisu_FV(0:NVisu_FV,0:PP_N))
CALL FV_Build_VisuVdm(PP_N,Vdm_GaussN_NVisu_FV)
#endif
ALLOCATE(U_NVisu(PP_nVar_loc,0:NVisu,0:NVisu,0:ZDIM(NVisu),1:(nElems-nFV_Elems)))
U_NVisu = 0.
ALLOCATE(Coords_NVisu(1:3,0:NVisu,0:NVisu,0:ZDIM(NVisu),1:(nElems-nFV_Elems)))

DG_iElem = 0
FV_iElem = 0

DO iElem=1,nElems
#if FV_ENABLED
  ! DG Element
  IF (FV_Elems(iElem).EQ.0) THEN
#endif
    DG_iElem = DG_iElem+1
    ! Create coordinates of visualization points
    CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,DG_iElem))
    ! Interpolate solution onto visu grid
    CALL ChangeBasisVolume(PP_nVar,PP_N,NVisu,Vdm_GaussN_NVisu,U(1:PP_nVar,:,:,:,iElem),U_NVisu(1:PP_nVar,:,:,:,DG_iElem))
#if FV_ENABLED
    U_NVisu(PP_nVar_loc-1,:,:,:,DG_iElem) = IndValue(iElem)
    U_NVisu(PP_nVar_loc,:,:,:,DG_iElem) = FV_Elems(iElem)
  ! FV Element
  ELSE
    FV_iElem = FV_iElem+1

    CALL ChangeBasisVolume(3,PP_N,NVisu_FV,Vdm_GaussN_NVisu_FV,Elem_xGP(1:3,:,:,:,iElem),FV_Coords_NVisu(1:3,:,:,:,FV_iElem))
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      CALL ConsToPrim(UPrim ,U(:,i,j,k,iElem))
      DO kk=0,PP_dim-2; DO jj=0,1; DO ii=0,1
        kkk = k*2+kk
        jjj = j*2+jj
        iii = i*2+ii
#if FV_RECONSTRUCT
        dx = MERGE(  -FV_dx_XI_L(j,k,i,iElem),  FV_dx_XI_R(j,k,i,iElem),ii.EQ.0)
        dy = MERGE( -FV_dx_ETA_L(i,k,j,iElem), FV_dx_ETA_R(i,k,j,iElem),jj.EQ.0)
        UPrim2 = UPrim + gradUxi(:,j,k,i,iElem) * dx + gradUeta(:,i,k,j,iElem) * dy
#if PP_dim == 3
        dz = MERGE(-FV_dx_ZETA_L(i,j,k,iElem),FV_dx_ZETA_R(i,j,k,iElem),kk.EQ.0)
        UPrim2 = UPrim2 + gradUzeta(:,i,j,k,iElem) * dz
#endif
#else
        UPrim2 = UPrim
#endif
        CALL PrimToCons(UPrim2(1:PP_nVarPrim), FV_U_NVisu(:,iii,jjj,kkk,FV_iElem))
      END DO; END DO; END DO
    END DO; END DO; END DO
    FV_U_NVisu(PP_nVar_loc-1,:,:,:,FV_iElem) = IndValue(iElem)
    FV_U_NVisu(PP_nVar_loc,:,:,:,FV_iElem) = FV_Elems(iElem)
  END IF
#endif
END DO !iElem

ALLOCATE(StrVarNames_loc(PP_nVar_loc))
DO iVar=1,PP_nVar
  StrVarNames_loc(iVar) = StrVarNames(iVar)
END DO ! iVar=1,PP_nVar
#if FV_ENABLED
  StrVarNames_loc(PP_nVar_loc-1) = "IndValue"
  StrVarNames_loc(PP_nVar_loc  ) = "FV_Elems"
#endif

! Visualize data
SELECT CASE(OutputFormat)
CASE(OUTPUTFORMAT_TECPLOT)
    CALL CollectiveStop(__STAMP__,'Tecplot output removed due to license issues (possible GPL incompatibility).')
CASE(OUTPUTFORMAT_TECPLOTASCII)
    CALL CollectiveStop(__STAMP__,'Tecplot output removed due to license issues (possible GPL incompatibility).')
CASE(OUTPUTFORMAT_PARAVIEW)
#if FV_ENABLED
  FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DG',OutputTime))
#else
  FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))
#endif
  Coords_NVisu_p => Coords_NVisu
  U_NVisu_p => U_NVisu
  CALL WriteDataToVTK(PP_nVar_loc,NVisu,nElems-nFV_Elems,StrVarNames_loc,Coords_NVisu_p,U_NVisu_p,TRIM(FileString_DG),dim=PP_dim,DGFV=0)
#if FV_ENABLED
  FileString_FV=TRIM(TIMESTAMP(TRIM(ProjectName)//'_FV',OutputTime))
  FV_Coords_NVisu_p => FV_Coords_NVisu
  FV_U_NVisu_p => FV_U_NVisu
  CALL WriteDataToVTK(PP_nVar_loc,NVisu_FV,nFV_Elems,StrVarNames_loc,FV_Coords_NVisu_p,FV_U_NVisu_p,TRIM(FileString_FV),dim=PP_dim,DGFV=1)

  IF (MPIRoot) THEN
    ! write multiblock file
    FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))
    CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_DG,FileString_FV)
  ENDIF
#endif
END SELECT

DEALLOCATE(U_NVisu)
DEALLOCATE(Coords_NVisu)
#if FV_ENABLED
DEALLOCATE(FV_U_NVisu)
DEALLOCATE(FV_Coords_NVisu)
DEALLOCATE(Vdm_GaussN_NVisu_FV)
#endif

END SUBROUTINE Visualize


!==================================================================================================================================
!> Creates or opens file for structured output of data time series in comma seperated value or Tecplot ASCII format.
!> Default is comma seperated value.
!> Searches the file for a dataset at restart time and tries to resume at this point.
!> Otherwise a new file is created.
!==================================================================================================================================
SUBROUTINE InitOutputToFile(Filename,ZoneName,nVar,VarNames,lastLine)
! MODULES
USE MOD_Globals
USE MOD_Output_Vars,  ONLY: ProjectName,ASCIIOutputFormat
USE MOD_Restart_Vars, ONLY: RestartTime
USE MOD_StringTools,  ONLY: INTTOSTR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileName         !< file to be written, without data type extension
CHARACTER(LEN=*),INTENT(IN)   :: ZoneName         !< name of zone (e.g. names of boundary conditions), used for tecplot
INTEGER,INTENT(IN)            :: nVar             !< number of variables
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVar)   !< variable names to be written
REAL,INTENT(OUT),OPTIONAL     :: lastLine(nVar+1) !< last written line to search for, when appending to the file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: stat            !< File IO status
INTEGER                        :: ioUnit          !< IO Unit
INTEGER                        :: i,iMax          !< Counter for header lines
REAL                           :: dummytime       !< Simulation time read from file
LOGICAL                        :: file_exists     !< marker if file exists and is valid
CHARACTER(LEN=255)             :: FileName_loc    !< FileName with data type extension
CHARACTER(LEN=255)             :: tmpStr          !< FileName with data type extension (iterator)
CHARACTER(LEN=4)               :: tmpExt          !< Data type extension
INTEGER                        :: counter         !< Incrementing counter for file backups
!==================================================================================================================================
IF(.NOT.MPIRoot) RETURN

IF(PRESENT(lastLine)) lastLine=-HUGE(1.)

! Append data type extension to FileName
IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
  FileName_loc = TRIM(FileName)//'.csv'
  tmpExt       = '.csv'
ELSE
  FileName_loc = TRIM(FileName)//'.dat'
  tmpExt       = '.dat'
END IF

! Check for file
file_exists = FILEEXISTS(FileName_loc)
IF (RestartTime.LT.0) file_exists = .FALSE.
!! File processing starts here open old and extract information or create new file.
ioUnit = 0

IF(file_exists)THEN ! File exists and append data
  OPEN(NEWUNIT  = ioUnit             , &
       FILE     = TRIM(FileName_loc) , &
       FORM     = 'FORMATTED'        , &
       STATUS   = 'OLD'              , &
       POSITION = 'APPEND'           , &
       RECL     = 50000              , &
       IOSTAT = stat                 )
  IF(stat.NE.0)THEN
    WRITE(UNIT_stdOut,'(A)') ' File '//TRIM(FileName_loc)// ' is invalid. Rewriting file...'
    file_exists=.FALSE.
  END IF
END IF

IF(file_exists)THEN
  ! If we have a restart we need to find the position from where to move on.
  ! Read the values from the previous analyse interval, get the CPUtime
  WRITE(UNIT_stdOut,'(A)')              ' | Opening file '//TRIM(FileName_loc)
  WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' Searching file for time stamp...'

  REWIND(ioUnit)
  ! Loop over header and try to read the first data line. Header size depends on output format.
  iMax =MERGE(2,4,ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV)
  DO i=1,iMax
    READ(ioUnit,*,IOSTAT=stat)
    IF(stat.NE.0)THEN
      ! file is broken, rewrite
      file_exists=.FALSE.
      WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' failed. Writing new file.'
      EXIT
    END IF
  END DO
END IF

IF(file_exists)THEN
  ! Loop until we have found the position
  Dummytime = 0.
  stat=0
  DO WHILE (Dummytime.LT.RestartTime .AND. stat.EQ.0)
    READ(ioUnit,*,IOSTAT=stat) Dummytime
  END DO
  IF(stat.EQ.0)THEN
    ! read final dataset
    IF(PRESENT(lastLine))THEN
      BACKSPACE(ioUnit)
      READ(ioUnit,*,IOSTAT=stat) lastLine
    END IF

    ! Check if there is more than one line left, i.e. lines that will be deleted
    stat    = 0
    BACKSPACE(ioUnit)
    DO WHILE (stat.EQ.0)
      READ(ioUnit,*,IOSTAT=stat) Dummytime
    END DO

    ! Perform a backup of the old file, last line is counted double
    IF (Dummytime.GT.RestartTime) THEN
      counter = 1
      ! Find the first free slot
      DO WHILE(FILEEXISTS(TRIM(FileName)//'.'//TRIM(ADJUSTL(INTTOSTR(counter)))//TRIM(tmpExt)))
        counter = counter + 1
      END DO
      tmpStr  = TRIM(FileName)//'.'//TRIM(ADJUSTL(INTTOSTR(counter)))//TRIM(tmpExt)
      ! Move the file to the free slot
      WRITE(Unit_StdOut,'(A,A,A,A)') ' |> Copying existing file ',TRIM(FileName_loc),' to ',TRIM(tmpStr)
      CALL EXECUTE_COMMAND_LINE('cp '//TRIM(FileName_loc)//' '//TRIM(tmpStr), WAIT=.FALSE.)
    END IF

    ! Rewind back to the beginning of the file
    REWIND(ioUnit)
    ! Loop over header and try to read the first data line. Header size depends on output format.
    iMax = MERGE(2,4,ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV)
    DO i = 1,iMax
      READ(ioUnit,*,IOSTAT=stat)
    END DO

    ! Restore previous position
    Dummytime = 0.
    stat      = 0
    DO WHILE (Dummytime.LT.RestartTime .AND. stat.EQ.0)
      READ(ioUnit,*,IOSTAT=stat) Dummytime
    END DO

    ! Delete from here to end of file
    BACKSPACE(ioUnit)
    ENDFILE(ioUnit)
    WRITE(UNIT_stdOut,'(A,ES15.5)',ADVANCE='YES')' Searching file for time stamp successfull. Resuming file at time ',Dummytime
  ELSE
    WRITE(UNIT_stdOut,'(A)'       ,ADVANCE='YES')' Searching file time for stamp failed. Appending data to end of file.'
  END IF
END IF
CLOSE(ioUnit)

IF(.NOT.file_exists)THEN ! No restart create new file
  OPEN(NEWUNIT= ioUnit             ,&
       FILE   = TRIM(FileName_loc) ,&
       STATUS = 'UNKNOWN'          ,&
       ACCESS = 'SEQUENTIAL'       ,&
       IOSTAT = stat               )
  IF (stat.NE.0) CALL Abort(__STAMP__,'ERROR: cannot open '//TRIM(FileName_loc))

  ! Create a new file with the CSV or Tecplot header
  IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
    WRITE(ioUnit,'(A)',ADVANCE='NO') 'Time'
    DO i=1,nVar
      WRITE(ioUnit,'(A)',ADVANCE='NO') ','//TRIM(VarNames(i))
    END DO
  ELSE
    WRITE(ioUnit,*)'TITLE="'//TRIM(ZoneName)//','//TRIM(ProjectName)//'"'
    WRITE(ioUnit,'(A)',ADVANCE='NO')'VARIABLES = "Time"'
    DO i=1,nVar
      WRITE(ioUnit,'(A)',ADVANCE='NO') ',"'//TRIM(VarNames(i))//'"'
    END DO
    WRITE(ioUnit,'(A)',ADVANCE='YES')
    WRITE(ioUnit,'(A)') 'ZONE T="'//TRIM(ZoneName)//','//TRIM(ProjectName)//'"'
  END IF
  CLOSE(ioUnit) ! outputfile
END IF
END SUBROUTINE InitOutputToFile


!==================================================================================================================================
!> Outputs formatted data into a text file.
!> Format is either comma seperated value or tecplot.
!==================================================================================================================================
SUBROUTINE OutputToFile(FileName,time,nVar,output)
! MODULES
USE MOD_Globals
USE MOD_Output_Vars,  ONLY: ASCIIOutputFormat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileName                 !< name of file to be written, without data type extension
INTEGER,INTENT(IN)            :: nVar(2)                  !< 1: number of variables, 2: number of time samples
REAL,INTENT(IN)               :: time(nVar(2))            !< array of output times
REAL,INTENT(IN)               :: output(nVar(1)*nVar(2))  !< array containing one dataset vector per output time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: openStat                !< File IO status
CHARACTER(LEN=50)              :: formatStr               !< format string for the output and Tecplot header
INTEGER                        :: ioUnit,i
CHARACTER(LEN=255)             :: FileName_loc            ! FileName with data type extension
!==================================================================================================================================
! Append data type extension to FileName
IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
  FileName_loc = TRIM(FileName)//'.csv'
ELSE
  FileName_loc = TRIM(FileName)//'.dat'
END IF

OPEN(NEWUNIT  = ioUnit             , &
     FILE     = TRIM(FileName_loc) , &
     FORM     = 'FORMATTED'        , &
     STATUS   = 'OLD'              , &
     POSITION = 'APPEND'           , &
     RECL     = 50000              , &
     IOSTAT = openStat             )
IF (openStat.NE.0) CALL Abort(__STAMP__,'ERROR: cannot open '//TRIM(FileName_loc))

! Choose between CSV and tecplot output format
IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
  ! Create format string for the variable output: WITH COMMA SEPARATION
  WRITE(formatStr,'(A10,I2,A18)')'(E23.14E5,',nVar(1),'(",",1X,E23.14E5))'
ELSE
  ! Create format string for the variable output: WITH BLANK SEPARATION
  WRITE(formatStr,'(A10,I2,A14)')'(E23.14E5,',nVar(1),'(1X,E23.14E5))'
END IF

DO i=1,nVar(2)
  WRITE(ioUnit,formatstr) time(i),output(nVar(1)*(i-1)+1:nVar(1)*i)
END DO

CLOSE(ioUnit) ! outputfile

END SUBROUTINE OutputToFile


!==================================================================================================================================
!> Deallocate arrays of output routines
!==================================================================================================================================
SUBROUTINE FinalizeOutput()
! MODULES
USE MOD_Output_Vars,ONLY:Vdm_GaussN_NVisu,Vdm_N_NOut,OutputInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(Vdm_GaussN_NVisu)
SDEALLOCATE(Vdm_N_NOut)
OutputInitIsDone = .FALSE.

END SUBROUTINE FinalizeOutput

END MODULE MOD_Output
