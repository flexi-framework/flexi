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
#include "eos.h"

!==================================================================================================================================
!> Contains analyze routines specific to the Navierstokes equations
!==================================================================================================================================
MODULE MOD_AnalyzeEquation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyzeEquation
  MODULE PROCEDURE InitAnalyzeEquation
END INTERFACE

INTERFACE AnalyzeEquation
  MODULE PROCEDURE AnalyzeEquation
END INTERFACE

INTERFACE FinalizeAnalyzeEquation
  MODULE PROCEDURE FinalizeAnalyzeEquation
END INTERFACE


PUBLIC:: AnalyzeEquation, InitAnalyzeEquation, FinalizeAnalyzeEquation
!==================================================================================================================================

PUBLIC::DefineParametersAnalyzeEquation
CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("AnalyzeEquation")
CALL prms%CreateLogicalOption('CalcBodyForces'   , "Set true to compute body forces at walls"         , '.FALSE.')
CALL prms%CreateLogicalOption('CalcBulkState'    , "Set true to compute the flows bulk quantities"    , '.FALSE.')
CALL prms%CreateLogicalOption('CalcMeanFlux'     , "Set true to compute mean flux through boundaries" , '.FALSE.')
CALL prms%CreateLogicalOption('CalcWallVelocity' , "Set true to compute velocities at wall boundaries", '.FALSE.')
CALL prms%CreateLogicalOption('CalcTotalStates'  , "Set true to compute total states (e.g. Tt,pt)"    , '.FALSE.')
CALL prms%CreateLogicalOption('CalcTimeAverage'  , "Set true to compute time averages"                , '.FALSE.')
CALL prms%CreateLogicalOption('WriteBodyForces'  , "Set true to write bodyforces to file"             , '.TRUE.')
CALL prms%CreateLogicalOption('WriteBulkState'   , "Set true to write bulk state to file"             , '.TRUE.')
CALL prms%CreateLogicalOption('WriteMeanFlux'    , "Set true to write mean flux to file"              , '.TRUE.')
CALL prms%CreateLogicalOption('WriteWallVelocity', "Set true to write wall velolcities file"          , '.TRUE.')
CALL prms%CreateLogicalOption('WriteTotalStates' , "Set true to write total states to file"           , '.TRUE.')
CALL prms%CreateStringOption( 'VarNameAvg'       , "Names of variables to be time-averaged"           , multiple=.TRUE.)
CALL prms%CreateStringOption( 'VarNameFluc'      , "Names of variables for which Flucs (time-averaged&
                                                   & square of the variable) should be computed.&
                                                   & Required for computing actual fluctuations."      , multiple=.TRUE.)
END SUBROUTINE DefineParametersAnalyzeEquation


!==================================================================================================================================
!> Initializes variables necessary for NavierStokes specific analyze subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_Equation_Vars,      ONLY: StrVarNamesPrim,StrVarNames
USE MOD_ReadInTools,        ONLY: GETLOGICAL
USE MOD_Mesh_Vars,          ONLY: nBCs,BoundaryType,BoundaryName
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_TimeAverage,        ONLY: InitCalcTimeAverage
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i
!==================================================================================================================================
! Get the various analysis/output variables
doCalcBodyForces    = GETLOGICAL('CalcBodyForces')
doCalcBulkState     = GETLOGICAL('CalcBulkState')
doCalcMeanFlux      = GETLOGICAL('CalcMeanFlux')
doCalcWallVelocity  = GETLOGICAL('CalcWallVelocity')
doCalcTotalStates   = GETLOGICAL('CalcTotalStates')
doWriteBodyForces   = GETLOGICAL('WriteBodyForces')
doWriteBulkState    = GETLOGICAL('WriteBulkState')
doWriteMeanFlux     = GETLOGICAL('WriteMeanFlux')
doWriteWallVelocity = GETLOGICAL('WriteWallVelocity')
doWriteTotalStates  = GETLOGICAL('WriteTotalStates')
doCalcTimeAverage   = GETLOGICAL('CalcTimeAverage')

! Generate wallmap
ALLOCATE(isWall(nBCs))
DO i=1,nBCs
  SELECT CASE(BoundaryType(i,BC_TYPE))
  CASE(3,4,9)
    isWall(i)=.TRUE.
  CASE DEFAULT
    isWall(i)=.FALSE.
  END SELECT
END DO
maxlen=MAX(MAXVAL(LEN_TRIM(BoundaryName))+1,14)

IF(.NOT.ANY(isWall))THEN
  doCalcBodyForces=.FALSE.
  doCalcWallVelocity=.FALSE.
END IF

! Initialize eval routines
IF(MPIRoot)THEN
  IF(doCalcBodyForces.AND.doWriteBodyForces)THEN
    ALLOCATE(Filename_BodyForce(nBCs))
    DO i=1,nBCs
      IF(.NOT.isWall(i)) CYCLE
      FileName_BodyForce(i) = TRIM(ProjectName)//'_BodyForces_'//TRIM(BoundaryName(i))
      CALL InitOutputToFile(FileName_BodyForce(i),TRIM(BoundaryName(i)),9,&
           [CHARACTER(7) :: "x-Force","y-Force","z-Force","Fp_x","Fp_y","Fp_z","Fv_x","Fv_y","Fv_z"])
    END DO
  END IF
  IF(doCalcWallVelocity.AND.doWriteWallVelocity)THEN
    ALLOCATE(Filename_WallVel(nBCs))
    DO i=1,nBCs
      IF(.NOT.isWall(i)) CYCLE
      FileName_WallVel(i) = TRIM(ProjectName)//'_WallVel_'//TRIM(BoundaryName(i))
      CALL InitOutputToFile(FileName_WallVel(i),TRIM(BoundaryName(i)),3,&
           [CHARACTER(7) :: "MeanVel","MinVel","MaxVel"])! gfortran hates mixed length arrays
    END DO
  END IF
  IF(doCalcTotalStates.AND.doWriteTotalStates)THEN
    ALLOCATE(Filename_TotalStates(nBCs))
    DO i=1,nBCs
      IF(BoundaryType(i,BC_TYPE).EQ.1) CYCLE
      FileName_TotalStates(i) = TRIM(ProjectName)//'_TotalStates_'//TRIM(BoundaryName(i))
      CALL InitOutputToFile(FileName_TotalStates(i),TRIM(BoundaryName(i)),4,&
           [CHARACTER(4) :: "pt","p","Tt","Mach"])
    END DO
  END IF
  IF(doCalcBulkState.AND.doWriteBulkState)THEN
    FileName_Bulk  = TRIM(ProjectName)//'_Bulk'
    CALL InitOutputToFile(FileName_Bulk,'Bulk',2*PP_nVar-1,[StrVarNamesPrim,StrVarNames(2:PP_nVar)])
  END IF
  IF(doCalcMeanFlux.AND.doWriteMeanFlux)THEN
    ALLOCATE(Filename_MeanFlux(nBCs))
    DO i=1,nBCs
      IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LE.0)) CYCLE
      FileName_MeanFlux(i) = TRIM(ProjectName)//'_MeanFlux_'//TRIM(BoundaryName(i))
      CALL InitOutputToFile(FileName_MeanFlux(i),TRIM(BoundaryName(i)),PP_nVar,StrVarNames)
    END DO
  END IF
END IF

IF(doCalcTimeAverage)  CALL InitCalcTimeAverage()

END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> Wrapper routine for the equation system specific analyze routines. Will call the specific subroutines to calculate the quantities
!> set in the parameter file and the respective output routines.
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_CalcBodyForces,     ONLY: CalcBodyForces
USE MOD_Mesh_Vars,          ONLY: BoundaryName,nBCs,BoundaryType
USE MOD_Output,             ONLY: OutputToFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                              !< Current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)               :: formatStr
REAL,DIMENSION(3,nBCs)          :: Fv,Fp,BodyForce ! viscous/pressure/resulting body force
REAL,DIMENSION(PP_nVar,nBCs)    :: MeanFlux
REAL,DIMENSION(4,nBCs)          :: meanTotals
REAL,DIMENSION(nBCs)            :: meanV,maxV,minV
REAL                            :: BulkPrim(PP_nVarPrim),BulkCons(PP_nVar)
INTEGER                         :: i
!==================================================================================================================================
! Calculate derived quantities
IF(doCalcBodyforces)   CALL CalcBodyforces(Bodyforce,Fp,Fv)
IF(doCalcWallVelocity) CALL CalcWallVelocity(maxV,minV,meanV)
IF(doCalcMeanFlux)     CALL CalcMeanFlux(MeanFlux)
IF(doCalcBulkState)    CALL CalcBulkState(bulkPrim,bulkCons)
IF(doCalcTotalStates)  CALL CalcKessel(meanTotals)

IF(MPIRoot)THEN
  IF(doCalcBodyforces)THEN
    WRITE(UNIT_stdOut,*)'BodyForces (Pressure, Friction) : '
    WRITE(formatStr,'(A,I2,A)')'(A',maxlen,',6ES18.9)'
    DO i=1,nBCs
      IF(.NOT.isWall(i)) CYCLE
      IF (doWriteBodyForces) &
        CALL OutputToFile(FileName_BodyForce(i),(/Time/),(/9,1/),(/BodyForce(:,i),Fp(:,i),Fv(:,i)/))
      WRITE(UNIT_stdOut,formatStr) ' '//TRIM(BoundaryName(i)),Fp(:,i),Fv(:,i)
    END DO
  END IF

  IF(doCalcWallVelocity)THEN
    WRITE(UNIT_stdOut,*)'Wall Velocities (mean/min/max)  : '
    WRITE(formatStr,'(A,I2,A)')'(A',maxlen,',3ES18.9)'
    DO i=1,nBCs
      IF(.NOT.isWall(i)) CYCLE
      IF (doWriteWallVelocity) &
        CALL OutputToFile(FileName_WallVel(i),(/Time/),(/3,1/),(/meanV(i),minV(i),maxV(i)/))
      WRITE(UNIT_stdOut,formatStr) ' '//TRIM(BoundaryName(i)),meanV(i),minV(i),maxV(i)
    END DO
  END IF

  IF(doCalcMeanFlux)THEN
    WRITE(formatStr,'(A,I2,A,I2,A)')'(A',maxlen,',',PP_nVar,'ES18.9)'
    WRITE(UNIT_stdOut,*)'MeanFlux through boundaries     : '
    DO i=1,nBCs
      IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LE.0)) CYCLE
      IF (doWriteMeanFlux) &
        CALL OutputToFile(FileName_MeanFlux(i),(/Time/),(/PP_nVar,1/),MeanFlux(:,i))
      WRITE(UNIT_stdOut,formatStr) ' '//TRIM(BoundaryName(i)),MeanFlux(:,i)
    END DO
  END IF  !(doCalcBodyforces)

  IF(doCalcBulkState)THEN
    IF (doWriteBulkState) &
      CALL OutputToFile(FileName_Bulk,(/Time/),(/PP_nVarPrim+PP_nVar-1,1/),(/BulkPrim,BulkCons(2:PP_nVar)/))
    WRITE(formatStr,'(A,I2,A)')'(A14,',PP_nVarPrim,'ES18.9)'
    WRITE(UNIT_stdOut,formatStr)' Bulk Prims : ',bulkPrim
    WRITE(formatStr,'(A,I2,A)')'(A14,',PP_nVar,'ES18.9)'
    WRITE(UNIT_stdOut,formatStr)' Bulk Cons  : ',bulkCons
  END IF

  IF(doCalcTotalStates)THEN
    WRITE(UNIT_stdOut,*)'Mean total states at boundaries : '
    WRITE(formatStr,'(A,I2,A)')'(A',maxlen,',4ES18.9)'
    DO i=1,nBCs
      IF(BoundaryType(i,BC_TYPE).EQ.1) CYCLE
      IF (doWriteTotalStates) &
        CALL OutputToFile(FileName_TotalStates(i),(/Time/),(/4,1/),meanTotals(:,i) )
      WRITE(UNIT_stdOut,formatStr) ' '//TRIM(BoundaryName(i)),MeanTotals(:,i)
    END DO
  END IF
END IF ! MPIRoot

END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Calculates bulk quantities over whole domain
!==================================================================================================================================
SUBROUTINE CalcBulkState(BulkPrim,BulkCons)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U,UPrim
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_w
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)                :: BulkPrim(PP_nVarPrim)                   !< Primitive bulk quantities
REAL,INTENT(OUT)                :: BulkCons(PP_nVar)                       !< Conservative bulk quantities
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: IntegrationWeight
INTEGER                         :: iElem,i,j,k
#if USE_MPI
REAL                            :: box(PP_nVar+PP_nVarPrim)
#endif
!==================================================================================================================================
BulkPrim=0.
BulkCons=0.
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) THEN ! FV Element
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      IntegrationWeight=FV_w(i)*FV_w(j)*FV_w(k)/sJ(i,j,k,iElem,1)
      BulkCons         =BulkCons+U(:,i,j,k,iElem)*IntegrationWeight
      BulkPrim         =BulkPrim+UPrim(:,i,j,k,iElem)*IntegrationWeight
    END DO; END DO; END DO !i,j,k
  ELSE ! DG element
#endif
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
      BulkCons         =BulkCons+U(:,i,j,k,iElem)*IntegrationWeight
      BulkPrim         =BulkPrim+UPrim(:,i,j,k,iElem)*IntegrationWeight
    END DO; END DO; END DO !i,j,k
#if FV_ENABLED
  END IF
#endif
END DO ! iElem

#if USE_MPI
Box(1:PP_nVarPrim)=BulkPrim; Box(PP_nVarPrim+1:PP_nVarPrim+PP_nVar)=BulkCons
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,box,PP_nVar+PP_nVarPrim,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  BulkPrim=Box(1:PP_nVarPrim); BulkCons=Box(PP_nVarPrim+1:PP_nVarPrim+PP_nVar)
ELSE
  CALL MPI_REDUCE(Box         ,0  ,PP_nVar+PP_nVarPrim,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif

BulkPrim=BulkPrim/Vol
BulkCons=BulkCons/Vol

END SUBROUTINE CalcBulkState


!===================================================================================================================================
!> Computes total quantities
!===================================================================================================================================
SUBROUTINE CalcKessel(meanTotals)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars,           ONLY: U_master
USE MOD_Mesh_Vars,         ONLY: SurfElem
USE MOD_Mesh_Vars,         ONLY: nBCSides,BC,BoundaryType,nBCs
USE MOD_Analyze_Vars,      ONLY: wGPSurf,Surf
USE MOD_EOS_Vars,          ONLY: Kappa,R,sKappaM1,KappaM1
#if FV_ENABLED
USE MOD_FV_Vars,           ONLY: FV_Elems_master,FV_w
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(OUT)             :: meanTotals(4,nBCs)           !< total and static pressure pt,p
                                                             !< total temperature Tt and Mach
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: dA,c,mach
REAL                           :: primvar(1:14),UE(PP_2Var)
INTEGER                        :: SideID,i,j,iBC
!===================================================================================================================================
meanTotals= 0.
DO SideID=1,nBCSides
  iBC=BC(SideID)
  IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE

  DO j=0,PP_NZ; DO i=0,PP_N
    ! TODO: ATTENTION: Temperature of UE not filled!!!
    UE(EXT_CONS)=U_master(CONS,i,j,SideID)
    UE(EXT_SRHO)=1./UE(EXT_DENS)
    UE(EXT_VELV)=VELOCITY_HE(UE)
    UE(EXT_PRES)=PRESSURE_HE(UE)

    PrimVar(1:3)=UE(EXT_VELV)

    ! VelocityMagnitude
    PrimVar(4)=SQRT(SUM(PrimVar(1:3)*PrimVar(1:3)))
    ! Pressure
    PrimVar(5)=UE(EXT_PRES)
    ! VelocitySound
    PrimVar(6)=SPEEDOFSOUND_HE(UE)
    ! Mach
    PrimVar(7)=PrimVar(4)/PrimVar(6)
    ! Temperature
    PrimVar(8)=TEMPERATURE_HE(UE)
    ! EnergyStagnation
    PrimVar(9)=UE(EXT_ENER)*UE(EXT_SRHO)
    ! EnthalpyStagnation
    PrimVar(10)=PrimVar(9)+PrimVar(5)*UE(EXT_SRHO)
    ! Entropy
    PrimVar(11)= ENTROPY_H(UE,PrimVar(8))
    ! Potential Temperature
!    PrimVar(12)=PrimVar(8)/(PrimVar(5)/P0)**(1.-sKappa)
    c=PrimVar(6)
    Mach=PrimVar(7)
    ! Total Temperature
    PrimVar(13)=TOTAL_TEMPERATURE_H(PrimVar(8),Mach)
    ! Total Pressure
    PrimVar(14)=TOTAL_PRESSURE_H(PrimVar(5),Mach)
    ! Calculate velocity magnitude
 !   locV=SQRT(Vel(1)*Vel(1)+Vel(2)*Vel(2)+Vel(3)*Vel(3))
!    maxV(iBC)=MAX(maxV(iBC),locV)
 !   minV(iBC)=MIN(minV(iBC),locV)
#if FV_ENABLED
    IF (FV_Elems_master(SideID).EQ.1) THEN ! FV element
      dA=FV_w(i)*FV_w(j)*SurfElem(i,j,1,SideID)
    ELSE
#endif
      dA=wGPSurf(i,j)*SurfElem(i,j,0,SideID)
#if FV_ENABLED
    END IF
#endif
    meanTotals(1,iBC)=meanTotals(1,iBC)+Primvar(14)*dA ! pt
    meanTotals(2,iBC)=meanTotals(2,iBC)+Primvar(5)*dA  ! p
    meanTotals(3,iBC)=meanTotals(3,iBC)+Primvar(13)*dA ! Tt
    meanTotals(4,iBC)=meanTotals(4,iBC)+mach*dA        ! Ma
  END DO; END DO
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,meanTotals,4*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(meanTotals  ,0         ,4*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif

IF(.NOT.MPIRoot) RETURN

DO iBC=1,nBCs
  IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE
  MeanTotals(:,iBC)=MeanTotals(:,iBC)/Surf(iBC)
END DO
END SUBROUTINE CalcKessel


!==================================================================================================================================
!> Calculate velocity at walls
!==================================================================================================================================
SUBROUTINE CalcWallVelocity(maxV,minV,meanV)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY: UPrim_master
USE MOD_Mesh_Vars,            ONLY: SurfElem
USE MOD_Mesh_Vars,            ONLY: nBCSides,BC,nBCs
USE MOD_Analyze_Vars,         ONLY: wGPSurf,Surf
USE MOD_AnalyzeEquation_Vars, ONLY: isWall
#if FV_ENABLED
USE MOD_FV_Vars,              ONLY: FV_Elems_master,FV_w
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)               :: maxV(nBCs)          !< Maximum of wall velocity per boundary
REAL,INTENT(OUT)               :: minV(nBCs)          !< Minimum of wall velocity per boundary
REAL,INTENT(OUT)               :: meanV(nBCs)         !< Mean of wall velocity per boundary
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: dA,Vel(3),locV
INTEGER                        :: iSide,i,j,iBC
!==================================================================================================================================
minV =  1.e14
maxV = -1.e14
meanV= 0.
DO iSide=1,nBCSides
  iBC=BC(iSide)
  IF(.NOT.isWall(iBC)) CYCLE
  DO j=0,PP_NZ; DO i=0,PP_N
    Vel=UPrim_master(2:4,i,j,iSide)
    ! Calculate velocity magnitude
    locV=SQRT(DOT_PRODUCT(vel,vel))
    maxV(iBC)=MAX(maxV(iBC),locV)
    minV(iBC)=MIN(minV(iBC),locV)
#if FV_ENABLED
    IF (FV_Elems_master(iSide).EQ.1) THEN ! FV element
      dA=FV_w(i)*FV_w(j)*SurfElem(i,j,1,iSide)
    ELSE
#endif
      dA=wGPSurf(i,j)*SurfElem(i,j,0,iSide)
#if FV_ENABLED
    END IF
#endif
    meanV(iBC)=meanV(iBC)+locV*dA
  END DO; END DO
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,maxV ,nBCs,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,minV ,nBCs,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,meanV,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(maxV        ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(minV        ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(meanV       ,0    ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif
DO iBC=1,nBCs
  IF(.NOT.isWall(iBC)) CYCLE
  MeanV(iBC)=MeanV(iBC)/Surf(iBC)
END DO

END SUBROUTINE CalcWallVelocity


!==================================================================================================================================
!> Calculate the mean fluxes on the boundaries
!==================================================================================================================================
SUBROUTINE CalcMeanFlux(MeanFlux)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars,           ONLY: Flux_master
USE MOD_Analyze_Vars,      ONLY: wGPSurf,Surf
USE MOD_Mesh_Vars,         ONLY: nSides,nMPISides_YOUR,AnalyzeSide,nBCs,BoundaryType
#if FV_ENABLED
USE MOD_FV_Vars,           ONLY: FV_Elems_master,FV_w
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)               :: MeanFlux(PP_nVar,nBCs)        !< Mean flux in each conservative variable per boundary
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSide,iSurf,i,j
!==================================================================================================================================
MeanFlux=0.
DO iSide=1,nSides-nMPISides_YOUR
  iSurf=AnalyzeSide(iSide)
  IF(iSurf.EQ.0) CYCLE
#if FV_ENABLED
  IF (FV_Elems_master(iSide).EQ.1) THEN ! FV element
    DO j=0,PP_NZ; DO i=0,PP_N
      ! Don't multiply with Surfelem, its already contained in the fluxes
      MeanFlux(:,iSurf)=MeanFlux(:,iSurf)+Flux_master(:,i,j,iSide)*FV_w(i)*FV_w(j)
    END DO; END DO
  ELSE ! DG element
#endif
    DO j=0,PP_NZ; DO i=0,PP_N
      ! Don't multiply with Surfelem, its already contained in the fluxes
      MeanFlux(:,iSurf)=MeanFlux(:,iSurf)+Flux_master(:,i,j,iSide)*wGPSurf(i,j)
    END DO; END DO
#if FV_ENABLED
  END IF
#endif
END DO

#if USE_MPI
i=PP_nVar*nBCs
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,MeanFlux,i,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(MeanFlux    ,0       ,i,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif

DO i=1,nBCs
  IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LT.0)) CYCLE
  MeanFlux(:,i)=MeanFlux(:,i)/Surf(i)
END DO

END SUBROUTINE CalcMeanFlux


!==================================================================================================================================
!> Finalizes variables necessary for analyze subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
USE MOD_AnalyzeEquation_Vars
USE MOD_TimeAverage,        ONLY: FinalizeTimeAverage
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(isWall)
SDEALLOCATE(FileName_BodyForce)
SDEALLOCATE(FileName_WallVel)
SDEALLOCATE(FileName_MeanFlux)
SDEALLOCATE(FileName_TotalStates)

IF (doCalcTimeAverage) CALL FinalizeTimeAverage()
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
