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
#include "eos.h"

#if FV_ENABLED
#error "This testcase is not tested with FV"
#endif

!==================================================================================================================================
!> Ercoftac periodic hill testcase
!> http://www.kbwiki.ercoftac.org/w/index.php/Abstr:2D_Periodic_Hill_Flow
!==================================================================================================================================
MODULE MOD_TestCase
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE DefineParametersTestcase
  MODULE PROCEDURE DefineParametersTestcase
End INTERFACE

INTERFACE InitTestcase
  MODULE PROCEDURE InitTestcase
END INTERFACE

INTERFACE FinalizeTestcase
  MODULE PROCEDURE FinalizeTestcase
END INTERFACE

INTERFACE ExactFuncTestcase
  MODULE PROCEDURE ExactFuncTestcase
END INTERFACE

INTERFACE CalcForcing
  MODULE PROCEDURE CalcForcing
END INTERFACE

INTERFACE TestcaseSource
  MODULE PROCEDURE TestcaseSource
END INTERFACE

INTERFACE AnalyzeTestcase
  MODULE PROCEDURE DO_NOTHING_LOG
END INTERFACE

INTERFACE GetBoundaryFluxTestcase
  MODULE PROCEDURE GetBoundaryFluxTestcase
END INTERFACE

INTERFACE GetBoundaryFVgradientTestcase
  MODULE PROCEDURE GetBoundaryFVgradientTestcase
END INTERFACE

INTERFACE Lifting_GetBoundaryFluxTestcase
  MODULE PROCEDURE Lifting_GetBoundaryFluxTestcase
END INTERFACE

PUBLIC:: DefineParametersTestcase
PUBLIC:: InitTestcase
PUBLIC:: FinalizeTestcase
PUBLIC:: ExactFuncTestcase
PUBLIC:: TestcaseSource
PUBLIC:: CalcForcing
PUBLIC:: AnalyzeTestcase
PUBLIC:: GetBoundaryFluxTestcase
PUBLIC:: GetBoundaryFVgradientTestcase
PUBLIC:: Lifting_GetBoundaryFluxTestcase

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersTestcase()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Testcase")
CALL prms%CreateRealOption('massFlowRef',      "Prescribed massflow for testcase."                                 , '1.')
CALL prms%CreateRealOption('Forcing_MaxMemory',"Maximum amount of memory to be used to buffer testcase forcing log data. "//&
                                               "If memory is exceeded before analyze level, log files are written.", '100.')
CALL prms%CreateStringOption('massFlowBCName', "Name of BC at which massflow is computed."                         , 'INFLOW')
END SUBROUTINE DefineParametersTestcase


!==================================================================================================================================
!> Initialize testcase specific variables
!==================================================================================================================================
SUBROUTINE InitTestcase()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TestCase_Vars
USE MOD_ReadInTools,    ONLY: GETREAL,GETSTR,GETINT
USE MOD_Mesh_Vars,      ONLY: nBCs,BoundaryName
USE MOD_Output_Vars,    ONLY: ProjectName
USE MOD_Output,         ONLY: InitOutputToFile
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)       :: massFlowBCName
INTEGER                  :: ioUnit,openStat,i
REAL                     :: maxMemory
CHARACTER(LEN=20)        :: varnames(4)
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE PERIODIC HILL...'

#if FV_ENABLED
CALL CollectiveStop(__STAMP__, &
  'The testcase has not been implemented for FV yet!')
#endif

massFlowRef    = GETREAL('massFlowRef')
massFlowBCName = GETSTR( 'massFlowBCName')
maxMemory      = GETREAL('Forcing_MaxMemory')             ! Max buffer (100MB)
MaxBuffer      = MaxMemory*131072/5                       != size in bytes/nLogVars
massFlowPrev   = massFlowRef

massFlowBC=-1
DO i=1,nBCs
  IF(TRIM(BoundaryName(i)).EQ.TRIM(massFlowBCName)) massFlowBC=i
END DO
IF(massFlowBC.EQ.-1) CALL Abort(__STAMP__,'No inflow BC found.')

IF(.NOT.MPIRoot) RETURN

Filename = TRIM(ProjectName)//'_Stats'
varnames(1) = 'dpdx'
varnames(2) = 'bulkVel'
varnames(3) = 'massFlowRateGlobal'
varnames(4) = 'massFlowRatePeriodic'
CALL InitOutputToFile(Filename,'Statistics',4,varnames)

SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE PERIODIC HILL DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitTestcase



!==================================================================================================================================
!> Specifies all the initial conditions for the periodic hill testcase.
!==================================================================================================================================
SUBROUTINE ExactFuncTestcase(tIn,x,Resu,Resu_t,Resu_tt)
! MODULES
USE MOD_Globals,      ONLY: Abort
USE MOD_Equation_Vars,ONLY: RefStatePrim,IniRefState
USE MOD_EOS,          ONLY: PrimToCons
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: x(3)        !< position in physical coordinates
REAL,INTENT(IN)                 :: tIn         !< current simulation time
REAL,INTENT(OUT)                :: Resu(5)     !< exact fuction evaluated at tIn, returning state in conservative variables
REAL,INTENT(OUT)                :: Resu_t(5)   !< first time deriv of exact fuction
REAL,INTENT(OUT)                :: Resu_tt(5)  !< second time deriv of exact fuction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: xloc,h,prim(PP_nVarPrim)
!==================================================================================================================================

! periodic hill initial condition (geometry scaled down to 1 by a factor of 28)
xloc= x(1)*28.
IF(xloc.GT.54) xloc=28.*9.-xloc
IF((xloc.GT.-1.e-10).AND.(xloc.LE.9.))THEN
!Between x=0. and x=9.
h=min(28.,&
      2.800000000000E+01          +0.000000000000E+00*xloc     &
     +6.775070969851E-03*xloc**2  -2.124527775800E-03*xloc**3)
ELSEIF((xloc.GT.9.).AND.(xloc.LE.14.))THEN
!Between x=9. and x=14.
h=  2.507355893131E+01          +9.754803562315E-01*xloc       &
   -1.016116352781E-01*xloc**2  +1.889794677828E-03*xloc**3

ELSEIF((xloc.GT.14.).AND.(xloc.LE.20.))THEN
!Between x=14. and x=20.
h=  2.579601052357E+01          +8.206693007457E-01*xloc       &
   -9.055370274339E-02*xloc**2  +1.626510569859E-03*xloc**3

ELSEIF((xloc.GT.20.).AND.(xloc.LE.30.))THEN
!Between x=20. and x=30.
h=  4.046435022819E+01          -1.379581654948E+00*xloc       &
   +1.945884504128E-02*xloc**2  -2.070318932190E-04*xloc**3

ELSEIF((xloc.GT.30.).AND.(xloc.LE.40.))THEN
!Between x=30. and x=40.
h=  1.792461334664E+01          +8.743920332081E-01*xloc       &
   -5.567361123058E-02*xloc**2  +6.277731764683E-04*xloc**3

ELSEIF((xloc.GT.40.).AND.(xloc.LE.54.))THEN
!Between x=40. and x=54.
h=max(0.,&
      5.639011190988E+01          -2.010520359035E+00*xloc     &
     +1.644919857549E-02*xloc**2  +2.674976141766E-05*xloc**3)
ELSEIF(xloc.GT.54.)THEN
h=  0.
ELSE
  CALL Abort(__STAMP__,&
             'Wrong hill geometry')
END IF
h=h/28.
Prim    = RefStatePrim(:,IniRefState)
Prim(2) = Prim(2)*2.025/(3.025-h)
Prim(6) = 0. ! T does not matter for prim to cons
CALL PrimToCons(Prim,Resu)
Resu_t =0.
Resu_tt=0.
END SUBROUTINE ExactFuncTestcase



!==================================================================================================================================
!> Specifies the forcing term for the periodic hill testcase
!==================================================================================================================================
SUBROUTINE CalcForcing(t,dt)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TestCase_Vars
USE MOD_DG_Vars,        ONLY: U,U_master,U_slave
USE MOD_Analyze_Vars,   ONLY: wGPSurf,wGPVol,Surf,Vol,writeData_dt,tWriteData
USE MOD_Mesh_Vars,      ONLY: SurfElem,AnalyzeSide,nSides,nMPISides_YOUR,sJ
USE MOD_Mesh_Vars,      ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< current solution time
REAL,INTENT(IN)                 :: dt                     !< current time step
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem,SideID
REAL                            :: massFlow,massFlowGlobal,massFlowPeriodic,tmp
#if USE_MPI
REAL                            :: box(3)
#endif
!==================================================================================================================================
massFlowGlobal=0.
massFlowPeriodic=0.
BulkVel =0.
! Periodic hill testcase
! 1. Get Massflux at the periodic inlet (hill crest)
DO SideID=1,nSides-nMPISides_YOUR
  IF(AnalyzeSide(SideID).NE.massFlowBC) CYCLE
  DO j=0,PP_N; DO i=0,PP_N
    tmp     =0.5*wGPSurf(i,j)*SurfElem(i,j,0,SideID)
    massFlowPeriodic=massFlowPeriodic+(U_master(2,i,j,SideID)+U_slave(2,i,j,SideID))*tmp
  END DO; END DO
END DO

DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    massFlowGlobal = massFlowGlobal+U(MOM1,i,j,k,iElem)*wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
    BulkVel = BulkVel+U(MOM1,i,j,k,iElem)/U(DENS,i,j,k,iElem)*wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
END DO

#if USE_MPI
box(1) = massFlowGlobal; box(2) = massFlowPeriodic; box(3) = BulkVel
CALL MPI_ALLREDUCE(MPI_IN_PLACE,box,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
massFlowGlobal = box(1); massFlowPeriodic = box(2); BulkVel = box(3)
#endif

BulkVel = BulkVel/Vol
massFlowGlobal = massFlowGlobal/9. ! averaged mass flow over channel length
IF(firstTimestep)THEN
  dtPrev=dt
  massFlowPrev = massFlowRef
  firstTimeStep=.FALSE.
END IF

massFlow = massFlowGlobal ! use global

! 2. compute forcing term dp/dx
! dpdx_n+1 = dpdx_n - (mRef - 2*m_n + m_n-1) / (surf*dt)
dpdx = dpdxPrev - 0.3*(massFlowRef - 2*massFlow + massFlowPrev) / (Surf(massFlowBC)*dt)

!! proposed by Lenormand:
!!dpdx = dpdxPrev - (alpha*(massFlow-massFlowRef) + beta*(massFlowPrev - massFlowRef)) / (Surf(massFlowBC))
!alpha=2.
!beta =-0.2
!massFlowPredictor = massFlow+dt/dtPrev*(massFlow-massFlowPrev)
!dpdx = dpdxPrev + (alpha*(massFlowPredictor-massFlowRef) + beta*(massFlow - massFlowRef)) / (Surf(massFlowBC))

massFlowPrev = massFlow
dpdxPrev     = dpdx
dtPrev       = dt
IF(MPIRoot)THEN
  IF(.NOT.ALLOCATED(writeBuf))THEN
    Buffer = MIN(CEILING((2.0*WriteData_dt)/dt),MaxBuffer) ! use 2x security
    ALLOCATE(writeBuf(5,buffer))
    SWRITE(*,*) 'Buffer for massflow logging is ', buffer
  END IF
  ioCounter=ioCounter+1
  writeBuf(:,ioCounter) = (/t, dpdx, BulkVel, MassFlowGlobal, MassFlowPeriodic /)
  IF((ioCounter.GE.buffer).OR.((tWriteData-t-dt).LT.1e-10))THEN
    CALL WriteStats()
  END IF
END IF

END SUBROUTINE CalcForcing


!==================================================================================================================================
!> Apply forcing term for periodic hill testcase.
!> - Streamwise pressure gradient computed by predefined massflow
!==================================================================================================================================
SUBROUTINE TestcaseSource(Ut)
! MODULES
USE MOD_PreProc
USE MOD_TestCase_Vars,ONLY:dpdx,BulkVel
USE MOD_Mesh_Vars,    ONLY:sJ,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================
! Periodic hill testcase
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    Ut(2,i,j,k,iElem)=Ut(2,i,j,k,iElem)-dpdx/sJ(i,j,k,iElem,0)
    Ut(5,i,j,k,iElem)=Ut(5,i,j,k,iElem)-dpdx*BulkVel/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
END DO
END SUBROUTINE TestcaseSource


!==================================================================================================================================
!> Writes statistics of testcase to file
!==================================================================================================================================
SUBROUTINE WriteStats()
! MODULES
USE MOD_Globals      ,ONLY:Abort
USE MOD_TestCase_Vars,ONLY:writeBuf,FileName,ioCounter
USE MOD_Output,       ONLY:OutputToFile
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: ioUnit,openStat,i
!==================================================================================================================================
CALL OutputToFile(FileName,writeBuf(1,1:ioCounter),(/4,ioCounter/),RESHAPE(writeBuf(2:5,1:ioCounter),(/4*ioCounter/)))
ioCounter=0

END SUBROUTINE WriteStats


!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE FinalizeTestcase()
! MODULES
USE MOD_Globals      ,ONLY:MPIRoot
USE MOD_TestCase_Vars,ONLY:writeBuf
IMPLICIT NONE
!==================================================================================================================================
IF(MPIRoot) DEALLOCATE(writeBuf)
END SUBROUTINE


! SUBROUTINE DO_NOTHING(optionalREAL,optionalREAL2)
! IMPLICIT NONE
! REAL,OPTIONAL,INTENT(IN)    :: optionalREAL,optionalREAL2
! END SUBROUTINE DO_NOTHING


SUBROUTINE DO_NOTHING_LOG(optionalREAL,optionalLOG)
IMPLICIT NONE
REAL,OPTIONAL,INTENT(IN)    :: optionalREAL
LOGICAL,OPTIONAL,INTENT(IN) :: optionalLOG
END SUBROUTINE DO_NOTHING_LOG


SUBROUTINE GetBoundaryFluxTestcase(SideID,t,Nloc,Flux,UPrim_master,                   &
#if PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: SideID  !< ID of current side
REAL,INTENT(IN)      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)   :: Nloc    !< polynomial degree
REAL,INTENT(IN)      :: UPrim_master( PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))    !< inner surface solution
#if PARABOLIC
REAL,INTENT(IN)      :: gradUx_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !> inner surface solution gradients in x-direction
REAL,INTENT(IN)      :: gradUy_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !> inner surface solution gradients in y-direction
REAL,INTENT(IN)      :: gradUz_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !> inner surface solution gradients in z-direction
#endif /*PARABOLIC*/
REAL,INTENT(IN)      :: NormVec (  3,0:Nloc,0:ZDIM(Nloc))  !< normal vectors on surfaces
REAL,INTENT(IN)      :: TangVec1(  3,0:Nloc,0:ZDIM(Nloc))  !< tangential1 vectors on surfaces
REAL,INTENT(IN)      :: TangVec2(  3,0:Nloc,0:ZDIM(Nloc))  !< tangential2 vectors on surfaces
REAL,INTENT(IN)      :: Face_xGP(  3,0:Nloc,0:ZDIM(Nloc))  !< positions of surface flux points
REAL,INTENT(OUT)     :: Flux(PP_nVar,0:Nloc,0:ZDIM(Nloc))  !< resulting boundary fluxes
!==================================================================================================================================
END SUBROUTINE GetBoundaryFluxTestcase


SUBROUTINE GetBoundaryFVgradientTestcase(SideID,t,gradU,UPrim_master)
USE MOD_PreProc
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                                   !< ID of current side
REAL,INTENT(IN)    :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)   :: gradU       (PP_nVarPrim,0:PP_N,0:PP_NZ) !< FV boundary gradient
!==================================================================================================================================
END SUBROUTINE GetBoundaryFVgradientTestcase


SUBROUTINE Lifting_GetBoundaryFluxTestcase(SideID,t,UPrim_master,Flux)
USE MOD_PreProc
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                                   !< ID of current side
REAL,INTENT(IN)    :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)   :: Flux(     PP_nVarLifting,0:PP_N,0:PP_NZ) !< lifting boundary flux
!==================================================================================================================================
END SUBROUTINE Lifting_GetBoundaryFluxTestcase

END MODULE MOD_TestCase
