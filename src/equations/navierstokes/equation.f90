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

!==================================================================================================================================
!> Soubroutines necessary for calculating Navier-Stokes equations
!==================================================================================================================================
MODULE MOD_Equation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE

INTERFACE GetPrimitiveStateSurface
  MODULE PROCEDURE GetPrimitiveStateSurface
END INTERFACE

INTERFACE GetConservativeStateSurface
  MODULE PROCEDURE GetConservativeStateSurface
END INTERFACE

INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC:: DefineParametersEquation,InitEquation,FinalizeEquation
PUBLIC:: GetPrimitiveStateSurface,GetConservativeStateSurface
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_ReadInTools,ONLY: prms
USE MOD_Riemann    ,ONLY: DefineParametersRiemann
#ifdef SPLIT_DG
USE MOD_SplitFlux  ,ONLY: DefineParametersSplitDG
#endif /*SPLIT_DG*/
#if EDDYVISCOSITY
USE MOD_EddyVisc,   ONLY: DefineParametersEddyVisc
#endif
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(      'IniRefState',  "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('RefState',     "State(s) in primitive variables (density, velx, vely, velz, pressure).",&
                                                multiple=.TRUE.)
CALL prms%CreateStringOption(   'BCStateFile',  "File containing the reference solution on the boundary to be used as BC.")

CALL DefineParametersRiemann()
#ifdef SPLIT_DG
CALL DefineParametersSplitDG()
#endif /*SPLIT_DG*/
#if EDDYVISCOSITY
CALL DefineParametersEddyVisc()
#endif /*EDDYVISCOSITY*/
END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Set parameters needed by equation modules and initialize equations as well as boundary conditions and testcases
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars
USE MOD_Eos               ,ONLY: InitEos,PrimToCons
USE MOD_EOS_Vars          ,ONLY: R
USE MOD_Exactfunc         ,ONLY: InitExactFunc
USE MOD_ReadInTools       ,ONLY: CountOption,GETREALARRAY,GETSTR
USE MOD_TestCase          ,ONLY: InitTestcase
USE MOD_Riemann           ,ONLY: InitRiemann
USE MOD_GetBoundaryFlux,   ONLY: InitBC
USE MOD_CalcTimeStep      ,ONLY: InitCalctimestep
#if EDDYVISCOSITY
USE MOD_EddyVisc          ,ONLY: InitEddyVisc
#endif
#ifdef SPLIT_DG
USE MOD_SplitFlux         ,ONLY: InitSplitDG
#endif /*SPLIT_DG*/
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
REAL    :: UE(PP_2Var)
!==================================================================================================================================
IF(EquationInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitEquation not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT NAVIER-STOKES...'

s43=4./3.
s23=2./3.

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

! Read in boundary parameters
IniRefState  = 0

CALL InitExactFunc()
CALL InitEOS()

! Read Boundary information / RefStates / perform sanity check
nRefState=CountOption('RefState')
IF(IniRefState.GT.nRefState)THEN
  CALL CollectiveStop(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState):',IniRefState,REAL(nRefState))
END IF

IF(nRefState .GT. 0)THEN
  ALLOCATE(RefStatePrim(PP_nVarPrim,nRefState))
  ALLOCATE(RefStateCons(PP_nVar    ,nRefState))
  DO i=1,nRefState
    RefStatePrim(1:5,i)  = GETREALARRAY('RefState',5)
#if PP_dim==2
    IF(RefStatePrim(VEL3,i).NE.0.) THEN
      SWRITE(UNIT_stdOut,'(A)')' You are computing in 2D! RefStatePrim(4) will be set to zero!'
      RefStatePrim(VEL3,i)=0.
    END IF
#endif
    ! TODO: ATTENTION only sRho and Pressure of UE filled!!!
    UE(EXT_SRHO) = 1./RefStatePrim(DENS,i)
    UE(EXT_PRES) = RefStatePrim(PRES,i)
    RefStatePrim(TEMP,i) = TEMPERATURE_HE(UE)
    CALL PrimToCons(RefStatePrim(:,i),RefStateCons(:,i))
  END DO
END IF

! boundary state filename if present
BCStateFile=GETSTR('BCStateFile','nonexistingfile')

! Initialize Riemann solvers to be in volume and on BCs
CALL InitRiemann()

! Initialize timestep calculation
CALL InitCalctimestep()

#if EDDYVISCOSITY
! Initialize eddyViscosity
CALL InitEddyVisc()
#endif

#ifdef SPLIT_DG
! Initialize SplitDG
CALL InitSplitDG()
#endif /*SPLIT_DG*/
CALL InitBC()

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT NAVIER-STOKES DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

! Initialize current testcase
CALL InitTestcase()

END SUBROUTINE InitEquation


!==================================================================================================================================
!> Converts conservative solution vector to primitive variables
!>
!> Two possibilities for sides if using non-Lobatto node sets:
!> 1. Convert U_master/slave to prims (used):
!>    prims consistent to cons, but inconsistent to prim volume
!>    cheap and simple, no communication and mortars required.
!>    Using this version the primitive solution is no longer a polynomial.
!> 2. Compute UPrim_master/slave from volume UPrim
!>    UPrim_master/slave consistent to UPrim, but inconsistent to U_master/slave
!>    more expensive, communication and mortars required.
!>    This version gives thermodynamically inconsistant states at sides.
!>
!> TODO: Provide switch for these two versions.
!==================================================================================================================================
SUBROUTINE GetPrimitiveStateSurface(U_master,U_slave,UPrim_master,UPrim_slave)
! MODULES
USE MOD_Preproc
USE MOD_EOS,      ONLY: ConsToPrim
USE MOD_Mesh_Vars,ONLY: firstInnerSide,firstMPISide_YOUR,lastMPISide_YOUR,nSides
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_master(    CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(IN)  :: U_slave(     CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
REAL,INTENT(OUT) :: UPrim_master(PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(OUT) :: UPrim_slave( PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i,j,iSide
!==================================================================================================================================
DO iSide=1,nSides
  IF(iSide.GE.firstMPISide_YOUR.AND.iSide.LE.lastMPISide_YOUR) CYCLE
  DO j=0,PP_NZ; DO i=0,PP_N
    CALL ConsToPrim(UPrim_master(:,i,j,iSide),U_master(:,i,j,iSide))
  END DO; END DO
END DO
DO iSide=firstInnerSide,lastMPISide_YOUR
  DO j=0,PP_NZ; DO i=0,PP_N
    CALL ConsToPrim(UPrim_slave(:,i,j,iSide),U_slave(:,i,j,iSide))
  END DO; END DO
END DO

!! Version 2: Compute UPrim_master/slave from volume UPrim
!
!#if USE_MPI
!! Prolong to face for MPI sides - send direction
!CALL StartReceiveMPIData(UPrim_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE
!CALL ProlongToFaceCons(PP_N,UPrim,UPrim_master,UPrim_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
!CALL U_Mortar(UPrim_master,UPrim_slave,doMPISides=.TRUE.)
!CALL StartSendMPIData(   UPrim_slave,DataSizeSide,firstSlaveSide,lastSlaveSide,MPIRequest_U(:,RECV),SendID=2) ! Send YOUR
!#endif /*USE_MPI*/
!
!CALL ProlongToFaceCons(PP_N,UPrim,UPrim_master,UPrim_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
!CALL U_Mortar(UPrim_master,UPrim_slave,doMPISides=.FALSE.)
!
!#if USE_MPI
!! Complete send / receive
!CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U) !Send YOUR - receive MINE
!#endif /*USE_MPI*/
END SUBROUTINE GetPrimitiveStateSurface

!==================================================================================================================================
!> Converts primitive variables to conservative solution vector at surfaces.
!> Routine requires mask so that conversion is only done on masked sides.
!==================================================================================================================================
SUBROUTINE GetConservativeStateSurface(UPrim_master,UPrim_slave,U_master,U_slave, mask_master, mask_slave, mask_ref)
! MODULES
USE MOD_Preproc
USE MOD_EOS,      ONLY: PrimToCons
USE MOD_Mesh_Vars,ONLY: firstInnerSide,firstMPISide_YOUR,lastMPISide_YOUR,nSides
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim_master(PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(IN)    :: UPrim_slave( PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
REAL,INTENT(OUT)   :: U_master(    CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(OUT)   :: U_slave(     CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
INTEGER,INTENT(IN) :: mask_master( 1:nSides)                     !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_slave ( 1:nSides)                     !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_ref                                   !< reference value for mask comparison
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i,j,SideID
!==================================================================================================================================
DO SideID=1,nSides
  IF ((firstMPISide_YOUR.LE.SideID).AND.(SideID.LE.lastMPISide_YOUR)) CYCLE
  IF (mask_master(SideID).EQ.mask_ref) THEN
    DO j=0,PP_NZ; DO i=0,PP_N
      CALL PrimToCons(UPrim_master(:,i,j,SideID),U_master(:,i,j,SideID))
    END DO; END DO
  END IF
END DO
DO SideID=firstInnerSide,lastMPISide_YOUR
  IF (mask_slave(SideID).EQ.mask_ref) THEN
    DO j=0,PP_NZ; DO i=0,PP_N
      CALL PrimToCons(UPrim_slave(:,i,j,SideID),U_slave(:,i,j,SideID))
    END DO; END DO
  END IF
END DO
END SUBROUTINE


!==================================================================================================================================
!> Finalizes equation, calls finalize for testcase and Riemann
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars
USE MOD_TestCase        ,ONLY: FinalizeTestcase
USE MOD_Riemann         ,ONLY: FinalizeRiemann
USE MOD_CalcTimeStep    ,ONLY: FinalizeCalctimestep
#if EDDYVISCOSITY
USE MOD_EddyVisc        ,ONLY: FinalizeEddyVisc
#endif /*EDDYVISCOSITY*/
USE MOD_GetBoundaryFlux, ONLY: FinalizeBC
IMPLICIT NONE
!==================================================================================================================================
CALL FinalizeTestcase()
CALL FinalizeRiemann()
CALL FinalizeCalctimestep()
#if EDDYVISCOSITY
CALL FinalizeEddyVisc()
#endif /*EDDYVISCOSITY*/
CALL FinalizeBC()
SDEALLOCATE(RefStatePrim)
SDEALLOCATE(RefStateCons)
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation
