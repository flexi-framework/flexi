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
!> General routines for RANS equations with Spalart-Allmaras turbulence model
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
USE MOD_ReadInTools,ONLY: prms,addStrListEntry
USE MOD_Riemann    ,ONLY: DefineParametersRiemann
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(      'IniRefState',  "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('RefState',     "State(s) in primitive variables (density, velx, vely, velz, pressure).",&
                                                multiple=.TRUE.)
CALL prms%CreateStringOption(   'BCStateFile',  "File containing the reference solution on the boundary to be used as BC.")
CALL prms%CreateRealOption(     'PrTurb',       "Prandtl number"                                       , '0.9')
CALL prms%CreateLogicalOption(  'includeTrip',  "Switch on to include a trip term in the SA equations.", '.FALSE.')
CALL prms%CreateLogicalOption(  'DebugSA',      "Switch on to include debug output for SA equation."   , '.FALSE.')

CALL DefineParametersRiemann()
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
USE MOD_ReadInTools       ,ONLY: CountOption,GETREALARRAY,GETSTR,GETREAL,GETLOGICAL
USE MOD_TestCase          ,ONLY: InitTestcase
USE MOD_Riemann           ,ONLY: InitRiemann
USE MOD_GetBoundaryFlux,   ONLY: InitBC
USE MOD_CalcTimeStep      ,ONLY: InitCalctimestep
USE MOD_Mesh_Vars         ,ONLY: nElems,Elem_xGP,offsetElem,MeshFile,ElemToSide,Face_xGP
USE MOD_HDF5_Input        ,ONLY: ReadArray,OpenDataFile,CloseDataFile,GetDataSize,ReadAttribute
USE MOD_IO_HDF5
#if PP_dim == 3
USE MOD_2D                ,ONLY: ExpandArrayTo3D
#else
USE MOD_2D                ,ONLY: to2D_rank4
#endif
#if FV_ENABLED
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisVolume
USE MOD_FV_Vars          ,ONLY: FV_Vdm
#endif
#if USE_MPI
USE MOD_Mesh_Readin      ,ONLY: ELEMIPROC
USE MOD_MPI_Vars
USE MOD_MPI
#endif
USE MOD_IO_HDF5          ,ONLY:AddToFieldData,FieldOut
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i
REAL               :: UE(PP_2Var)
INTEGER            :: iElem,j,k
INTEGER            :: HSize_proc(4)
REAL               :: RefStatePrimTmp(6)
CHARACTER(LEN=255) :: FileName
REAL,ALLOCATABLE   :: SAd_local(:,:,:,:)
INTEGER            :: tripElem,tripLocSide
LOGICAL            :: file_exists
!==================================================================================================================================
IF(EquationInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitEquation not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RANS WITH SA...'

s43=4./3.
s23=2./3.

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

! Read in boundary parameters
IniRefState  = 0

CALL InitExactFunc()
CALL InitEOS()

! SA-specific parameters
includeTrip = GETLOGICAL('includeTrip')
PrTurb      = GETREAL(   'PrTurb')
ALLOCATE(SAd(0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE,nElems))
! We choose a large number as our default for the walldistance, since it basically means we calculate free turbulence away from a
! wall. The square-root is taken since the value get's squared in the auxilliary functions, and this prevents errounus arithmetic
! operations to take place.
SAd = SQRT(HUGE(1.))
! Read-in of walldistance
FileName = MeshFile(1:INDEX(MeshFile,'_mesh.h5')-1)//'_walldistance.h5'
file_exists = FILEEXISTS(FileName)
IF (file_exists) THEN
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetDataSize(File_ID,'walldistance',nDims,HSize)
  IF (HSize(1)-1.NE.PP_N) CALL Abort(__STAMP__,"Polynomial degree of walldistance file does not match!")
  ALLOCATE(SAd_local(0:HSize(1)-1,0:HSize(2)-1,0:HSize(3)-1,nElems))
  HSize_proc = INT(HSize)
  HSize_proc(4) = nElems
  CALL ReadArray('walldistance',4,&
                 HSize_proc,&
                 offsetElem,4,RealArray=SAd_local)
#if PP_dim == 3
  IF (HSize(3).EQ.1) THEN
    ! Walldistance was created by 2D tool, expand in third dimension
    CALL ExpandArrayTo3D(4,(/PP_N+1,PP_N+1,1,nElems/),3,PP_N+1,SAd_local,SAd(:,:,:,0,:))
  ELSE
    ! 3D walldistance tool and 3D Flexi
    SAd(:,:,:,0,:) = SAd_local
  END IF
#else
  IF (HSize(3).EQ.1) THEN
    ! 2D Walldistace tool and 2D Flexi
    SAd(:,:,:,0,:) = SAd_local
  ELSE
    ! Walldistance was created by 3D tool => reduce third space dimension
    CALL to2D_rank4((/0,0,0,1/),(/PP_N,PP_N,PP_N,nElems/),3,SAd_local)
    SAd(:,:,:,0,:) = SAd_local
  END IF
#endif
  DEALLOCATE(HSize)
  IF (includeTrip) THEN
    CALL ReadAttribute(File_ID,'TripX',2,RealArray=TripX)
    CALL ReadAttribute(File_ID,'TripElem',1,IntScalar=tripElem)
#if USE_MPI
    tripOnProc = ((tripElem.GT.offsetElem+1).AND.(tripElem.LT.(offsetElem+1+nElems)))
#else
    tripOnProc = .TRUE.
#endif
#if USE_MPI
    tripRoot = ELEMIPROC(TripElem)
#endif
    IF (tripOnProc) THEN
      CALL ReadAttribute(File_ID,'TripPQ',2,IntArray=TripPQ)
      CALL ReadAttribute(File_ID,'TripLocSide',1,IntScalar=tripLocSide)
      TripElem = TripElem - offsetElem ! From global to local elem index
      tripSideID = ElemToSide(E2S_SIDE_ID,triplocSide,tripElem)
    END IF
  END IF
  CALL CloseDataFile()
  DEALLOCATE(SAd_local)
ELSE
  includeTrip = .FALSE.
  SWRITE(UNIT_stdOut, *) "WARNING: No walldistance file found! Scaling with walldistance deactivated!"
END IF

IF (includeTrip) THEN
  ALLOCATE(SAdt(0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE,nElems))
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      SAdt(i,j,k,0,iElem) = NORM2(Elem_xGP(1:2,i,j,k,iElem)-TripX)
    END DO; END DO; END DO! i,j,k=0,PP_N
  END DO ! iElem
  IF (tripOnProc) dXt = NORM2(Face_xGP(1:2,0,0,0,tripSideID)-Face_xGP(1:2,PP_N,0,0,tripSideID))/(PP_N+1)
#if USE_MPI
  CALL MPI_BCAST(dXt,1,MPI_DOUBLE_PRECISION,tripRoot,MPI_COMM_FLEXI,iError)
#endif
END IF

doSADebug = GETLOGICAL('DebugSA')
IF (doSADebug) THEN
  ALLOCATE(SADebug(4,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  CALL AddToFieldData(FieldOut,(/4,PP_N+1,PP_N+1,PP_NZ+1,nElems/),'SADebug',(/'Prod','Dest','Trip','Diff'/),RealArray=SADebug)
END IF

#if FV_ENABLED
! Calculate wall distance at sub cell nodes. This assumes that the walldistance can be adequately represented as a polynomial!
! TODO: Replace with seperately calculated wall distance, directly on FV points
DO iElem=1,nElems
  CALL ChangeBasisVolume(PP_N,PP_N,FV_Vdm,SAd(:,:,:,0,iElem),SAd(:,:,:,1,iElem))
END DO ! iElem
#endif

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
    RefStatePrimTmp = GETREALARRAY('RefState',6)
    RefStatePrim(1:5,i)  = RefStatePrimTmp(1:5)
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
    RefStatePrim(NUSA,i) = RefStatePrimTmp(6)
    CALL PrimToCons(RefStatePrim(:,i),RefStateCons(:,i))
  END DO
END IF

! boundary state filename if present
BCStateFile=GETSTR('BCStateFile','nonexistingfile')

! Initialize Riemann solvers to be in volume and on BCs
CALL InitRiemann()

! Initialize timestep calculation
CALL InitCalctimestep()


CALL InitBC()

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RANS WITH SA DONE!'
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: U_master(    CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(IN)    :: U_slave(     CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
REAL,INTENT(INOUT) :: UPrim_master(PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(INOUT) :: UPrim_slave( PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim_master(PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(IN)    :: UPrim_slave( PRIM,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
REAL,INTENT(INOUT) :: U_master(    CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(INOUT) :: U_slave(     CONS,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
INTEGER,INTENT(IN) :: mask_master(1:nSides)                      !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_slave (1:nSides)                      !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_ref                                   !< reference value for mask comparison
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,SideID
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
USE MOD_GetBoundaryFlux, ONLY: FinalizeBC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL FinalizeTestcase()
CALL FinalizeRiemann()
CALL FinalizeCalctimestep()
CALL FinalizeBC()
SDEALLOCATE(RefStatePrim)
SDEALLOCATE(RefStateCons)
SDEALLOCATE(SAd)
SDEALLOCATE(SAdt)
SDEALLOCATE(SADebug)
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation
