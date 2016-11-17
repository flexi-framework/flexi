#include "flexi.h"
#if EQNSYSNR == 2 /* NAVIER-STOKES */ 
#include "eos.h"
#endif

!==================================================================================================================================
!> This module contains the all indicators useable e.g. for Shock-Capturing/Limiting
!>
!> Each indicator function must have the following interface:
!>   Indicator_func(U)
!> where:
!>   REAL,INTENT(IN)  :: U(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!> The indicator is stored in the array IndValue
!==================================================================================================================================
MODULE MOD_Indicator
! MODULES
IMPLICIT NONE

PRIVATE

LOGICAL :: doCalcIndicator=.FALSE. !< switch whether to compute indicator

INTEGER,PARAMETER :: INDTYPE_DG           = 0
INTEGER,PARAMETER :: INDTYPE_FV           = 1
INTEGER,PARAMETER :: INDTYPE_PERSSON      = 2
INTEGER,PARAMETER :: INDTYPE_JAMESON      = 8
INTEGER,PARAMETER :: INDTYPE_DUCROS       = 9
INTEGER,PARAMETER :: INDTYPE_HALFHALF     = 3
INTEGER,PARAMETER :: INDTYPE_CHECKERBOARD = 33

INTERFACE InitIndicator
  MODULE PROCEDURE InitIndicator
END INTERFACE

INTERFACE CalcIndicator
  MODULE PROCEDURE CalcIndicator
END INTERFACE

INTERFACE IndPersson
  MODULE PROCEDURE IndPersson
END INTERFACE

#if EQNSYSNR == 2 /* NAVIER-STOKES */ 
#if PARABOLIC
INTERFACE DucrosIndicator
  MODULE PROCEDURE DucrosIndicator
END INTERFACE
#endif /* PARABOLIC */

#if FV_ENABLED
INTERFACE JamesonIndicator
  MODULE PROCEDURE JamesonIndicator
END INTERFACE
#endif /* FV_ENABLED */
#endif /* EQNSYSNR == 2 */


INTERFACE FinalizeIndicator
  MODULE PROCEDURE FinalizeIndicator
END INTERFACE

PUBLIC::doCalcIndicator
PUBLIC::InitIndicator
PUBLIC::CalcIndicator
PUBLIC::IndPersson
PUBLIC::FinalizeIndicator
!==================================================================================================================================

PUBLIC::DefineParametersIndicator
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersIndicator()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Indicator")
CALL prms%CreateIntFromStringOption('IndicatorType',"Specify type of indicator to be used: DG, FV, Persson,"//&
                                             "Ducros, halfhalf, checkerboard",&
                                             'DG')
CALL addStrListEntry('IndicatorType','dg',          INDTYPE_DG)
CALL addStrListEntry('IndicatorType','fv',          INDTYPE_FV)
CALL addStrListEntry('IndicatorType','persson',     INDTYPE_PERSSON)
CALL addStrListEntry('IndicatorType','halfhalf',    INDTYPE_HALFHALF)
CALL addStrListEntry('IndicatorType','checkerboard',INDTYPE_CHECKERBOARD)
CALL addStrListEntry('IndicatorType','jameson',     INDTYPE_JAMESON)
CALL addStrListEntry('IndicatorType','ducros',      INDTYPE_DUCROS)
CALL prms%CreateIntOption('IndVar',        "Specify variable upon which indicator is applied, for general indicators.",&
                                           '1')
CALL prms%CreateRealOption('IndStartTime', "Specify physical time when indicator evalution starts. Before this time"//&
                                           "a high indicator value is returned from indicator calculation."//&
                                           "(Idea: FV everywhere at begin of computation to smooth solution)", '0.0')
END SUBROUTINE DefineParametersIndicator


!==================================================================================================================================
!> Initialize indicators
!==================================================================================================================================
SUBROUTINE InitIndicator()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Indicator_Vars
USE MOD_ReadInTools    ,ONLY: GETINT,GETREAL,GETINTFROMSTR
USE MOD_Mesh_Vars      ,ONLY: nElems
USE MOD_IO_HDF5        ,ONLY: AddToElemData
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF(IndicatorInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitIndicator not ready to be called or already called.")
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT INDICATORS...'

! Read in  parameters
IndicatorType = GETINTFROMSTR('IndicatorType')
SELECT CASE(IndicatorType)
CASE(INDTYPE_JAMESON)
#if !(FV_ENABLED)
  CALL Abort(__STAMP__, &
      "Jameson indicator only works with FV_ENABLED.")
#endif
CASE(INDTYPE_DUCROS)
#if !(PARABOLIC)
  CALL Abort(__STAMP__, &
      "Ducros indicator not available without PARABOLIC!")
#endif
CASE(-1) ! legacy
  IndicatorType=INDTYPE_DG
END SELECT

IndStartTime = GETREAL('IndStartTime')
ALLOCATE(IndValue(nElems))
IndValue=0.
CALL AddToElemData('Indicator',RealArray=IndValue)

IndVar = GETINT('IndVar','1')

IndicatorInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT INDICATOR DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitIndicator

SUBROUTINE CalcIndicator(U,t)
!==================================================================================================================================
! Perform calculation of the limiter
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Indicator_Vars ,ONLY: IndicatorType,IndVar,IndValue,IndStartTime
USE MOD_Mesh_Vars      ,ONLY: offsetElem,Elem_xGP,nElems
#if FV_ENABLED
USE MOD_FV_Vars        ,ONLY: FV_Elems,FV_sVdm
#endif /* FV_ENABLED */
#if PARABOLIC
USE MOD_Lifting_Vars   ,ONLY: gradUx,gradUy,gradUz
#endif
USE MOD_ChangeBasis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
! U
REAL,INTENT(INOUT),TARGET :: U(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL,INTENT(IN)           :: t
! U
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iElem
#if FV_ENABLED
REAL,TARGET               :: U_DG(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
#endif
REAL,POINTER              :: U_P(:,:,:,:)
!==================================================================================================================================
! if time is before IndStartTime return high Indicator value (FV)
IF (t.LT.IndStartTime) THEN
  IndValue = 1.E16
  RETURN
END IF

SELECT CASE (IndicatorType)
CASE(INDTYPE_DG) ! no indicator, just a high value to trigger filtering
  IndValue=-100
CASE(INDTYPE_FV) ! indicator everywhere
  IndValue = 100
CASE(INDTYPE_PERSSON) ! Modal Persson indicator
  DO iElem=1,nElems
#if FV_ENABLED
    IF (FV_Elems(iElem).EQ.0) THEN ! DG Element 
#endif      
      U_P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N) => U(:,:,:,:,iElem)
#if FV_ENABLED
    ELSE
      CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_sVdm,U(:,:,:,:,iElem),U_DG)
      U_P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N) => U_DG
    END IF
#endif
    IndValue(iElem) = IndPersson(U_P(IndVar,:,:,:))
  END DO ! iElem
#if EQNSYSNR == 2 /* NAVIER-STOKES */ 
#if FV_ENABLED
CASE(INDTYPE_JAMESON) 
  IndValue = JamesonIndicator(U)
#endif
#if PARABOLIC
CASE(INDTYPE_DUCROS) 
  IndValue = DucrosIndicator(gradUx,gradUy,gradUz)
#endif
#endif /* NAVIER-STOKES */
CASE(INDTYPE_HALFHALF)  ! half/half
  DO iElem=1,nElems
    IF (Elem_xGP(1,0,0,0,iElem).GT.0.0) THEN
      IndValue(iElem) = 100
    ELSE
      IndValue(iElem) = -100
    END IF
  END DO ! iElem

CASE(INDTYPE_CHECKERBOARD) ! every second element (checkerboard like)
   DO iElem = 1, nElems
    IF (MOD(iElem+offsetElem,2).EQ.0) THEN 
      IndValue(iElem) = -100
    ELSE
      IndValue(iElem) =  100
    END IF
  END DO ! iElem = 1, nElems
CASE DEFAULT ! unknown Indicator Type
  CALL abort(__STAMP__,&
    "Unknown IndicatorType!")
END SELECT

END SUBROUTINE CalcIndicator


FUNCTION IndPersson(U) RESULT(IndValue)
!==================================================================================================================================
!> Determine, if given a modal representation solution "U_Modal" is oscillating
!> Indicator value is scaled to \f$\sigma=0 \ldots 1\f$
!> Suggested by Persson et al.
!==================================================================================================================================
USE MOD_PreProc
USE MOD_Interpolation_Vars, ONLY:sVdm_Leg
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: U(0:PP_N,0:PP_N,0:PP_N)
REAL               :: IndValue
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: nDeg,iDeg,i,j,k,l
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N) :: U_Xi
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N) :: U_Eta
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N) :: U_Modal
!==================================================================================================================================

! Transform nodal solution to a modal representation
U_Xi   = 0.
U_Eta  = 0.
U_Modal= 0.
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO l=0,PP_N
  U_Xi(i,j,k)    = U_Xi(i,j,k)    + sVdm_Leg(i,l)*U(l,j,k)
END DO ; END DO ; END DO ; END DO 
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO l=0,PP_N
  U_Eta(i,j,k)   = U_Eta(i,j,k)   + sVdm_Leg(j,l)*U_Xi(i,l,k) 
END DO ; END DO ; END DO ; END DO 
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO l=0,PP_N
  U_Modal(i,j,k) = U_Modal(i,j,k) + sVdm_Leg(k,l)*U_Eta(i,j,l) 
END DO ; END DO ; END DO ; END DO 

! Adapted Persson indicator
IndValue=TINY(0.)
nDeg=MIN(PP_N-1,1)
!nDeg=0
DO iDeg=0,nDeg
  ! Build maximum of 1D indicators
  ! Xi
  IndValue=MAX(IndValue,SUM(U_Modal(PP_N-iDeg:PP_N-iDeg,:,:))**2 /  &
                        (SUM(U_Modal(0:PP_N-iDeg,:,:))**2+EPSILON(0.)))
  ! Eta
  IndValue=MAX(IndValue,SUM(U_Modal(:,PP_N-iDeg:PP_N-iDeg,:))**2 /  &
                        (SUM(U_Modal(:,0:PP_N-iDeg,:))**2+EPSILON(0.)))
  ! Zeta
  IndValue=MAX(IndValue,SUM(U_Modal(:,:,PP_N-iDeg:PP_N-iDeg))**2 /  &
                        (SUM(U_Modal(:,:,0:PP_N-iDeg))**2+EPSILON(0.)))
END DO
! Normalize indicator value
IndValue=LOG10(IndValue)

END FUNCTION IndPersson

#if EQNSYSNR == 2 /* NAVIER-STOKES */ 
#if PARABOLIC
FUNCTION DucrosIndicator(gradUx, gradUy, gradUz) RESULT(IndValue) 
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: nElems,sJ
USE MOD_Analyze_Vars       ,ONLY: wGPVol
#if FV_ENABLED
USE MOD_FV_Vars            ,ONLY: FV_Elems
#endif
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: gradUx(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL,INTENT(IN)    :: gradUy(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL,INTENT(IN)    :: gradUz(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL               :: IndValue(1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: VorticityLoc(3),Vorticity2,IntegrationWeight
REAL    :: divV2
!==================================================================================================================================
DO iElem=1,nElems
  IndValue = 0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      VorticityLoc(1)=gradUy(4,i,j,k,iElem)-gradUz(3,i,j,k,iElem)  ! dw/dy-dv/dz
      VorticityLoc(2)=gradUz(2,i,j,k,iElem)-gradUx(4,i,j,k,iElem)  ! du/dz-dw/dx
      VorticityLoc(3)=gradUx(3,i,j,k,iElem)-gradUy(2,i,j,k,iElem)  ! dv/dx-du/dy
      Vorticity2=SUM(VorticityLoc(:)**2)
      
      divV2 = (gradUx(2,i,j,k,iElem) + gradUy(3,i,j,k,iElem) + gradUz(4,i,j,k,iElem))**2

      IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem,FV_Elems(iElem))
      IF (Vorticity2.LT.100) Vorticity2 = 0.
      IF (divV2.LT.100) divV2 = 0.
      IndValue = IndValue + divV2 /(divV2 + Vorticity2 + 1e-15)* IntegrationWeight
      !IndValue = IndValue + divV2 * IntegrationWeight
      !IndValue = IndValue + Vorticity2 * IntegrationWeight
  END DO; END DO; END DO
  !IndValue = (EXP(IndValue/ElemVol)-EXP(0.))/(EXP(1.)-EXP(0.))
END DO ! iElem
END FUNCTION DucrosIndicator
#endif /* PARABOLIC */

#if FV_ENABLED
FUNCTION JamesonIndicator(U) RESULT(IndValue)
USE MOD_PreProc
USE MOD_Globals
USE MOD_Indicator_Vars     ,ONLY: IndVar
USE MOD_EOS_Vars           ,ONLY: KappaM1
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus
USE MOD_Mesh_Vars          ,ONLY: nElems,nSides
USE MOD_Mesh_Vars          ,ONLY: firstMortarInnerSide,lastMortarInnerSide,firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars          ,ONLY: firstBCSide,lastBCSide
USE MOD_Mesh_Vars          ,ONLY: MortarType,ElemToSide
USE MOD_Mesh_Vars          ,ONLY: sJ
USE MOD_Mappings           ,ONLY: SideToVol
USE MOD_Analyze_Vars       ,ONLY: wGPVol
USE MOD_ProlongToFace1     ,ONLY: ProlongToFace1
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY: MPIRequest_U,MPIRequest_Flux,nNbProcs
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_FillMortar1        ,ONLY: U_Mortar1,Flux_Mortar1
USE MOD_FV_Vars            ,ONLY: FV_Elems,FV_Elems_master,FV_Elems_slave
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)           :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL                      :: IndValue(1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: ElemVol,IntegrationWeight
INTEGER                   :: i,j,k,flip,p,q,SideID,ijk(3),iSide,iElem
REAL                      :: v(-1:PP_N+1,-1:PP_N+1,-1:PP_N+1),vmin,vmax
REAL                      :: UJameson(1:1,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL                      :: UJameson_master(1:1,0:PP_N,0:PP_N,1:nSides)
REAL                      :: UJameson_slave( 1:1,0:PP_N,0:PP_N,1:nSides)
REAL                      :: UE(1:PP_2Var)
INTEGER                   :: TMP(1:nElems)
INTEGER                   :: TMP_master(1:nSides)
INTEGER                   :: TMP_slave( 1:nSides)
INTEGER                   :: DataSizeSide_loc
INTEGER                   :: firstMortarSideID,lastMortarSideID
INTEGER                   :: MortarSideID,tf
!==================================================================================================================================
! Fill UJameson with conservative variable or pressure
SELECT CASE(IndVar)
CASE(1:PP_nVar)
  UJameson(1,:,:,:,:) = U(IndVar,:,:,:,:)
CASE(6)
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      UE(CONS)=U(:,i,j,k,iElem)
      UE(SRHO)=1./UE(DENS)
      UE(VELV)=VELOCITY_HE(UE)
      UJameson(1,i,j,k,iElem)=PRESSURE_HE(UE)
    END DO; END DO; END DO! i,j,k=0,PP_N
  END DO ! iElem
END SELECT

! dummy parameters to store FV_Elems,FV_Elems_master/slave before overwriting them 
! to force ProlongToFace to use FV version everywhere
TMP        = FV_Elems
TMP_master = FV_Elems_master
TMP_slave  = FV_Elems_slave
FV_Elems        = 1
FV_Elems_master = 1
FV_Elems_slave  = 1

! prolongate UJameson to the faces (FV everywhere) 
! and bring it from the big to the small mortar faces
#if USE_MPI
CALL ProlongToFace1(PP_N,UJameson,UJameson_master,UJameson_slave,L_Minus,L_Plus,doMPiSides=.TRUE.)
#endif
CALL ProlongToFace1(PP_N,UJameson,UJameson_master,UJameson_slave,L_Minus,L_Plus,doMPiSides=.FALSE.)
#if USE_MPI
! revert the temporal forcing to use FV everywhere in the ProlongToFace
FV_Elems        = TMP       
FV_Elems_master = TMP_master
FV_Elems_slave  = TMP_slave 
CALL U_Mortar1(UJameson_master,UJameson_slave,doMPiSides=.TRUE.)
#endif
CALL U_Mortar1(UJameson_master,UJameson_slave,doMPiSides=.FALSE.)

! communicate UJameson_master from master to slave
! communicate UJameson_slave  from slave  to master
#if USE_MPI
DataSizeSide_loc = (PP_N+1)**2
CALL StartReceiveMPIData(UJameson_slave ,DataSizeSide_loc,1,nSides,MPIRequest_U(   :,SEND),SendID=2) !  U_slave: slave -> master
CALL StartSendMPIData(   UJameson_slave ,DataSizeSide_loc,1,nSides,MPIRequest_U(   :,RECV),SendID=2) !  U_slave: slave -> master
CALL StartReceiveMPIData(UJameson_master,DataSizeSide_loc,1,nSides,MPIRequest_Flux(:,SEND),SendID=1) !  U_master: master -> slave
CALL StartSendMPIData(   UJameson_master,DataSizeSide_loc,1,nSides,MPIRequest_Flux(:,RECV),SendID=1) !  U_master: master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)   
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux)
#endif

! bring UJameson from small to big mortar faces
DO tf=0,1
  ! ATTENTION: call Flux_Mortar1 with swapped slave/master arguments, since we use it to bring the small
  !            side UJameson to the big Mortar UJameson!
  CALL Flux_Mortar1(UJameson_slave,UJameson_master,doMPISides=(tf.EQ.0),weak=.FALSE.)
  firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,tf.EQ.0) 
   lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,tf.EQ.0) 
  DO MortarSideID=firstMortarSideID,lastMortarSideID
    SELECT CASE(MortarType(1,MortarSideID))
    CASE(1) !1->4
      UJameson_slave(:,:,:,MortarSideID) = 0.25 * UJameson_slave(:,:,:,MortarSideID)
    CASE(2) !1->2 in eta
      UJameson_slave(:,:,:,MortarSideID) = 0.5  * UJameson_slave(:,:,:,MortarSideID)
    CASE(3) !1->2 in xi
      UJameson_slave(:,:,:,MortarSideID) = 0.5  * UJameson_slave(:,:,:,MortarSideID)
    END SELECT
  END DO
END DO

! evaluate the Jameson indicator for each element
DO iElem=1,nElems
  v(0:PP_N,0:PP_N,0:PP_N) = UJameson(1,0:PP_N,0:PP_N,0:PP_N,iElem)

  DO iSide=1,6
    Flip   = ElemToSide(E2S_FLIP,   iSide,iElem)
    SideID = ElemToSide(E2S_SIDE_ID,iSide,iElem)
    IF ((firstBCSide.LE.SideID.AND.SideID.LE.lastBCSide)) THEN ! BC side
      ! if we are at a BC side, then we have to use the local data, which is prolongated into UJameson_master 
      DO q=0,PP_N; DO p=0,PP_N
        ijk = SideToVol(PP_N,-1,p,q,Flip,iSide)
        v(ijk(1),ijk(2),ijk(3)) = UJameson_master(1,p,q,SideID)
      END DO; END DO ! p,q=0,PP_N
    ELSE
      IF (Flip.EQ.0) THEN ! non BC side
        ! Master side => use data from slave side
        DO q=0,PP_N; DO p=0,PP_N
          ijk = SideToVol(PP_N,-1,p,q,Flip,iSide)
          v(ijk(1),ijk(2),ijk(3)) = UJameson_slave(1,p,q,SideID)
        END DO; END DO ! p,q=0,PP_N
      ELSE
        ! Slave side => use data from master side 
        DO q=0,PP_N; DO p=0,PP_N
          ijk = SideToVol(PP_N,-1,p,q,Flip,iSide)
          v(ijk(1),ijk(2),ijk(3)) = UJameson_master(1,p,q,SideID)
        END DO; END DO ! p,q=0,PP_N
      END IF
    END IF
  END DO

  ElemVol = 0.0
  IndValue(iElem) = 0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem,FV_Elems(iElem))
    ElemVol = ElemVol + IntegrationWeight
    vmin = MIN(v(i,j,k),v(i-1,j,k),v(i+1,j,k))
    vmin = MIN(vmin,    v(i,j-1,k),v(i,j+1,k))
    vmin = MIN(vmin,    v(i,j,k-1),v(i,j,k+1))
    vmax = MAX(v(i,j,k),v(i-1,j,k),v(i+1,j,k))
    vmax = MAX(vmax,    v(i,j-1,k),v(i,j+1,k))
    vmax = MAX(vmax,    v(i,j,k-1),v(i,j,k+1))
    IndValue(iElem) = IndValue(iElem) + ABS(vmin-2.*v(i,j,k)+vmax)/ABS(vmin+2.*v(i,j,k)+vmax) * IntegrationWeight
  END DO; END DO; END DO! i,j,k=0,PP_N
  IndValue(iElem) = IndValue(iElem) / ElemVol
END DO ! iElem
END FUNCTION JamesonIndicator
#endif

#endif /* EQNSYSNR == 2 */


!==================================================================================================================================
!> Deallocate indicator variables
!==================================================================================================================================
SUBROUTINE FinalizeIndicator()
! MODULES
USE MOD_Indicator_Vars
IMPLICIT NONE
!==================================================================================================================================
IndicatorInitIsDone=.FALSE.
SDEALLOCATE(IndValue)
END SUBROUTINE FinalizeIndicator

END MODULE MOD_Indicator
