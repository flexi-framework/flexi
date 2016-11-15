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
!> Routine performing time averaging of variables and the preparation to computing fluctuations
!> For the fluctuations we rely on the fact that \f$ U^{'} U^{'} = \overline{U}^2 - \overline{U^2} \f$
!> The terms computed in this routine are therefore the TimeAvg: \f$ \overline{U} \f$ and 
!> the squared solution denoted by Fluc: \f$ \overline{U^2} \f$ 
!==================================================================================================================================
MODULE MOD_TimeAverage
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER                        :: nMaxVarAvg,nMaxVarFluc
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitCalcTimeAverage
  MODULE PROCEDURE InitCalcTimeAverage
END INTERFACE

INTERFACE FinalizeTimeAverage
  MODULE PROCEDURE FinalizeTimeAverage
END INTERFACE

INTERFACE CalcTimeAverage
  MODULE PROCEDURE CalcTimeAverage
END INTERFACE

PUBLIC::InitCalcTimeAverage, FinalizeTimeAverage, CalcTimeAverage
!==================================================================================================================================
CONTAINS



!==================================================================================================================================
!> Initializes the time averaging variables and builds map from fluctuation quantities to required time averaged variables
!> (e.g. if VelocityMagnitude fluctuations are to be computed the time averages of the velocity components u,v,w are computed)
!> - only variables specified in the variable list can be averaged
!==================================================================================================================================
SUBROUTINE InitCalcTimeAverage()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools, ONLY: CountOption,GETSTR,GETLOGICAL,GETINT
USE MOD_Mesh_Vars,   ONLY: nElems
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVar,iVar2
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAvgIni(:), VarNamesAvgList(:), VarNamesFlucList(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesFlucIni(:)
LOGICAL,ALLOCATABLE            :: hasAvgVars(:)
!==================================================================================================================================
nVarAvg  = CountOption('VarNameAvg')
nVarFluc = CountOption('VarNameFluc')
IF((nVarAvg.EQ.0).AND.(nVarFluc.EQ.0))THEN
  CALL CollectiveStop(__STAMP__, &
    'No quantities for time averaging have been specified. Please specify quantities or disable time averaging!')
#if FV_ENABLED
ELSE
  CALL CollectiveStop(__STAMP__, &
    'Timeaveraging has not been implemented for FV yet!')
#endif
END IF

! Define variables to be averaged
nMaxVarAvg=15
ALLOCATE(VarNamesAvgList(nMaxVarAvg))
VarNamesAvgList(1)  ='Density'
VarNamesAvgList(2)  ='MomentumX'
VarNamesAvgList(3)  ='MomentumY'
VarNamesAvgList(4)  ='MomentumZ'
VarNamesAvgList(5)  ='EnergyStagnationDensity'
VarNamesAvgList(6)  ='VelocityX'
VarNamesAvgList(7)  ='VelocityY'
VarNamesAvgList(8)  ='VelocityZ'
VarNamesAvgList(9)  ='VelocityMagnitude'
VarNamesAvgList(10) ='Pressure'
VarNamesAvgList(11) ='VelocitySound'
VarNamesAvgList(12) ='Mach'
VarNamesAvgList(13) ='Temperature'
VarNamesAvgList(14) ='TotalTemperature'
VarNamesAvgList(15) ='TotalPressure'

nMaxVarFluc=19
ALLOCATE(VarNamesFlucList(nMaxVarFluc),hasAvgVars(nMaxVarFluc))
hasAvgVars=.TRUE.
!define fluctuation variables
VarNamesFlucList(1)  ='Density'
VarNamesFlucList(2)  ='MomentumX'
VarNamesFlucList(3)  ='MomentumY'
VarNamesFlucList(4)  ='MomentumZ'
VarNamesFlucList(5)  ='EnergyStagnationDensity'
VarNamesFlucList(6)  ='VelocityX'
VarNamesFlucList(7)  ='VelocityY'
VarNamesFlucList(8)  ='VelocityZ'
VarNamesFlucList(9)  ='VelocityMagnitude'
VarNamesFlucList(10) ='Pressure'
VarNamesFlucList(11) ='VelocitySound'
VarNamesFlucList(12) ='Mach'
VarNamesFlucList(13) ='Temperature'
!VarNamesFlucList(14) ='EnthalpyStagnation'
!VarNamesFlucList(15) ='Entropy'
!VarNamesFlucList(16) ='VorticityX'
!VarNamesFlucList(17) ='VorticityY'
!VarNamesFlucList(18) ='VorticityZ'
!VarNamesFlucList(19) ='VorticityMagnitude'
VarNamesFlucList(14) = 'uv'
VarNamesFlucList(15) = 'uw'
VarNamesFlucList(16) = 'vw'
VarNamesFlucList(17) = 'DR_u'; hasAvgVars(17)=.FALSE.
VarNamesFlucList(18) = 'DR_S'; hasAvgVars(18)=.FALSE.
VarNamesFlucList(19) = 'TKE';  hasAvgVars(19)=.FALSE.

! Read VarNames from ini file
ALLOCATE(VarNamesAvgIni(nVarAvg),VarNamesFlucIni(nVarFluc))
DO iVar=1,nVarAvg
  VarNamesAvgIni(iVar)=GETSTR('VarNameAvg')
END DO
DO iVar=1,nVarFluc
  VarNamesFlucIni(iVar)=GETSTR('VarNameFluc')
END DO

! Check which variables have to be calculated and create mappings to global variable index (1:nVarout)
! CalcAvgTmp(1,:) for normal variables, CalcAvgTmp(2,:) for fluctuations
ALLOCATE(CalcAvg(nMaxVarAvg),CalcFluc(nMaxVarFluc))
CalcAvg=.FALSE.
CalcFluc=.FALSE.

! check each average from ini file
DO iVar=1,nVarAvg
  ! check if avg from ini file is in avg list
  iVar2 = GETMAPBYNAME(VarNamesAvgIni(iVar),VarNamesAvgList,nMaxVarAvg)
  IF(iVar2.NE.-1)THEN
    CalcAvg(iVar2) = .TRUE.
  ELSE
    CALL CollectiveStop(__STAMP__, &
    'Specified varname does not exist: ' // VarNamesAvgIni(iVar))
  END IF
END DO

! check each fluctuation from ini file
DO iVar=1,nVarFluc
  ! check if fluc from ini file is in fluc list
  iVar2 = GETMAPBYNAME(VarNamesFlucIni(iVar),VarNamesFlucList,nMaxVarFluc)
  IF(iVar2.NE.-1)THEN
    CalcFluc(iVar2) = .TRUE.
  ELSE
    CALL CollectiveStop(__STAMP__, &
    'Specified varname does not exist: ' // VarNamesFlucIni(iVar))
  END IF

  ! if fluctuation is set also compute base variable
  iVar2 = GETMAPBYNAME(VarNamesFlucIni(iVar),VarNamesAvgList,nMaxVarAvg)
  IF(iVar2.NE.-1) CalcAvg(iVar2) = .TRUE.
END DO

! For fluctuations with mixed base vars
IF(CalcFluc(GETMAPBYNAME('uv',VarNamesFlucList,nMaxVarFluc)))THEN !uv
   CalcAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF
IF(CalcFluc(GETMAPBYNAME('uw',VarNamesFlucList,nMaxVarFluc)))THEN !uw
   CalcAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF
IF(CalcFluc(GETMAPBYNAME('vw',VarNamesFlucList,nMaxVarFluc)))THEN !vw
   CalcAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF
IF(CalcFluc(GETMAPBYNAME('DR_u',VarNamesFlucList,nMaxVarFluc)).OR.&
        CalcFluc(GETMAPBYNAME('DR_S',VarNamesFlucList,nMaxVarFluc)).OR.&
        CalcFluc(GETMAPBYNAME('TKE',VarNamesFlucList,nMaxVarFluc)))THEN !Dissipation,TKE
   CalcAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME  ('Density',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF
nVarAvg=0 ! recount nVarAvg
DO iVar=1,nMaxVarAvg
  IF(CalcAvg(iVar)) nVarAvg=nVarAvg+1
END DO

! Set indices (iAvg,iFluc) and Varnames
ALLOCATE(VarNamesFlucOut(nVarFluc),VarNamesAvgOut(nVarAvg))
ALLOCATE(iAvg(nMaxVarAvg),iFluc(nMaxVarFluc))
! iAvg     -> Mapping from VariableList to index in UAvg array
! iFluc    -> Mapping from index in UFluc array to index in UAvg array
!             (e.g. for mixed term uv: iFluc(1,1) -> u iFluc(2,1) -> v)

VarNamesFlucOut(:)=''
VarNamesAvgOut(:)=''
nVarAvg=0
nVarFluc=0
iAvg=0
iFluc=0
! Build map for avg
DO iVar=1,nMaxVarAvg
  IF(CalcAvg(iVar))THEN
    nVarAvg=nVarAvg+1
    iAvg(iVar)=nVarAvg
    VarNamesAvgOut(nVarAvg) = TRIM(VarNamesAvgList(iVar))
  END IF
END DO
! Build map from fluclist to calcfluc
DO iVar=1,nMaxVarFluc
  IF(CalcFluc(iVar).AND.hasAvgVars(iVar))THEN
    nVarFluc=nVarFluc+1
    iFluc(iVar)=nVarFluc
    VarNamesFlucOut(nVarFluc) = TRIM(VarNamesFlucList(iVar))
  END IF
END DO
nVarFlucHasAvg=nVarFluc
ALLOCATE(FlucAvgMap(2,nVarFlucHasAvg))
FlucAvgMap=0
DO iVar=1,nMaxVarFluc
  IF(CalcFluc(iVar).AND.(.NOT.hasAvgVars(iVar)))THEN
    nVarFluc=nVarFluc+1
    iFluc(iVar)=nVarFluc
    VarNamesFlucOut(nVarFluc) = TRIM(VarNamesFlucList(iVar))
  END IF
END DO

! set map from fluc array to avg array needed to compute fluc
DO iVar=1,nMaxVarFluc
  IF((iFluc(iVar).NE.0).AND.hasAvgVars(iVar))THEN
    iVar2 = GETMAPBYNAME(VarNamesFlucList(iVar),VarNamesAvgList,nMaxVarAvg)
    IF(iVar2.GT.0) FlucAvgMap(:,iFluc(iVar))=iAvg(iVar2)
    IF(iVar.EQ.GETMAPBYNAME('uv',VarNamesFlucList,nMaxVarFluc)) THEN !uv
      FlucAvgMap(1,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg))
      FlucAvgMap(2,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg))
    END IF
    IF(iVar.EQ.GETMAPBYNAME('vw',VarNamesFlucList,nMaxVarFluc)) THEN !vw
      FlucAvgMap(1,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg))
      FlucAvgMap(2,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg))
    END IF
    IF(iVar.EQ.GETMAPBYNAME('uw',VarNamesFlucList,nMaxVarFluc)) THEN !uw
      FlucAvgMap(1,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg))
      FlucAvgMap(2,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg))
    END IF
  END IF
END DO

#if !(PARABOLIC)
IF(CalcFluc(17).OR.CalcFluc(18))THEN
  CALL CollectiveStop(__STAMP__,&
    'Cannot compute dissipation. Not compiled with parabolic flag.')
END IF
#endif /* PARABOLIC */

! Allocate arrays
ALLOCATE(UAvg(nVarAvg,0:PP_N,0:PP_N,0:PP_N,nElems),UFluc(nVarFluc,0:PP_N,0:PP_N,0:PP_N,nElems))
UAvg = 0.
UFluc = 0.
dtOld=0.
dtAvg=0.

DEALLOCATE(VarNamesAvgList,VarNamesAvgIni,VarNamesFlucIni)
END SUBROUTINE InitCalcTimeAverage


!==================================================================================================================================
!> Return index of string VarName in array VarNameList
!==================================================================================================================================
FUNCTION GETMAPBYNAME(VarName,VarNameList,nVarList)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: VarName                 !< string to be compared
CHARACTER(LEN=*),INTENT(IN)    :: VarNameList(nVarList)   !< list of strings to be searched
INTEGER,INTENT(IN)             :: nVarList                !< length of list
INTEGER                        :: GETMAPBYNAME            !< index of VarName in VarNameList
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i
!==================================================================================================================================
GETMAPBYNAME=-1
DO i=1,nVarList
  IF(TRIM(VarName).EQ.TRIM(VarNameList(i)))THEN
    GETMAPBYNAME=i
    RETURN
  END IF
END DO
END FUNCTION



!==================================================================================================================================
!> Compute time averages by trapezoidal rule
!> TODO: extend description
!==================================================================================================================================
SUBROUTINE CalcTimeAverage(Finalize,dt,t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: U
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: GradUx,GradUy,GradUz
#endif
USE MOD_Mesh_Vars    ,ONLY: MeshFile,nElems
USE MOD_HDF5_Output  ,ONLY: WriteTimeAverage
USE MOD_EOS          ,ONLY: ConsToPrim
USE MOD_EOS_Vars     ,ONLY: Kappa
USE MOD_Analyze_Vars ,ONLY: WriteData_dt
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN)              :: Finalize               !< finalized trapezoidal rule and output file
REAL,INTENT(IN)                 :: dt                     !< current timestep for averaging
REAL,INTENT(IN),OPTIONAL        :: t                      !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
REAL                            :: tFuture
REAL                            :: dtStep
REAL                            :: vel(3), Mach
REAL                            :: tmpVars(nVarAvg,0:PP_N,0:PP_N,0:PP_N)
LOGICAL                         :: getPrims=.FALSE.
REAL                            :: prim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N),UE(PP_2Var)
#if PARABOLIC
INTEGER                         :: p,q
REAL                            :: GradVel(1:3,1:3), Shear(1:3,1:3)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
dtStep = (dtOld+dt)*0.5
IF(Finalize) dtStep = dt*0.5
dtAvg  = dtAvg+dtStep
dtOld  = dt
IF(ANY(CalcAvg(6:nMaxVarAvg))) getPrims=.TRUE.

DO iElem=1,nElems
  IF(getPrims)THEN
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      CALL ConsToPrim(prim(:,i,j,k),U(:,i,j,k,iElem))
    END DO; END DO; END DO
  END IF

  ! Compute time averaged variables and fluctuations of these variables
  IF(CalcAvg(1)) &  !'Density'
    tmpVars(iAvg(1),:,:,:) = U(1,:,:,:,iElem)

  IF(CalcAvg(2)) &  !'MomentumX'
    tmpVars(iAvg(2),:,:,:) = U(2,:,:,:,iElem)

  IF(CalcAvg(3)) &  !'MomentumY'
    tmpVars(iAvg(3),:,:,:) = U(3,:,:,:,iElem)

  IF(CalcAvg(4)) &  !'MomentumZ'
    tmpVars(iAvg(4),:,:,:) = U(4,:,:,:,iElem)

  IF(CalcAvg(5)) &  !'EnergyStagnationDensity'
    tmpVars(iAvg(5),:,:,:) = U(5,:,:,:,iElem)

  IF(CalcAvg(6)) &  !'VelocityX'
    tmpVars(iAvg(6),:,:,:)=prim(2,:,:,:)

  IF(CalcAvg(7)) &  !'VelocityY'
    tmpVars(iAvg(7),:,:,:)=prim(3,:,:,:)

  IF(CalcAvg(8)) &  !'VelocityZ'
    tmpVars(iAvg(8),:,:,:)=prim(4,:,:,:)

  IF(CalcAvg(9))THEN  !'VelocityMagnitude'
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      tmpVars(iAvg(9),i,j,k)=SQRT(SUM(prim(2:4,i,j,k)**2))
    END DO; END DO; END DO
  END IF

  IF(CalcAvg(10)) & !'Pressure'
    tmpVars(iAvg(10),:,:,:)=prim(5,:,:,:)

  IF(CalcAvg(11))THEN !'VelocitySound'
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      UE(CONS)=U(:,i,j,k,iElem)
      UE(PRIM)=prim(:,i,j,k)
      UE(SRHO)=1./prim(1,i,j,k)
      tmpVars(iAvg(11),i,j,k)=SPEEDOFSOUND_HE(UE)
    END DO; END DO; END DO
  END IF

  IF(CalcAvg(12))THEN !'Mach'
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      UE(CONS)=U(:,i,j,k,iElem)
      UE(PRIM)=prim(:,i,j,k)
      UE(SRHO)=1./prim(1,i,j,k)
      tmpVars(iAvg(12),i,j,k)=SQRT(SUM(prim(2:4,i,j,k)**2)/SPEEDOFSOUND_HE(UE))
    END DO; END DO; END DO
  END IF

  IF(CalcAvg(13))THEN !'Temperature'
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      tmpVars(iAvg(13),i,j,k)=prim(6,i,j,k)
    END DO; END DO; END DO
  END IF

  IF(CalcAvg(14))THEN !'TotalTemperature'
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      UE(CONS)=U(:,i,j,k,iElem)
      UE(PRIM)=prim(:,i,j,k)
      UE(SRHO)=1./prim(1,i,j,k)
      Mach=SQRT(SUM(prim(2:4,i,j,k)**2)/SPEEDOFSOUND_HE(UE))
      tmpVars(iAvg(14),i,j,k)=TOTAL_TEMPERATURE_H(UE(TEMP),Mach)
    END DO; END DO; END DO
  END IF

  IF(CalcAvg(15))THEN !'TotalPressure
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      UE(CONS)=U(:,i,j,k,iElem)
      UE(PRIM)=prim(:,i,j,k)
      UE(SRHO)=1./prim(1,i,j,k)
      Mach=SQRT(SUM(prim(2:4,i,j,k)**2)/SPEEDOFSOUND_HE(UE))
      tmpVars(iAvg(15),:,:,:)=TOTAL_PRESSURE_H(UE(PRES),Mach)
    END DO; END DO; END DO
  END IF

  UAvg(:,:,:,:,iElem)= UAvg (:,:,:,:,iElem) + tmpVars(1:nVarAvg,:,:,:)*dtStep
  IF(nVarFluc.GT.0)&
    UFluc(1:nVarFlucHasAvg,:,:,:,iElem) = UFluc(1:nVarFlucHasAvg,:,:,:,iElem) + &
                                 tmpVars(FlucAvgMap(1,1:nVarFlucHasAvg),:,:,:)*tmpVars(FlucAvgMap(2,1:nVarFlucHasAvg),:,:,:)*dtStep

END DO ! iElem

#if PARABOLIC
IF(CalcFluc(17).OR.CalcFluc(18))THEN  !'Dissipation via vel gradients'
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      GradVel(:,1)=GradUx(2:4,i,j,k,iElem)
      GradVel(:,2)=GradUy(2:4,i,j,k,iElem)
      GradVel(:,3)=GradUz(2:4,i,j,k,iElem)
      IF(CalcFluc(17))THEN
        Shear=0.5*(Gradvel+TRANSPOSE(GradVel))
        DO p=1,3; DO q=1,3
          UFluc(iFluc(17),i,j,k,iElem) =UFluc(iFluc(17),i,j,k,iElem) + Shear(p,q)*Shear(p,q)*dtStep
        END DO; END DO
      END IF
      IF(CalcFluc(18))THEN
        DO p=1,3; DO q=1,3
          UFluc(iFluc(18),i,j,k,iElem) =UFluc(iFluc(18),i,j,k,iElem) +  GradVel(p,q)*GradVel(p,q)*dtStep
        END DO; END DO
      END IF
    END DO; END DO; END DO
  END DO ! iElem
END IF
#endif /* PARABOLIC */

IF(CalcFluc(19))THEN  !'TKE'
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      vel(1:3) =U(2:4,i,j,k,iElem)/U(1,i,j,k,iElem)
      UFluc(iFluc(19),i,j,k,iElem)=UFluc(iFluc(19),i,j,k,iElem)+SUM(vel(:)**2)*dtStep
    END DO; END DO; END DO
  END DO ! iElem
END IF

! Calc time average and write solution to file
IF(Finalize)THEN
  UAvg =UAvg /dtAvg
  UFluc=UFluc/dtAvg
  tFuture=t+WriteData_dt
  CALL WriteTimeAverage(TRIM(MeshFile),t,tFuture,VarNamesAvgOut,VarNamesFlucOut,UAvg,UFluc,dtAvg,nVarAvg,nVarFluc)
  UAvg=0.
  UFluc=0.
  dtAvg=0.
  dtOld=0.
END IF

END SUBROUTINE CalcTimeAverage



!==================================================================================================================================
!> Finalizes the time averaging routines
!==================================================================================================================================
SUBROUTINE FinalizeTimeAverage()
! MODULES
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(CalcAvg)
SDEALLOCATE(iAvg)
SDEALLOCATE(iFluc)
SDEALLOCATE(UAvg)
SDEALLOCATE(UFluc)
SDEALLOCATE(VarNamesAvgOut)
SDEALLOCATE(VarNamesFlucOut)
END SUBROUTINE FinalizeTimeAverage

END MODULE MOD_TimeAverage
