!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2026 Prof. Andrea Beck
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
!----------------------------------------------------------------------------------------------------------------------------------

LOGICAL           :: doIndicatorBaseFlow    = .FALSE. !< switch whether to compute indicator on baseflow

INTEGER,PARAMETER :: INDTYPE_DG             = 0
INTEGER,PARAMETER :: INDTYPE_FV             = 1
INTEGER,PARAMETER :: INDTYPE_PERSSON        = 2
INTEGER,PARAMETER :: INDTYPE_JAMESON        = 8
INTEGER,PARAMETER :: INDTYPE_DUCROS         = 9
INTEGER,PARAMETER :: INDTYPE_DUCROSTIMESJST = 10
INTEGER,PARAMETER :: INDTYPE_HALFHALF       = 3
INTEGER,PARAMETER :: INDTYPE_CHECKERBOARD   = 33
INTEGER,PARAMETER :: INDTYPE_BLEND          = 42
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: doIndicatorBaseFlow
PUBLIC:: DefineParametersIndicator
PUBLIC:: InitIndicator
PUBLIC:: CalcIndicator
PUBLIC:: FinalizeIndicator
!==================================================================================================================================

CONTAINS

#include "indicator_func.t90"

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersIndicator()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Indicator")
CALL prms%CreateIntFromStringOption('IndicatorType',"Specify type of indicator to be used: DG, FV, Persson,"//&
                                             "Ducros, DucrosTimesJST, halfhalf, checkerboard",&
                                             'DG')
CALL addStrListEntry('IndicatorType','dg',            INDTYPE_DG)
CALL addStrListEntry('IndicatorType','fv',            INDTYPE_FV)
CALL addStrListEntry('IndicatorType','blend'         ,INDTYPE_BLEND)
CALL addStrListEntry('IndicatorType','persson',       INDTYPE_PERSSON)
CALL addStrListEntry('IndicatorType','halfhalf',      INDTYPE_HALFHALF)
CALL addStrListEntry('IndicatorType','checkerboard',  INDTYPE_CHECKERBOARD)
CALL addStrListEntry('IndicatorType','jameson',       INDTYPE_JAMESON)
CALL addStrListEntry('IndicatorType','ducros',        INDTYPE_DUCROS)
CALL addStrListEntry('IndicatorType','ducrostimesjst',INDTYPE_DUCROSTIMESJST)
CALL prms%CreateIntOption('IndVar',        "Specify variable upon which indicator is applied, for general indicators.",&
                                           multiple=.TRUE.)
CALL prms%CreateRealOption('IndStartTime', "Specify physical time when indicator evalution starts. Before this time"//&
                                           "a high indicator value is returned from indicator calculation."//&
                                           "(Idea: FV everywhere at begin of computation to smooth solution)", '0.0')
CALL prms%CreateIntOption('nModes',        "Number of highest modes to be checked for Persson modal indicator.",'2')
CALL prms%CreateLogicalOption('FVBoundaries',  "Use FV discretization in element that contains a side of a certain BC_TYPE", '.FALSE.')
CALL prms%CreateIntOption    ('FVBoundaryType',"BC_TYPE that should be discretized with FV."//&
                                               "Set it to BC_TYPE, setting 0 will apply FV to all BC Sides",multiple=.TRUE.)
END SUBROUTINE DefineParametersIndicator


!==================================================================================================================================
!> Initialize indicators
!==================================================================================================================================
SUBROUTINE InitIndicator()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Indicator_Vars
USE MOD_ReadInTools         ,ONLY: GETINT,GETREAL,GETINTFROMSTR,GETLOGICAL,CountOption
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_IO_HDF5             ,ONLY: AddToElemData,ElementOut
USE MOD_Overintegration_Vars,ONLY: NUnder
USE MOD_Filter_Vars         ,ONLY: NFilter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: nModes_In
INTEGER                                  :: iBC,nFVBoundaryType
INTEGER                                  :: iVar
!==================================================================================================================================
IF(IndicatorInitIsDone) &
  CALL CollectiveStop(__STAMP__, "InitIndicator not ready to be called or already called.")

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT INDICATOR...'

! Read in  parameters
IndicatorType = GETINTFROMSTR('IndicatorType')


SELECT CASE(IndicatorType)
CASE(INDTYPE_BLEND)
#if FV_ENABLED != 2
  CALL ABORT(__STAMP__, &
        "Fixed blending factor only supported for FV_BLENDING.")
#endif /* FV_ENABLED == 2 */
CASE(INDTYPE_JAMESON)
#if EQNSYSNR != 2 /* NOT NAVIER-STOKES */
  CALL ABORT(__STAMP__, &
      "Jameson indicator only works with Navier-Stokes equations.")
#endif /* EQNSYSNR != 2 */
#if FV_ENABLED >= 2
  CALL ABORT(__STAMP__, &
        "Jameson indicator only works with FV switching.")
#endif /* FV_ENABLED == 2 */
CASE(INDTYPE_DUCROS)
#if !(PARABOLIC)
  CALL ABORT(__STAMP__, &
      "Ducros indicator not available without PARABOLIC!")
#endif
#if EQNSYSNR != 2 /* NOT NAVIER-STOKES */
  CALL ABORT(__STAMP__, &
      "Ducros indicator only works with Navier-Stokes equations.")
#endif /* EQNSYSNR != 2 */
#if FV_ENABLED >= 2
  CALL ABORT(__STAMP__, &
      "Ducros indicator only works with FV switching.")
#endif /* FV_ENABLED == 2 */
CASE(INDTYPE_DUCROSTIMESJST)
#if !(PARABOLIC)
  CALL ABORT(__STAMP__, &
      "Ducros*JST indicator not available without PARABOLIC!")
#endif
#if EQNSYSNR != 2 /* NOT NAVIER-STOKES */
  CALL ABORT(__STAMP__, &
      "Ducros*JST indicator only works with Navier-Stokes equations.")
#endif /* EQNSYSNR != 2 */
#if FV_ENABLED >= 2
  CALL ABORT(__STAMP__, &
      "Ducros*JST indicator only works with FV switching.")
#endif /* FV_ENABLED == 2 */
CASE(INDTYPE_PERSSON)
  ! number of modes to be checked by Persson indicator
  nModes_In = GETINT('nModes')
  ! For overintegration, the last PP_N-Nunder modes are empty. Add them to nModes, so we check non-empty ones
  nModes_In = nModes_In+PP_N-MIN(NUnder,NFilter)
  ! Safety checks: At least one mode must be left and only values >0 make sense
  nModes = MAX(1,MIN(PP_N-1,nModes_In))
  IF (nModes.NE.nModes_In) THEN
    SWRITE(UNIT_stdOut,'(A,I0)') 'WARNING: nModes set by user not within range [1,PP_N-1]. Was instead set to nModes=', nModes
  END IF
#if FV_ENABLED >= 2
  T_FV   = 0.5*10**(-1.8*(PP_N+1)**.25) ! Eq.(42) in: S. Hennemann et al., J.Comp.Phy., 2021
  sdT_FV = s_FV/T_FV
#if EQNSYSNR != 2 /* NOT NAVIER-STOKES */
  CALL ABORT(__STAMP__, &
      "Persson indicator for FV-Blending only works with Navier-Stokes equations.")
#endif /* EQNSYSNR != 2 */
#endif /*FV_ENABLED*/
CASE(-1) ! legacy
  IndicatorType=INDTYPE_DG
END SELECT

IndStartTime = GETREAL('IndStartTime')
ALLOCATE(IndValue(nElems))
IndValue=0.
CALL AddToElemData(ElementOut,'IndValue',RealArray=IndValue)

nIndVars = CountOption('IndVar')
! minimum number of indvars = 1
nIndVars = MAX(nIndVars,1)
ALLOCATE(IndVar(nIndVars))
DO iVar=1,nIndVars
  IndVar(iVar) = GETINT('IndVar','1')
END DO

! FV element at boundaries
FVBoundaries    = GETLOGICAL('FVBoundaries')
nFVBoundaryType = CountOption('FVBoundaryType')
ALLOCATE(FVBoundaryType(nFVBoundaryType))
DO iBC=1,nFVBoundaryType
  FVBoundaryType(iBC) = GETINT('FVBoundaryType','0')! which BCType should be at an FV element? Default value means every BC will be FV
END DO

IndicatorInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT INDICATOR DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitIndicator


!==================================================================================================================================
!> Perform calculation of the indicator.
!==================================================================================================================================
SUBROUTINE CalcIndicator(U,t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Indicator_Vars   ,ONLY: IndicatorType,IndValue,IndStartTime
USE MOD_Mesh_Vars        ,ONLY: offsetElem,Elem_xGP,nElems
#if PARABOLIC && EQNSYSNR == 2
USE MOD_Lifting_Vars     ,ONLY: gradUx,gradUy,gradUz
#endif
#if FV_ENABLED == 2
USE MOD_FV_Blending      ,ONLY: FV_ExtendAlpha
#if PP_NodeType == 1
USE MOD_FV_Blending      ,ONLY: FV_CommAlpha
#endif
USE MOD_FV_Vars          ,ONLY: FV_alpha,FV_alpha_min,FV_alpha_max,FV_alpha_fix,FV_doExtendAlpha
USE MOD_Indicator_Vars   ,ONLY: sdT_FV,T_FV
#elif FV_ENABLED == 3
USE MOD_FV_Vars          ,ONLY: FV_alpha,FV_alpha_min,FV_alpha_max,FV_alpha_fix
USE MOD_Indicator_Vars   ,ONLY: sdT_FV,T_FV
#else
USE MOD_FV_Vars          ,ONLY: FV_Elems,FV_sVdm
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisVolume
#endif /*FV_ENABLED==2*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT),TARGET :: U(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)   !< Solution
REAL,INTENT(IN)           :: t                                             !< Simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iElem
#if FV_ENABLED == 1
REAL,POINTER              :: U_P(:,:,:,:)
REAL,TARGET               :: U_DG(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================

! if time is before IndStartTime return high Indicator value (FV)
IF (t.LT.IndStartTime) THEN
#if FV_ENABLED == 1
  IndValue = HUGE(1.)
#elif FV_ENABLED >= 2
  FV_alpha = FV_alpha_max
#endif /*FV_ENABLED*/
  RETURN
END IF

SELECT CASE (IndicatorType)
CASE(INDTYPE_DG) ! no indicator, just a high value to trigger filtering
#if FV_ENABLED == 1
  IndValue=-100
#else
  FV_alpha = 0.
#endif /*FV_ENABLED == 1*/
CASE(INDTYPE_FV) ! indicator everywhere
#if FV_ENABLED == 1
  IndValue = 100
#else
  FV_alpha = FV_alpha_max
#endif /*FV_ENABLED == 1*/
CASE(INDTYPE_BLEND) ! fixed blending factor
#if FV_ENABLED >= 2
  FV_alpha = FV_alpha_fix
#endif
CASE(INDTYPE_PERSSON) ! Modal Persson indicator
#if FV_ENABLED == 2
  DO iElem=1,nElems
    IndValue(iElem) = IndPerssonBlend(U(:,:,:,:,iElem))
    FV_alpha(iElem)  = 1. / (1. + EXP(-sdT_FV * (IndValue(iElem) - T_FV)))
    ! Limit to alpha_max
    FV_alpha(iElem) = MIN(FV_alpha(iElem),FV_alpha_max)
  END DO ! iElem
  CALL FV_ExtendAlpha(FV_alpha)
    ! Do not compute FV contribution for elements below threshold
  DO iElem=1,nElems
    IF (FV_alpha(iElem) .LT. FV_alpha_min) FV_alpha(iElem) = 0.
  END DO ! iElem
#if PP_NodeType == 1
  IF (.NOT.FV_doExtendAlpha) CALL FV_CommAlpha(FV_alpha)
#endif

#elif FV_ENABLED == 3
  DO iElem=1,nElems
    IndValue(iElem) = IndPerssonBlend(U(:,:,:,:,iElem))
    IndValue(iElem) = 1. / (1. + EXP(-sdT_FV * (IndValue(iElem) - T_FV)))
    ! Limit to alpha_max
    FV_alpha(:,:,:,:,iElem) = MIN(IndValue(iElem),FV_alpha_max)
  END DO ! iElem
  ! Do not compute FV contribution for elements below threshold
  DO iElem=1,nElems
    IF (MAXVAL(FV_alpha(:,:,:,:,iElem)) .LT. FV_alpha_min) FV_alpha(:,:,:,:,iElem) = 0.
  END DO ! iElem
#else
  DO iElem=1,nElems
    IF (FV_Elems(iElem).EQ.0) THEN ! DG Element
      U_P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) => U(:,:,:,:,iElem)
    ELSE
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_sVdm,U(:,:,:,:,iElem),U_DG)
      U_P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) => U_DG
    END IF
    IndValue(iElem) = IndPersson(U_P)
  END DO ! iElem
#endif /*FV_ENABLED==2*/
#if EQNSYSNR == 2 /* NAVIER-STOKES */
#if FV_ENABLED
CASE(INDTYPE_JAMESON)
  IndValue = JamesonIndicator(U)
#endif
#if PARABOLIC
CASE(INDTYPE_DUCROS)
  IndValue = DucrosIndicator(gradUx,gradUy,gradUz)
CASE(INDTYPE_DUCROSTIMESJST)
  IndValue = JamesonIndicator(U) * DucrosIndicator(gradUx,gradUy,gradUz)
#endif /*PARABOLIC*/
#endif /* NAVIER-STOKES */
CASE(INDTYPE_HALFHALF)  ! half/half
  DO iElem=1,nElems
    IF (Elem_xGP(1,0,0,0,iElem).GT.0.0) THEN
#if FV_ENABLED == 1
      IndValue(iElem) = 100
#elif FV_ENABLED == 2
      FV_alpha(iElem) = FV_alpha_max
#elif FV_ENABLED == 3
      FV_alpha(:,:,:,:,iElem) = FV_alpha_max
#endif /*FV_ENABLED*/
    ELSE
#if FV_ENABLED == 1
      IndValue(iElem) = -100
#elif FV_ENABLED == 2
      FV_alpha(iElem) = 0.
#elif FV_ENABLED == 3
      FV_alpha(:,:,:,:,iElem) = 0.
#endif /*FV_ENABLED*/
    END IF
  END DO ! iElem
CASE(INDTYPE_CHECKERBOARD) ! every second element (checkerboard like)
   DO iElem = 1, nElems
    IF (MOD(iElem+offsetElem,2).EQ.0) THEN
#if FV_ENABLED == 1
      IndValue(iElem) = -100
#elif FV_ENABLED == 2
      FV_alpha(iElem) = 0.
#elif FV_ENABLED == 3
      FV_alpha(:,:,:,:,iElem) = 0.
#endif /*FV_ENABLED*/
    ELSE
#if FV_ENABLED == 1
      IndValue(iElem) = 100
#elif FV_ENABLED == 2
      FV_alpha(iElem) = FV_alpha_max
#elif FV_ENABLED == 3
      FV_alpha(:,:,:,:,iElem) = FV_alpha_max
#endif /*FV_ENABLED*/
    END IF
  END DO ! iElem = 1, nElems
CASE DEFAULT ! unknown Indicator Type
  CALL ABORT(__STAMP__,&
    "Unknown IndicatorType!")
END SELECT

! obtain indicator value for elements that contain domain boundaries
CALL IndFVBoundaries()

END SUBROUTINE CalcIndicator


!==================================================================================================================================
!> Calculate IndValue for domain boundaries. There are two options:
!> 1.) All boundaries are treated by FV, i.e. FVBoundaryType=0. In this case, all elements that have a side that is a BCSide
!>     will permanently be  FV elements
!> 2.) FVBoundaryType=BC_TYPE, i.e. only elements that contain a BCSide of BC_TYPE will be FV elements
!==================================================================================================================================
SUBROUTINE IndFVBoundaries()
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if FV_ENABLED == 1
USE MOD_Indicator_Vars,  ONLY: IndValue
#else
USE MOD_FV_Vars,         ONLY: FV_alpha,FV_alpha_max
#endif
USE MOD_Indicator_Vars,  ONLY: FVBoundaries,FVBoundaryType
USE MOD_Equation_Vars,   ONLY: nBCByType,BCSideID
USE MOD_Mesh_Vars,       ONLY: SideToElem
USE MOD_Mesh_Vars,       ONLY: nBCs,BoundaryType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,iBC,nBCLoc,BCType,iSide,ElemID
!==================================================================================================================================
IF (FVBoundaries) THEN
  DO iBC=1,nBCs
    BCType = BoundaryType(iBC,BC_TYPE)
    nBCLoc = nBCByType(iBC)
    IF (BCType.EQ.1) CYCLE ! no FV Boundaries at periodic BC
    IF (ANY(FVBoundaryType.EQ.BCType) .OR. ANY(FVBoundaryType.EQ.0)) THEN
      DO iSide=1,nBCLoc
        SideID=BCSideID(iBC,iSide)
        ElemID = SideToElem(S2E_ELEM_ID,SideID)
#if FV_ENABLED == 1
        IndValue(ElemID) = 100.E3
#elif FV_ENABLED == 2
        FV_alpha(ElemID) = FV_alpha_max
#elif FV_ENABLED == 3
        FV_alpha(:,:,:,:,ElemID) = FV_alpha_max
#endif
      END DO !iSide
    ENDIF
  END DO
ELSE
  RETURN
END IF
END SUBROUTINE IndFVBoundaries


!==================================================================================================================================
!> Deallocate indicator variables
!==================================================================================================================================
SUBROUTINE FinalizeIndicator()
! MODULES
USE MOD_Indicator_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(IndVar)
SDEALLOCATE(IndValue)
SDEALLOCATE(FVBoundaryType)

IndicatorInitIsDone=.FALSE.

END SUBROUTINE FinalizeIndicator

END MODULE MOD_Indicator
