!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
#if FV_ENABLED
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Module for the Finite Volume sub-cells shock capturing.
!>
!> DG elements, that are detected to contain a shock/high gradients/oscillations/..., can be switched to a Finite Volume scheme.
!> A DG element of polynomial degree N is subdivided into (N+1)^dim sub-cells (to each Gauss Point/DOF one FV sub-cell).
!> The FV sub-cells of such an element are updated using FV method with 2nd order TVD reconstruction (slope limiters).
!==================================================================================================================================
MODULE MOD_FV
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DefineParametersFV
  MODULE PROCEDURE DefineParametersFV
END INTERFACE

INTERFACE InitFV
  MODULE PROCEDURE InitFV
END INTERFACE

INTERFACE FV_DGtoFV
  MODULE PROCEDURE FV_DGtoFV
END INTERFACE

INTERFACE FV_PrimToCons
  MODULE PROCEDURE FV_PrimToCons
END INTERFACE

INTERFACE FV_ConsToPrim
  MODULE PROCEDURE FV_ConsToPrim
END INTERFACE

INTERFACE FinalizeFV
  MODULE PROCEDURE FinalizeFV
END INTERFACE

PUBLIC::DefineParametersFV
PUBLIC::InitFV
PUBLIC::FV_DGtoFV
PUBLIC::FV_PrimToCons
PUBLIC::FV_ConsToPrim
PUBLIC::FinalizeFV
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FV
!==================================================================================================================================
SUBROUTINE DefineParametersFV()
! MODULES
USE MOD_Globals
USE MOD_FV_Basis    ,ONLY: DefineParametersFV_Basis
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
#if FV_RECONSTRUCT
USE MOD_FV_Limiter  ,ONLY: DefineParametersFV_Limiter
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('FV')
! FV Switching
CALL prms%CreateLogicalOption('FV_SwitchConservative',"Perform FV/DG switch in reference element"                                 &
                                                     ,'.TRUE.')
CALL prms%CreateLogicalOption('FV_IniSupersample'    ,"Supersample initial solution inside each sub-cell and take mean value \n"//&
                                                      " as average sub-cell value."                                               &
                                                     ,'.TRUE.')
CALL prms%CreateLogicalOption('FV_IniSharp'          ,"Maintain a sharp interface in the initial solution in the FV region"       &
                                                     ,'.FALSE.')
CALL prms%CreateRealOption(   'FV_IndUpperThreshold' ,"Upper threshold: Element is switched from DG to FV if indicator \n"      //&
                                                      "rises above this value"                                                    &
                                                     ,'99.')
CALL prms%CreateRealOption(   'FV_IndLowerThreshold' ,"Lower threshold: Element is switched from FV to DG if indicator \n"      //&
                                                      "falls below this value"                                                    &
                                                     ,'-99.')
CALL prms%CreateLogicalOption('FV_toDG_indicator'    ,"Apply additional Persson indicator to check if DG solution after \n"     //&
                                                      " switch from FV to DG is valid."                                           &
                                                     ,'.FALSE.')
CALL prms%CreateRealOption   ('FV_toDG_limit'        ,"Threshold for FV_toDG_indicator")
CALL prms%CreateLogicalOption('FV_toDGinRK'          ,"Allow switching of FV elements to DG during Runge Kutta stages. \n"      //&
                                                      "This may violated the DG timestep restriction of the element."             &
                                                     ,'.FALSE.')
! FV Blending
CALL prms%CreateRealOption(   'FV_alpha_min'          ,"Lower bound for alpha (all elements below threshold are treated as pure DG)"&
                                                      ,'0.01')
CALL prms%CreateRealOption(   'FV_alpha_max'          ,"Maximum value for alpha",'0.5' )
CALL prms%CreateRealOption(   'FV_alpha_fix'          , "Specify a fixed Blending factor for IndicatorType blend.", '0.0')
CALL prms%CreateRealOption(   'FV_alpha_ExtScale'     ,"Scaling factor for elpha if extended into neighboring elements",'0.5' )
CALL prms%CreateIntOption(    'FV_nExtendAlpha'       ,"Number of times alpha should be passed to neighbor elements per timestep",&
                                                       '1' )
CALL prms%CreateLogicalOption('FV_doExtendAlpha'      ,"Blending factor is prolongated into neighboring elements", '.FALSE.')

#if FV_RECONSTRUCT
CALL DefineParametersFV_Limiter()
#endif
CALL DefineParametersFV_Basis()
END SUBROUTINE DefineParametersFV

!==================================================================================================================================
!> Read in parameters needed for FV sub-cells (indicator min/max and type of limiter) and allocate several arrays.
!> Build metrics for FV sub-cells and performe initial switch from DG to FV sub-cells for all troubled cells.
!==================================================================================================================================
SUBROUTINE InitFV()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis               ,ONLY: InitializeVandermonde
USE MOD_FV_Vars
USE MOD_FV_Basis
#if FV_ENABLED == 1
USE MOD_Filter_Vars         ,ONLY: NFilter
USE MOD_Indicator_Vars      ,ONLY: nModes,IndicatorType
USE MOD_Overintegration_Vars,ONLY: NUnder
#endif
USE MOD_IO_HDF5             ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars           ,ONLY: nElems,nSides
USE MOD_ReadInTools
#if FV_RECONSTRUCT
USE MOD_FV_Limiter
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF(.NOT.FVInitBasisIsDone)THEN
   CALL CollectiveStop(__STAMP__,&
     'InitFV not ready to be called or already called.')
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FV...'

! Read flag, which allows to perform the switching from FV to DG in the reference element
switchConservative = GETLOGICAL("FV_SwitchConservative")

#if FV_ENABLED == 1
! Read minimal and maximal threshold for the indicator
FV_IndLowerThreshold = GETREAL('FV_IndLowerThreshold')
FV_IndUpperThreshold = GETREAL('FV_IndUpperThreshold')

! Read flag indicating, if an additional Persson indicator should check if a FV sub-cells element really contains no oscillations
! anymore.
FV_toDG_indicator = GETLOGICAL('FV_toDG_indicator')
IF (FV_toDG_indicator) THEN
  FV_toDG_limit = GETREAL('FV_toDG_limit')
  ! If the main indicator is not already the persson indicator, then we need to read in the parameters
  IF (IndicatorType .NE. 2) THEN
    nModes = GETINT('nModes')
    nModes = MAX(1,MIN(PP_N-1,nModes+PP_N-MIN(NUnder,NFilter)))
    SWRITE(UNIT_stdOut,'(A,I0)') ' | nModes = ', nModes
  END IF
END IF

! Read flag, which allows switching from FV to DG between the stages of a Runge-Kutta time step
! (this might lead to instabilities, since the time step for a DG element is smaller)
FV_toDGinRK = GETLOGICAL("FV_toDGinRK")

! Options for initial solution
FV_IniSharp       = GETLOGICAL("FV_IniSharp")
IF (.NOT.FV_IniSharp) FV_IniSupersample = GETLOGICAL("FV_IniSupersample")

#elif FV_ENABLED == 2
! Initialize parameters for FV Blending
FV_alpha_min = GETREAL('FV_alpha_min')
FV_alpha_max = GETREAL('FV_alpha_max')
FV_alpha_fix = GETREAL('FV_alpha_fix')
FV_doExtendAlpha = GETLOGICAL('FV_doExtendAlpha')
IF (FV_doExtendAlpha) THEN
  FV_nExtendAlpha = GETINT('FV_nExtendAlpha')
  FV_alpha_extScale = GETREAL('FV_alpha_extScale')
  IF ((FV_alpha_extScale.GT.1.) .OR. (FV_alpha_extScale.LT.0.)) CALL ABORT(__STAMP__,&
                                      'The parameter FV_alpha_extScale has to be between 0. and 1.!')
ENDIF

ALLOCATE(FV_alpha(1:nElems))
ALLOCATE(FV_alpha_master(nSides))
ALLOCATE(FV_alpha_slave( nSides))
CALL AddToElemData(ElementOut,'FV_alpha',FV_alpha)
#endif /*FV_ENABLED*/

#if FV_RECONSTRUCT
CALL InitFV_Limiter()
#endif

ALLOCATE(FV_Elems(nElems)) ! holds information if element is DG (0) or FV (1)
! All cells are initially DG cells
FV_Elems = 0
CALL AddToElemData(ElementOut,'FV_Elems',IntArray=FV_Elems) ! append this array to HDF5 output files

! The elementwise information of 'FV_Elems' is also needed at the faces and therefore
! is 'prolongated' to the faces into the arrays 'FV_Elems_master/slave'.
! The additional 'FV_Elems_Sum' array sums up these two arrays in the following way:
!     FV_Elems_Sum = FV_Elems_master + 2 * FV_Elems_slave
! This leads to the following information stored in 'FV_Elems_Sum' per face:
!             FV_Elems_Sum  |  0 |  1 |  2 |  3 |
!   master side element is  | DG | FV | DG | FV |
!    slave side element is  | DG | DG | FV | FV |
!ALLOCATE(FV_Elems_master(1:nSides)) ! moved to InitFV_Metrics, since needed there for U_Mortar routine
ALLOCATE(FV_Elems_slave(1:nSides))
ALLOCATE(FV_Elems_Sum(1:nSides))
FV_Elems_master = 0
FV_Elems_slave = 0
FV_Elems_Sum = 0

! arrays for FV/DG statistics
ALLOCATE(FV_Elems_counter(nElems))
ALLOCATE(FV_Elems_Amount(nElems))
FV_Elems_counter  = 0
FV_Switch_counter = 0
FV_Elems_Amount = 0
CALL AddToElemData(ElementOut,'FV_Elems_Amount',RealArray=FV_Elems_Amount)

#if FV_RECONSTRUCT
! Allocate array for multi purposes:
! - For FV elements it stores the slope between the nodes next and second next to the interface.
!    |  x    x    x    x  |                          | = face, x = node
!                  <-->   ^ to this interface
!                    ^ the slope between those nodes
! - For DG elements it stores the solution at the nodes next to the interface.
!    |  x    x    x    x  |                          | = face, x = node
!                         ^ to this interface
!                      ^ the solution at this node
ALLOCATE(FV_multi_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(FV_multi_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
FV_multi_slave = 0.0
FV_multi_master = 0.0

! Allocate array for FD-gradient over faces
!    | x  x  x  x | x  x  x  x |                     | = face, x = node
!               <--->
!                  ^ the slope over the face
ALLOCATE(FV_surf_gradU(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))

! The gradients of the primitive variables are stored at each volume integration point and
! are computed by limiting the slopes to the two adjacent points in the respective direction.
! These are physical gradients, but they are labeled ...xi/eta/zeta, since they are the slopes
! along the xi-/eta-/zeta-lines in physical space. These slopes are required to reconstruct
! the solution at the sub-cell boundaries.
ALLOCATE(gradUxi  (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems))
ALLOCATE(gradUeta (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems))
ALLOCATE(gradUzeta(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems))
gradUxi=0.
gradUeta=0.
gradUzeta=0.
#if PARABOLIC
! Same as gradUxi/eta/zeta, but instead of a TVD-limiter the mean value of the slopes to the
! adjacent points is used. These slopes are used to calculate the physical gradients in
! x-/y-/z-direction, which are required for the parabolic/viscous flux.
! Therefore the array size is adjusted to the number of lifting variables instead of the primitives.
! The gradients in x-/y-/z-direction are stored in the gradUx/y/z arrays of the lifting.
ALLOCATE(gradUxi_central  (PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUeta_central (PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUzeta_central(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems))
gradUxi_central  =0.
gradUeta_central =0.
gradUzeta_central=0.
#endif /* PARABOLIC */
#endif /* FV_RECONSTRUCT */

FVInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT FV DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitFV

!==================================================================================================================================
!> Interpolate face solution from DG representation to FV subcells.
!> Interpolation is done either conservatively in reference space or non-conservatively in phyiscal space.
!==================================================================================================================================
PPURE SUBROUTINE FV_InterpolateDG2FV_Face(nVar,U_In,sJ_In)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars          ,ONLY: switchConservative,FV_Vdm
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVar                               !< number of variables
REAL,INTENT(INOUT) :: U_In(nVar,0:PP_N,0:PP_NZ)          !< state vector to be switched from DG to FV representation
REAL,INTENT(IN)    :: sJ_In(0:PP_N,0:PP_NZ,0:FV_SIZE)    !< inverse of Jacobian determinant at each Gauss point
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
IF (switchConservative) THEN
  ! Transform the DG solution into the reference element
  DO q=0,PP_NZ; DO p=0,PP_N
    U_In(:,p,q)=U_In(:,p,q)/sJ_In(p,q,0)
  END DO; END DO
  ! Perform interpolation from DG to FV
  CALL ChangeBasisSurf(nVar,PP_N,PP_N,FV_Vdm,U_In)
  ! Transform back to physical space
  DO q=0,PP_NZ; DO p=0,PP_N
    U_In(:,p,q)=U_In(:,p,q)*sJ_In(p,q,1)
  END DO; END DO
ELSE
  CALL ChangeBasisSurf(nVar,PP_N,PP_N,FV_Vdm,U_In)
END IF
END SUBROUTINE FV_InterpolateDG2FV_Face

!==================================================================================================================================
!> Switch DG solution at faces between a DG element and an FV sub-cells element to Finite Volume.
!==================================================================================================================================
PPURE SUBROUTINE FV_DGtoFV(nVar,U_master,U_slave)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars          ,ONLY: FV_Elems_Sum
USE MOD_Mesh_Vars        ,ONLY: firstInnerSide,lastMPISide_MINE,nSides
USE MOD_Mesh_Vars        ,ONLY: sJ_master,sJ_slave
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVar                                   !< number of solution variables
REAL,INTENT(INOUT) :: U_master(nVar,0:PP_N,0:PP_NZ,1:nSides) !< Solution on master side
REAL,INTENT(INOUT) :: U_slave (nVar,0:PP_N,0:PP_NZ,1:nSides) !< Solution on slave side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: firstSideID,lastSideID,SideID
!==================================================================================================================================
firstSideID = firstInnerSide
lastSideID  = lastMPISide_MINE

DO SideID=firstSideID,lastSideID
  IF (FV_Elems_Sum(SideID).EQ.2) THEN
    CALL FV_InterpolateDG2FV_Face(nVar,U_master(:,:,:,SideID),sJ_master(1,:,:,SideID,:))
  ELSE IF (FV_Elems_Sum(SideID).EQ.1) THEN
    CALL FV_InterpolateDG2FV_Face(nVar,U_slave( :,:,:,SideID),sJ_slave( 1,:,:,SideID,:))
  END IF
END DO

END SUBROUTINE FV_DGtoFV

!==================================================================================================================================
!> Prim to cons for FV
!==================================================================================================================================
PPURE SUBROUTINE FV_PrimToCons(nVarPrim,nVar,UPrim_master,UPrim_slave,U_master,U_slave)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars          ,ONLY: FV_Elems_Sum
USE MOD_Mesh_Vars        ,ONLY: firstInnerSide,lastMPISide_MINE,nSides
USE MOD_EOS              ,ONLY: PrimToCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarPrim                                       !< number of solution variables
INTEGER,INTENT(IN) :: nVar                                           !< number of solution variables
REAL,INTENT(IN)    :: UPrim_master(nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< Solution on master side
REAL,INTENT(IN)    :: UPrim_slave (nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< Solution on slave side
REAL,INTENT(INOUT) :: U_master    (nVar    ,0:PP_N,0:PP_NZ,1:nSides) !< Solution on master side
REAL,INTENT(INOUT) :: U_slave     (nVar    ,0:PP_N,0:PP_NZ,1:nSides) !< Solution on slave side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: firstSideID,lastSideID,SideID
!==================================================================================================================================
firstSideID = firstInnerSide
lastSideID  = lastMPISide_MINE

DO SideID=firstSideID,lastSideID
  IF (FV_Elems_Sum(SideID).EQ.2) THEN
    CALL PrimToCons(PP_N,UPrim_master(:,:,:,SideID),U_master(:,:,:,SideID))
  ELSE IF (FV_Elems_Sum(SideID).EQ.1) THEN
    CALL PrimToCons(PP_N,UPrim_slave(:,:,:,SideID), U_slave(:,:,:,SideID))
  END IF
END DO

END SUBROUTINE FV_PrimToCons

!==================================================================================================================================
!> Cons to prim for FV
!==================================================================================================================================
PPURE SUBROUTINE FV_ConsToPrim(nVarPrim,nVar,UPrim_master,UPrim_slave,U_master,U_slave)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars          ,ONLY: FV_Elems_Sum
USE MOD_Mesh_Vars        ,ONLY: firstInnerSide,lastMPISide_MINE,nSides
USE MOD_EOS              ,ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarPrim                                       !< number of solution variables
INTEGER,INTENT(IN) :: nVar                                           !< number of solution variables
REAL,INTENT(INOUT) :: UPrim_master(nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< Solution on master side
REAL,INTENT(INOUT) :: UPrim_slave (nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< Solution on slave side
REAL,INTENT(IN)    :: U_master    (nVar    ,0:PP_N,0:PP_NZ,1:nSides) !< Solution on master side
REAL,INTENT(IN)    :: U_slave     (nVar    ,0:PP_N,0:PP_NZ,1:nSides) !< Solution on slave side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: firstSideID,lastSideID,SideID
!==================================================================================================================================
firstSideID = firstInnerSide
lastSideID  = lastMPISide_MINE

DO SideID=firstSideID,lastSideID
  IF (FV_Elems_Sum(SideID).EQ.2) THEN
    CALL ConsToPrim(PP_N,UPrim_master(:,:,:,SideID),U_master(:,:,:,SideID))
  ELSE IF (FV_Elems_Sum(SideID).EQ.1) THEN
    CALL ConsToPrim(PP_N,UPrim_slave(:,:,:,SideID), U_slave(:,:,:,SideID))
  END IF
END DO

END SUBROUTINE FV_ConsToPrim

!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeFV()
! MODULES
USE MOD_FV_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(FV_Elems)
!SDEALLOCATE(FV_Elems_master) ! moved to mesh.f90
SDEALLOCATE(FV_Elems_slave)
SDEALLOCATE(FV_Elems_Counter)
SDEALLOCATE(FV_Elems_Amount)
SDEALLOCATE(FV_Elems_Sum)
#if FV_RECONSTRUCT
SDEALLOCATE(FV_surf_gradU)
SDEALLOCATE(FV_multi_master)
SDEALLOCATE(FV_multi_slave)
SDEALLOCATE(gradUxi)
SDEALLOCATE(gradUeta)
SDEALLOCATE(gradUzeta)
#if PARABOLIC
SDEALLOCATE(gradUxi_central)
SDEALLOCATE(gradUeta_central)
SDEALLOCATE(gradUzeta_central)
#endif
#endif
#if FV_ENABLED == 2
SDEALLOCATE(FV_alpha)
SDEALLOCATE(FV_alpha_slave )
SDEALLOCATE(FV_alpha_master)
#endif

FVInitIsDone=.FALSE.
END SUBROUTINE FinalizeFV

END MODULE MOD_FV
#endif /* FV_ENABLED */
