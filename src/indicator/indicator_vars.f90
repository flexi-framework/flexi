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

!==================================================================================================================================
!> Contains variables relevant for indicators
!==================================================================================================================================
MODULE MOD_Indicator_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                :: IndicatorInitIsDone=.FALSE.
INTEGER                :: IndicatorType               !< Type of indicator used
INTEGER                :: nIndVars                    !< number of variables to be checked by the indicator
INTEGER,ALLOCATABLE    :: IndVar(:)                   !< variable on which indicator is applied (only cons)
INTEGER                :: nModes                      !< number of modes to check for Persson modal indicator
REAL,ALLOCATABLE       :: IndValue(:)                 !< indicator output
REAL                   :: IndStartTime                !< specify starting time of indicator evaluation
LOGICAL                :: FVBoundaries = .FALSE.      !< specify if FV element is set at boundaries
INTEGER,ALLOCATABLE    :: FVBoundaryType(:)           !< select at which kind of BCs FV elements should be used
#if FV_ENABLED == 2 || FV_ENABLED == 3
REAL                   :: T_FV                        !< Threshold for FV blending as function of PP_N
REAL,PARAMETER         :: s_FV = 9.21024              !< "Sharpness factor" for FV blending function according to
                                                      !< Eq.(45) in: S. Hennemann et al., J.Comp.Phy., 2021
REAL                   :: sdT_FV                      !< Precomputed variable for FV blending
#endif
#if FV_ENABLED == 3
REAL                   :: SobelFilterMatrix(1:2,1:3,1:3)
REAL,ALLOCATABLE       :: VdmGaussEqui(:,:)
#endif
!==================================================================================================================================
END MODULE MOD_Indicator_Vars
