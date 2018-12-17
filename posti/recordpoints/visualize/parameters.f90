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
!===================================================================================================================================
!> Contains global variables provided by the output routines
!===================================================================================================================================
MODULE MOD_ParametersVisu
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                             :: justVisualizeState
CHARACTER(len=255)                  :: ProjectName
CHARACTER(len=255)                  :: RP_DefFile 
CHARACTER(len=255),ALLOCATABLE      :: GroupNames_visu(:)
INTEGER                             :: nGroups_visu
LOGICAL                             :: OutputTimeData
LOGICAL                             :: doFluctuations
LOGICAL                             :: equiTimeSpacing
LOGICAL                             :: OutputTimeAverage
LOGICAL                             :: OutputLines,OutputPlanes,OutputPoints
LOGICAL                             :: Line_GlobalCoords
LOGICAL                             :: Line_LocalCoords
LOGICAL                             :: Line_LocalVel
LOGICAL                             :: Plane_LocalCoords
LOGICAL                             :: Plane_LocalVel
LOGICAL                             :: Plane_doBLProps
INTEGER                             :: Plane_BLvelScaling
LOGICAL                             :: usePrims
LOGICAL                             :: RP_SET_defined
!--------------------------------------------------
! Filter
LOGICAL                             :: doFilter
INTEGER                             :: FilterMode
REAL                                :: FilterWidth
!--------------------------------------------------
! Spectral Analysis, FFT, PSD
LOGICAL                             :: doSpec
LOGICAL                             :: doPSD
LOGICAL                             :: doFFT
INTEGER                             :: nBlocks 
INTEGER                             :: BlockSize
REAL                                :: cutoffFreq,samplingFreq
LOGICAL                             :: doHanning
LOGICAL                             :: fourthDeriv,ThirdOct
REAL                                :: u_infPhys,chordPhys
REAL                                :: Line_LocalVel_vec(3)
REAL                                :: Mu0
!--------------------------------------------------
! Turbulence
LOGICAL                             :: doTurb
INTEGER                             :: nVarVisu  
CHARACTER(len=255),ALLOCATABLE      :: VarNamevisu(:)
INTEGER                             :: OutputFormat
INTEGER                             :: Skip ! nur jeder skipte RP sample wird eingelesen



INTEGER,ALLOCATABLE               :: DepTable(:,:)
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarNamesAll(:)
INTEGER                           :: nVarDep                 ! 
INTEGER                           :: nVarCalc
INTEGER                           :: nVarVisuTotal
INTEGER,ALLOCATABLE               :: mapCalc(:)
INTEGER,ALLOCATABLE               :: mapVisu(:)
!===================================================================================================================================
END MODULE MOD_ParametersVisu
