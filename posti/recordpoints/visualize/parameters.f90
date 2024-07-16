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
!===================================================================================================================================
!> Contains global variables provided by the visualize recordpoints tool
!===================================================================================================================================
MODULE MOD_ParametersVisu
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                             :: justVisualizeState      !< If no output variables have been specified, visualize everything
                                                               !< that is found in the HDF5 file
! General input parameters
CHARACTER(len=255)                  :: ProjectName             !<  Name of the project
CHARACTER(len=255)                  :: RP_DefFile              !< Path to the *RPset.h5 file
CHARACTER(len=255),ALLOCATABLE      :: GroupNames_visu(:)      !< Name(s) of the group(s) to visualize
INTEGER                             :: nGroups_visu            !< Number of groups to visualize
LOGICAL                             :: OutputTimeData          !< Should the time series be written?
LOGICAL                             :: doFluctuations          !< Should the fluctuations be computed and written?
LOGICAL                             :: equiTimeSpacing         !< Interpolate the temporal data to equidistant time steps (needed for
                                                               !< FFT)
LOGICAL                             :: OutputTimeAverage       !< Should the time average be computed and written?
LOGICAL                             :: OutputLines             !< General option to turn off the output of lines
LOGICAL                             :: OutputPlanes            !< General option to turn off the output of planes
LOGICAL                             :: OutputPoints            !< General option to turn off the output of points
LOGICAL                             :: OutputBoxes             !< General option to turn off the output of boxes
LOGICAL                             :: Line_LocalCoords        !< Set to use local instead of global coordinates along planes
LOGICAL                             :: Line_LocalVel           !< Set to use local instead of global velocities along planes
LOGICAL                             :: Plane_LocalCoords       !< Set to use local instead of global coordinates along planes
LOGICAL                             :: Plane_LocalVel          !< Set to use local instead of global velocities along planes
LOGICAL                             :: Plane_doBLProps         !< Set to use local instead of global velocities along planes
INTEGER                             :: Plane_BLvelScaling      !< 0: no scaling, 1: laminar, 2: turbulent
LOGICAL                             :: Box_LocalCoords         !< Set to use local instead of global coordinates along boxes
LOGICAL                             :: Box_LocalVel            !< Set to use local instead of global velocities along boxes
LOGICAL                             :: Box_doBLProps           !< Set to use local instead of global velocities along boxes
INTEGER                             :: Box_BLvelScaling        !< 0: no scaling, 1: laminar, 2: turbulent
LOGICAL                             :: usePrims                !< State directly gives the primitive variables
LOGICAL                             :: RP_SET_defined          !< True if definition file has been specified
!--------------------------------------------------
! Filter
LOGICAL                             :: doFilter                !< Set to perform temporal filtering for each RP
INTEGER                             :: FilterMode              !< Set to 0 for low pass filter and to 1 for high pass filter
REAL                                :: FilterWidth             !< Width of the temporal filter
!--------------------------------------------------
! Ensemble Averaging
LOGICAL                             :: doEnsemble              !< Set to perform ensemble averaging for each RP
REAL                                :: EnsemblePeriod          !< Period used for ensemble averaging
REAL                                :: EnsembleFreq            !< Frequency used for ensemble averaging
REAL                                :: Kappa                   !< heat capacity ratio / isentropic exponent
!--------------------------------------------------
! Spectral Analysis, FFT, PSD
LOGICAL                             :: doSpec                  !< Set if any spectral analysis (FFT,PSD) is performed
LOGICAL                             :: doPSD                   !< Calculate the power spectral density of the time signal
LOGICAL                             :: doFFT                   !< Calculate a fast Fourier transform of the time signal
INTEGER                             :: nBlocks                 !< Number of blocks for spectral averaging
INTEGER                             :: BlockSize               !< Size of blocks in samples
REAL                                :: cutoffFreq              !< Cutoff frequency of spectral analysis
REAL                                :: samplingFreq            !< Alternative of specfying number of blocks
LOGICAL                             :: doHanning               !< Perform windowing using Hann function
LOGICAL                             :: fourthDeriv             !< Calculate fourth temporal derivative
LOGICAL                             :: ThirdOct                !<
REAL                                :: u_infPhys               !<
REAL                                :: chordPhys               !<
REAL                                :: Line_LocalVel_vec(3)    !< Vector used in transformation to local velocity on line
!--------------------------------------------------
! Turbulence
LOGICAL                             :: doTurb                  !< Compute temporal FFT and calculate turbulent quantities
REAL                                :: Mu0                     !< Viscosity (for turbulent quantitites)
!--------------------------------------------------
! Output
INTEGER                             :: OutputFormat            !< 0: ParaView (vts), 2: HDF5
INTEGER                             :: Skip                    !< Skip every nth recorded time step
!--------------------------------------------------
! Variable mappings
INTEGER,ALLOCATABLE                   :: DepTable(:,:)         !< Array containing all available variables and the dependecies
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarNamesAll(:)        !< Name of all availabe variables
CHARACTER(len=255),ALLOCATABLE        :: VarNameVisu(:)        !< Name of variables to visualize
INTEGER                               :: nVarVisu              !< Number of quantities to visualize
INTEGER                               :: nVarDep               !< Number of dependent variables
INTEGER                               :: nVarCalc              !< Number of variables that must be calculated
INTEGER,ALLOCATABLE                   :: mapCalc(:)            !< Mapping from dependent to calculated variables
INTEGER,ALLOCATABLE                   :: mapVisu(:)            !< Mapping from dependent to visualized variables
!===================================================================================================================================
END MODULE MOD_ParametersVisu
