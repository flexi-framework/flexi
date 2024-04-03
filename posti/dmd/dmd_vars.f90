!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz 
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
!> Contains global variables provided by the dmd routines
!===================================================================================================================================
MODULE MOD_DMD_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: nFiles                    !< Number of input files to perform dmd on
INTEGER             :: nVar_State                !< Number of variables in DG_Solution array
INTEGER             :: N_StateZ                  !< Polynomial degree of input state in 3rd Dimension, is set to 0 for 2D-Computation
INTEGER             :: N_StateZ_Out              !< Polynomial degree of output state in 3rd Dimension, is set to 0 for 2D-Output
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames_State(:)           !< List of varnames in state
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames_TimeAvg(:) !< List of varnames in TimeAvg-File 
INTEGER             :: N_State                   !< Polynomial degree of input state
INTEGER             :: nElems_State              !< Number of elements in state mesh
INTEGER             :: nDoFs                     !< Number of degrees of freedom of input state file
CHARACTER(LEN=255)  :: NodeType_State            !< NodeType of the input state (Gauss/Gauss-Lobatto)
CHARACTER(LEN=255)  :: MeshFile_State            !< Mesh file name of input states
REAL                :: Time_State                !< Time of the first State
REAL                :: TimeEnd_State             !< Time of the last State

!-----------------------------------------------------------------------------------------------------------------------------------
! DMD User Input Vars
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255), ALLOCATABLE  :: VarNameDMD(:)             !< Name of Variable to visualize
INTEGER                          :: nVarDMD                   !< Name of Variable to visualize
REAL                             :: SvdThreshold              !< Defines relative lower bound of singular values
INTEGER                          :: nModes                    !< Number of Modes to be visualized
LOGICAL                          :: sortFreq                  !< Decide if modes are sorted by frequency or amplitude
LOGICAL                          :: PlotSingleMode            !< Decide if a single mode is plotted
REAL                             :: ModeFreq                  !< Specify the mode frequency.
LOGICAL                          :: useBaseFlow               !< Using Basflow or not 
CHARACTER(LEN=255)               :: BaseFlow                  !< Name of the Baseflow-File
LOGICAL                          :: use2D                     !< Set T to compute the DMD on 2D data

!-----------------------------------------------------------------------------------------------------------------------------------
! DMD Vars
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                             :: dt                        !> Equidistant time interval 
REAL,ALLOCATABLE                 :: freq(:)                   !> Frequency of the DMD-Modes
COMPLEX ,ALLOCATABLE             :: Phi(:,:)                  !> State of the DMD-Modes
COMPLEX ,ALLOCATABLE             :: lambda(:)                 !> Logarithmic mapping of the eigenvales 
COMPLEX ,ALLOCATABLE             :: sigmaSort(:)              !> Eigenvalues of the DMD-Modes
COMPLEX ,ALLOCATABLE             :: alpha(:)                  !> Coefficients of the linearcombination of the data sequence
DOUBLE PRECISION,ALLOCATABLE     :: K(:,:)                    !> Snapshot-Matrix


!===================================================================================================================================
END MODULE MOD_DMD_Vars
