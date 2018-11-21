!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames_State(:)           !< List of varnames in state
INTEGER             :: N_State                   !< Polynomial degree of input state
INTEGER             :: nElems_State              !< Number of elements in state mesh
INTEGER             :: nDoFs                     !< Number of degrees of freedom of input state file
CHARACTER(LEN=255)  :: NodeType_State            !< NodeType of the input state (Gauss/Gauss-Lobatto)
CHARACTER(LEN=255)  :: MeshFile_State            !< Mesh file name of input states
REAL                :: Time_State
REAL                :: TimeEnd_State

!-----------------------------------------------------------------------------------------------------------------------------------
! DMD User Input Vars
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)  :: VarNameDMD                !< Name of Variable to visualize
REAL                :: SvdThreshold              !< Defines relative lower bound of singular values
INTEGER             :: nModes                    !< Number of Modes to be visualized

!-----------------------------------------------------------------------------------------------------------------------------------
! DMD Vars
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                    :: dt                    !>
REAL,ALLOCATABLE        :: freq(:)               !> 
COMPLEX ,ALLOCATABLE    :: Phi(:,:)              !> 
COMPLEX ,ALLOCATABLE    :: lambda(:)             !> 
COMPLEX ,ALLOCATABLE    :: sigmaSort(:)        !> 
COMPLEX ,ALLOCATABLE    :: alpha(:)        !>
DOUBLE PRECISION,ALLOCATABLE    :: K(:,:)        !> 


!===================================================================================================================================
END MODULE MOD_DMD_Vars
