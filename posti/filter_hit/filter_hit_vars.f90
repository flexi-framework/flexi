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
!> Contains global variables provided by the AnalyzeHIT routines
!===================================================================================================================================
MODULE MOD_Filter_HIT_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER               :: N_Visu                    !<
INTEGER               :: N_Filter

! FFT vars
INTEGER               :: N_FFT                     !<
INTEGER               :: Nc                        !<
REAL,ALLOCATABLE      :: LocalXYZ(:,:,:,:)         !<
REAL,ALLOCATABLE      :: LocalK(:,:,:,:)           !<
INTEGER(KIND=8)       :: plan                      !<
INTEGER               :: startijk(3)
INTEGER               :: endijk(3)
INTEGER               :: startw(3)
INTEGER               :: endw(3)                   !<
COMPLEX,PARAMETER     :: II = CMPLX(0.,1.0)

! state file variables
LOGICAL               :: FieldDataExists=.FALSE.
REAL                  :: Time_HDF5                 !<
INTEGER               :: nVar_HDF5
INTEGER               :: nVarField_HDF5
INTEGER               :: N_HDF5
INTEGER               :: nElems_HDF5
CHARACTER(LEN=255)    :: NodeType_HDF5
CHARACTER(LEN=255)    :: MeshFile_HDF5
CHARACTER(LEN=255)    :: ProjectName_HDF5
!===================================================================================================================================
END MODULE MOD_Filter_HIT_Vars
