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

CHARACTER(LEN=255)    :: MeshFile                  !< mesh file
INTEGER               :: InitSpec                  !<
INTEGER               :: N_FFT                     !<
INTEGER               :: N_Visu                    !<
INTEGER               :: Nc                        !<
REAL,ALLOCATABLE      :: U_HDF5(:,:,:,:,:)
REAL,ALLOCATABLE      :: Uloc(:,:,:,:)             !<
COMPLEX,ALLOCATABLE   :: Uloc_c(:,:,:,:)           !<
COMPLEX,ALLOCATABLE   :: U_j(:,:,:,:)              !<
COMPLEX,ALLOCATABLE   :: U_k(:,:,:,:)              !<
COMPLEX,ALLOCATABLE   :: U_FFT(:,:,:,:)            !<
REAL,ALLOCATABLE      :: LocalXYZ(:,:,:,:)         !<
REAL,ALLOCATABLE      :: LocalK(:,:,:,:)           !<
COMPLEX, ALLOCATABLE  :: F_vv(:,:,:,:,:)           !<
COMPLEX, ALLOCATABLE  :: fhat(:,:,:,:)             !<
COMPLEX, ALLOCATABLE  :: phat(:,:,:)               !<
COMPLEX               :: II
REAL                  :: TwoPi                     !<
REAL                  :: Pi
REAL                  :: scalefactor               !<
REAL                  :: maxTol                    !<
REAL                  :: abortTol                  !<
INTEGER               :: EndW(3)                   !<
INTEGER(KIND=8)       :: plan                      !<

INTEGER               :: nVar_HDF5
INTEGER               :: N_HDF5
INTEGER               :: nElems_HDF5
CHARACTER(LEN=255)    :: NodeType_HDF5
CHARACTER(LEN=255)    :: MeshFile_HDF5
CHARACTER(LEN=255)    :: Project_Name

REAL                  :: Time_HDF5                     !<
INTEGER               :: N_FA(3)
INTEGER               :: Nunder

LOGICAL               :: validHDF5

INTEGER               :: startijk(3)
INTEGER               :: startw(3)
INTEGER               :: endijk(3)

!analyze variables
LOGICAL   :: AnalyzeInitIsDone=.FALSE.
LOGICAL   :: CalcErrorNorms=.FALSE.
INTEGER   :: kmax
REAL, ALLOCATABLE :: E_k(:),T_k(:),MuSGS_K(:),k_shell_glob(:),numDiss_k(:),eps_k(:),eMean_k(:)
REAL,ALLOCATABLE        :: MuSGS_Avg(:)
INTEGER   :: FileUnit_HIT
CHARACTER(LEN=255)::Filename_HIT
INTEGER   :: N_Filter

REAL,ALLOCATABLE        :: U_Global(:,:,:,:)
!===================================================================================================================================
END MODULE MOD_Filter_HIT_Vars
