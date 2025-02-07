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
!> Variables needed for the evaluation of the record points
!==================================================================================================================================
MODULE MOD_RecordPoints_Vars
! MODULES
#if USE_MPI
USE mpi_f08
#endif

IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255) :: RPDefFile               !< File containing element-local parametric recordpoint coordinates and structure
LOGICAL            :: RecordPointsInitIsDone = .FALSE. !< mark wheter recordpoints init routine is finished
LOGICAL            :: RP_inUse  = .FALSE.     !< mark whether recordpoints should be evaluated during computation
LOGICAL            :: RP_onProc = .FALSE.     !< marks wheter current proc has RPs
LOGICAL            :: RP_fileExists = .FALSE. !< flag if RP file for analyze level has been created
INTEGER            :: RP_Buffersize           !< no. of time samples (size of RP_Data)
INTEGER            :: RP_MaxBuffersize        !< max. allowed no. of time samples
INTEGER            :: RP_SamplingOffset       !< sampling rate (each .. iterations)
INTEGER            :: nRP                     !< no. of RP on proc
INTEGER            :: nGlobalRP               !< total no. of RP
INTEGER            :: offsetRP                !< offset for each proc in global RP list
INTEGER            :: iSample=0               !< no of samples in array
INTEGER            :: nSamples=0              !< total no. samples in case of multiple io steps
INTEGER            :: chunkSamples=0          !< time samples per chunk for IO (first iSample in file)
INTEGER,ALLOCATABLE:: RP_ElemID(:)            !< mapping from RP->Elem (nRP)
#if FV_ENABLED
INTEGER,ALLOCATABLE:: FV_RP_ijk(:,:)          !< ijk-index of FV subcell nearest to record point [1:3,nRP]
#endif
REAL,ALLOCATABLE   :: L_xi_RP(:,:)            !< Lagrange basis evaluated at RP coords (xi-dir)
REAL,ALLOCATABLE   :: L_eta_RP(:,:)           !< Lagrange basis evaluated at RP coords (eta-dir)
REAL,ALLOCATABLE   :: L_zeta_RP(:,:)          !< Lagrange basis evaluated at RP coords (zeta-dir)
REAL,ALLOCATABLE   :: RP_Data(:,:,:)          !< solution evaluated at RPs (nvar,nRP,nSamples)
REAL,ALLOCATABLE   :: lastSample(:,:)         !< solution evaluated at RPs (nvar,nRP,nSamples)
CHARACTER(LEN=255) :: StrVarNames(PP_nVar)    !< RP variables names for output

!----------------------------------------------------------------------------------------------------------------------------------
! MPI Communicator for RPs
!----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
INTEGER            :: myRPrank                !< rank within RP communicator
INTEGER            :: nRP_Procs               !< number of procs with RPs
TYPE(MPI_Comm)     :: RP_COMM=MPI_COMM_NULL   !< MPI RP communicator
#endif /* USE_MPI */

END MODULE MOD_recordPoints_Vars
