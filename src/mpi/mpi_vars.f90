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
!> Contains the variables that are used to control non-blocking communication
!==================================================================================================================================
MODULE MOD_MPI_Vars
#if USE_MPI
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: MPIRequest_U(:,:)        !< communication handle for the surface solution
INTEGER,ALLOCATABLE :: MPIRequest_Flux(:,:)     !< communication handle for the surface flux
INTEGER,ALLOCATABLE :: MPIRequest_Avg2dSend(:)  !< communication handle for averaging in 2D
INTEGER,ALLOCATABLE :: MPIRequest_Avg2DRecv(:)  !< communication handle for averaging in 2D
#if FV_ENABLED
INTEGER,ALLOCATABLE :: MPIRequest_FV_Elems(:,:) !< communication handle for the FV_Elems array
INTEGER,ALLOCATABLE :: MPIRequest_FV_gradU(:,:) !< communication handle for the slopes of the FV reconstruction
#if ((FV_ENABLED == 2) && (PP_NodeType == 1))
INTEGER,ALLOCATABLE :: MPIRequest_FV_U(:,:)     !< communication handle for the surface solution (pureFV)
INTEGER,ALLOCATABLE :: MPIRequest_FV_Flux(:,:)  !< communication handle for the surface flux     (pureFV)
#endif
#if FV_RECONSTRUCT
INTEGER,ALLOCATABLE :: MPIRequest_Rec_SM(:,:)   !< communication handle for extended primitive solution for preconditioner
INTEGER,ALLOCATABLE :: MPIRequest_Rec_MS(:,:)   !< communication handle for extended primitive solution for preconditioner
#endif
#endif
#if EDDYVISCOSITY
INTEGER,ALLOCATABLE :: MPIRequest_SGS(:,:)      !< communication handle for the SGS Model Indicator
#endif
#if PARABOLIC
INTEGER,ALLOCATABLE :: MPIRequest_gradU(:,:,:)  !< communication handle for the surface gradients
#endif /*PARABOLIC*/
INTEGER             :: nSendVal                 !< number of values to be sent
INTEGER             :: nRecVal                  !< number of values to be received
INTEGER             :: DataSizeSide             !< datasize of one face, =PP_nVar*(PP_N+1)**2
INTEGER             :: DataSizeSidePrim         !< datasize of one face for (primitive) gradients, =PP_nVarPrim*(PP_N+1)**2
INTEGER             :: DataSizeSideGrad         !< datasize of one face for one value, =1*(PP_N+1)**2 (in case of less liftet)
INTEGER             :: DataSizeSideSGS          !< datasize of one face for one value, =1*(PP_N+1)**2

INTEGER             :: SideID_start,SideID_end
INTEGER             :: nNbProcs                 !< number of neighbor procs, is set in ReadMesh
INTEGER,ALLOCATABLE :: NbProc(:)                !< list of neighbor procs; allocated from 1:nNbProcs, is set in ReadMesh
INTEGER,ALLOCATABLE :: nMPISides_Proc(:)        !< number of mpisides for all neighbor procs, is set in ReadMesh
INTEGER,ALLOCATABLE :: nMPISides_MINE_Proc(:)   !< number of MINE mpisides for all neighbor procs, is set in setLocalSideIDs
                                                !< (prepare_mesh.f90)
INTEGER,ALLOCATABLE :: nMPISides_YOUR_Proc(:)   !< number of YOUR mpisides for all neighbor procs, is set in setLocalSideIDs
                                                !< (prepare_mesh.f90)

INTEGER,ALLOCATABLE :: offsetMPISides_MINE(:)   !< gives position of send/recv block in *_MINE arrays,allocated from 0:nNbProcs, is
                                                !< set in setLocalSideIDs (prepare_mesh.f90)
INTEGER,ALLOCATABLE :: offsetMPISides_YOUR(:)   !< gives position of send/recv block in *_YOUR arrays,allocated from 0:nNbProcs,
                                                !< is set in setLocalSideIDs (prepare_mesh.f90)
INTEGER,ALLOCATABLE :: offsetElemMPI(:)         !< gives offsetposotion of elements of all procs, allocated from 0:nProcessors set
                                                !< in ReadMesh
INTEGER,ALLOCATABLE :: nMPISides_send(:,:)      !< number of sides to send, (1:nNbProcs,1:2), last index: 1: SEND MINE, 2: SEND YOUR
INTEGER,ALLOCATABLE :: nMPISides_rec(:,:)       !< number of sides to receive, (1:nNbProcs,1:2), last index: 1: RECEIVE YOUR,
                                                !< 2: RECEIVE MINE
INTEGER,ALLOCATABLE :: OffsetMPISides_send(:,:) !< offset of sides to send,(1:nNbProcs,1:2), last index: 1: SEND MINE, 2: SEND YOUR
INTEGER,ALLOCATABLE :: OffsetMPISides_rec(:,:)  !< offset of sides to receive, (1:nNbProcs,1:2), last index: 1: RECEIVE YOUR,
                                                !< 2: RECEIVE MINE
!==================================================================================================================================
#endif
END MODULE MOD_MPI_Vars
