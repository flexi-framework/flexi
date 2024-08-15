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
!==================================================================================================================================
!> Contains global variables used by the Avg2D module
!==================================================================================================================================
MODULE MOD_Avg2D_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                :: doAvg2D                    !< Do average data in spatial domain
REAL,ALLOCATABLE       :: RecvBufferAvg2D(:,:,:,:,:) !< Reciving Buffers for MPI communication
REAL,ALLOCATABLE       :: SendBufferAvg2D(:,:,:,:,:) !< Send Buffers for MPI communication
REAL,ALLOCATABLE       :: UAvg2D(:,:,:,:,:)          !< 2D averaged solution
REAL,ALLOCATABLE       :: UAvg2DLocal(:,:,:,:)       !< Locally 2D averaged U
REAL,ALLOCATABLE       :: UAvg2DGlobal(:,:,:,:)      !< Globally 2D averaged U
INTEGER                :: minIJ(2)                   !< Min i and j global at this processor
INTEGER                :: nVarsAvg2D                 !< Number of variables to be averaged in 2D
INTEGER,ALLOCATABLE    :: RecvProcs(:)               !< List of Master Processors that receive elements from this processor (Slave)
INTEGER,ALLOCATABLE    :: SendProcs(:)               !< List of Slave Processors that send elements to this Processor (Master)
INTEGER,ALLOCATABLE    :: nElemsToRecv(:)            !< List of number of elements to be received from different slave procs
INTEGER,ALLOCATABLE    :: nElemsToSend(:)            !< List of number of elements to be send to different master procs
INTEGER                :: nSendProcs,nRecvProcs      !< Number of sending and receiving procs
INTEGER,ALLOCATABLE    :: ElemsToSendIJSorted(:,:)   !< List of Elements to be sent sorted i,j ascending order (RecvProc,nElems2S)
INTEGER,ALLOCATABLE    :: ElemsToRecvIJSorted(:,:)   !< List of Elements to be receive sorted i,j ascending order (RecvProc,nElems2S)
REAL                   :: SpanWidth                  !< Global span width
INTEGER,ALLOCATABLE    :: MyElemsIJ(:,:,:)           !< List of elems at each i,j to sort back
INTEGER,ALLOCATABLE    :: nElemsIJ(:,:)              !< Number of elmes in that are on a specific i,j
INTEGER,ALLOCATABLE    :: iDiffElem(:)               !< Mapping of iElem to different Elements in i,j plane

INTEGER                :: Avg2DDir                   !< Direction in which the averaging is performed
INTEGER                :: IJK_Mask(2)
!==================================================================================================================================
END MODULE MOD_Avg2D_Vars
