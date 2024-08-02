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
!> Contains global variables used by the baseflow.
!==================================================================================================================================
MODULE MOD_Baseflow_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                 :: doBaseflow       = .FALSE.      !< Switch on to compute baseflow of conservative mean solution
LOGICAL                 :: doBaseflowRMS    = .FALSE.      !< Switch on to compute baseflow of Reynolds stresses
LOGICAL                 :: initBaseflowDone = .FALSE.      !< Logical indicating if initbaseflow() has been called
CHARACTER(LEN=255)      :: BaseFlowFile                    !< Filename of baseflow file
REAL                    :: BaseFlowTime                    !< Time stamp of baseflow file

REAL,ALLOCATABLE,TARGET :: BaseFlowPrim(:,:,:,:,:)         !< Local tmp array to compute baseflow in primitive variabels
REAL,ALLOCATABLE,TARGET :: BaseFlowPrimFiltered(:,:,:,:,:) !< Local tmp array to compute baseflow in primitive variabels
REAL,ALLOCATABLE,TARGET :: BaseFlow(:,:,:,:,:)             !< Local baseflow depending on local filter width
REAL,ALLOCATABLE,TARGET :: BaseFlowFiltered(:,:,:,:,:)     !< Local baseflow depending on local filter width and Filter
REAL,ALLOCATABLE,TARGET :: BaseFlowRMS(:,:,:,:,:)          !< Local baseflow of Reynolds stresses

REAL,ALLOCATABLE        :: TimeFilterWidthBaseflow(:)      !< Array that contains the element local temporal filter width
REAL,ALLOCATABLE        :: fac(:)                          !< Array that contains the element local exponential time filter factor

! Selective filter of baseflow
LOGICAL,DIMENSION(3)    :: doSelectiveFilter=.FALSE.       !< Logical for Filter
INTEGER,DIMENSION(3)    :: SelectiveFilter                 !< Polynomial degree of selective filter in each direction
REAL,ALLOCATABLE        :: SelectiveFilterMatrix(:,:,:)    !< Filter Matrix of the selective Filter
!==================================================================================================================================
END MODULE MOD_Baseflow_Vars
