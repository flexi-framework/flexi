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
MODULE MOD_OutputRPVisu_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: nSamples_out                !< number of visualisation points is NVisu+1
INTEGER                       :: nCoords                     !< number of visualisation points is NVisu+1
CHARACTER(LEN=255),ALLOCATABLE:: CoordNames(:)               !< including space and time coordinates
REAL,ALLOCATABLE              :: RPData_out(:,:,:)           !< output state
REAL,ALLOCATABLE              :: RPDataTimeAvg_out(:,:)      !< time average state
REAL,ALLOCATABLE              :: RPDataRMS_out(:,:)          !< RMS of fluctuations
LOGICAL                       :: OutputInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_OutputRPVisu_Vars
