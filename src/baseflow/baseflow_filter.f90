!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
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
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Subroutines needed for the general base flow based on a moving time average of the instationary flow field, also known as Pruett
!> damping. See "The temporally filtered Navier–Stokes equations: Properties of the residual stress" for details.
!==================================================================================================================================
MODULE MOD_BaseFlow_Filter
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC :: InitBaseFlowFilter
PUBLIC :: BaseFlowFilter
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Perform initialization of baseflow selective filter routine
!==================================================================================================================================
SUBROUTINE InitBaseFlowFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_BaseFlow_Vars,      ONLY: doSelectiveFilter,SelectiveFilter,SelectiveFilterMatrix
USE MOD_Interpolation_Vars, ONLY: Vdm_Leg,sVdm_Leg
USE MOD_ReadInTools,        ONLY: GETINTARRAY
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j
!==================================================================================================================================
!Setup a FilterMatrix if necessary
SelectiveFilter(:) = GETINTARRAY('SelectiveFilter',3)

#if PP_dim==2
SWRITE(UNIT_stdOut,'(A)') "Ignoring third dimension of SelectiveFilter vector as simulation is twodimensional!"
#endif

DO i=1,PP_dim
  IF (SelectiveFilter(i) .EQ. -999) THEN
    SelectiveFilter(i) = PP_N
  END IF
  IF (SelectiveFilter(i) .NE. PP_N) THEN
    doSelectiveFilter(i)=.TRUE.
  END IF
END DO

IF (ANY(doSelectiveFilter)) THEN
  ALLOCATE(SelectiveFilterMatrix(PP_dim,0:PP_N,0:PP_N))
  SelectiveFilterMatrix=0.
  DO i=1,PP_dim
    DO j=0,SelectiveFilter(i)
      SelectiveFilterMatrix(i,j,j)=1.
    END DO
    SelectiveFilterMatrix(i,:,:) = MATMUL(MATMUL(Vdm_Leg,SelectiveFilterMatrix(i,:,:)),sVdm_Leg)
  END DO
END IF

END SUBROUTINE InitBaseFlowFilter

!==================================================================================================================================
!> Perform selective filter of baseflow
!==================================================================================================================================
SUBROUTINE BaseFlowFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_BaseFlow_Vars,      ONLY: doSelectiveFilter,SelectiveFilterMatrix,BaseFlowFiltered
USE MOD_Filter,             ONLY: Filter_Selective
USE MOD_Mesh_Vars,          ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iElem
!==================================================================================================================================
IF(ANY(doSelectiveFilter)) THEN
    DO i=1,PP_dim
      DO iElem=1,nElems
        CALL Filter_Selective(PP_nVar,SelectiveFilterMatrix(i,:,:),BaseFlowFiltered(:,:,:,:,iElem),doSelectiveFilter)
      END DO ! iElem
    END DO ! PP_dim
END IF

END SUBROUTINE BaseFlowFilter

END MODULE MOD_BaseFlow_Filter
