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
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Subroutines needed for the general base flow based on a moving time average of the instationary flow field, also known as Pruett
!> damping. See "The temporally filtered Navier–Stokes equations: Properties of the residual stress" for details.
!==================================================================================================================================
MODULE MOD_Baseflow_Filter
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitBaseflowFilter
  MODULE PROCEDURE InitBaseflowFilter
END INTERFACE

INTERFACE BaseflowFilter
  MODULE PROCEDURE BaseflowFilter
END INTERFACE

PUBLIC :: InitBaseflowFilter
PUBLIC :: BaseflowFilter
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Perform init of baseflow selective filter routine
!==================================================================================================================================
SUBROUTINE InitBaseflowFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Baseflow_Vars,      ONLY: doSelectiveFilter,SelectiveFilter,SelectiveFilterMatrix
USE MOD_ReadInTools,        ONLY: GETINTARRAY
USE MOD_Interpolation_Vars, ONLY: Vdm_Leg,sVdm_Leg
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j
!==================================================================================================================================
!Setup a FilterMatrix if necessary
SelectiveFilter(:) = GETINTARRAY('SelectiveFilter',3)
DO i=1,3
  IF (SelectiveFilter(i) .EQ. -999) THEN
    SelectiveFilter(i) = PP_N
  END IF
  IF (SelectiveFilter(i) .NE. PP_N) THEN
    doSelectiveFilter(i)=.TRUE.
  END IF
END DO
IF (ANY(doSelectiveFilter)) THEN
  ALLOCATE(SelectiveFilterMatrix(3,0:PP_N,0:PP_N))
  SelectiveFilterMatrix=0.
  DO i=1,3
    DO j=0,SelectiveFilter(i)
      SelectiveFilterMatrix(i,j,j)=1.
    END DO
    SelectiveFilterMatrix(i,:,:) = MATMUL(MATMUL(Vdm_Leg,SelectiveFilterMatrix(i,:,:)),sVdm_Leg)
  END DO
END IF
END SUBROUTINE InitBaseflowFilter

!==================================================================================================================================
!> Perform selective filter of baseflow
!==================================================================================================================================
SUBROUTINE BaseflowFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Baseflow_Vars,      ONLY: doSelectiveFilter,SelectiveFilterMatrix,BaseflowFiltered
#if EQNSYSNR == 2 /* NAVIER-STOKES */
USE MOD_Baseflow_Vars,      ONLY: doBaseFlowRMS
#endif /* NAVIER-STOKES */
USE MOD_ReadInTools,        ONLY: GETINTARRAY
USE MOD_Filter,             ONLY: Filter_Selective
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iElem
!==================================================================================================================================
IF(ANY(doSelectiveFilter)) THEN
#if EQNSYSNR == 2 /* NAVIER-STOKES */
  ! Compute Baseflow of RMS values
  IF(doBaseFlowRMS) THEN
    DO i=1,3
      DO iElem=1,nElems
        CALL Filter_Selective(PP_nVar+PP_nVarRMS,SelectiveFilterMatrix(i,:,:),BaseFlowFiltered(:,:,:,:,iElem),doSelectiveFilter)
      END DO ! iElem
    END DO
  ELSE
#endif /* NAVIER-STOKES */
    DO i=1,3
      DO iElem=1,nElems
        CALL Filter_Selective(PP_nVar,SelectiveFilterMatrix(i,:,:),BaseFlowFiltered(:,:,:,:,iElem),doSelectiveFilter)
      END DO ! iElem
    END DO
#if EQNSYSNR == 2 /* NAVIER-STOKES */
  END IF
#endif /* NAVIER-STOKES */
END IF
END SUBROUTINE BaseflowFilter

END MODULE MOD_Baseflow_Filter
