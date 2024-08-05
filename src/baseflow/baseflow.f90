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
MODULE MOD_Baseflow
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersBaseflow
  MODULE PROCEDURE DefineParametersBaseflow
END INTERFACE

INTERFACE InitBaseflow
  MODULE PROCEDURE InitBaseflow
END INTERFACE

INTERFACE UpdateBaseflow
  MODULE PROCEDURE UpdateBaseflow
END INTERFACE

INTERFACE FinalizeBaseflow
  MODULE PROCEDURE FinalizeBaseflow
END INTERFACE

PUBLIC :: DefineParametersBaseflow
PUBLIC :: InitBaseflow
PUBLIC :: UpdateBaseflow
PUBLIC :: FinalizeBaseflow
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters of basflow to compute a moving time-average during simulation runtime
!==================================================================================================================================
SUBROUTINE DefineParametersBaseflow()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Baseflow")
CALL prms%CreateLogicalOption( 'doBaseflow'             ,"Switch on to calculate a baseflow."                                     &
                                                        ,'.FALSE.')
CALL prms%CreateStringOption(  'BaseFlowFile'           ,"FLEXI file (e.g. baseflow, TimeAvg) from which baseflow is read."       &
                                                        ,'none')
CALL prms%CreateIntOption(     'BaseFlowRefState'       ,"Specify which refstate should be used in no baseflowfile is given.")
CALL prms%CreateIntArrayOption('SelectiveFilter'        ,"Filter Mean to another polynomial degree.",'(/-999,-999,-999/)')
CALL prms%CreateRealOption(    'TimeFilterWidthBaseflow',"Temporal filter width of exponential, explicit time filter.")
END SUBROUTINE DefineParametersBaseflow


!==================================================================================================================================
!> Perform init of baseflow
!==================================================================================================================================
SUBROUTINE InitBaseflow()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Baseflow_Vars
USE MOD_Baseflow_Filter,    ONLY: InitBaseflowFilter,BaseflowFilter
USE MOD_Baseflow_Readin,    ONLY: ReadBaseFlow
USE MOD_Equation_Vars,      ONLY: RefStateCons
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_ReadInTools,        ONLY: GETSTR,GETINT,GETREAL,GETINTARRAY,GETLOGICAL
USE MOD_Restart_Vars,       ONLY: doRestart,RestartTime
USE MOD_StringTools,        ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iElem,BaseFlowRefState
CHARACTER(LEN=255) :: FileName
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT BASEFLOW ...'

! Check if baseflow is switched on, otherwise throw a notification
IF (.NOT. doBaseflow)    doBaseflow    = GETLOGICAL('doBaseflow')

IF(.NOT. doBaseflow) THEN
  SWRITE(UNIT_stdOut,'(A)') ' | Baseflow is not computed. Switch on if desired ...'
  SWRITE(UNIT_stdOut,'(A)')' INIT BASEFLOW DONE!'
  SWRITE(UNIT_stdOut,'(132("-"))')
  RETURN
END IF

ALLOCATE(TimeFilterWidthBaseflow(nElems))
ALLOCATE(fac(nElems))
fac                     = 0.
TimeFilterWidthBaseflow = 1./GETREAL("TimeFilterWidthBaseflow")
BaseFlowFile            = GETSTR('BaseFlowFile')

! Check if the previous baseflow file is available
IF (doRestart .AND. STRICMP(TRIM(BaseflowFile),'none')) THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BaseFlow',RestartTime))//'.h5'
  SWRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO') ' | Searching for baseflow file "',TRIM(FileName),'" ...'
  IF (FILEEXISTS(FileName)) THEN
    BaseFlowFile = FileName
    SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' SUCCESS!'
  ELSE
    SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' FAILED!'
  END IF
END IF

ALLOCATE(BaseFlow(        PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(BaseFlowFiltered(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

BaseFlow         = 0.
BaseFlowFiltered = 0.

IF (.NOT. (STRICMP(TRIM(BaseflowFile),'none'))) THEN
  CALL ReadBaseFlow(BaseFlowFile)
ELSE
  BaseFlowRefState = GETINT('BaseflowRefState')
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      BaseFlow(:,i,j,k,iElem) = RefStateCons(:,BaseFlowRefState)
    END DO; END DO; END DO
  END DO
END IF

! Filtering of Baseflow
BaseFlowFiltered = BaseFlow
CALL InitBaseflowFilter()
CALL BaseflowFilter()

InitBaseflowDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT BASEFLOW DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitBaseflow


!==================================================================================================================================
!>Integrate the Pruett baseflow (time-filtered solution) in time using a simple Euler forward approach:
!>\f$ \frac{d}{dt} \bar{u} \approx \frac{\bar{u}^{n+1} - \bar{u}^{n}}{\Delta t}= \frac{u^n-\bar{u}^n}{\Delta} \f$
!==================================================================================================================================
SUBROUTINE UpdateBaseflow(dt)
! MODULES
USE MOD_PreProc
USE MOD_Baseflow_Filter,ONLY: BaseflowFilter
USE MOD_Baseflow_Vars
USE MOD_DG_Vars,        ONLY: U
USE MOD_Mesh_Vars,      ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: dt                                        !< Current timestep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iElem
!==================================================================================================================================
DO iElem=1,nElems
  fac(iElem) = MIN(dt*TimeFilterWidthBaseflow(iElem),1.0)
  BaseFlow(:,:,:,:,iElem) = BaseFlow(:,:,:,:,iElem)  + (U(:,:,:,:,iElem) - BaseFlow(:,:,:,:,iElem))*fac(iElem)
END DO ! iElem

! Selective Filter
BaseflowFiltered = BaseFlow
CALL BaseflowFilter()

END SUBROUTINE UpdateBaseflow

!==================================================================================================================================
!> Finalizes variables necessary for baseflow.
!==================================================================================================================================
SUBROUTINE FinalizeBaseflow()
! MODULES
USE MOD_Baseflow_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(Baseflow)
SDEALLOCATE(BaseflowFiltered)
SDEALLOCATE(TimeFilterWidthBaseflow)
SDEALLOCATE(fac)
SDEALLOCATE(SelectiveFilterMatrix)
END SUBROUTINE FinalizeBaseflow

END MODULE MOD_Baseflow
