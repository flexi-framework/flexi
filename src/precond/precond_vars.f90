#include "flexi.h"
MODULE MOD_Precond_Vars
!===================================================================================================================================
! Contains global variables used by the Timedisc modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE,TARGET  :: invP(:,:,:)  !inverse of block Jacobian for each element (1:nDOF_elem,1:nDOFelem,1:nElems)
INTEGER               :: PrecondType
INTEGER(KIND=8)       :: PrecondIter               !< Defines how often preconditioner is built
INTEGER               :: DebugMatrix
LOGICAL               :: doVol,doSurf
LOGICAL               :: PrecondInitIsDone
LOGICAL               :: EulerPrecond
INTEGER               :: SolveSystem
LOGICAL               :: NoFillIn
LOGICAL               :: DoDisplayPrecond
!===================================================================================================================================
END MODULE MOD_Precond_Vars
