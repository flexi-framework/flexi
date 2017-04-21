!===================================================================================================================================
!> Contains global variables for the Posti tool
!===================================================================================================================================
MODULE MOD_Parameters
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: NSuper
REAL                          :: maxTol                   !< max overshoot in param coords (1+maxTol)
LOGICAL                       :: doVisuRP

INTEGER                       :: nCoords                  !< number of visualisation points is NVisu+1
CHARACTER(LEN=255),ALLOCATABLE:: CoordNames(:)            !< including space and time coordinates
LOGICAL                       :: OutputInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Parameters
