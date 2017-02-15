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
LOGICAL                       :: OutputInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_OutputRPVisu_Vars
