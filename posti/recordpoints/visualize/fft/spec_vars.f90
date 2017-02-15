MODULE MOD_spec_Vars
!===================================================================================================================================
! Contains global variables provided by the output routines
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: nSamples_spec               ! number of visualisation points is NVisu+1
REAL,ALLOCATABLE              :: RPData_spec(:,:,:)          ! output state 
REAL,ALLOCATABLE              :: RPData_freq(:)              ! time average state 
INTEGER                       :: nSamples_oct               ! number of visualisation points is NVisu+1
REAL,ALLOCATABLE              :: RPData_Oct(:,:,:)          ! output state 
REAL,ALLOCATABLE              :: RPData_freqOct(:)              ! time average state 
!===================================================================================================================================
END MODULE MOD_spec_Vars
