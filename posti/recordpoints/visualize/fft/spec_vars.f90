!===================================================================================================================================
! Contains global variables provided by spectral analysis routines
!===================================================================================================================================
MODULE MOD_Spec_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: nSamples_spec               !> Number of RP samples used for spectral analysis
REAL,ALLOCATABLE              :: RPData_spec(:,:,:)          !> Spectral data
REAL,ALLOCATABLE              :: RPData_freq(:)              !> List of frequencies
INTEGER                       :: nSamples_oct                !> Number of RP samples used for 1/3 octave analysis
REAL,ALLOCATABLE              :: RPData_Oct(:,:,:)           !> Spectral data for 1/3 octave analysis
REAL,ALLOCATABLE              :: RPData_freqOct(:)           !> List of frequencies for 1/3 octave analysis
!===================================================================================================================================
END MODULE MOD_Spec_Vars
