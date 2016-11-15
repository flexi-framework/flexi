MODULE MOD_Indicator_Vars
!==================================================================================================================================
! Contains variables relevant for indicators
!==================================================================================================================================
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                :: IndicatorInitIsDone=.FALSE.
INTEGER                :: IndicatorType  ! Type of indicator used: 0 = no indicator, 1 = Persson
INTEGER                :: IndVar         ! variable on which indicator is applied (only cons)
REAL,ALLOCATABLE       :: IndValue(:)    ! indicator output
REAL                   :: IndStartTime
!==================================================================================================================================
END MODULE MOD_Indicator_Vars
