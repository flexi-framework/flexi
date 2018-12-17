!===================================================================================================================================
!> Contains global variables used for/by the Interpolation Module
!===================================================================================================================================
MODULE MOD_RPInterpolation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                         :: CalcTimeAverage  !< Logical set if time averages should be calculated
REAL,ALLOCATABLE                :: dt(:)            !< Time step between each RP sample
REAL                            :: dt_out           !< Equidistant time step
REAL                            :: TEnd             !< Time of final RP sample
!===================================================================================================================================
END MODULE MOD_RPInterpolation_Vars
