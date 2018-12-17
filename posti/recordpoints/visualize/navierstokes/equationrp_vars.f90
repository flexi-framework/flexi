!===================================================================================================================================
!> Contains global variables provided by the visualize recordpoints navier stokes module
!===================================================================================================================================
MODULE MOD_EquationRP_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

! LOCAL TRANSFORMATION -------------------------------------------------------------------------------------------------------------
INTEGER                            :: nVecTrans                    !< Number of vector quantities that should be transformed
INTEGER,ALLOCATABLE                :: TransMap(:,:)                !< Mapping to those vector quantities
LOGICAL,ALLOCATABLE                :: is2D(:)                      !< Indicating if one of those quantities is two dimensional
! BOUNDARY LAYER PROPERTIES --------------------------------------------------------------------------------------------------------
INTEGER                            :: nBLProps                     !< Number of avariables for boundary layer properties
CHARACTER(LEN=255),ALLOCATABLE     :: VarNames_BLProps(:)          !< Variable names of boundary layer properties

LOGICAL                            :: EquationRPInitIsDone=.FALSE. !< Switch to signal that init is done
!===================================================================================================================================
END MODULE MOD_EquationRP_Vars
