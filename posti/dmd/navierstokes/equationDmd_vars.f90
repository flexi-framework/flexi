MODULE MOD_EquationDMD_Vars
!===================================================================================================================================
! Contains global variables provided by the post processing routines
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
! OUTPUT ---------------------------------------------------------------------------------------------------------------------------
LOGICAL                            :: EquationDMInitIsDone      =.FALSE.
LOGICAL                            :: ConsAvail     =.TRUE.
INTEGER                            :: PrimMap(5)
INTEGER                            :: nVecTrans
INTEGER,ALLOCATABLE                :: TransMap(:,:)
LOGICAL,ALLOCATABLE                :: is2D(:)

INTEGER,ALLOCATABLE               :: DepTable(:,:)
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarNamesAll(:)
INTEGER                           :: nVarDep                 !
INTEGER                           :: nVarCalc
INTEGER                           :: nVarVisuTotal
INTEGER,ALLOCATABLE               :: mapCalc(:)
INTEGER,ALLOCATABLE               :: mapVisu(:)
!===================================================================================================================================
END MODULE MOD_EquationDMD_Vars
