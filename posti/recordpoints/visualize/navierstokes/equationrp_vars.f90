MODULE MOD_EquationRP_Vars
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
CHARACTER(LEN=255)                 :: strOutputFile, strMeshFile, strLinkFile
LOGICAL                            :: InitOutputDone=.FALSE.
LOGICAL                            :: EquationRPInitIsDone      =.FALSE.
LOGICAL                            :: ConsAvail     =.TRUE.
INTEGER                            :: PrimMap(5)
INTEGER                            :: nVecTrans
INTEGER,ALLOCATABLE                :: TransMap(:,:)
LOGICAL,ALLOCATABLE                :: is2D(:)
INTEGER                            :: nBLProps
CHARACTER(LEN=255),ALLOCATABLE     :: VarNames_BLProps(:)
!===================================================================================================================================
END MODULE MOD_EquationRP_Vars
