#include "flexi.h"

MODULE MOD_SparseILU_Vars
!===================================================================================================================================
! Contains global variables used for the linear operations
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                     :: nMTriangle
REAL                                        :: epsZero
INTEGER,ALLOCATABLE,DIMENSION(:)            :: nUNonZeros,nLNonZeros
REAL,ALLOCATABLE,DIMENSION(:,:)             :: Dinv

TYPE tILU                                                                     ! ILU for each element
 REAL,ALLOCATABLE,DIMENSION(:)              :: Entry
 INTEGER,ALLOCATABLE,DIMENSION(:)           :: IEntry,JEntry
END TYPE                                                                     
TYPE(tILU), ALLOCATABLE                     :: IL(:)                          ! Incomplete Lower matrix
TYPE(tILU), ALLOCATABLE                     :: IU(:)                          ! Incomplete Upper matrix 
!===================================================================================================================================
END MODULE MOD_SparseILU_Vars
