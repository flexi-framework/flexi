!===================================================================================================================================
!> Contains global variables used for/by the RPSet
!===================================================================================================================================
MODULE MOD_RPSetVisuVisu_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                         :: RPSetInitIsDone = .FALSE.
INTEGER                         :: nRP_global,nRP_HDF5
INTEGER                         :: nGroups,nLines,nPoints,nPlanes
CHARACTER(LEN=255),ALLOCATABLE  :: GroupNames(:)
LOGICAL,ALLOCATABLE             :: OutputGroup(:)
INTEGER,ALLOCATABLE             :: Points_IDlist(:)
INTEGER,ALLOCATABLE             :: Points_GroupIDlist(:)
REAL,ALLOCATABLE                :: x_RP(:,:)
REAL,ALLOCATABLE                :: xF_RP(:,:)
INTEGER,ALLOCATABLE             :: RPOutMap(:)

TYPE tLine         
  CHARACTER(LEN=255)            :: Name
  INTEGER                       :: GroupID
  INTEGER                       :: nRP
  INTEGER,ALLOCATABLE           :: IDlist(:)
  REAL,ALLOCATABLE              :: LocalCoord(:)
  REAL,ALLOCATABLE              :: Tmat(:,:)
END TYPE tLine

TYPE tPlane
  CHARACTER(LEN=255)            :: Name
  INTEGER                       :: GroupID
  INTEGER                       :: Type=0 ! 0 - standard, 1 - sphere, 2 - BLPlane
  INTEGER                       :: nRP(2) ! RP resolution in i,j direction i: P1->P2, 
  INTEGER,ALLOCATABLE           :: IDlist(:,:)
  REAL,ALLOCATABLE              :: NormVec(:,:),TangVec(:,:),LocalCoord(:,:,:)
  REAL,ALLOCATABLE              :: BLProps(:,:)
END TYPE tPlane

TYPE(tLine),POINTER             :: Lines(:)
TYPE(tPlane),POINTER            :: Planes(:)

!===================================================================================================================================
END MODULE MOD_RPSetVisuVisu_Vars
