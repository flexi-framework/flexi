!===================================================================================================================================
!> Contains global variables used for/by the RP data module
!===================================================================================================================================
MODULE MOD_RPData_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: nSamples_global !> Total number of samples (all RP data files)
INTEGER                         :: nVar_HDF5       !> Number of variable in the HDF5 file
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames_HDF5(:)!> Name of the variables in the HDF5 file
REAL,ALLOCATABLE                :: RPData(:,:,:)   !> Global array containing the data of all samples
REAL,ALLOCATABLE                :: RPTime(:)       !> The time value of all samples
!-----------------------------------------------------------------------------------------------------------------------------------
! Output Buffer
!-----------------------------------------------------------------------------------------------------------------------------------
!> Type that is used to collect all the record point data from the different files before merging them in the RPData array,
!> organized in a linked list
TYPE tRPDataSet
  REAL,ALLOCATABLE              :: data(:,:,:) !> Actual data array
  INTEGER                       :: nSamples    !> Number of (local) samples in the current data set
  TYPE(tRPDataSet),POINTER      :: nextset     !> Pointer to the next set in the linked list
END TYPE tRPDataSet

TYPE(tRPDataSet),POINTER        :: firstset, actualset !> Pointers to first and current data set in the linked list

!===================================================================================================================================

INTERFACE getNewRPDataSet
  MODULE PROCEDURE getNewRPDataSet
END INTERFACE 

PUBLIC :: getNewRPDataSet

CONTAINS 

!===================================================================================================================================
!> Routine to create a new entry in the linked list of RPdata
!===================================================================================================================================
SUBROUTINE getNewRPDataSet(RPDataSet,nSamples_in)
! MODULES
USE MOD_RPSetVisuVisu_Vars            ,ONLY: nRP_global
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nSamples_in
TYPE(tRPDataSet),POINTER      :: RPDataSet
!===================================================================================================================================
ALLOCATE(RPDataSet)
ALLOCATE(RPDataSet%data(0:nVar_HDF5,1:nRP_global,1:nSamples_in))
RPDataSet%data=0.
NULLIFY(RPDataSet%nextset)
RPDataSet%nSamples=nSamples_in
END SUBROUTINE getNewRPDataSet

END MODULE MOD_RPData_Vars
