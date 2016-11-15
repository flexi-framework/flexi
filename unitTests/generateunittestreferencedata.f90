#include "flexi.h"

!==================================================================================================================================
!> Module for generating a reference binary file for the unittests. 
!==================================================================================================================================
MODULE MOD_GenerateUnittestReferenceData
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE GenerateUnittestReferenceData
  MODULE PROCEDURE GenerateUnittestReferenceData
END INTERFACE

PUBLIC::GenerateUnittestReferenceData
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Writes several arrays to a binary file, that are used for the unittest. Therewith the unittest do not have to compute 
!> geometrical quantities, index arrays, ...
!==================================================================================================================================
SUBROUTINE GenerateUnittestReferenceData()
! MODULES
USE MOD_Mesh_Vars
USE MOD_Interpolation_Vars
USE MOD_DG_Vars
! Local Variables
CHARACTER(LEN=255)             :: Filename 
!==================================================================================================================================
Filename = "UnittestElementData.bin"
! Save the calculated solution to a binary file for later comparison
OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(Filename),FORM='unformatted')  ! replace an existing file or create a new one
WRITE(10) nElems,SideToElem,firstMPISide_YOUR,lastMPISide_MINE,nSides,S2V3,CS2V2,V2S2,L_Minus,L_Plus,L_HatPlus,L_HatMinus,sJ
CLOSE(10) ! close the file
WRITE(*,*) "Generated Unittest reference data into: ", TRIM(Filename)
END SUBROUTINE GenerateUnittestReferenceData

END MODULE MOD_GenerateUnittestReferenceData
