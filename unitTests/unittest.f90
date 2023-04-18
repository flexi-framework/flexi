#include "flexi.h"

!==================================================================================================================================
!> Module containing general routines used in Unittests
!==================================================================================================================================
MODULE MOD_Unittest
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE GenerateUnittestReferenceData
  MODULE PROCEDURE GenerateUnittestReferenceData
END INTERFACE

INTERFACE ReadInReferenceElementData
  MODULE PROCEDURE ReadInReferenceElementData
END INTERFACE

PUBLIC::GenerateUnittestReferenceData,ReadInReferenceElementData
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
!----------------------------------------------------------------------------------------------------------------------------------
! Local Variables
CHARACTER(LEN=255)             :: Filename 
!==================================================================================================================================
#if FV_ENABLED == 0
WRITE(*,*) 'To generate reference data for the unit tests, please compile FLEXI with FV activated.'
STOP
#endif

#if PP_dim == 3
Filename = "UnittestElementData3D.bin"
#else
Filename = "UnittestElementData2D.bin"
#endif
! Save the calculated solution to a binary file for later comparison
OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(Filename),FORM='unformatted')  ! replace an existing file or create a new one
WRITE(10) nElems,SideToElem,firstMPISide_YOUR,lastMPISide_MINE,nSides,S2V2,L_Minus,L_Plus,L_HatPlus,L_HatMinus,sJ
CLOSE(10) ! close the file
WRITE(*,*) "Generated Unittest reference data into: ", TRIM(Filename)
END SUBROUTINE GenerateUnittestReferenceData

!==================================================================================================================================
!> Read in the data for the curved reference element, allocate the necessary arrays beforehand.
!==================================================================================================================================
SUBROUTINE ReadInReferenceElementData()
! MODULES
USE MOD_Unittest_Vars
USE MOD_Mesh_Vars,             ONLY: SideToElem,S2V2,nElems,nSides,firstMPISide_YOUR,lastMPISide_MINE,sJ
USE MOD_Interpolation_Vars,    ONLY: L_Minus,L_Plus
USE MOD_DG_Vars,               ONLY: L_HatPlus,L_HatMinus
!----------------------------------------------------------------------------------------------------------------------------------
! Local Variables
CHARACTER(LEN=255)             :: Filename 
INTEGER                        :: Flip_lower,Flip_upper,locSide_lower,locSide_upper
!==================================================================================================================================
! Dimensions for mappings
#if PP_dim == 3
  Flip_lower = 0
  Flip_upper = 4
  locSide_lower = 1
  locSide_upper = 6
#else
  Flip_lower = 0
  Flip_upper = 1
  locSide_lower = 2
  locSide_upper = 5
#endif

#if PP_dim == 3
Filename = "UnittestElementData3D.bin"
#else
Filename = "UnittestElementData2D.bin"
#endif

! Read in data from single curved element
ALLOCATE(SideToElem(1:5,1:nSidesRef))
ALLOCATE(S2V2(1:2,0:NRef,0:NRefZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper))
ALLOCATE(L_Minus(0:NRef))
ALLOCATE(L_Plus(0:NRef))
ALLOCATE(L_HatMinus(0:NRef))
ALLOCATE(L_HatPlus(0:NRef))
ALLOCATE(sJ(0:NRef,0:NRef,0:NRefZ,0:FV_SIZE,1:nElemsRef))
OPEN(UNIT = 10, STATUS='old',FILE=TRIM(Filename),FORM='unformatted')  ! open an existing file
READ(10) nElems,SideToElem,firstMPISide_YOUR,lastMPISide_MINE,nSides,S2V2,L_Minus,L_Plus,L_HatPlus,L_HatMinus,sJ
CLOSE(10) ! close the file

END SUBROUTINE ReadInReferenceElementData

END MODULE MOD_Unittest
