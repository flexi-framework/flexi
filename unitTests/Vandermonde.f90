#include "flexi.h"

!==================================================================================================================================
!> Unit test 'Vandermonde'
!> Test the routine: 'GetVandermonde', from module: 'Interpolation'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM Vandermonde
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_Interpolation_Vars,  ONLY: NodeTypeG,NodeTypeGL,NodeTypeCL,NodeTypeVISU
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Vdm_In_Out(0:6,0:5,1:512),Vdm_Out_In(0:5,0:6,1:512)
REAL                           :: Vdm_In_Out_ref(0:6,0:5,1:512),Vdm_Out_In_ref(0:5,0:6,1:512)
INTEGER                        :: iIn,iOut,iNodeTypeIn,iNodeTypeOut,iModal,nArgs,k
CHARACTER(LEN=*),PARAMETER     :: BinaryString='Vandermonde.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,modal,equal
CHARACTER(LEN=255)             :: argument,NodeTypeIn,NodeTypeOut
!==================================================================================================================================
! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) THEN
  CALL GET_COMMAND_ARGUMENT(1,argument)
  IF (argument.EQ.TRIM('--generate-reference')) THEN
    doGenerateReference = .TRUE.
  ELSE
    WRITE(*,*) 'ERROR - Unknown command line argument.'
    STOP -1
  END IF
END IF

! Generate Vandermondes from source
Vdm_In_Out = 0.
Vdm_Out_In = 0.
! All combinations of input and output polynomial degrees
k = 1
DO iIn=2,5; DO iOut=3,6
  DO iNodeTypeIn=1,4; DO iNodeTypeOut=1,4
    SELECT CASE(iNodeTypeIn)
    CASE(1)
      NodeTypeIn = NodeTypeG
    CASE(2)
      NodeTypeIn = NodeTypeGL
    CASE(3)
      NodeTypeIn = NodeTypeCL
    CASE(4)
      NodeTypeIn = NodeTypeVISU
    END SELECT
    SELECT CASE(iNodeTypeOut)
    CASE(1)
      NodeTypeOut = NodeTypeG
    CASE(2)
      NodeTypeOut = NodeTypeGL
    CASE(3)
      NodeTypeOut = NodeTypeCL
    CASE(4)
      NodeTypeOut = NodeTypeVISU
    END SELECT
    DO iModal = 1,2
      SELECT CASE(iModal)
      CASE(1)
        modal = .FALSE.
      CASE(2)
        modal = .TRUE.
      END SELECT
      CALL GetVandermonde(iIn,NodeTypeIn,iOut,NodeTypeOut,Vdm_In_Out(0:iOut,0:iIn,k),Vdm_Out_In(0:iIn,0:iOut,k),modal)
      k=k+1
    END DO
  END DO; END DO
END DO; END DO


IF (doGenerateReference) THEN
  ! Save the calculated solution to a binary file for later comparison
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Vdm_Out_In,Vdm_In_Out
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)
  
  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) Vdm_Out_In_ref,Vdm_In_Out_ref
    CLOSE(10) ! close the file
    ! Check if the computed and the reference solutions are within a given tolerance
    equal = .TRUE.
    DO k = 1,512
      DO iIn=0,5; DO iOut=0,6
        equal = ALMOSTEQUALABSORREL(Vdm_Out_In(iIn,iOut,k),Vdm_Out_In_ref(iIn,iOut,k),100.*PP_RealTolerance) .AND. equal
        equal = ALMOSTEQUALABSORREL(Vdm_In_Out(iOut,iIn,k),Vdm_In_Out_ref(iOut,iIn,k),100.*PP_RealTolerance) .AND. equal
      END DO; END DO
    END DO
    IF (.NOT.equal) THEN
      WRITE(*,*) 'ERROR - Calculated Vandermondes deviate from reference.'
      STOP -1
    ELSE
      WRITE(*,*) 'Checked Vandermondes against stored data -- SUCCESSFUL.'
    END IF
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM Vandermonde


