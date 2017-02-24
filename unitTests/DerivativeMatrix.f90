#include "flexi.h"

!==================================================================================================================================
!> Unit test 'DerivativeMatrix'
!> Test the routine: 'GetDerivativeMatrix', from module: 'Interpolation'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM DerivativeMatrix
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation,       ONLY: GetDerivativeMatrix
USE MOD_Interpolation_Vars,  ONLY: NodeTypeG,NodeTypeGL,NodeTypeCL,NodeTypeVISU
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: D(0:10,0:10,1:40),D_ref(0:10,0:10,1:40)
INTEGER                        :: i,j,k,iNodeType,nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryString='DerivativeMatrix.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal
CHARACTER(LEN=255)             :: argument,NodeTypeIn
!==================================================================================================================================
! Initalize derivative matrix
D = 0.

! Generate derivative matrix from current source
k = 1 ! Counter for global matrix storage
DO iNodeType=1,4 ! Loop over all node tyoes
  SELECT CASE(iNodeType)
  CASE(1)
    NodeTypeIn = NodeTypeG
  CASE(2)
    NodeTypeIn = NodeTypeGL
  CASE(3)
    NodeTypeIn = NodeTypeCL
  CASE(4)
    NodeTypeIn = NodeTypeVISU
  END SELECT
  DO i=1,10 ! Loop over polynomial degrees
    CALL GetDerivativeMatrix(i,NodeTypeIn,D(0:i,0:i,k))
    k = k + 1
  END DO
END DO


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

IF (doGenerateReference) THEN
  ! Save the calculated solution to a binary file for later comparison
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) D
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)
  
  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) D_ref
    CLOSE(10) ! close the file
    ! Check if the computed and the reference solutions are within a given tolerance
    equal =  .TRUE.
    DO i=0,10; DO j=0,10; DO k=1,40
      equal = ALMOSTEQUALABSORREL(D(i,j,k),D_ref(i,j,k),100.*PP_RealTolerance) .AND. equal
    END DO; END DO; END DO
    ! Check against real tolerance times 100
    IF (.NOT.equal) THEN
      WRITE(*,*) 'ERROR - Calculated derivative matrizes deviate from reference.'
      STOP -1
    ELSE
      WRITE(*,*) 'Checked derivative matrizes against stored data -- SUCCESSFUL.'
    END IF
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM DerivativeMatrix
