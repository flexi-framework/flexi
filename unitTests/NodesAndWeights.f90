#include "flexi.h"

!==================================================================================================================================
!> Unit test 'NodesAndWeights'
!> Test the routine: 'GetNodesAndWeights', from module: 'Interpolation'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM NodesAndWeights
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation,       ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars,  ONLY: NodeTypeG,NodeTypeGL,NodeTypeCL,NodeTypeVISU
USE MOD_Basis,               ONLY: EQUALTOTOLERANCE

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: xi(0:10,1:10,1:4),xi_ref(0:10,1:10,1:4)
REAL                           :: w(0:10,1:10,1:4),w_ref(0:10,1:10,1:4)
REAL                           :: wBary(0:10,1:10,1:4),wBary_ref(0:10,1:10,1:4)
INTEGER                        :: i,j,k,nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryString='NodesAndWeights.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal
CHARACTER(LEN=255)             :: argument
!==================================================================================================================================

! Generate nodes and weights from current source
xi = 0.
w = 0.
wBary = 0.
DO i=1,10
  CALL GetNodesAndWeights(i,NodeTypeG,   xi(0:i,i,1),w(0:i,i,1),wBary(0:i,i,1))
  CALL GetNodesAndWeights(i,NodeTypeGL,  xi(0:i,i,2),w(0:i,i,2),wBary(0:i,i,2))
  CALL GetNodesAndWeights(i,NodeTypeCL,  xi(0:i,i,3),w(0:i,i,3),wBary(0:i,i,3))
  CALL GetNodesAndWeights(i,NodeTypeVISU,xi(0:i,i,4),w(0:i,i,4),wBary(0:i,i,4))
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
  WRITE(10) xi,w,wBary
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)
  
  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) xi_ref,w_ref,wBary_ref
    CLOSE(10) ! close the file
    ! Check if the computed and the reference solutions are within a given tolerance
    equal =  .TRUE.
    DO i=0,10; DO j=1,10; DO k=1,4
      equal = EQUALTOTOLERANCE(xi(i,j,k),xi_ref(i,j,k),50.*PP_RealTolerance) .AND. equal
      equal = EQUALTOTOLERANCE(w(i,j,k),w_ref(i,j,k),50.*PP_RealTolerance) .AND. equal
      equal = EQUALTOTOLERANCE(wBary(i,j,k),wBary_ref(i,j,k),50.*PP_RealTolerance) .AND. equal
    END DO; END DO; END DO
    IF (.NOT.equal) THEN
      WRITE(*,*) 'ERROR - Calculated nodes and weights deviate from reference.'
      STOP -1
    END IF
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM NodesAndWeights
