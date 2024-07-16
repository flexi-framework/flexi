!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
!> Unit test 'ReadInToolsUnitTest'
!> Test the module: MOD_ReadInTools
!==================================================================================================================================
PROGRAM ReadInToolsUnitTest
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_MPI       ,ONLY: InitMPI
#ifdef VDM_ANALYTICAL
USE MOD_Mathtools ,ONLY: INVERSE_LU
#else
USE MOD_Mathtools ,ONLY: getSPDInverse
#endif
USE MOD_Mathtools ,ONLY: INVERSE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: nArgs,i,j
INTEGER,PARAMETER         :: nDim=3,nDim2=4
REAL,DIMENSION(nDim,nDim) :: A,AInv
REAL,DIMENSION(nDim2,nDim2) :: B,BInv
LOGICAL                   :: debug
!==================================================================================================================================
debug=.TRUE.
CALL InitMPI()
! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) CALL abort(__STAMP__,'ERROR - Unknown command line argument.')

WRITE (*,*) "\nTest 1 of 3 ---------------------------------------------------\n"

! Invert A
A(1,1) = 16.0
A(1,2) = -8.0
A(1,3) = -4.0

A(2,1) = -8.0
A(2,2) = 29.0
A(2,3) = 12.0

A(3,1) = -4.0
A(3,2) = 12.0
A(3,3) = 41.0

! Debugging
IF(debug)THEN
  WRITE (*,*) "A ="
  DO i = 1, nDim
    DO j = 1, nDim
      write(unit=*, FMT="(E24.12)", ADVANCE="NO") A(i,j)
    END DO ! j = 1, nDim
    print*,""
  END DO ! i = 1, nDim
END IF ! debug
#ifdef VDM_ANALYTICAL
! Computes AInv via analytical expression (only works for Lagrange polynomials, hence the "analytical"
! pre-processor flag) when Lapack fails
! For Bezier (Bernstein basis) polynomial: use INVERSE_LU function
WRITE (*,*) "Testing AInv = INVERSE_LU(A)"
AInv = INVERSE_LU(A)
#else
WRITE (*,*) "Testnig AInv = getSPDInverse(nDim,A)"
AInv = getSPDInverse(nDim,A)
#endif /*VDM_ANALYTICAL*/
! Debugging
IF(debug)THEN
  WRITE (*,*) "AInv ="
  DO i = 1, nDim
    DO j = 1, nDim
      write(unit=*, FMT="(E24.12)", ADVANCE="NO") AInv(i,j)
    END DO ! j = 1, nDim
    print*,""
  END DO ! i = 1, nDim
  A=MATMUL(A,AInv)
  WRITE (*,*) "A*AInv ="
  DO i = 1, nDim
    DO j = 1, nDim
      write(unit=*, FMT="(E24.12)", ADVANCE="NO") A(i,j)
    END DO ! j = 1, nDim
    print*,""
  END DO ! i = 1, nDim
END IF ! debug

IF(ANY(ISNAN(MATMUL(A,AInv))))THEN
  CALL abort(__STAMP__,'MATMUL(A,AInv) has NaNs')
END IF

WRITE (*,*) "\nTest 2 of 3 ---------------------------------------------------\n"

! Invert A
A(1,1) = 16.0
A(1,2) = -8.0
A(1,3) = -4.0

A(2,1) = -8.0
A(2,2) = 29.0
A(2,3) = 12.0

A(3,1) = -4.0
A(3,2) = 12.0
A(3,3) = 41.0

! Debugging
IF(debug)THEN
  WRITE (*,*) "A ="
  DO i = 1, nDim
    DO j = 1, nDim
      WRITE (unit=*, FMT="(E24.12)", ADVANCE="NO") A(i,j)
    END DO ! j = 1, nDim
    print*,""
  END DO ! i = 1, nDim
END IF ! debug

! inverse
WRITE (*,*) "Testing AInv(:,:) = INVERSE(A(:,:))"
AInv(:,:) = INVERSE(A(:,:))

! Debugging
IF(debug)THEN
  WRITE (*,*) "AInv ="
  DO i = 1, nDim
    DO j = 1, nDim
      write(unit=*, FMT="(E24.12)", ADVANCE="NO") AInv(i,j)
    END DO ! j = 1, nDim
    print*,""
  END DO ! i = 1, nDim
  A=MATMUL(A,AInv)
  WRITE (*,*) "A*AInv ="
  DO i = 1, nDim
    DO j = 1, nDim
      WRITE (unit=*, FMT="(E24.12)", ADVANCE="NO") A(i,j)
    END DO ! j = 1, nDim
    WRITE (*,*) ""
  END DO ! i = 1, nDim
END IF ! debug

IF(ANY(ISNAN(MATMUL(A,AInv))))THEN
  CALL abort(__STAMP__,'MATMUL(A,AInv) has NaNs')
END IF

WRITE (*,*) "\nTest 3 of 3 ---------------------------------------------------\n"

! Invert B
B(1,1) = 16.0
B(1,2) = -8.0
B(1,3) = -4.0
B(1,4) = 0.0

B(2,1) = -8.0
B(2,2) = 0.0
B(2,3) = 12.0
B(2,4) = 2.0

B(3,1) = -4.0
B(3,2) = 12.0
B(3,3) = 41.0
B(3,4) = 44.0

B(4,1) = 123.0
B(4,2) = 132.0
B(4,3) = 31.0
B(4,4) = 1.0

! Debugging
IF(debug)THEN
  WRITE (*,*) "B ="
  DO i = 1, nDim2
    DO j = 1, nDim2
      WRITE(unit=*, FMT="(E24.12)", ADVANCE="NO") B(i,j)
    END DO ! j = 1, nDim2
    WRITE (*,*) ""
  END DO ! i = 1, nDim2
END IF ! debug

! inverse
WRITE (*,*) "Testing BInv(:,:) = INVERSE(B(:,:))"
BInv(:,:) = INVERSE(B(:,:))

! Debugging
IF(debug)THEN
  WRITE (*,*) "BInv ="
  DO i = 1, nDim2
    DO j = 1, nDim2
      write(unit=*, FMT="(E24.12)", ADVANCE="NO") BInv(i,j)
    END DO ! j = 1, nDim2
    WRITE (*,*) ""
  END DO ! i = 1, nDim2
  B=MATMUL(B,BInv)
  WRITE (*,*) "B*BInv ="
  DO i = 1, nDim2
    DO j = 1, nDim2
      WRITE(unit=*, FMT="(E24.12)", ADVANCE="NO") B(i,j)
    END DO ! j = 1, nDim2
    WRITE (*,*) ""
  END DO ! i = 1, nDim2
END IF ! debug

IF(ANY(ISNAN(MATMUL(B,BInv))))THEN
  CALL abort(__STAMP__,'MATMUL(B,BInv) has NaNs')
END IF

#if USE_MPI
! we also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError.NE.MPI_SUCCESS) CALL abort(__STAMP__,'MPI finalize error')
#endif

END PROGRAM ReadInToolsUnitTest
