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

!==================================================================================================================================
!> Wrapper routines for LAPACK
!==================================================================================================================================
MODULE MOD_Lapack
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

! Explicit interfaces alone do not connect the subroutines to their implementations.
! > he EXTERNAL attribute is automatically implied when defining them in an INTERFACE block
INTERFACE
  ! Interface for DSYEV
  SUBROUTINE DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
    INTEGER, INTENT(IN) :: N, LDA, LWORK
    DOUBLE PRECISION, INTENT(INOUT) :: A(LDA, *)
    DOUBLE PRECISION, INTENT(OUT) :: W(N)
    DOUBLE PRECISION, INTENT(INOUT) :: WORK(LWORK)
    INTEGER, INTENT(OUT) :: INFO
  END SUBROUTINE DSYEV

  ! Interface for DGETRF
  SUBROUTINE DGETRF(M, N, A, LDA, IPIV, INFO)
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    INTEGER, INTENT(IN) :: M, N, LDA
    DOUBLE PRECISION, INTENT(INOUT) :: A(LDA, *)
    INTEGER, INTENT(OUT) :: IPIV(MIN(M, N))
    INTEGER, INTENT(OUT) :: INFO
  END SUBROUTINE DGETRF

  ! Interface for DGETRI
  SUBROUTINE DGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    INTEGER, INTENT(IN) :: N, LDA, LWORK
    DOUBLE PRECISION, INTENT(INOUT) :: A(LDA, *)
    INTEGER, INTENT(IN) :: IPIV(N)
    DOUBLE PRECISION, INTENT(INOUT) :: WORK(LWORK)
    INTEGER, INTENT(OUT) :: INFO
  END SUBROUTINE DGETRI
END INTERFACE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSYEV
PUBLIC :: DGETRF
PUBLIC :: DGETRI
!==================================================================================================================================

END MODULE MOD_Lapack
