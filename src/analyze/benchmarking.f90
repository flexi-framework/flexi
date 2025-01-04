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
!> Contains benchmarking routines measuring the amount of flops
!==================================================================================================================================
MODULE MOD_Benchmarking
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
#ifdef PAPI
LOGICAL :: doMeasureFlops
#endif
! Public Part ----------------------------------------------------------------------------------------------------------------------

PUBLIC:: InitBenchmarking
PUBLIC:: Benchmarking
PUBLIC:: FinalizeBenchmarking
!==================================================================================================================================

CONTAINS

SUBROUTINE InitBenchmarking()
!==================================================================================================================================
!> Initializes variables necessary for benchmarking subroutines
!==================================================================================================================================
! MODULES
USE MOD_ReadInTools,        ONLY: GETLOGICAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#ifdef PAPI
doMeasureFlops=GETLOGICAL('doMeasureFlops')
#endif
END SUBROUTINE InitBenchmarking


SUBROUTINE Benchmarking()
!==================================================================================================================================
!> Calls routine to measure amount of flops
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#ifdef PAPI
IF(doMeasureFlops) CALL MeasureFlops()
#endif
END SUBROUTINE Benchmarking


#ifdef PAPI
SUBROUTINE MeasureFlops()
!==================================================================================================================================
!> Caclulates and summarizes amount of flops
!==================================================================================================================================
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_Preproc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(C_FLOAT)        :: realtime,proctime,mflops
INTEGER(C_INT)       :: err
INTEGER(C_LONG_LONG) :: flpops
!==================================================================================================================================
CALL PAPIF_FLOPS(realtime,proctime,flpops,mflops,err)
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,mflops,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  WRITE(UNIT_stdOut,'(A14,ES18.9)')' Sim-MFLOPS : ',mflops
ELSE
  CALL MPI_REDUCE(mflops,0           ,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif
END SUBROUTINE MeasureFlops
#endif


SUBROUTINE FinalizeBenchmarking()
!==================================================================================================================================
!> Finalizes variables necessary for benchmarking subroutines
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeBenchmarking

END MODULE MOD_Benchmarking
