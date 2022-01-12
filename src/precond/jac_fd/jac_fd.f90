!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
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

!===================================================================================================================================
!> Module that contains routines to compute the DG residual dRdU using a finite difference approach. Call the DG operator for each
!> degree of freedom with a slightly pertubated state. Then use a first order FD to compute the derivative w.r.t that DOF.
!> This is a very accurate but super slow preconditioner, mainly used for debugging (compare against approximate precond).
!> ATTENTION: For MPI computations, this only works if all processor have the same number of elements (otherwise, some processors
!> will call the DG operator more often than other).
!===================================================================================================================================
MODULE MOD_Jac_FD
! MODULES
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Jac_FD
  MODULE PROCEDURE Jac_FD
END INTERFACE

PUBLIC::Jac_FD
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Computes the Finite Difference-Derivative dRdU for the Calculation of Preconditioner P
!> Attention: dRdU = 0 (in precond.f90)
!===================================================================================================================================
SUBROUTINE Jac_FD(t,iElem,dRdU)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars, ONLY: nDOFVarElem,R_Xk,Xk,sreps0_O1,reps0_O1
USE MOD_DG,            ONLY: DGTimeDerivative_WeakForm
USE MOD_DG_Vars,       ONLY: U,Ut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                    :: t        !< current stage time
INTEGER,INTENT(IN)                 :: iElem    !< element counter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: dRdU(1:nDOFVarElem,1:nDOFVarElem) !< Jacobian of the DG operator, computed using FD approach
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iVar,r,s,ii,jj,kk,iiVar,i,j,k
INTEGER                            :: iProc
!===================================================================================================================================
Xk=U ! Store non-pertubated state
CALL DGTimeDerivative_WeakForm(t)
R_xk=Ut ! Store the DG operator for the non-pertubated state

! For each degree of freedom, subsequently pertubate that state and compute the new DG operator. Then, calculate the first order FD
! approximation of dRdU.
#if USE_MPI
DO iProc=0,nProcessors-1
#else
iProc = 0
#endif
  s=1
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          IF (iProc.EQ.myRank) THEN
            U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k,iElem) + reps0_O1
          END IF
          CALL DGTimeDerivative_WeakForm(t)
          IF (iProc.EQ.myRank) THEN
            U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k,iElem)
          END IF
          IF (iProc.EQ.myRank) THEN
            r=1
            DO kk=0,PP_NZ
              DO jj=0,PP_N
                DO ii=0,PP_N
                  DO iiVar=1,PP_nVar
                    dRdU(r,s) = dRdU(r,s)+(Ut(iiVar,ii,jj,kk,iElem)-R_xk(iiVar,ii,jj,kk,iElem))*sreps0_O1
                    r=r+1
                  END DO !iiVar
                END DO !ii
              END DO !jj
            END DO !kk
          END IF
          s=s+1
        END DO !PP_nVar
      END DO !i
    END DO !j
  END DO !k
#if USE_MPI
END DO !iProc
#endif

! cleanup
CALL DGTimeDerivative_WeakForm(t)

END SUBROUTINE Jac_FD
END MODULE MOD_Jac_FD
