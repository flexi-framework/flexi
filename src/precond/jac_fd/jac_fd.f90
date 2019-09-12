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

MODULE MOD_Jac_FD
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U 
!===================================================================================================================================
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
USE MOD_Implicit_Vars, ONLY: nDOFVarElem,R_Xk,Xk,sreps0,reps0
USE MOD_DG,            ONLY: DGTimeDerivative_WeakForm
USE MOD_DG_Vars,       ONLY: U,Ut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                    :: t
INTEGER,INTENT(IN)                 :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: dRdU(1:nDOFVarElem,1:nDOFVarElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                            :: iVar,r,s,ii,jj,kk,iiVar,i,j,k
#if USE_MPI
INTEGER                            :: iProc
#endif
!===================================================================================================================================
Xk=U
CALL DGTimeDerivative_WeakForm(t)
R_xk=Ut !linearization Ut of Xk for FD

dRdU=0.

#if USE_MPI
DO iProc=0,nProcessors-1
  s=1
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          IF (iProc.EQ.myRank) THEN
            U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k,iElem) + reps0
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
                    dRdU(r,s) = dRdU(r,s)+(Ut(iiVar,ii,jj,kk,iElem)-R_xk(iiVar,ii,jj,kk,iElem))*sreps0
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
END DO !iProc
#else
s=1
DO k=0,PP_NZ
  DO j=0,PP_N
    DO i=0,PP_N
      DO iVar=1,PP_nVar
        U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k,iElem) + reps0
        CALL DGTimeDerivative_WeakForm(t)
        U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k,iElem) 
        r=1
        DO kk=0,PP_NZ
          DO jj=0,PP_N
            DO ii=0,PP_N
              DO iiVar=1,PP_nVar
                dRdU(r,s) = dRdU(r,s)+(Ut(iiVar,ii,jj,kk,iElem)-R_xk(iiVar,ii,jj,kk,iElem))*sreps0
                r=r+1
              END DO !iiVar
            END DO !ii
          END DO !jj
        END DO !kk
        s=s+1
      END DO !PP_nVar
    END DO !i
  END DO !j
END DO !k
#endif

!cleanup
CALL DGTimeDerivative_WeakForm(t)

END SUBROUTINE Jac_FD
END MODULE MOD_Jac_FD
