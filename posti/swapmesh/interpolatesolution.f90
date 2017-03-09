!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
!> Evaluates the old solution at some given parametric coordinates of the new mesh
!===================================================================================================================================
MODULE MOD_InterpolateSolution
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InterpolateSolution
  MODULE PROCEDURE InterpolateSolution
END INTERFACE

PUBLIC :: InterpolateSolution
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Evaluates the old solution at some given parametric coordinates of the new mesh
!===================================================================================================================================
SUBROUTINE InterpolateSolution()
! MODULES
USE MOD_Globals
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys,BarycentricWeights
USE MOD_Interpolation,     ONLY: GetNodesAndWeights
USE MOD_SwapMesh_Vars,     ONLY: equalElem,Vdm_GPNState_GPNNew
USE MOD_SwapMesh_Vars,     ONLY: NState,NInter,NNew,NodeTypeState,RefState,nVar_State
USE MOD_SwapMesh_Vars,     ONLY: Vdm_CLNInter_GPNNew
USE MOD_SwapMesh_Vars,     ONLY: UOld,xiInter,InterToElem,nElemsNew,IPDone
USE MOD_DG_Vars,           ONLY: U
USE MOD_ChangeBasis,       ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: i,j,k    ! GPs in old element
INTEGER                       :: ii,jj,kk ! GPs in new element
INTEGER                       :: iElemOld,iElemNew
REAL                          :: Utmp(1:nVar_State,0:NInter,0:NInter,0:NInter)
REAL                          :: L_xi(  0:NState,0:NInter,0:NInter,0:NInter)
REAL                          :: L_eta( 0:NState,0:NInter,0:NInter,0:NInter)
REAL                          :: L_zeta(0:NState,0:NInter,0:NInter,0:NInter)
REAL                          :: L_eta_zeta
REAL                          :: xGP(0:NState),wBaryGP(0:NState)
REAL                          :: Time
!===================================================================================================================================
! GPs and Barycentric weights for solution
CALL GetNodesAndWeights(NState,NodeTypeState,xGP)
CALL BarycentricWeights(NState,xGP,wBaryGP)

SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' INTERPOLATE STATE TO NEW MESH...'
Time=FLEXITIME()
U=0.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iElemNew,iElemOld,L_xi,ii,jj,kk,L_eta,L_zeta,Utmp,i,j,k)
!$OMP DO
DO iElemNew=1,nElemsNew
  ! equal elements
  IF(equalElem(iElemNew).GT.0) THEN
    iElemOld=equalElem(iElemNew)
    IF(NState.EQ.Nnew)THEN
      U(:,:,:,:,iElemNew)=UOld(:,:,:,:,iElemOld)
    ELSE
      CALL ChangeBasis3D(nVar_State,NState,NNew,Vdm_GPNState_GPNNew,UOld(:,:,:,:,iElemOld),U(:,:,:,:,iElemNew))
    END IF
    CYCLE
  END IF

  ! rest
  ! Calculate basis function values for xiInter
  DO kk=0,NInter; DO jj=0,NInter; DO ii=0,NInter
    IF(.NOT.IPDone(ii,jj,kk,iElemNew)) CYCLE
    CALL LagrangeInterpolationPolys(xiInter(1,ii,jj,kk,iElemNew),NState,xGP,wBaryGP,L_xi(  :,ii,jj,kk))
    CALL LagrangeInterpolationPolys(xiInter(2,ii,jj,kk,iElemNew),NState,xGP,wBaryGP,L_eta( :,ii,jj,kk))
    CALL LagrangeInterpolationPolys(xiInter(3,ii,jj,kk,iElemNew),NState,xGP,wBaryGP,L_zeta(:,ii,jj,kk))
  END DO; END DO; END DO

  ! Evaluate old solution at interpolation points and transform to new state polynomial degree and node type
  DO kk=0,NInter; DO jj=0,NInter; DO ii=0,NInter
    IF(.NOT.IPDone(ii,jj,kk,iElemNew))THEN
      ! If the interpolation point was not found, overwrite with reference state
      ! If no reference state was given, the program already aborted
      Utmp(:,ii,jj,kk)=RefState
      CYCLE
    ELSE
      Utmp(:,ii,jj,kk)=0.
    END IF
    iElemOld = InterToElem(ii,jj,kk,iElemNew)
    DO k=0,NState
      DO j=0,NState
        L_eta_zeta=L_eta(j,ii,jj,kk)*L_zeta(k,ii,jj,kk)
        DO i=0,NState
          Utmp(:,ii,jj,kk) = Utmp(:,ii,jj,kk) + UOld(:,i,j,k,iElemOld)*L_xi(i,ii,jj,kk)*L_eta_zeta
        END DO !i
      END DO !j
    END DO !k
  END DO; END DO; END DO
  CALL ChangeBasis3D(nVar_State,NInter,NNew,Vdm_CLNInter_GPNNew,Utmp,U(:,:,:,:,iElemNew))
END DO
!$OMP END DO
!$OMP END PARALLEL
Time=FLEXITIME() -Time
SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',Time,'s]'

END SUBROUTINE InterpolateSolution

END MODULE MOD_InterpolateSolution
