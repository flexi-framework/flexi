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
USE MOD_SwapMesh_Vars,     ONLY: NState,NInter,NNew,NodeTypeState,RefState,nVar_State,NodeTypeOut
USE MOD_SwapMesh_Vars,     ONLY: Vdm_CLNInter_GPNNew
USE MOD_SwapMesh_Vars,     ONLY: UOld,xiInter,InterToElem,nElemsNew,IPDone
USE MOD_SwapMesh_Vars,     ONLY: Elem_IJK,ExtrudeTo3D,ExtrudeK
USE MOD_SwapMesh_Vars,     ONLY: ExtrudePeriodic,nElemsOld_IJK
USE MOD_DG_Vars,           ONLY: U
USE MOD_ChangeBasisByDim,  ONLY: ChangeBasisVolume
#if USE_OPENMP
USE OMP_Lib,               ONLY: OMP_GET_WTIME
#endif
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
REAL                          :: Utmp(1:nVar_State,0:NInter,0:NInter,0:ZDIM(NInter))
REAL                          :: L_xi(    0:NState,0:NInter,0:NInter,0:ZDIM(NInter))
REAL                          :: L_eta(   0:NState,0:NInter,0:NInter,0:ZDIM(NInter))
REAL                          :: L_zeta(  0:NState,0:NInter,0:NInter,0:NInter)
REAL                          :: L_eta_zeta
REAL                          :: xGP(0:NState),wBaryGP(0:NState)
REAL                          :: StartT,EndT
INTEGER                       :: jElemNew,jElemOld,iElemExtrusion,iIter
!===================================================================================================================================
! GPs and Barycentric weights for solution
CALL GetNodesAndWeights(NState,NodeTypeState,xGP)
CALL BarycentricWeights(NState,xGP,wBaryGP)

SWRITE(UNIT_stdOut,'(A,A)',ADVANCE='NO')' INTERPOLATE STATE TO NEW MESH ...',ACHAR(13)
StartT = FLEXITIME()
U      = 0.

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iElemNew,iElemOld,ii,jj,kk,L_xi,L_eta,L_zeta,L_eta_zeta,Utmp,i,j,k)
!$OMP DO
DO iElemNew=1,nElemsNew
  IF (ExtrudePeriodic) THEN
    ! Only the choosen layers are interpolated if an periodic extrusion is performed
    IF (Elem_IJK(3,iElemNew).GT.nElemsOld_IJK(3)) CYCLE
  END IF
  IF (ExtrudeTo3D) THEN
    ! Only the choosen layer is interpolated if an extrusion is performed
    IF (Elem_IJK(3,iElemNew).NE.ExtrudeK) CYCLE
  END IF

  ! Equal elements
  IF(equalElem(iElemNew).GT.0) THEN
    iElemOld=equalElem(iElemNew)
    IF((NState.EQ.Nnew).AND.(NodeTypeOut.EQ.NodetypeState))THEN
      U(:,:,:,:,iElemNew)=UOld(:,:,:,:,iElemOld)
    ELSE
      CALL ChangeBasisVolume(nVar_State,NState,NNew,Vdm_GPNState_GPNNew,UOld(:,:,:,:,iElemOld),U(:,:,:,:,iElemNew))
    END IF
    CYCLE
  END IF

  ! -- Remaining elements
  ! Calculate basis function values for xiInter
  DO kk=0,ZDIM(NInter); DO jj=0,NInter; DO ii=0,NInter
    IF(.NOT.IPDone(ii,jj,kk,iElemNew)) CYCLE
    CALL LagrangeInterpolationPolys(xiInter(1,ii,jj,kk,iElemNew),NState,xGP,wBaryGP,L_xi(  :,ii,jj,kk))
    CALL LagrangeInterpolationPolys(xiInter(2,ii,jj,kk,iElemNew),NState,xGP,wBaryGP,L_eta( :,ii,jj,kk))
#if PP_dim == 3
    CALL LagrangeInterpolationPolys(xiInter(3,ii,jj,kk,iElemNew),NState,xGP,wBaryGP,L_zeta(:,ii,jj,kk))
#endif
  END DO; END DO; END DO

  ! Evaluate old solution at interpolation points and transform to new state polynomial degree and node type
  DO kk=0,ZDIM(NInter); DO jj=0,NInter; DO ii=0,NInter
    IF(.NOT.IPDone(ii,jj,kk,iElemNew))THEN
      ! If the interpolation point was not found, overwrite with reference state
      ! If no reference state was given, the program already aborted
      Utmp(:,ii,jj,kk)=RefState
      CYCLE
    ELSE
      Utmp(:,ii,jj,kk)=0.
    END IF

    iElemOld = InterToElem(ii,jj,kk,iElemNew)
    DO k=0,ZDIM(NState)
      DO j=0,NState
#if PP_dim == 3
        L_eta_zeta=L_eta(j,ii,jj,kk)*L_zeta(k,ii,jj,kk)
#else
        L_eta_zeta=L_eta(j,ii,jj,kk)
#endif
        DO i=0,NState
          Utmp(:,ii,jj,kk) = Utmp(:,ii,jj,kk) + UOld(:,i,j,k,iElemOld)*L_xi(i,ii,jj,kk)*L_eta_zeta
        END DO !i
      END DO !j
    END DO !k
  END DO; END DO; END DO
  CALL ChangeBasisVolume(nVar_State,NInter,NNew,Vdm_CLNInter_GPNNew,Utmp,U(:,:,:,:,iElemNew))
END DO
!$OMP END DO
!$OMP END PARALLEL

! Copy the data from the extrusion layer to all other elements in the z direction
IF (ExtrudeTo3D) THEN
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iElemNew,iElemOld,L_xi,ii,jj,kk,L_eta,L_zeta,L_eta_zeta,Utmp,i,j,k)
!$OMP DO
  DO iElemNew=1,nElemsNew
    ! Skip the extrusion layer, already done
    IF (Elem_IJK(3,iElemNew).EQ.ExtrudeK) CYCLE
    ! Search for the corresponding element in the extrusion layer
    DO jElemNew = 1, nElemsNew
      IF (ALL(Elem_IJK(:,jElemNew).EQ.(/Elem_IJK(1,iElemNew),Elem_IJK(2,iElemNew),ExtrudeK/))) THEN
        iElemExtrusion = jElemNew
        EXIT
      END IF
    END DO ! jElemNew = 1, nElemsNew

    ! Copy data
    U(:,:,:,:,iElemNew) = U(:,:,:,:,iElemExtrusion)
  END DO
!$OMP END DO
!$OMP END PARALLEL
END IF

! Copy the data from the extrusion mesh to all other elements in the z direction
IF (ExtrudePeriodic) THEN
  DO iElemNew=1,nElemsNew
    IF (Elem_IJK(3,iElemNew).LE.nElemsOld_IJK(3)) CYCLE ! Skip the extrusion layer, already done
    ! Search for the corresponding element in the extrusion layer
    iIter = REAL(Elem_IJK(3,iElemNew)-1)/REAL(nElemsOld_IJK(3))
    DO jElemNew = 1, nElemsNew
      IF (ALL(Elem_IJK(:,jElemNew).EQ.(/Elem_IJK(1,iElemNew),Elem_IJK(2,iElemNew), Elem_IJK(3,iElemNew) -iIter*nElemsOld_IJK(3)/))) THEN
        iElemExtrusion = jElemNew
        EXIT
      END IF
    END DO ! jElemNew = 1, nElemsNew
    ! Copy data
    U(:,:,:,:,iElemNew) = U(:,:,:,:,iElemExtrusion)
  END DO
END IF

EndT = FLEXITIME()
CALL DisplayMessageAndTime(EndT-StartT, 'INTERPOLATE STATE TO NEW MESH DONE', DisplayLine=.FALSE.)

END SUBROUTINE InterpolateSolution

END MODULE MOD_InterpolateSolution
