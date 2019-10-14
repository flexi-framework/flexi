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
#include "eos.h"

!===================================================================================================================================
!> Main module for the analytical block-Jacobian preconditioner. Contains initialization and the main routine, which calls the
!> surface and volume contributions from their own modules.
!===================================================================================================================================
MODULE MOD_Jac_Ex
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitJac_Ex
  MODULE PROCEDURE InitJac_Ex
END INTERFACE

INTERFACE Jac_Ex
  MODULE PROCEDURE Jac_Ex
END INTERFACE

INTERFACE FinalizeJac_Ex
  MODULE PROCEDURE FinalizeJac_Ex
END INTERFACE


PUBLIC::InitJac_Ex,Jac_Ex,FinalizeJac_Ex
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Inititializes the variables required for the analytical block-Jacobi preconditioner 
!===================================================================================================================================
SUBROUTINE InitJac_Ex()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Jac_Ex_Vars
USE MOD_ReadInTools        ,ONLY: GETLOGICAL
USE MOD_Interpolation_Vars ,ONLY: L_minus,L_plus
USE MOD_DG_Vars            ,ONLY: L_Hatminus,L_Hatplus
USE MOD_Mesh_Vars          ,ONLY: nElems
#if PARABOLIC
USE MOD_Jac_br2            ,ONLY: Build_BR2_SurfTerms
USE MOD_Precond_Vars       ,ONLY: HyperbolicPrecond
#endif
#if FV_ENABLED && FV_RECONSTRUCT
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                      :: i,j,iLocSide
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT ANALYTICAL BLOCK JACOBIAN...'

! Pre-computed variables needed in the surface integral jacobians
ALLOCATE(LL_minus(0:PP_N,0:PP_N), LL_plus(0:PP_N,0:PP_N))
DO j=0,PP_N
  DO i=0,PP_N
    LL_minus(i,j) = L_Hatminus(i)*L_minus(j)
    LL_plus(i,j)  = L_Hatplus(i) *L_plus(j)
  END DO
END DO 

#if PARABOLIC
ALLOCATE(L_mp(0:PP_N,6))
#if PP_dim==3
DO iLocSide=1,6
#else
DO iLocSide=2,5
#endif 
  SELECT CASE(iLocSide)
  CASE(XI_MINUS,ETA_MINUS,ZETA_MINUS)
    DO i=0,PP_N
      L_mp (i,iLocSide) = L_minus(i)
    END DO !i
  CASE(XI_PLUS,ETA_PLUS,ZETA_PLUS)
    DO i=0,PP_N
      L_mp (i,iLocSide) = L_plus(i)
    END DO !i
  END SELECT
END DO !iLocSide
#endif /*PARABOLIC*/

! Extended arrays for FV reconstruction: Lenghts and primitive variables for each element extended by one layer of the neighbouring
! elements
#if FV_ENABLED && FV_RECONSTRUCT
ALLOCATE(FV_sdx_XI_extended  (0:PP_N,0:PP_NZ,0:PP_N+1,nElems))          ! 1. / FV_dx_XI   Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_sdx_ETA_extended (0:PP_N,0:PP_NZ,0:PP_N+1,nElems))          ! 1. / FV_dx_ETA  Attention: storage order is (i,k,j,iElem)
#if PP_dim == 3
ALLOCATE(FV_sdx_ZETA_extended(0:PP_N,0:PP_N ,0:PP_N+1,nElems))          ! 1. / FV_dx_ZETA Attention: storage order is (i,j,k,iElem)
ALLOCATE(UPrim_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,1:nElems))! extended primitive solution vector
#else
ALLOCATE(UPrim_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1, 0:PP_NZ ,1:nElems))! extended primitive solution vector
#endif
#endif

#if PARABOLIC
IF(HyperbolicPrecond.EQV..FALSE.) THEN
  ALLOCATE(JacLiftingFlux(PP_nVarPrim,PP_nVarPrim,0:PP_N,0:PP_NZ,6))
  CALL Build_BR2_SurfTerms()
END IF
#endif /*PARABOLIC*/

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT BLOCK JACOBIAN DONE!'
END SUBROUTINE InitJac_Ex

!===================================================================================================================================
!> Main routine to compute the analytical Jacobi. Calls the corresponding volume and surface routines and applies the Jacobian of
!> the mapping at the end. 
!===================================================================================================================================
SUBROUTINE Jac_Ex(t,iElem,BJ,doVol,doSurf)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars             ,ONLY:nDOFVarElem
USE MOD_JacSurfInt                ,ONLY:DGJacSurfInt
#if PARABOLIC
USE MOD_Precond_Vars              ,ONLY:HyperbolicPrecond
USE MOD_Jac_br2                   ,ONLY:FillJacLiftingFlux
USE MOD_Jac_Ex_Vol                ,ONLY:DGVolIntGradJac
#endif
USE MOD_Jac_Ex_Vol                ,ONLY:DGVolIntJac
#if FV_ENABLED
USE MOD_FV_Vars                   ,ONLY:FV_Elems
USE MOD_Jac_Ex_Vol                ,ONLY:FVVolIntJac
#if PARABOLIC
USE MOD_Jac_Ex_Vol                ,ONLY:FVVolIntGradJac
#endif
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t                                         !< current simulation time
INTEGER,INTENT(IN) :: iElem                                     !< index of current element
LOGICAL,INTENT(IN) :: dovol, dosurf                             !< logicals indicating to do derivative of volume/surface-integral 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: BJ(1:nDOFVarElem,1:nDOFVarElem)             !< block-Jacobian of current element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: FVEM
!===================================================================================================================================
#if FV_ENABLED
FVEM = FV_Elems(iElem)
#else
FVEM = 0
#endif
#if PARABOLIC
IF (.NOT.(HyperbolicPrecond)) THEN
  CALL FillJacLiftingFlux(t,iElem)
END IF
#endif /*PARABOLIC*/
IF(doVol)THEN
  IF(FVEM.EQ.0)THEN
    CALL DGVolIntJac(BJ,iElem) !without sJ!      !d(F^a+F^v)/dU partial
#if FV_ENABLED
  ELSE
    CALL FVVolIntJac(BJ,iElem)
#endif
  END IF
#if PARABOLIC
  IF (.NOT.(HyperbolicPrecond)) THEN
    IF(FVEM.EQ.0)THEN
      CALL DGVolIntGradJac(BJ,iElem)               !d(F^v)/dQ*dQ/dU
#if FV_ENABLED
    ELSE
      CALL FVVolIntGradJac(BJ,iElem)               !d(F^v)/dQ*dQ/dU
#endif
    END IF
  END IF
#endif /*PARABOLIC*/
END IF!doVol
IF(doSurf) THEN
  CALL DGJacSurfInt(t,BJ,iElem) 
END IF
CALL Apply_sJ(BJ,iElem)

END SUBROUTINE Jac_Ex


!===================================================================================================================================
!> Transformation of the block-Jacobian to physical coordinates (multiplication with sJ)
!===================================================================================================================================
SUBROUTINE Apply_sJ(BJ,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Implicit_Vars ,ONLY:nDOFVarElem
USE MOD_Mesh_Vars     ,ONLY:sJ
#if FV_ENABLED
USE MOD_FV_Vars       ,ONLY:FV_Elems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: BJ(nDOFVarElem,nDOFVarElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: r,s,i,j,k
INTEGER :: FVEM
!===================================================================================================================================
#if FV_ENABLED
FVEM = FV_Elems(iElem)
#else
FVEM = 0
#endif
DO s=0,nDOFVarElem-1,PP_nVar
  r=0
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = -BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)*sJ(i,j,k,iElem,FVEM)
        r=r+PP_nVar
      END DO !i
    END DO !j
  END DO !k
END DO ! s

END SUBROUTINE Apply_sJ

!===================================================================================================================================
!> Deallocate global variables
!===================================================================================================================================
SUBROUTINE FinalizeJac_Ex()
! MODULES
USE MOD_Jac_ex_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(LL_minus)
SDEALLOCATE(LL_plus)
SDEALLOCATE(l_mp)
#if PARABOLIC
SDEALLOCATE(R_Minus)
SDEALLOCATE(R_Plus)
SDEALLOCATE(JacLiftingFlux)
#endif /*PARABOLIC*/
#if FV_ENABLED && FV_RECONSTRUCT
SDEALLOCATE(UPrim_extended)
SDEALLOCATE(FV_sdx_XI_extended)
SDEALLOCATE(FV_sdx_ETA_extended)
#if PP_dim == 3
SDEALLOCATE(FV_sdx_ZETA_extended)
#endif
#endif
END SUBROUTINE FinalizeJac_Ex

END MODULE MOD_Jac_Ex
