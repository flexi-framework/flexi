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

!==================================================================================================================================
!> Compute and integrate force from fluid onto wall surfaces (e.g. adiabatic, isothermal, Euler walls)
!==================================================================================================================================
MODULE MOD_CalcBodyForces
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CalcBodyForces
  MODULE PROCEDURE CalcBodyForces
END INTERFACE

PUBLIC :: CalcBodyForces
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Control routine for CalcBodyforces
!==================================================================================================================================
SUBROUTINE CalcBodyForces(BodyForce,Fp,Fv)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars,         ONLY:UPrim_master
#if PARABOLIC
USE MOD_Lifting_Vars,    ONLY:gradUx_master,gradUy_master,gradUz_master
#endif
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem,nBCSides,BC,nBCs
USE MOD_AnalyzeEquation_Vars,ONLY:isWall
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)               :: Fp(3,nBCs)              !< integrated pressure force per wall BC
REAL,INTENT(OUT)               :: Fv(3,nBCs)              !< integrated friction force per wall BC
REAL,INTENT(OUT)               :: BodyForce(3,nBCs)       !< Sum of pressure/friction force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Fp_loc(3)
#if PARABOLIC
REAL                           :: Fv_loc(3)
#endif
INTEGER                        :: SideID,iBC
#if USE_MPI
REAL                           :: Box(6,nBCs)
#endif /*USE_MPI*/
!==================================================================================================================================
! Calculate body forces  ! Attention: during the initialization phase no face data / gradients available!

Fp=0.
Fv=0.
BodyForce=0.
DO SideID=1,nBCSides
  iBC=BC(SideID)
  IF(.NOT.isWall(iBC)) CYCLE
  ! Calculate pressure force (Euler wall / Navier-Stokes wall)
  CALL CalcPressureForce(Fp_loc,UPrim_master(5,:,:,SideID),SurfElem(:,:,0,SideID),NormVec(:,:,:,0,SideID))
  Fp(:,iBC)=Fp(:,iBC)+Fp_loc
#if PARABOLIC
  ! Calculate viscous force (Navier-Stokes wall)
  CALL CalcViscousForce(Fv_loc,                      &
                        UPrim_master(:,:,:,SideID),  &
                        gradUx_master(:,:,:,SideID), &
                        gradUy_master(:,:,:,SideID), &
                        gradUz_master(:,:,:,SideID), &
                        SurfElem(:,:,0,SideID),      &
                        NormVec(:,:,:,0,SideID))
  Fv(:,iBC)=Fv(:,iBC)+Fv_loc
#endif
END DO

#if USE_MPI
Box(1:3,1:nBCs)=Fv; Box(4:6,1:nBCs)=Fp
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Box,6*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  Fv=Box(1:3,1:nBCs); Fp=Box(4:6,1:nBCs)
  BodyForce=Fv+Fp
ELSE
  CALL MPI_REDUCE(Box         ,0  ,6*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif

END SUBROUTINE CalcBodyForces



!==================================================================================================================================
!> Compute integral pressure force per face
!==================================================================================================================================
SUBROUTINE CalcPressureForce(Fp,p_Face,SurfElem,NormVec)
! MODULES
USE MOD_PreProc
USE MOD_Analyze_Vars,      ONLY:wGPSurf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)               :: p_Face(0:PP_N,0:PP_NZ)        !< (IN) pressure on face
REAL, INTENT(IN)               :: SurfElem(0:PP_N,0:PP_NZ)      !< (IN) face surface
REAL, INTENT(IN)               :: NormVec(3,0:PP_N,0:PP_NZ)     !< (IN) face normal vectors
REAL, INTENT(OUT)              :: Fp(3)                        !< (OUT) integrated pressure force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: dA
INTEGER                        :: i, j
!==================================================================================================================================
Fp=0.
DO j=0,PP_NZ; DO i=0,PP_N
  dA=wGPSurf(i,j)*SurfElem(i,j)
  Fp=Fp+p_Face(i,j)*NormVec(:,i,j)*dA
END DO; END DO
END SUBROUTINE CalcPressureForce


#if PARABOLIC
!==================================================================================================================================
!> Compute integral viscous force per face (only if compiled with parabolic terms)
!==================================================================================================================================
SUBROUTINE CalcViscousForce(Fv,UPrim_Face,gradUx_Face,gradUy_Face,gradUz_Face,SurfElem,NormVec)
! MODULES
USE MOD_PreProc
USE MOD_Viscosity
USE MOD_Analyze_Vars, ONLY:wGPSurf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)               :: UPrim_Face( PP_nVarPrim,0:PP_N,0:PP_NZ) !< (IN) primitive solution on face
REAL, INTENT(IN)               :: gradUx_Face(PP_nVarPrim,0:PP_N,0:PP_NZ) !< (IN) sln. gradients x-dir on face
REAL, INTENT(IN)               :: gradUy_Face(PP_nVarPrim,0:PP_N,0:PP_NZ) !< (IN) sln. gradients y-dir on face
REAL, INTENT(IN)               :: gradUz_Face(PP_nVarPrim,0:PP_N,0:PP_NZ) !< (IN) sln. gradients z-dir on face
REAL, INTENT(IN)               :: SurfElem(0:PP_N,0:PP_NZ)                !< (IN) face surface
REAL, INTENT(IN)               :: NormVec(3,0:PP_N,0:PP_NZ)               !< (IN) face normal vectors
REAL, INTENT(OUT)              :: Fv(3)                                  !< (OUT) integrated pressure force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: tau(3,3)                  ! Viscous stress tensor
REAL                           :: muS
REAL                           :: GradV(3,3),DivV,prim(PP_nVarPrim)
INTEGER                        :: i, j
!==================================================================================================================================
Fv       =0.

DO j=0,PP_NZ; DO i=0,PP_N
  ! calculate viscosity
  prim = UPrim_Face(:,i,j)
  muS=VISCOSITY_PRIM(prim)

  ! velocity gradients
  GradV(:,1)=gradUx_Face(LIFT_VELV,i,j)
  GradV(:,2)=gradUy_Face(LIFT_VELV,i,j)
#if PP_dim==3
  GradV(:,3)=gradUz_Face(LIFT_VELV,i,j)
#else
  GradV(:,3)=0.
#endif

  ! Velocity divergence
  DivV=GradV(1,1)+GradV(2,2)+GradV(3,3)
  ! Calculate shear stress tensor
  tau=muS*(GradV + TRANSPOSE(GradV))
  tau(1,1)=tau(1,1)-2./3.*muS*DivV
  tau(2,2)=tau(2,2)-2./3.*muS*DivV
#if PP_dim==3
  tau(3,3)=tau(3,3)-2./3.*muS*DivV
#endif
  ! Calculate viscous force vector
  Fv=Fv+MATMUL(tau,NormVec(:,i,j))*wGPSurf(i,j)*SurfElem(i,j)
END DO; END DO

Fv=-Fv  ! Change direction to get the force acting on the wall

END SUBROUTINE CalcViscousForce
#endif /*PARABOLIC*/

END MODULE MOD_CalcBodyForces
