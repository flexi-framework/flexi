!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz 
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
!> Contains the routines for the calculation of the analytical flux jacobians of the different equation systems
!===================================================================================================================================
MODULE MOD_Jacobian
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE dConsdPrim
  MODULE PROCEDURE dConsdPrim
END INTERFACE

INTERFACE dPrimdCons
  MODULE PROCEDURE dPrimdCons
END INTERFACE

INTERFACE dConsdPrimTemp
  MODULE PROCEDURE dConsdPrim
END INTERFACE

INTERFACE dPrimTempdCons
  MODULE PROCEDURE dPrimdCons
END INTERFACE

PUBLIC::EvalAdvFluxJacobian
#if PARABOLIC
PUBLIC::EvalDiffFluxJacobian
PUBLIC::EvalFluxGradJacobian
#endif
PUBLIC::dConsdPrim,dPrimdCons,dConsdPrimTemp,dPrimTempdCons
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Linear Scalar Advection Difussion:
!> The Jacobian of the analytical advective Flux with respect to the Variable U
!===================================================================================================================================
SUBROUTINE EvalAdvFluxJacobian(U,UPrim,fJac,gJac,hJac)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY:AdvVel
USE MOD_DG_Vars           ,ONLY:nDOFElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar,nDOFElem),INTENT(IN)              :: U
REAL,DIMENSION(PP_nVarPrim,nDOFElem),INTENT(IN)          :: UPrim
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar,nDOFElem),INTENT(OUT) :: fJac,gJac,hJac             ! Cartesian fluxes (iVar,i,j,k)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
fJac = AdvVel(1)
gJac = AdvVel(2)
#if PP_dim==3
hJac = AdvVel(3)
#endif
END SUBROUTINE EvalAdvFluxJacobian


#if PARABOLIC
!===================================================================================================================================
!> The Jacobian of the diffusion flux with respect to the conservative variables U
!===================================================================================================================================
SUBROUTINE EvalDiffFluxJacobian(nDOF_loc,U,UPrim,gradUx,gradUy,gradUz,fJac,gJac,hJac &
#if EDDYVISCOSITY
                                ,muSGS &
#endif
                                )
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------  
!----------------------------------------------------------------------------------------------------------------------------------  
INTEGER,INTENT(IN)                                   :: nDOF_loc             !< number of degrees of freedom
REAL,DIMENSION(PP_nVar        ,nDOF_loc),INTENT(IN)  :: U                    !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim    ,nDOF_loc),INTENT(IN)  :: UPrim                !< solution in primitive variables
REAL,DIMENSION(PP_nVarPrim    ,nDOF_loc),INTENT(IN)  :: gradUx,gradUy,gradUz !< primitive gradients
REAL,DIMENSION(PP_nVar,PP_nVar,nDOF_loc),INTENT(OUT) :: fJac,gJac,hJac       !< Derivative of the Cartesian fluxes (iVar,i,j,k)
#if EDDYVISCOSITY
REAL,DIMENSION(1              ,nDOF_loc),INTENT(IN)  :: muSGS                !< eddyviscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
fJac = 0.
gJac = 0.
hJac = 0.
END SUBROUTINE EvalDiffFluxJacobian

!===================================================================================================================================
!> Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
SUBROUTINE EvalFluxGradJacobian(nDOF_loc,U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz &
#if EDDYVISCOSITY
                               ,muSGS &
#endif
                               )
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                              :: nDOF_loc !< number of degrees of freedom
REAL,DIMENSION(PP_nVar    ,nDOF_loc),INTENT(IN) :: U        !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim,nDOF_loc),INTENT(IN) :: UPrim    !< solution in primitive variables
#if EDDYVISCOSITY
REAL,DIMENSION(1          ,nDOF_loc),INTENT(IN) :: muSGS    !< eddyviscosity
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar,nDOF_loc),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz !<
                                                        !> Jacobian of the diffusive Cartesian fluxes (iVar,i,j,k) w.r.t gradients
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
fJacQx = -DiffC
fJacQy = 0.
fJacQz = 0.

gJacQx = 0.
gJacQy = -DiffC
gJacQz = 0.

#if PP_dim==3
hJacQx = 0.
hJacQy = 0.
hJacQz = -DiffC
#endif
END SUBROUTINE EvalFluxGradJacobian
#endif /*PARABOLIC*/

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive variables
!===================================================================================================================================
SUBROUTINE dConsdPrim(UPrim,Jac)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: Jac      !< cons to prim Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Jac = 1.
END SUBROUTINE dConsdPrim

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive variables
!===================================================================================================================================
SUBROUTINE dPrimdCons(UPrim,Jac)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: Jac      !< prim to cons Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Jac = 1.
END SUBROUTINE dPrimdCons

END MODULE MOD_Jacobian
