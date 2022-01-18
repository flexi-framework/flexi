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

!==================================================================================================================================
!> Routines providing initialization and initial solutions for the linear advection-diffusion equation
!==================================================================================================================================
MODULE MOD_Equation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE

INTERFACE GetPrimitiveStateSurface
  MODULE PROCEDURE GetPrimitiveStateSurface
END INTERFACE

INTERFACE GetConservativeStateSurface
  MODULE PROCEDURE GetConservativeStateSurface
END INTERFACE

INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE


PUBLIC::InitEquation,FinalizeEquation
PUBLIC:: GetPrimitiveStateSurface,GetConservativeStateSurface
PUBLIC::DefineParametersEquation
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateRealArrayOption('AdvVel',       "Advection velocity for advection part of LinAdv-Diff.")
CALL prms%CreateRealOption(     'DiffC',        "Diffusion constant for diffusion part of LinAdv-Diff.")
END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Read equation parameters (advection velocity, diffusion coeff, exact function)  from the ini file
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,        ONLY:GETREALARRAY,GETREAL
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Exactfunc          ,ONLY:InitExactFunc
USE MOD_Equation_Vars
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitLinearScalarAdvection not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SCALAR LINADV...'

! Read the velocity vector from ini file
AdvVel = GETREALARRAY('AdvVel',3)
#if PP_dim==2
! Make sure advection velocity is 0 in third dimension for two-dimensional computations,
! computing wave speeds etc. will get easier.
IF(AdvVel(3).NE.0.) THEN
  SWRITE(UNIT_stdOut,'(A)')' You are computing in 2D! AdvVel(3) will be set to zero!'
  AdvVel(3) = 0.
END IF
#endif
! Read the diffusion constant from ini file
DiffC  = GETREAL('DiffC','0.')

! Call initialization of exactfunc
CALL InitExactFunc()

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LINADV DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEquation


!==================================================================================================================================
!> Converts conservative solution vector to primitive variables
!==================================================================================================================================
SUBROUTINE GetPrimitiveStateSurface(U_master,U_slave,UPrim_master,UPrim_slave)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY:nSides
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_master(        PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(IN)  :: U_slave(         PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
REAL,INTENT(OUT) :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(OUT) :: UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Copy the coservative state to the primitive arrays
UPrim_slave = U_slave
UPrim_master = U_master
END SUBROUTINE GetPrimitiveStateSurface


!==================================================================================================================================
!> Converts primite solution vector to conservative variables
!==================================================================================================================================
SUBROUTINE GetConservativeStateSurface(UPrim_master,UPrim_slave,U_master,U_slave, mask_master, mask_slave, mask_ref)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY: nSides
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(IN)    :: UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
REAL,INTENT(OUT)   :: U_master(        PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(OUT)   :: U_slave(         PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
INTEGER,INTENT(IN) :: mask_master(1:nSides)                             !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_slave (1:nSides)                             !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_ref                                          !< reference value for mask comparison
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Copy the coservative state to the primitive arrays
U_master = UPrim_master
U_slave  = UPrim_slave
END SUBROUTINE GetConservativeStateSurface

!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars,ONLY:EquationInitIsDone
IMPLICIT NONE
!==================================================================================================================================
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation

