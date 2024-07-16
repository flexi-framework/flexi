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
!> This module contains routines necessary for the Equation of State. For the linear scalar advection diffusion equation,
!> no EOS is needed. This module provides some dummy routines due to compatibility issues with the Navier-Stokes equation system.
!==================================================================================================================================
MODULE MOD_EOS
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
  MODULE PROCEDURE ConsToPrim_Side
  MODULE PROCEDURE ConsToPrim_Elem
  MODULE PROCEDURE ConsToPrim_Volume
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
  MODULE PROCEDURE PrimToCons_Side
  MODULE PROCEDURE PrimToCons_Elem
  MODULE PROCEDURE PrimToCons_Volume
END INTERFACE

INTERFACE ConsToEntropy
  MODULE PROCEDURE ConsToEntropy
  MODULE PROCEDURE ConsToEntropy_Volume
END INTERFACE

INTERFACE EntropyToCons
  MODULE PROCEDURE EntropyToCons
  MODULE PROCEDURE EntropyToCons_Side
END INTERFACE


PUBLIC::ConsToPrim,PrimToCons
PUBLIC::ConsToEntropy,EntropyToCons
PUBLIC::DefineParametersEos,InitEos
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters (dummy)
!==================================================================================================================================
PURE SUBROUTINE DefineParametersEos()
! MODULES
IMPLICIT NONE
!==================================================================================================================================
! dummy routine
END SUBROUTINE DefineParametersEos

!==================================================================================================================================
!> Initialize EOS (dummy)
!==================================================================================================================================
PURE SUBROUTINE InitEos()
! MODULES
IMPLICIT NONE
!==================================================================================================================================
! dummy routine
END SUBROUTINE InitEos


!==================================================================================================================================
!> Dummy routine, necessary for output routines due to compatibility to Navier-Stokes equation system
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim(prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(PP_nVar)     !< conservative variables
REAL,INTENT(OUT) :: prim(PP_nVarPrim) !< primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! copy cons to prim (PP_nVar = PP_nVarPrim)
prim = cons
END SUBROUTINE ConsToPrim

!==================================================================================================================================
!> Dummy routine, necessary for output routines due to compatibility to Navier-Stokes equation system
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim_Side(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                 !< polynomial degree
REAL,INTENT(IN)    :: cons(    PP_nVar,0:Nloc,0:ZDIM(Nloc))  !< conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))  !< primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! copy cons to prim (PP_nVar = PP_nVarPrim)
prim = cons
END SUBROUTINE ConsToPrim_Side

!==================================================================================================================================
!> Dummy routine, necessary for output routines due to compatibility to Navier-Stokes equation system
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim_Elem(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                           !< polynomial degree
REAL,INTENT(IN)    :: cons(    PP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc))   !< conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc))   !< primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! copy cons to prim (PP_nVar = PP_nVarPrim)
prim = cons
END SUBROUTINE ConsToPrim_Elem

!==================================================================================================================================
!> Dummy routine, necessary for output routines due to compatibility to Navier-Stokes equation system
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                                  !< polynomial degree
REAL,INTENT(IN)    :: cons(    PP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! copy cons to prim (PP_nVar = PP_nVarPrim)
prim = cons
END SUBROUTINE ConsToPrim_Volume

!==================================================================================================================================
!> Dummy routine, necessary for output routines due to compatibility to Navier-Stokes equation system
!==================================================================================================================================
PURE SUBROUTINE PrimToCons(prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: prim(PP_nVarPrim) !< vector of primitive variables
REAL,INTENT(OUT) :: cons(PP_nVar)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! copy prim to cons (PP_nVar = PP_nVarPrim)
cons = prim
END SUBROUTINE PrimToCons

!==================================================================================================================================
!> Dummy routine, necessary for output routines due to compatibility to Navier-Stokes equation system
!==================================================================================================================================
PURE SUBROUTINE PrimToCons_Side(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: prim(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! copy prim to cons (PP_nVar = PP_nVarPrim)
cons = prim
END SUBROUTINE PrimToCons_Side

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PURE SUBROUTINE PrimToCons_Elem(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                                  !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc))          !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc))          !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
cons = prim
END SUBROUTINE PrimToCons_Elem

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PURE SUBROUTINE PrimToCons_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                                      !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
cons = prim
END SUBROUTINE PrimToCons_Volume

!==================================================================================================================================
!> Transformation from conservative variables U to entropy vector, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PPURE SUBROUTINE ConsToEntropy(entropy,cons)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: cons    !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: entropy !< vector of entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
entropy = cons
END SUBROUTINE ConsToEntropy

!==================================================================================================================================
!> Transformation from entropy to conservative variables U, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PPURE SUBROUTINE EntropyToCons(entropy,cons)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)   :: entropy !< vector of entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT)  :: cons    !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: s,entropy2(PP_nVar),rhoe
!==================================================================================================================================
cons = entropy
END SUBROUTINE EntropyToCons

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PPURE SUBROUTINE ConsToEntropy_Volume(Nloc,entropy,cons)
! MODULES
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                                  !< local polynomial degree of solution representation
REAL,INTENT(OUT)   :: entropy(PP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)  !< vector of entropy variables
REAL,INTENT(IN)    :: cons(PP_nVar   ,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)  !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iElem
!==================================================================================================================================
entropy = cons
END SUBROUTINE ConsToEntropy_Volume

!> Transformation from primitive to conservative variables on a single side
!==================================================================================================================================
PPURE SUBROUTINE EntropyToCons_Side(Nloc,entropy,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                  !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: entropy(PP_nVar,0:Nloc,0:ZDIM(Nloc))  !< vector of entropy variables
REAL,INTENT(OUT)   :: cons(PP_nVar   ,0:Nloc,0:ZDIM(Nloc))  !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
cons = entropy
END SUBROUTINE EntropyToCons_Side

END MODULE MOD_EOS
