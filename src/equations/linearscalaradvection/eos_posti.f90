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

!==================================================================================================================================
!> Contains all the routines to calculate the (equation system and EOS dependent) conservative/primitive/derived quantities. 
!> Dependency table will be filled in here.
!==================================================================================================================================
MODULE MOD_EOS_Posti
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE GetMaskCons
  MODULE PROCEDURE GetMaskCons
END INTERFACE

INTERFACE GetMaskPrim
  MODULE PROCEDURE GetMaskPrim
END INTERFACE

INTERFACE GetMaskGrad
  MODULE PROCEDURE GetMaskGrad
END INTERFACE

#if FV_ENABLED && FV_RECONSTRUCT
INTERFACE AppendNeededPrims
  MODULE PROCEDURE AppendNeededPrims
END INTERFACE
PUBLIC :: AppendNeededPrims
#endif

PUBLIC :: GetMaskCons
PUBLIC :: GetMaskPrim
PUBLIC :: GetMaskGrad
PUBLIC :: CalcQuantities
!==================================================================================================================================

CONTAINS


#if FV_ENABLED && FV_RECONSTRUCT
SUBROUTINE AppendNeededPrims(mapCalc,mapCalc_FV,nVarCalc) 
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)      :: mapCalc(nVarTotalEOS)
INTEGER,INTENT(OUT)     :: mapCalc_FV(nVarTotalEOS)
INTEGER,INTENT(OUT)     :: nVarCalc
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar
!===================================================================================================================================
mapCalc_FV = mapCalc

! renumerate mapCalc_FV
nVarCalc   = 0
DO iVar=1,nVarTotalEOS
  IF (mapCalc_FV(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapCalc_FV(iVar) = nVarCalc
  END IF
END DO
END SUBROUTINE AppendNeededPrims
#endif

FUNCTION GetMaskCons() 
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars
USE MOD_Equation_Vars ,ONLY: StrVarNames
USE MOD_StringTools   ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER :: GetMaskCons(nVarTotalEOS)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar,iVar2
!===================================================================================================================================
GetMaskCons = 0
DO iVar=1,nVarTotalEOS
  DO iVar2=1,PP_nVar
    IF (STRICMP(StrVarNames(iVar2),DepNames(iVar))) THEN
      GetMaskCons(iVar) = 1
    END IF
  END DO
END DO
END FUNCTION GetMaskCons


FUNCTION GetMaskPrim() 
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER :: GetMaskPrim(nVarTotalEOS)
!===================================================================================================================================
GetMaskPrim = 0
END FUNCTION GetMaskPrim


FUNCTION GetMaskGrad()
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER :: GetMaskGrad(nVarTotalEOS)
!===================================================================================================================================
GetMaskGrad = DepTableEOS(:,0)
END FUNCTION GetMaskGrad


SUBROUTINE CalcQuantities(nVarCalc,nVal,iElems,mapCalc,UCalc,maskCalc,gradUx,gradUy,gradUz)
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarCalc
INTEGER,INTENT(IN) :: nVal(:)
INTEGER,INTENT(IN) :: iElems(:)
INTEGER,INTENT(IN) :: mapCalc(nVarTotalEOS)
INTEGER,INTENT(IN) :: maskCalc(nVarTotalEOS)
REAL,INTENT(OUT)   :: UCalc(PRODUCT(nVal),1:nVarCalc)
REAL,DIMENSION(1:PP_nVarPrim,PRODUCT(nVal)),INTENT(IN),OPTIONAL :: gradUx,gradUy,gradUz
!===================================================================================================================================
STOP 'Not available for linear scalar advection.'
END SUBROUTINE CalcQuantities


SUBROUTINE CalcDerivedQuantity(iVarCalc,DepName,nVarCalc,nVal,iElems,mapCalc,UCalc,gradUx,gradUy,gradUz)
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)            :: iVarCalc
CHARACTER(LEN=255),INTENT(IN) :: DepName
INTEGER,INTENT(IN)            :: nVarCalc
INTEGER,INTENT(IN)            :: nVal(:)
INTEGER,INTENT(IN)            :: iElems(:)
INTEGER,INTENT(IN)            :: mapCalc(nVarTotalEOS)
REAL,INTENT(INOUT)            :: UCalc(PRODUCT(nVal),1:nVarCalc)
REAL,DIMENSION(1:PP_nVarPrim,PRODUCT(nVal)),INTENT(IN),OPTIONAL :: gradUx,gradUy,gradUz
!===================================================================================================================================
STOP 'Not available for linear scalar advection.'
END SUBROUTINE CalcDerivedQuantity

END MODULE MOD_EOS_Posti
