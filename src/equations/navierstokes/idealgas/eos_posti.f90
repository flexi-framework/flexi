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
#include "eos.h"

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
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE FillDepTable
  MODULE PROCEDURE FillDepTable
END INTERFACE

#if FV_ENABLED && FV_RECONSTRUCT
INTERFACE AppendNeededPrims
  MODULE PROCEDURE AppendNeededPrims
END INTERFACE
PUBLIC :: AppendNeededPrims
#endif

INTERFACE GetMaskCons
  MODULE PROCEDURE GetMaskCons
END INTERFACE

INTERFACE GetMaskPrim
  MODULE PROCEDURE GetMaskPrim
END INTERFACE

INTERFACE GetMaskGrad
  MODULE PROCEDURE GetMaskGrad
END INTERFACE

INTERFACE FillPressure
  MODULE PROCEDURE FillPressure
END INTERFACE

INTERFACE FillTemperature
  MODULE PROCEDURE FillTemperature
END INTERFACE

INTERFACE CalcQuantities
  MODULE PROCEDURE CalcQuantities
END INTERFACE

PUBLIC :: FillDepTable
PUBLIC :: GetMaskCons
PUBLIC :: GetMaskPrim
PUBLIC :: GetMaskGrad
PUBLIC :: CalcQuantities
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Build the actual dependencies. If we call the DGTimeDerivative_weakForm, the primitive variables will be calculated
!> in there and don't need to be computed from the conservative variables.
!> Also, recursively add the dependencies of variabels that a variable depends on to the dependency of the first variable.
!==================================================================================================================================
SUBROUTINE FillDepTable(DepTable_In,withGradients)
USE MOD_EOS_Posti_Vars
USE MOD_Equation_Vars ,ONLY: StrVarNames,StrVarNamesPrim
USE MOD_StringTools   ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(INOUT) :: DepTable_In(1:nVarTotalEOS,0:nVarTotalEOS)
LOGICAL,INTENT(IN)    :: withGradients
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVar,ivar2,iVarPrim
!===================================================================================================================================


! for each quantity copy from all quantities that this quantity depends on the dependencies.
DO iVar=1,nVarTotalEOS
  DepTable_In(iVar,iVar) = 1
  DO iVar2=1,iVar-1
    IF (DepTable_In(iVar,iVar2).EQ.1) THEN
      DepTable_In(iVar,:) = MAX(DepTable_In(iVar,:), DepTable_In(iVar2,:))
    END IF
  END DO
END DO
END SUBROUTINE FillDepTable

#if FV_ENABLED && FV_RECONSTRUCT
SUBROUTINE AppendNeededPrims(mapCalc,mapCalc_FV,nVarCalc) 
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)      :: mapCalc(nVarTotalEOS)
INTEGER,INTENT(OUT)     :: mapCalc_FV(nVarTotalEOS)
INTEGER,INTENT(OUT)     :: nVarCalc
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar
!===================================================================================================================================
mapCalc_FV = mapCalc
DO iVar=1,PP_nVar
  IF (mapCalc(iVar).GT.0) THEN
    mapCalc_FV = MAX(mapCalc_FV, DepTablePrimToCons(iVar,1:nVarTotalEOS))
  END IF
END DO

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
USE MOD_EOS_Posti_Vars
USE MOD_Equation_Vars ,ONLY: StrVarNames
USE MOD_StringTools   ,ONLY: STRICMP
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
USE MOD_EOS_Posti_Vars
USE MOD_Equation_Vars ,ONLY: StrVarNamesPrim
USE MOD_StringTools   ,ONLY: STRICMP
! INPUT / OUTPUT VARIABLES 
INTEGER :: GetMaskPrim(nVarTotalEOS)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar,iVar2
!===================================================================================================================================
GetMaskPrim = 0
DO iVar=1,nVarTotalEOS
  DO iVar2=1,PP_nVarPrim
    IF (STRICMP(StrVarNamesPrim(iVar2),DepNames(iVar))) THEN
      GetMaskPrim(iVar) = 1
    END IF
  END DO
END DO
END FUNCTION GetMaskPrim


FUNCTION GetMaskGrad()
USE MOD_EOS_Posti_Vars
INTEGER :: GetMaskGrad(nVarTotalEOS)
!===================================================================================================================================
GetMaskGrad = DepTableEOS(:,0)
END FUNCTION GetMaskGrad



SUBROUTINE FillPressure(nElems_calc,indices,Nloc,nElems,U,Pressure)
!==================================================================================================================================
! MODULES
USE MOD_Eos_Vars  ,ONLY: KappaM1
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
INTEGER,INTENT(IN) :: nElems
REAL,INTENT(IN)    :: U(PP_nVar,0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(OUT)   :: Pressure( 0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER         :: i,j,k,iElem,iElem_calc
REAL            :: UE(PP_2Var)
!==================================================================================================================================
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    UE(SRHO) = 1./U(1,i,j,k,iElem)
    UE(MOMV) = U(2:4,i,j,k,iElem)
    UE(ENER) = U(5,i,j,k,iElem)
    UE(VELV) = UE(SRHO)*U(2:4,i,j,k,iElem)
    Pressure(i,j,k,iElem_calc) = PRESSURE_HE(UE) 
  END DO; END DO; END DO! i,j,k=0,Nloc
END DO ! iElem
END SUBROUTINE FillPressure

PURE SUBROUTINE FillTemperature(nElems,Nloc,Density,Pressure,Temperature)
!==================================================================================================================================
! MODULES
USE MOD_Eos_Vars  ,ONLY: R
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: Density    (0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(IN)    :: Pressure   (0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(OUT)   :: Temperature(0:Nloc,0:Nloc,0:Nloc,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER         :: i,j,k,iElem
REAL            :: UE(PP_2Var)
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    UE(SRHO) = 1./Density(i,j,k,iElem)
    UE(PRES) = Pressure(i,j,k,iElem)
    Temperature(i,j,k,iElem) = TEMPERATURE_HE(UE) 
  END DO; END DO; END DO! i,j,k=0,Nloc
END DO 
END SUBROUTINE FillTemperature

SUBROUTINE CalcQuantities(nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc,maskCalc) 
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EOS_Posti_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN) :: nVarCalc
INTEGER,INTENT(IN) :: Nloc
INTEGER,INTENT(IN) :: nElems_loc
INTEGER,INTENT(IN) :: iElems(nElems_loc)
INTEGER,INTENT(IN) :: mapCalc(nVarTotalEOS)
REAL,INTENT(OUT)   :: UCalc(0:Nloc,0:Nloc,0:Nloc,1:nElems_loc,1:nVarCalc)
INTEGER,INTENT(IN) :: maskCalc(nVarTotalEOS)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVar,iVarCalc
!===================================================================================================================================
DO iVar=1,nVarTotalEOS
  iVarCalc = mapCalc(iVar)
  IF (iVarCalc.GT.0 .AND. maskCalc(iVar).GT.0) THEN
    SWRITE(*,*) "  ",TRIM(DepNames(iVar))
    CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc)
  END IF
END DO
END SUBROUTINE CalcQuantities


SUBROUTINE CalcDerivedQuantity(iVarCalc,DepName,nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Posti_Vars
USE MOD_EOS_Vars        ,ONLY: cp,kappa,R,sKappaM1
USE MOD_StringTools     ,ONLY: LowCase,KEYVALUE
USE MOD_Mesh_Vars       ,ONLY: nElems
USE MOD_DG_Vars         ,ONLY: U
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)            :: iVarCalc
CHARACTER(LEN=255),INTENT(IN) :: DepName
INTEGER,INTENT(IN)            :: nVarCalc
INTEGER,INTENT(IN)            :: Nloc
INTEGER,INTENT(IN)            :: nElems_loc
INTEGER,INTENT(IN)            :: iElems(nElems_loc)
INTEGER,INTENT(IN)            :: mapCalc(nVarTotalEOS)
REAL,INTENT(OUT)              :: UCalc(0:Nloc,0:Nloc,0:Nloc,1:nElems_loc,1:nVarCalc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iMom,iDens,iPres,iVel1,iVel2,iVel3,iVelM,iVelS,iEner,iTemp
CHARACTER(LEN=255) :: DepName_low
#if PARABOLIC
INTEGER            :: iVor1,iVor2,iVor3,iVorM
INTEGER            :: i,j,k,iElem
REAL               :: Denominator
#endif
!===================================================================================================================================
CALL LowCase(DepName,DepName_low)
SELECT CASE(DepName_low)
  CASE("momentumx")
    iDens = KEYVALUE(DepNames,mapCalc,"density"    )
    iVel1 = KEYVALUE(DepNames,mapCalc,"velocityx"  )
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVel1)
  CASE("momentumy")
    iDens = KEYVALUE(DepNames,mapCalc,"density"    )
    iVel2 = KEYVALUE(DepNames,mapCalc,"velocityy"  )
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVel2)
  CASE("momentumz")
    iDens = KEYVALUE(DepNames,mapCalc,"density"    )
    iVel3 = KEYVALUE(DepNames,mapCalc,"velocityz"  )
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVel3)
  CASE("energystagnationdensity")
    iDens = KEYVALUE(DepNames,mapCalc,"density"    )
    iVel1 = KEYVALUE(DepNames,mapCalc,"velocityx"  )
    iVel2 = KEYVALUE(DepNames,mapCalc,"velocityy"  )
    iVel3 = KEYVALUE(DepNames,mapCalc,"velocityz"  )
    iPres = KEYVALUE(DepNames,mapCalc,"pressure"   )
    ! TODO: use a function from eos.f90
    UCalc(:,:,:,:,iVarCalc) = sKappaM1*UCalc(:,:,:,:,iPres) + 0.5*UCalc(:,:,:,:,iDens)* &
        (UCalc(:,:,:,:,iVel1)**2 + UCalc(:,:,:,:,iVel2)**2 + UCalc(:,:,:,:,iVel3)**2)
  CASE("velocityx")
    iDens = KEYVALUE(DepNames,mapCalc,'density'  )
    iMom  = KEYVALUE(DepNames,mapCalc,'momentumx')
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iMom) / UCalc(:,:,:,:,iDens)
  CASE("velocityy")
    iDens = KEYVALUE(DepNames,mapCalc,'density'  )
    iMom  = KEYVALUE(DepNames,mapCalc,'momentumy')
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iMom) / UCalc(:,:,:,:,iDens)
  CASE("velocityz")
    iDens = KEYVALUE(DepNames,mapCalc,'density'  )
    iMom  = KEYVALUE(DepNames,mapCalc,'momentumz')
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iMom) / UCalc(:,:,:,:,iDens)
  CASE("pressure")
    CALL FillPressure(nElems_loc,iElems,PP_N,nElems,U,UCalc(:,:,:,:,iVarCalc))
  CASE("temperature")
    iDens = KEYVALUE(DepNames,mapCalc,"density" )
    iPres = KEYVALUE(DepNames,mapCalc,"pressure")
    CALL FillTemperature(nElems_loc,PP_N,UCalc(:,:,:,:,iDens),UCalc(:,:,:,:,iPres),UCalc(:,:,:,:,iVarCalc))
  CASE("velocitymagnitude")
    iVel1 = KEYVALUE(DepNames,mapCalc,"velocityx")
    iVel2 = KEYVALUE(DepNames,mapCalc,"velocityy")
    iVel3 = KEYVALUE(DepNames,mapCalc,"velocityz")
    UCalc(:,:,:,:,iVarCalc) = SQRT(UCalc(:,:,:,:,iVel1)**2 + UCalc(:,:,:,:,iVel2)**2 + UCalc(:,:,:,:,iVel3)**2)
  CASE("velocitysound")
    iDens = KEYVALUE(DepNames,mapCalc,"density" )
    iPres = KEYVALUE(DepNames,mapCalc,"pressure")
    UCalc(:,:,:,:,iVarCalc) = SQRT(Kappa*UCalc(:,:,:,:,iPres)/UCalc(:,:,:,:,iDens))
  CASE("mach")
    iVelM = KEYVALUE(DepNames,mapCalc,"velocitymagnitude")
    iVelS = KEYVALUE(DepNames,mapCalc,"velocitysound"    )
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iVelM)/UCalc(:,:,:,:,iVelS)
  CASE("energystagnation")
    iDens = KEYVALUE(DepNames,mapCalc,"density"                )
    iEner = KEYVALUE(DepNames,mapCalc,"energystagnationdensity")
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iEner)/UCalc(:,:,:,:,iDens)
  CASE("enthalpystagnation")
    iDens = KEYVALUE(DepNames,mapCalc,"density"         )
    iPres = KEYVALUE(DepNames,mapCalc,"pressure"        )
    iEner = KEYVALUE(DepNames,mapCalc,"energystagnation")
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iEner) + UCalc(:,:,:,:,iPres)/UCalc(:,:,:,:,iDens)
  CASE("entropy")
    iDens = KEYVALUE(DepNames,mapCalc,"density"    )
    iTemp = KEYVALUE(DepNames,mapCalc,"temperature")
    UCalc(:,:,:,:,iVarCalc) = R*(sKappaM1*LOG(UCalc(:,:,:,:,iTemp))) - LOG(UCalc(:,:,:,:,iDens))
  CASE("totaltemperature")
    iTemp = KEYVALUE(DepNames,mapCalc,"temperature"      )
    iVelM = KEYVALUE(DepNames,mapCalc,"velocitymagnitude")
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iTemp)+UCalc(:,:,:,:,iVelM)**2/(2*cp)
  CASE("totalpressure")
    iDens = KEYVALUE(DepNames,mapCalc,"density"          )
    iPres = KEYVALUE(DepNames,mapCalc,"pressure"         )
    iVelM = KEYVALUE(DepNames,mapCalc,"velocitymagnitude")
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iPres)+0.5*UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVelM)**2
#if PARABOLIC      
  CASE("vorticityx")
    CALL FillVorticity(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc),1)
  CASE("vorticityy")
    CALL FillVorticity(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc),2)
  CASE("vorticityz")
    CALL FillVorticity(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc),3)
  CASE("vorticitymagnitude")
    iVor1 = KEYVALUE(DepNames,mapCalc,"vorticityx")
    iVor2 = KEYVALUE(DepNames,mapCalc,"vorticityy")
    iVor3 = KEYVALUE(DepNames,mapCalc,"vorticityy")
    UCalc(:,:,:,:,iVarCalc) = SQRT(UCalc(:,:,:,:,iVor1)**2 + UCalc(:,:,:,:,iVor2)**2 + UCalc(:,:,:,:,iVor3)**2)
  CASE("helicity")
    iVel1 = KEYVALUE(DepNames,mapCalc,"velocityx" )
    iVel2 = KEYVALUE(DepNames,mapCalc,"velocityy" )
    iVel3 = KEYVALUE(DepNames,mapCalc,"velocityz" )
    iVor1 = KEYVALUE(DepNames,mapCalc,"vorticityx")
    iVor2 = KEYVALUE(DepNames,mapCalc,"vorticityy")
    iVor3 = KEYVALUE(DepNames,mapCalc,"vorticityy")
    iVelM = KEYVALUE(DepNames,mapCalc,"velocitymagnitude" )
    iVorM = KEYVALUE(DepNames,mapCalc,"vorticitymagnitude")
    UCalc(:,:,:,:,iVarCalc) = SQRT( UCalc(:,:,:,:,iVor1)*UCalc(:,:,:,:,iVel1) &
        +UCalc(:,:,:,:,iVor2)*UCalc(:,:,:,:,iVel2) &
        +UCalc(:,:,:,:,iVor3)*UCalc(:,:,:,:,iVel3))
    DO iElem=1,nElems_loc
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Denominator = UCalc(i,j,k,iElem,iVelM) * UCalc(i,j,k,iElem,iVorM)
        UCalc(i,j,k,iElem,iVarCalc) = MERGE(0., UCalc(i,j,k,iElem,iVarCalc)/Denominator ,ABS(Denominator).LE.1e-15)
      END DO; END DO; END DO! i,j,k=0,PP_N
    END DO ! iElem
  CASE("lambda2")
    CALL FillLambda2(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc))
  CASE("dilatation")
    CALL FillDilatation(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc))
  CASE("qcriterion")
    CALL FillQcriterion(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc))
  CASE("schlieren")
    CALL FillSchlieren(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc))
#endif
END SELECT
END SUBROUTINE CalcDerivedQuantity

#if PARABOLIC
PURE SUBROUTINE FillVorticity(nElems_calc,indices,Nloc,Vorticity,dir)
!==================================================================================================================================
! MODULES
USE MOD_Lifting_Vars, ONLY: gradUx,gradUy,gradUz
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
INTEGER,INTENT(IN) :: dir
REAL,INTENT(OUT)   :: Vorticity(0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER         :: iElem,iElem_calc
!==================================================================================================================================
SELECT CASE (dir)
CASE(1) ! VorticityX = dw/dy-dv/dz
  DO iElem_calc=1,nElems_calc
    iElem = indices(iElem_calc)
    Vorticity(:,:,:,iElem_calc) = gradUy(4,:,:,:,iElem) - gradUz(3,:,:,:,iElem)
  END DO ! iElem
CASE(2) ! VorticityY = du/dz-dw/dx
  DO iElem_calc=1,nElems_calc
    iElem = indices(iElem_calc)
    Vorticity(:,:,:,iElem_calc) = gradUz(2,:,:,:,iElem) - gradUx(4,:,:,:,iElem)
  END DO ! iElem
CASE(3) ! VorticityZ = dv/dx-du/dy
  DO iElem_calc=1,nElems_calc
    iElem = indices(iElem_calc)
    Vorticity(:,:,:,iElem_calc) = gradUx(3,:,:,:,iElem) - gradUy(2,:,:,:,iElem)
  END DO ! iElem
END SELECT
END SUBROUTINE FillVorticity

SUBROUTINE FillLambda2(nElems_calc,indices,Nloc,Lambda2)
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Lifting_Vars, ONLY: gradUx,gradUy,gradUz
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(OUT)   :: Lambda2(0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: i,j,k,iElem,iElem_calc
REAL               :: gradUmat(3,3)
INTEGER            :: INFO
REAL,DIMENSION(3)  :: Lambda
REAL,DIMENSION(16) :: WORK
!==================================================================================================================================
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    gradUmat(:,1)= gradUx(2:4,i,j,k,iElem)
    gradUmat(:,2)= gradUy(2:4,i,j,k,iElem)
    gradUmat(:,3)= gradUz(2:4,i,j,k,iElem)
    gradUmat=MATMUL(0.5*(gradUmat+TRANSPOSE(gradUmat)),0.5*(gradUmat+TRANSPOSE(gradUmat))) & ! S^2
            +MATMUL(0.5*(gradUmat-TRANSPOSE(gradUmat)),0.5*(gradUmat-TRANSPOSE(gradUmat)))   !Omega^2
    ! Jacobi-Subroutine used from LAPACK 
    CALL DSYEV('N', 'U', 3, gradUmat, 3, Lambda, WORK, 16, INFO )
    Lambda2(i,j,k,iElem_calc) = Lambda(2)
  END DO; END DO; END DO! i,j,k=0,PP_N
END DO ! iElem
END SUBROUTINE FillLambda2

PURE SUBROUTINE FillDilatation(nElems_calc,indices,Nloc,Dilatation)
!==================================================================================================================================
! MODULES
USE MOD_Lifting_Vars, ONLY: gradUx,gradUy,gradUz
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(OUT)   :: Dilatation(0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iElem,iElem_calc
!==================================================================================================================================
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
  Dilatation(:,:,:,iElem_calc) = gradUx(2,:,:,:,iElem) + gradUy(3,:,:,:,iElem) + gradUz(4,:,:,:,iElem)  
END DO ! iElem
END SUBROUTINE FillDilatation

PURE SUBROUTINE FillQcriterion(nElems_calc,indices,Nloc,Qcriterion)
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Lifting_Vars, ONLY: gradUx,gradUy,gradUz
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(OUT)   :: Qcriterion(0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: i,j,k,l,m,iElem,iElem_calc
REAL               :: gradUmat(3,3),Q_loc,S,Rot
!==================================================================================================================================
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    gradUmat(:,1)= gradUx(2:4,i,j,k,iElem)
    gradUmat(:,2)= gradUy(2:4,i,j,k,iElem)
    gradUmat(:,3)= gradUz(2:4,i,j,k,iElem)
    Q_loc=0.
    DO m=1,3
      DO l=1,3
        ! shear rate
        S=0.5*(gradUmat(l,m)+gradUmat(m,l))
        ! rotation rate
        Rot=0.5*(gradUmat(l,m)-gradUmat(m,l))
        Q_loc=Q_loc+Rot*Rot-S*S
      END DO  ! l
    END DO  ! m
    Qcriterion(i,j,k,iElem_calc)=0.5*Q_loc
  END DO; END DO; END DO! i,j,k=0,PP_N
END DO ! iElem
END SUBROUTINE FillQcriterion

PURE SUBROUTINE FillSchlieren(nElems_calc,indices,Nloc,Schlieren)
!==================================================================================================================================
! MODULES
USE MOD_Lifting_Vars, ONLY: gradUx,gradUy,gradUz
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(OUT)   :: Schlieren(0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iElem,iElem_calc
!==================================================================================================================================
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
  Schlieren(:,:,:,iElem_calc)=LOG10(SQRT(gradUx(1,:,:,:,iElem)**2 + gradUy(1,:,:,:,iElem)**2 + gradUz(1,:,:,:,iElem)**2)+1.0)
END DO ! iElem
END SUBROUTINE FillSchlieren
#endif



END MODULE MOD_EOS_Posti
