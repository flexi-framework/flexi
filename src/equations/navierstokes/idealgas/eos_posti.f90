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
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars,ONLY: nVarTotalEOS,DepTableEOS,DepNames
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
USE MOD_EOS_Posti_Vars,ONLY: nVarTotalEOS,DepTableEOS,DepNames
USE MOD_Equation_Vars ,ONLY: StrVarNamesPrim
USE MOD_StringTools   ,ONLY: STRICMP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
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
!==================================================================================================================================
! MODULES
USE MOD_EOS_Posti_Vars,ONLY: nVarTotalEOS,DepTableEOS
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
USE MOD_Globals,ONLY: MPIRoot
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
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: withGradients
INTEGER            :: iVar,iVarCalc
!===================================================================================================================================
withGradients=(PRESENT(gradUx).AND.PRESENT(gradUy).AND.PRESENT(gradUz))
print*, 'bin da',withGradients; CALL FLUSH()

DO iVar=1,nVarTotalEOS
  iVarCalc = mapCalc(iVar)
  IF (iVarCalc.GT.0 .AND. maskCalc(iVar).GT.0) THEN
    SWRITE(*,*) "  ",TRIM(DepNames(iVar))
    IF(withGradients)THEN
      CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,nVal,iElems,mapCalc,UCalc,gradUx,gradUy,gradUz)
    ELSE
      CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,nVal,iElems,mapCalc,UCalc)
    END IF
  END IF
END DO
END SUBROUTINE CalcQuantities


SUBROUTINE CalcDerivedQuantity(iVarCalc,DepName,nVarCalc,nVal,iElems,mapCalc,UCalc,gradUx,gradUy,gradUz)
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EOS_Posti_Vars
USE MOD_EOS_Vars
USE MOD_StringTools     ,ONLY: LowCase,KEYVALUE
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
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iMom1,iMom2,iMom3,iDens,iPres,iVel1,iVel2,iVel3,iVelM,iVelS,iEner,iEnst,iTemp
CHARACTER(LEN=255) :: DepName_low
REAL               :: UE(PP_2Var)
INTEGER            :: nElems_loc,Nloc,nDOF,nDims
#if PARABOLIC
INTEGER            :: iVorM
#endif
LOGICAL            :: withGradients
!===================================================================================================================================
nDims=UBOUND(nVal,1)
nElems_loc=nVal(nDims)
nDOF=PRODUCT(nVal(1:nDims-1))
Nloc=nVal(1)-1

!IF(.NOT.withGradients.AND.KEYVALUE(DepNames,DepTableEOS(:,0),DepName).EQ.1) &
!  STOP 'The selected variable requires the computation of gradients.'

!IF(dataSetDim.LT.KEYVALUE(DepNames,requiredDim,DepName))THEN
!  ! TODO: STOP or RETURN with warning
!  STOP 'Derived quantity cannot be computed for this dataset.'
!  !RETURN
!END IF

withGradients=(PRESENT(gradUx).AND.PRESENT(gradUy).AND.PRESENT(gradUz))

! Already retrieve mappings for conservatives, as theses are required for many quantities
iDens = KEYVALUE(DepNames,mapCalc,"density"    )
iMom1 = KEYVALUE(DepNames,mapCalc,'momentumx')
iMom2 = KEYVALUE(DepNames,mapCalc,'momentumy')
iMom3 = KEYVALUE(DepNames,mapCalc,'momentumz')
iVel1 = KEYVALUE(DepNames,mapCalc,"velocityx"  )
iVel2 = KEYVALUE(DepNames,mapCalc,"velocityy"  )
iVel3 = KEYVALUE(DepNames,mapCalc,"velocityz"  )
iPres = KEYVALUE(DepNames,mapCalc,"pressure"   )
iEner = KEYVALUE(DepNames,mapCalc,'energystagnationdensity')
iVorM = KEYVALUE(DepNames,mapCalc,"vorticitymagnitude")
iVelM = KEYVALUE(DepNames,mapCalc,"velocitymagnitude")
iTemp = KEYVALUE(DepNames,mapCalc,"temperature")
iVelS = KEYVALUE(DepNames,mapCalc,"velocitysound"    )
iEnst = KEYVALUE(DepNames,mapCalc,"energystagnation")

CALL LowCase(DepName,DepName_low)
SELECT CASE(DepName_low)
  CASE("momentumx")
    UCalc(:,iVarCalc) = UCalc(:,iDens)*UCalc(:,iVel1)
  CASE("momentumy")
    UCalc(:,iVarCalc) = UCalc(:,iDens)*UCalc(:,iVel2)
  CASE("momentumz")
    UCalc(:,iVarCalc) = UCalc(:,iDens)*UCalc(:,iVel3)
  CASE("energystagnationdensity")
    ! TODO: use a function from eos.f90
    UCalc(:,iVarCalc) = sKappaM1*UCalc(:,iPres) + 0.5*UCalc(:,iDens)* &
                        (UCalc(:,iVel1)**2 + UCalc(:,iVel2)**2 + UCalc(:,iVel3)**2)
  CASE("velocityx")
    UCalc(:,iVarCalc) = UCalc(:,iMom1) / UCalc(:,iDens)
  CASE("velocityy")
    UCalc(:,iVarCalc) = UCalc(:,iMom2) / UCalc(:,iDens)
  CASE("velocityz")
    UCalc(:,iVarCalc) = UCalc(:,iMom3) / UCalc(:,iDens)
  CASE("pressure")
    DO i=1,nDOF*nElems_loc
      UE(SRHO) = 1./UCalc(i,iDens)
      UE(MOM1) =    UCalc(i,iMom1)
      UE(MOM2) =    UCalc(i,iMom2)
      UE(MOM3) =    UCalc(i,iMom3)
      UE(ENER) =    UCalc(i,iEner)
      UE(VELV) = UE(SRHO)*UE(MOMV)
      UCalc(i,iVarCalc) = PRESSURE_HE(UE)
    END DO
  CASE("temperature")
    DO i=1,nDOF*nElems_loc
      UE(SRHO) = 1./UCalc(i,iDens)
      UE(PRES) =    UCalc(i,iPres)
      UCalc(i,iVarCalc) = TEMPERATURE_HE(UE)
    END DO
  CASE("velocitymagnitude")
    UCalc(:,iVarCalc) = SQRT((UCalc(:,iMom1)/UCalc(:,iDens))**2 &
                           + (UCalc(:,iMom2)/UCalc(:,iDens))**2 &
                           + (UCalc(:,iMom3)/UCalc(:,iDens))**2)
  CASE("velocitysound")
    UCalc(:,iVarCalc) = SQRT(Kappa*UCalc(:,iPres)/UCalc(:,iDens))
  CASE("mach")
    UCalc(:,iVarCalc) = UCalc(:,iVelM)/UCalc(:,iVelS)
  CASE("energystagnation")
    UCalc(:,iVarCalc) = UCalc(:,iEner)/UCalc(:,iDens)
  CASE("enthalpystagnation")
    UCalc(:,iVarCalc) = UCalc(:,iEnst) + UCalc(:,iPres)/UCalc(:,iDens)
  CASE("entropy")
    UCalc(:,iVarCalc) = R*(sKappaM1*LOG(UCalc(:,iTemp))) - LOG(UCalc(:,iDens))
  CASE("totaltemperature")
    UCalc(:,iVarCalc) = UCalc(:,iTemp)+UCalc(:,iVelM)**2/(2*cp)
  CASE("totalpressure")
    UCalc(:,iVarCalc) = UCalc(:,iPres)+0.5*UCalc(:,iDens)*UCalc(:,iVelM)**2
  CASE("pressuretimederiv")
     CALL FillPressureTimeDeriv(nElems_loc,iElems,Nloc,UCalc(:,iVarCalc))
#if PARABOLIC      
  CASE("vorticityx")
    UCalc(:,iVarCalc) = FillVorticity(1,nVal,gradUx,gradUy,gradUz)
  CASE("vorticityy")
    UCalc(:,iVarCalc) = FillVorticity(2,nVal,gradUx,gradUy,gradUz)
  CASE("vorticityz")
    UCalc(:,iVarCalc) = FillVorticity(3,nVal,gradUx,gradUy,gradUz)
  CASE("vorticitymagnitude")
    UCalc(:,iVarCalc) = SQRT(FillVorticity(1,nVal,gradUx,gradUy,gradUz)**2 &
                           + FillVorticity(2,nVal,gradUx,gradUy,gradUz)**2 &
                           + FillVorticity(3,nVal,gradUx,gradUy,gradUz)**2)
  CASE("helicity")
    UCalc(:,iVarCalc) = SQRT(FillVorticity(1,nVal,gradUx,gradUy,gradUz)*UCalc(:,iMom1)/UCalc(:,iDens) &
                            +FillVorticity(2,nVal,gradUx,gradUy,gradUz)*UCalc(:,iMom2)/UCalc(:,iDens) &
                            +FillVorticity(3,nVal,gradUx,gradUy,gradUz)*UCalc(:,iMom3)/UCalc(:,iDens))
    UCalc(:,iVarCalc) = MERGE(0., UCalc(:,iVarCalc)/(UCalc(:,iVelM)*UCalc(:,iVorM)),ABS(UCalc(:,iVelM)*UCalc(:,iVorM)).LE.1e-12)
  CASE("lambda2")
    UCalc(:,iVarCalc) = FillLambda2(nVal,gradUx,gradUy,gradUz)
  CASE("dilatation")
    UCalc(:,iVarCalc) = gradUx(2,:) + gradUy(3,:) + gradUz(4,:)  
  CASE("qcriterion")
    UCalc(:,iVarCalc) = FillQcriterion(nVal,gradUx,gradUy,gradUz)
  CASE("schlieren")
    UCalc(:,iVarCalc) = LOG10(SQRT(gradUx(1,:)**2 + gradUy(1,:)**2 + gradUz(1,:)**2)+1.0)
#endif
END SELECT
END SUBROUTINE CalcDerivedQuantity


SUBROUTINE FillPressureTimeDeriv(nElems_calc,indices,Nloc,PressureTDeriv)
!==================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY:nElems
USE MOD_Eos_Vars,ONLY: KappaM1
USE MOD_DG_Vars ,ONLY:Ut,U
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(OUT)   :: PressureTDeriv(0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iElem,iElem_calc
!==================================================================================================================================
print*,nElems,nElems_calc
IF(Nloc.NE.PP_N) STOP 'Not possible here'
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
IF(iElem.NE.iElem_calc) print*,'passt nicht'
  PressureTDeriv(:,:,:,iElem_calc)=KappaM1*(Ut(5,:,:,:,iElem)-1/U(1,:,:,:,iElem)*(  &
                                             U(2,:,:,:,iElem)*Ut(2,:,:,:,iElem)  &
                                           + U(3,:,:,:,iElem)*Ut(3,:,:,:,iElem)  &
                                           + U(4,:,:,:,iElem)*Ut(4,:,:,:,iElem)) &
                                       + 0.5/U(1,:,:,:,iElem)**2*Ut(1,:,:,:,iElem)*(  &
                                             U(2,:,:,:,iElem)*U(2,:,:,:,iElem)   &
                                           + U(3,:,:,:,iElem)*U(3,:,:,:,iElem)   &
                                           + U(4,:,:,:,iElem)*U(4,:,:,:,iElem)))
END DO ! iElem
END SUBROUTINE FillPressureTimeDeriv


#if PARABOLIC
FUNCTION FillVorticity(dir,nVal,gradUx,gradUy,gradUz) RESULT(Vorticity)
!==================================================================================================================================
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: dir
INTEGER,INTENT(IN) :: nVal(:)
REAL,DIMENSION(PP_nVarPrim,PRODUCT(nVal)),INTENT(IN) :: gradUx,gradUy,gradUz
REAL               :: Vorticity(PRODUCT(nVal))
!==================================================================================================================================
SELECT CASE (dir)
CASE(1) ! VorticityX = dw/dy-dv/dz
  Vorticity = gradUy(4,:) - gradUz(3,:)
CASE(2) ! VorticityY = du/dz-dw/dx
  Vorticity = gradUz(2,:) - gradUx(4,:)
CASE(3) ! VorticityZ = dv/dx-du/dy
  Vorticity = gradUx(3,:) - gradUy(2,:)
END SELECT
END FUNCTION FillVorticity



FUNCTION FillLambda2(nVal,gradUx,gradUy,gradUz) RESULT(Lambda2)
!==================================================================================================================================
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVal(:)
REAL,DIMENSION(PP_nVarPrim,PRODUCT(nVal)),INTENT(IN) :: gradUx,gradUy,gradUz
REAL               :: Lambda2(PRODUCT(nVal))
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: i
REAL               :: gradUmat(3,3)
INTEGER            :: INFO
REAL,DIMENSION(3)  :: Lambda
REAL,DIMENSION(16) :: WORK
!==================================================================================================================================
DO i=1,PRODUCT(nVal)
  gradUmat(:,1)= gradUx(2:4,i)
  gradUmat(:,2)= gradUy(2:4,i)
  gradUmat(:,3)= gradUz(2:4,i)
  gradUmat=MATMUL(0.5*(gradUmat+TRANSPOSE(gradUmat)),0.5*(gradUmat+TRANSPOSE(gradUmat))) & ! S^2
          +MATMUL(0.5*(gradUmat-TRANSPOSE(gradUmat)),0.5*(gradUmat-TRANSPOSE(gradUmat)))   !Omega^2
  ! Jacobi-Subroutine used from LAPACK 
  CALL DSYEV('N', 'U', 3, gradUmat, 3, Lambda, WORK, 16, INFO )
  Lambda2(i) = Lambda(2)
END DO
END FUNCTION FillLambda2


FUNCTION FillQcriterion(nVal,gradUx,gradUy,gradUz) RESULT(Qcriterion)
!==================================================================================================================================
! MODULES
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVal(:)
REAL,DIMENSION(PP_nVarPrim,PRODUCT(nVal)),INTENT(IN) :: gradUx,gradUy,gradUz
REAL               :: Qcriterion(PRODUCT(nVal))
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: i,m,l
REAL               :: gradUmat(3,3),Q_loc,S,Rot
!==================================================================================================================================
DO i=1,PRODUCT(nVal)
  gradUmat(:,1)= gradUx(2:4,i)
  gradUmat(:,2)= gradUy(2:4,i)
  gradUmat(:,3)= gradUz(2:4,i)
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
  Qcriterion(i)=0.5*Q_loc
END DO
END FUNCTION FillQcriterion
#endif



END MODULE MOD_EOS_Posti
