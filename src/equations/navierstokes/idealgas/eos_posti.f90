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
!==================================================================================================================================
!> For FV, the conservative variables must be calculated from the primitive ones. In this routine the primitive variables
!> needed to calculate all conservative variables that are needed are added to the dep table.
!==================================================================================================================================
SUBROUTINE AppendNeededPrims(mapDepToCalc,mapDepToCalc_FV,nVarCalc) 
! MODULES
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)      :: mapDepToCalc(nVarDepEOS)
INTEGER,INTENT(OUT)     :: mapDepToCalc_FV(nVarDepEOS)
INTEGER,INTENT(OUT)     :: nVarCalc
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar
!===================================================================================================================================
mapDepToCalc_FV = mapDepToCalc
DO iVar=1,PP_nVar
  IF (mapDepToCalc(iVar).GT.0) THEN
    mapDepToCalc_FV = MAX(mapDepToCalc_FV, DepTablePrimToCons(iVar,1:nVarDepEOS))
  END IF
END DO

! renumerate mapDepToCalc_FV
nVarCalc   = 0
DO iVar=1,nVarDepEOS
  IF (mapDepToCalc_FV(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapDepToCalc_FV(iVar) = nVarCalc
  END IF
END DO
END SUBROUTINE AppendNeededPrims
#endif


!==================================================================================================================================
!> Create a mask for the conservative variables. The mask has a length of nVarDepEOS and is 1 at the index of the
!> conservative variables and 0 everywhere else.
!==================================================================================================================================
FUNCTION GetMaskCons() 
! MODULES
USE MOD_EOS_Posti_Vars,ONLY: nVarDepEOS,DepTableEOS,DepNames
USE MOD_Equation_Vars ,ONLY: StrVarNames
USE MOD_StringTools   ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER :: GetMaskCons(nVarDepEOS)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar,iVar2
!===================================================================================================================================
GetMaskCons = 0
DO iVar=1,nVarDepEOS
  DO iVar2=1,PP_nVar
    IF (STRICMP(StrVarNames(iVar2),DepNames(iVar))) THEN
      GetMaskCons(iVar) = 1
    END IF
  END DO
END DO
END FUNCTION GetMaskCons


!==================================================================================================================================
!> Create a mask for the primitive variables. The mask has a length of nVarDepEOS and is 1 at the index of the
!> primitive variables and 0 everywhere else.
!==================================================================================================================================
FUNCTION GetMaskPrim() 
! MODULES
USE MOD_EOS_Posti_Vars,ONLY: nVarDepEOS,DepTableEOS,DepNames
USE MOD_Equation_Vars ,ONLY: StrVarNamesPrim
USE MOD_StringTools   ,ONLY: STRICMP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER :: GetMaskPrim(nVarDepEOS)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar,iVar2
!===================================================================================================================================
GetMaskPrim = 0
DO iVar=1,nVarDepEOS
  DO iVar2=1,PP_nVarPrim
    IF (STRICMP(StrVarNamesPrim(iVar2),DepNames(iVar))) THEN
      GetMaskPrim(iVar) = 1
    END IF
  END DO
END DO
END FUNCTION GetMaskPrim


!==================================================================================================================================
!> Create a mask for the variables that need gradients. The mask has a length of nVarDepEOS and is 1 at the index of the
!> variables that need gradients and 0 everywhere else.
!==================================================================================================================================
FUNCTION GetMaskGrad()
! MODULES
USE MOD_EOS_Posti_Vars,ONLY: nVarDepEOS,DepTableEOS
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER :: GetMaskGrad(nVarDepEOS)
!===================================================================================================================================
GetMaskGrad = DepTableEOS(:,0)
END FUNCTION GetMaskGrad


!==================================================================================================================================
!> Wrapper routine that is called when derived quantities should be calculated. For every variable that should be calculated,
!> the routine that does the actual calculation CalcDerivedQuantity is called with the appropriate arguments.
!==================================================================================================================================
SUBROUTINE CalcQuantities(nVarCalc,nVal,mapCalcMeshToGlobalMesh,mapDepToCalc,UCalc,maskCalc,gradUx,gradUy,gradUz,&
    NormVec,TangVec1,TangVec2)
! MODULES
USE MOD_Globals
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                              :: nVarCalc
INTEGER,INTENT(IN)                                              :: nVal(:)
INTEGER,INTENT(IN)                                              :: mapCalcMeshToGlobalMesh(:)
INTEGER,INTENT(IN)                                              :: mapDepToCalc(nVarDepEOS)
INTEGER,INTENT(IN)                                              :: maskCalc(nVarDepEOS)
REAL,INTENT(OUT)                                                :: UCalc(PRODUCT(nVal),1:nVarCalc)
REAL,DIMENSION(1:PP_nVarPrim,PRODUCT(nVal)),INTENT(IN),OPTIONAL :: gradUx,gradUy,gradUz
REAL,DIMENSION(1:3,PRODUCT(nVal)),INTENT(IN),OPTIONAL           :: NormVec,TangVec1,TangVec2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: withGradients
LOGICAL            :: withVectors
INTEGER            :: iVar,iVarCalc
!===================================================================================================================================
withGradients=(PRESENT(gradUx).AND.PRESENT(gradUy).AND.PRESENT(gradUz))
withVectors=(PRESENT(NormVec).AND.PRESENT(TangVec1).AND.PRESENT(TangVec2))

DO iVar=1,nVarDepEOS
  iVarCalc = mapDepToCalc(iVar)
  IF (iVarCalc.GT.0 .AND. maskCalc(iVar).GT.0) THEN
    SWRITE(*,*) "  ",TRIM(DepNames(iVar))
    IF(withGradients)THEN
      IF(withVectors)THEN
        CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,nVal,mapCalcMeshToGlobalMesh,mapDepToCalc,UCalc,gradUx,gradUy,gradUz,&
            NormVec,TangVec1,TangVec2)
      ELSE
        CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,nVal,mapCalcMeshToGlobalMesh,mapDepToCalc,UCalc,gradUx,gradUy,gradUz)
      END IF
    ELSE
      IF(withVectors)THEN
        CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,nVal,mapCalcMeshToGlobalMesh,mapDepToCalc,UCalc,&
            NormVec=NormVec,TangVec1=TangVec1,TangVec2=TangVec2)
      ELSE
        CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,nVal,mapCalcMeshToGlobalMesh,mapDepToCalc,UCalc)
      END IF
    END IF
  END IF
END DO
END SUBROUTINE CalcQuantities


!==================================================================================================================================
!> Routine that calculates a single derived quantity called DepName and stores them in UCalc(iVarCalc,:). 
!> The routine either does the calculation in itself if the quantity is simple to calculate or calls helper functions
!> for the more complex ones.
!==================================================================================================================================
SUBROUTINE CalcDerivedQuantity(iVarCalc,DepName,nVarCalc,nVal,mapCalcMeshToGlobalMesh,mapDepToCalc,UCalc,gradUx,gradUy,gradUz, &
    NormVec,TangVec1,TangVec2)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Posti_Vars
USE MOD_EOS_Vars
USE MOD_StringTools     ,ONLY: LowCase,KEYVALUE
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)                                              :: iVarCalc
CHARACTER(LEN=255),INTENT(IN)                                   :: DepName
INTEGER,INTENT(IN)                                              :: nVarCalc
INTEGER,INTENT(IN)                                              :: nVal(:)
INTEGER,INTENT(IN)                                              :: mapCalcMeshToGlobalMesh(:)
INTEGER,INTENT(IN)                                              :: mapDepToCalc(nVarDepEOS)
REAL,INTENT(INOUT)                                              :: UCalc(PRODUCT(nVal),1:nVarCalc)
REAL,DIMENSION(1:PP_nVarPrim,PRODUCT(nVal)),INTENT(IN),OPTIONAL :: gradUx,gradUy,gradUz
REAL,DIMENSION(1:3,PRODUCT(nVal)),INTENT(IN),OPTIONAL           :: NormVec,TangVec1,TangVec2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iMom1,iMom2,iMom3,iDens,iPres,iVel1,iVel2,iVel3,iVelM,iVelS,iEner,iEnst,iTemp
CHARACTER(LEN=255) :: DepName_low
REAL               :: UE(PP_2Var)
INTEGER            :: nElems_loc,Nloc,nDOF,nDims
#if PARABOLIC
INTEGER            :: iVorM,iWFriX,iWFriY,iWFriZ,iWFriMag
#endif
LOGICAL            :: withGradients
LOGICAL            :: withVectors
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
withVectors=(PRESENT(NormVec).AND.PRESENT(TangVec1).AND.PRESENT(TangVec2))

! Already retrieve mappings for conservatives, as theses are required for many quantities
iDens = KEYVALUE(DepNames,mapDepToCalc,"density"    )
iMom1 = KEYVALUE(DepNames,mapDepToCalc,'momentumx')
iMom2 = KEYVALUE(DepNames,mapDepToCalc,'momentumy')
iMom3 = KEYVALUE(DepNames,mapDepToCalc,'momentumz')
iVel1 = KEYVALUE(DepNames,mapDepToCalc,"velocityx"  )
iVel2 = KEYVALUE(DepNames,mapDepToCalc,"velocityy"  )
iVel3 = KEYVALUE(DepNames,mapDepToCalc,"velocityz"  )
iPres = KEYVALUE(DepNames,mapDepToCalc,"pressure"   )
iEner = KEYVALUE(DepNames,mapDepToCalc,'energystagnationdensity')
iVelM = KEYVALUE(DepNames,mapDepToCalc,"velocitymagnitude")
iTemp = KEYVALUE(DepNames,mapDepToCalc,"temperature")
iVelS = KEYVALUE(DepNames,mapDepToCalc,"velocitysound"    )
iEnst = KEYVALUE(DepNames,mapDepToCalc,"energystagnation")
#if PARABOLIC
iVorM = KEYVALUE(DepNames,mapDepToCalc,"vorticitymagnitude")
iWFriX = KEYVALUE(DepNames,mapDepToCalc,"wallfrictionx")
iWFriY = KEYVALUE(DepNames,mapDepToCalc,"wallfrictiony")
iWFriZ = KEYVALUE(DepNames,mapDepToCalc,"wallfrictionz")
iWFriMag = KEYVALUE(DepNames,mapDepToCalc,"wallfrictionmagnitude")
#endif

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
    DO i=1,PRODUCT(nVal)
      UE(SRHO) = 1./UCalc(i,iDens)
      UE(MOM1) =    UCalc(i,iMom1)
      UE(MOM2) =    UCalc(i,iMom2)
      UE(MOM3) =    UCalc(i,iMom3)
      UE(ENER) =    UCalc(i,iEner)
      UE(VELV) = UE(SRHO)*UE(MOMV)
      UCalc(i,iVarCalc) = PRESSURE_HE(UE)
    END DO
  CASE("temperature")
    DO i=1,PRODUCT(nVal)
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
     CALL FillPressureTimeDeriv(nElems_loc,mapCalcMeshToGlobalMesh,Nloc,UCalc(:,iVarCalc))
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
  CASE("normalizedhelicity")
    ! calculate normalized helicity:
    ! helicity = scalar product of velocity and vorticity
    ! normalized helicity = helicity divided by magnitude of velocity times magnitude of vorticity
    ! Will effectively return the cosine of the angle between velocity and vorticity and is thus
    ! normalized between -1 and 1.
    UCalc(:,iVarCalc) =  FillVorticity(1,nVal,gradUx,gradUy,gradUz)*UCalc(:,iMom1)/UCalc(:,iDens) &
                        +FillVorticity(2,nVal,gradUx,gradUy,gradUz)*UCalc(:,iMom2)/UCalc(:,iDens) &
                        +FillVorticity(3,nVal,gradUx,gradUy,gradUz)*UCalc(:,iMom3)/UCalc(:,iDens)
    ! Catch division by very small values of magnitude of velocity times magnitude of vorticity
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
IF (withVectors) THEN
  SELECT CASE(DepName_low)
#if PARABOLIC      
    CASE("wallfrictionx")
      UCalc(:,iVarCalc) = FillWallFriction(1,nVal,UCalc(:,iTemp),gradUx,gradUy,gradUz,NormVec)
    CASE("wallfrictiony")
      UCalc(:,iVarCalc) = FillWallFriction(2,nVal,UCalc(:,iTemp),gradUx,gradUy,gradUz,NormVec)
    CASE("wallfrictionz")
      UCalc(:,iVarCalc) = FillWallFriction(3,nVal,UCalc(:,iTemp),gradUx,gradUy,gradUz,NormVec)
    CASE("wallfrictionmagnitude")
      UCalc(:,iVarCalc) = SQRT(UCalc(:,iWFriX)**2+UCalc(:,iWFriY)**2+UCalc(:,iWFriZ)**2)
    CASE("wallheattransfer")
      UCalc(:,iVarCalc) = FillWallHeatTransfer(nVal,UCalc(:,iTemp),gradUx,gradUy,gradUz,NormVec)
    CASE("x+")
      CALL FillNonDimensionalGridSpacing(nElems_loc,mapCalcMeshToGlobalMesh,Nloc,1,UCalc(:,iTemp),UCalc(:,iWFriMag), &
                                     UCalc(:,iDens),UCalc(:,iVarCalc))
    CASE("y+")
      CALL FillNonDimensionalGridSpacing(nElems_loc,mapCalcMeshToGlobalMesh,Nloc,2,UCalc(:,iTemp),UCalc(:,iWFriMag), &
                                     UCalc(:,iDens),UCalc(:,iVarCalc))
    CASE("z+")
      CALL FillNonDimensionalGridSpacing(nElems_loc,mapCalcMeshToGlobalMesh,Nloc,3,UCalc(:,iTemp),UCalc(:,iWFriMag), &
                                     UCalc(:,iDens),UCalc(:,iVarCalc))
#endif
  END SELECT
END IF
END SUBROUTINE CalcDerivedQuantity


!==================================================================================================================================
!> Calculate the time derivative of the pressure from the time derivative of the conservative variables using the chain rule.
!==================================================================================================================================
SUBROUTINE FillPressureTimeDeriv(nElems_calc,indices,Nloc,PressureTDeriv)
! MODULES
USE MOD_Preproc
USE MOD_Eos_Vars,ONLY:KappaM1
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
IF(Nloc.NE.PP_N) STOP 'Not possible here'
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
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
!==================================================================================================================================
!> Calculate vorticity in direction dir.
!==================================================================================================================================
FUNCTION FillVorticity(dir,nVal,gradUx,gradUy,gradUz) RESULT(Vorticity)
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


!==================================================================================================================================
!> Calculate the lambda 2 criterion, see Jeong, Jinhee, and Fazle Hussain. "On the identification of a vortex." Journal of fluid
!> mechanics 285 (1995): 69-94. 
!> This criterion is the second eigenvalue of the symmetric tensor ${\bm {\cal S}}^2 + {\bm \Omega}^2$; 
!> here ${\bm {\cal S}}$ and ${\bm \Omega}$ are respectively the symmetric and antisymmetric parts of 
!> the velocity gradient tensor ${\bm \Delta}{\bm u}$.
!> Calculation of the eigenvalues is done using LAPACK.
!==================================================================================================================================
FUNCTION FillLambda2(nVal,gradUx,gradUy,gradUz) RESULT(Lambda2)
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
          +MATMUL(0.5*(gradUmat-TRANSPOSE(gradUmat)),0.5*(gradUmat-TRANSPOSE(gradUmat)))   ! Omega^2
  ! Jacobi-Subroutine used from LAPACK 
  CALL DSYEV('N', 'U', 3, gradUmat, 3, Lambda, WORK, 16, INFO )
  Lambda2(i) = Lambda(2)
END DO
END FUNCTION FillLambda2


!==================================================================================================================================
!> Calculate the Q criterion, which is the second invariant of the velocity gradient tensor.
!==================================================================================================================================
FUNCTION FillQcriterion(nVal,gradUx,gradUy,gradUz) RESULT(Qcriterion)
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

!==================================================================================================================================
!> Calculate the wall friction in direction dir.
!==================================================================================================================================
FUNCTION FillWallFriction(dir,nVal,Temperature,gradUx,gradUy,gradUz,NormVec) RESULT(WallFriction)
! MODULES
USE MOD_Eos_Vars
USE MOD_Viscosity
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                   :: dir,nVal(:)
REAL,DIMENSION(PRODUCT(nVal)),INTENT(IN)             :: Temperature
REAL,DIMENSION(PP_nVarPrim,PRODUCT(nVal)),INTENT(IN) :: gradUx,gradUy,gradUz
REAL,DIMENSION(1:3,PRODUCT(nVal)),INTENT(IN)         :: NormVec
REAL                                                 :: WallFriction(PRODUCT(nVal))
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL              :: tau(3,3)                  ! Viscous stress tensor
REAL              :: mu
REAL              :: GradV(3,3),DivV
REAL              :: WallFrictionLoc(3)
REAL              :: temp
INTEGER           :: i
!===================================================================================================================================
DO i=1,PRODUCT(nVal)
  temp=Temperature(i)
  mu=VISCOSITY_TEMPERATURE(temp)
  GradV(1:3,1)=gradUx(2:4,i)
  GradV(1:3,2)=gradUy(2:4,i)
  GradV(1:3,3)=gradUz(2:4,i)
  ! Velocity divergence
  DivV=GradV(1,1)+GradV(2,2)+GradV(3,3)
  ! Calculate shear stress tensor
  tau=mu*(GradV + TRANSPOSE(GradV))
  tau(1,1)=tau(1,1)-2./3.*mu*DivV
  tau(2,2)=tau(2,2)-2./3.*mu*DivV
  tau(3,3)=tau(3,3)-2./3.*mu*DivV
  ! Calculate viscous force vector
  WallFrictionLoc(1:3)=-1*MATMUL(tau,NormVec(:,i))
  WallFriction(i)=WallFrictionLoc(dir)
END DO
END FUNCTION FillWallFriction

!==================================================================================================================================
!> Calculate the wall heat transfer normal to the wall.
!==================================================================================================================================
FUNCTION FillWallHeatTransfer(nVal,Temperature,gradUx,gradUy,gradUz,NormVec) RESULT(WallHeatTransfer)
! MODULES
USE MOD_Eos_Vars
USE MOD_Viscosity
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                   :: nVal(:)
REAL,DIMENSION(PRODUCT(nVal)),INTENT(IN)             :: Temperature
REAL,DIMENSION(PP_nVarPrim,PRODUCT(nVal)),INTENT(IN) :: gradUx,gradUy,gradUz
REAL,DIMENSION(1:3,PRODUCT(nVal)),INTENT(IN)         :: NormVec
REAL                                                 :: WallHeatTransfer(PRODUCT(nVal))
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL              :: mu
REAL              :: GradTn
REAL              :: temp
INTEGER           :: i
!===================================================================================================================================
DO i=1,PRODUCT(nVal)
  ! Calculate the viscosity
  temp=Temperature(i)
  mu=VISCOSITY_TEMPERATURE(temp)
  ! Calculate temperature gradient in wall normal direction
  GradTn = gradUx(6,i)*NormVec(1,i) &
         + gradUy(6,i)*NormVec(2,i) &
         + gradUz(6,i)*NormVec(3,i)
  ! Calculate wall heat transfer
  WallHeatTransfer(i) = -1.*mu*Kappa*sKappaM1*R/Pr*gradTn
END DO
END FUNCTION FillWallHeatTransfer

!==================================================================================================================================
!> Calculate the non dimensional grid spacing. This is a special case since we need information about the mesh. Thus 
!> a lot of actual FLEXI vars and routines are used. Care has to be taken that those are all available!
!> The non-dimensional grid spacing is the length of the element expressed in wall unit length, which is defined as 
!> y+ = y *frictionVel / kinematicVisc where frictionVel = sqrt(WallFriction/rho). Additionaly, the length is normalized with
!> the number of solution points in each direction per cell.
!==================================================================================================================================
SUBROUTINE FillNonDimensionalGridSpacing(nSides_calc,mapBCSideToVisuSides,Nloc,dir,Temperature,WallFrictionMag,Density,wallDistance)
! MODULES
USE MOD_Eos_Vars
USE MOD_Preproc
USE MOD_Viscosity
USE MOD_Mesh_Vars           ,ONLY: NGeo,Elem_xGP,SideToElem,nBCSides
USE MOD_Interpolation       ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars  ,ONLY: NodeType,NodeTypeGL
USE MOD_ChangeBasis         ,ONLY: ChangeBasis3D
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                   :: nSides_calc
INTEGER,INTENT(IN)                                   :: mapBCSideToVisuSides(nBCSides)
INTEGER,INTENT(IN)                                   :: Nloc
INTEGER,INTENT(IN)                                   :: dir
REAL,DIMENSION(0:Nloc,0:Nloc,nSides_calc),INTENT(IN) :: WallFrictionMag
REAL,DIMENSION(0:Nloc,0:Nloc,nSides_calc),INTENT(IN) :: Density
REAL,DIMENSION(0:Nloc,0:Nloc,nSides_calc),INTENT(IN) :: Temperature
REAL,INTENT(OUT)                                     :: wallDistance(0:Nloc,0:Nloc,nSides_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL              :: mu,temp,fricVel
INTEGER           :: iSide,p,q,iSideVisu,i,ElemID,locSideID
REAL              :: refVec(3),scalProd,scalProdMax,xVec(3),yVec(3),zVec(3),yloc(3),tVec(3,2)
REAL              :: NodeCoords(3,0:PP_N,0:PP_N,0:PP_N)
REAL              :: Vdm_N_GLNloc(0:PP_N,0:PP_N)
!===================================================================================================================================
! We need a Vandermonde matrix to get the GP coordinates on GL points, i.e. including the surface of the grid cells
CALL GetVandermonde(PP_N,NodeType,PP_N,NodeTypeGL,Vdm_N_GLNloc)

! Loop over all boundary sides
DO iSide=1,nBCSides
  iSideVisu = mapBCSideToVisuSides(iSide) ! Index on visualization side
  IF (iSideVisu.GT.0) THEN ! Check if this side should be visualized
    ElemID        = SideToElem(S2E_ELEM_ID,iSide)
    locSideID     = SideToElem(S2E_LOC_SIDE_ID,iSide)

    ! Get element coordinates on GL points (they include the edges that we use to calculate the element length)
    CALL ChangeBasis3D(3,PP_N,Nloc,Vdm_N_GLNloc,Elem_xGP(:,:,:,:,ElemID),NodeCoords)

    ! Depending in the local sideID, get the edge vectors of the element in wall-normal direction (stored in yVec) and for the two
    ! wall-tangential directions (stored in tVec(:,1-2)). These vectors are the connection of the cell vertices, so they do not
    ! include any information about curvature!!!!
    SELECT CASE(locSideID)
    CASE(XI_MINUS,XI_PLUS)
     yVec(:)   = NodeCoords(:,NGeo,0,0)- NodeCoords(:,0,0,0)
     tVec(:,1) = NodeCoords(:,0,NGeo,0)- NodeCoords(:,0,0,0)
     tVec(:,2) = NodeCoords(:,0,0,NGeo)- NodeCoords(:,0,0,0)
    CASE(ETA_MINUS,ETA_PLUS)
     yVec(:)   = NodeCoords(:,0,NGeo,0)- NodeCoords(:,0,0,0)
     tVec(:,1) = NodeCoords(:,NGeo,0,0)- NodeCoords(:,0,0,0)
     tVec(:,2) = NodeCoords(:,0,0,NGeo)- NodeCoords(:,0,0,0)
    CASE(ZETA_MINUS,ZETA_PLUS)
     yVec(:)   = NodeCoords(:,0,0,NGeo)- NodeCoords(:,0,0,0)
     tVec(:,1) = NodeCoords(:,NGeo,0,0)- NodeCoords(:,0,0,0)
     tVec(:,2) = NodeCoords(:,0,NGeo,0)- NodeCoords(:,0,0,0)
    END SELECT
     
    ! For the two tangential vectors, we try to find out which one is pointing in the physical x-direction
    refVec=(/1.,0.,0./) ! Vector in x-direction
    scalProdMax=0.
    ! Loop over both tangential vectors
    DO i=1,2
     ! Calculate absolute value of scalar product between vector in x-direction and the tangential vectors. This gives us 
     ! an information about how much the current tangential vector is alinged with the x-direction.
     scalProd=ABS(SUM(tVec(:,i)*refVec))
     ! Check if this tangential vector is more aligned with the x-direction than the other one. Store the one
     ! that is most aligned in xVec.
     IF(scalProd .GT. scalProdMax) THEN
       scalProdMax=scalProd
       xVec=tVec(:,i)
       zVec=tVec(:,3-i)
     END IF
    END DO

    ! Calculate lengths of the cell = lengths of the edge vectors
    yloc(1)=NORM2(xVec)
    yloc(2)=NORM2(yVec)
    yloc(3)=NORM2(zVec)

    ! Normalize the distances with the number of points per direction per cell
    yloc=yloc/(PP_N+1)

    ! Non-dimensional wall distance is calculated using  y+ = y * frictionVel / kinematicVisc
    ! where frictionVel = sqrt(WallFriction/rho)
    DO q=0,Nloc; DO p=0,Nloc
      fricVel=SQRT(WallFrictionMag(p,q,iSideVisu))
      temp = Temperature(p,q,iSideVisu)
      mu   = VISCOSITY_TEMPERATURE(temp)
      wallDistance(p,q,iSideVisu) = yloc(dir)*fricVel*Density(p,q,iSideVisu)/mu
    END DO; END DO ! p,q=0,Nloc

  END IF ! iSideVisu.GT.0
END DO ! iSide = 1,nBCSides

END SUBROUTINE FillNonDimensionalGridSpacing
#endif


END MODULE MOD_EOS_Posti
