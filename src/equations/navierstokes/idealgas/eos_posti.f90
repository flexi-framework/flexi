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
!> Contains all the routines to calculate the (equation system and EOS dependant) conservative/primitive/derived quantities. 
!> Dependandy table will be filled in here.
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

INTERFACE InitDepTable
  MODULE PROCEDURE InitDepTable
END INTERFACE 

INTERFACE FillDepTable
  MODULE PROCEDURE FillDepTable
END INTERFACE

#if FV_ENABLED && FV_RECONSTRUCT
INTERFACE AppendNeededPrims
  MODULE PROCEDURE AppendNeededPrims
END INTERFACE
PUBLIC :: AppendNeededPrims
#endif

INTERFACE GetVarnames
  MODULE PROCEDURE GetVarnames
END INTERFACE

INTERFACE GetMaskCons
  MODULE PROCEDURE GetMaskCons
END INTERFACE

INTERFACE GetMaskPrim
  MODULE PROCEDURE GetMaskPrim
END INTERFACE

INTERFACE GetMaskGrad
  MODULE PROCEDURE GetMaskGrad
END INTERFACE

INTERFACE FillCopy
  MODULE PROCEDURE FillCopy
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

#if FV_ENABLED && FV_RECONSTRUCT
INTERFACE CalcConsFromPrim
  MODULE PROCEDURE CalcConsFromPrim
END INTERFACE
PUBLIC:: CalcConsFromPrim
#endif


PUBLIC :: InitDepTable
PUBLIC :: FillDepTable
PUBLIC :: GetVarnames
PUBLIC :: GetMaskCons
PUBLIC :: GetMaskPrim
PUBLIC :: GetMaskGrad
PUBLIC :: CalcQuantities
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Build a table containing the information on which other variabels all variables depent.
!==================================================================================================================================
SUBROUTINE InitDepTable() 
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
!================================================================================================================================
                                        

! ATTENTION: The first     5 variables must be the conservative ones
!            The following 5 variables must be the primitive ones
!                          E            
!                          n            
!                          e                 
!                          r                 
!                          g                 
!                          y                    E              V
!                          S            V       n              o
!                          t            e     E t   T          r
!                          a            l     n h   o          t
!                          g            o     e a   t          i
!                          n            c V   r l   a T        c
!                          a            i e   g p   l o        i
!                          t         T  t l   y y   T t        t
!                          i         e  y o   S S   e a  V V V y     D Q
!                G   M M M o V V V   m  M c   t t   m l  o o o M     i C S
!                r   o o o n e e e P p  a i   a a   p P  r r r a H   l r c 
!                a D m m m D l l l r e  g t   g g E e r  t t t g e L a i h  
!                d e e e e e o o o e r  n y   n n n r e  i i i n l a t t l  
!                i n n n n n c c c s a  i S   a a t a s  c c c i i m a e i  
!                e s t t t s i i i s t  t o M t t r t s  i i i t c b t r e  
!                n i u u u i t t t u u  u u a i i o u u  t t t u i d i i r  
!                t t m m m t y y y r r  d n c o o p r r  y y y d t a o o e  
!                s y X Y Z y X Y Z e e  e d h n n y e e  X Y Z e y 2 n n n  
DepTable(1 ,:)=(/0,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Density
DepTable(2 ,:)=(/0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! MomentumX
DepTable(3 ,:)=(/0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! MomentumY
DepTable(4 ,:)=(/0,0,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! MomentumZ
DepTable(5 ,:)=(/0,0,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! EnergyStagnationDensity
DepTable(6 ,:)=(/0,1,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VelocityX
DepTable(7 ,:)=(/0,1,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VelocityY
DepTable(8 ,:)=(/0,1,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VelocityZ
DepTable(9 ,:)=(/0,1,1,1,1,1,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Pressure
DepTable(10,:)=(/0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Temperature
                                                         
DepTable(11,:)=(/0,0,0,0,0,0,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VelocityMagnitude 
DepTable(12,:)=(/0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VelocitySound     
DepTable(13,:)=(/0,0,0,0,0,0,0,0,0,0,0, 1,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Mach              
DepTable(14,:)=(/0,1,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! EnergyStagnation  
DepTable(15,:)=(/0,1,0,0,0,0,0,0,0,1,0, 0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! EnthalpyStagnation
DepTable(16,:)=(/0,1,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Entropy           
DepTable(17,:)=(/0,0,0,0,0,0,0,0,0,0,1, 1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! TotalTemperature  
DepTable(18,:)=(/0,1,0,0,0,0,0,0,0,1,0, 1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! TotalPressure     
#if PARABOLIC                                                         
DepTable(19,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VorticityX
DepTable(20,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VorticityY
DepTable(21,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! VorticityZ
DepTable(22,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 1,1,1,0,0,0,0,0,0 /) ! VorticityMagnitude
DepTable(23,:)=(/1,0,0,0,0,0,1,1,1,0,0, 1,0,0,0,0,0,0,0, 1,1,1,1,0,0,0,0,0 /) ! Helicity
DepTable(24,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Lambda2
DepTable(25,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Dilatation
DepTable(26,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! QCriterion
DepTable(27,:)=(/1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 /) ! Schlieren
#endif

CALL FillDepNames()

#if FV_ENABLED && FV_RECONSTRUCT
!                                    E            
!                                    n            
!                                    e                 
!                                    r                 
!                                    g                 
!                                    y          
!                                    S          
!                                    t          
!                                    a          
!                                    g          
!                                    n          
!                                    a          
!                                    t         T
!                                    i         e
!                          G   M M M o V V V   m
!                          r   o o o n e e e P p 
!                          a D m m m D l l l r e  
!                          d e e e e e o o o e r  
!                          i n n n n n c c c s a  
!                          e s t t t s i i i s t  
!                          n i u u u i t t t u u  
!                          t t m m m t y y y r r  
!                          s y X Y Z y X Y Z e e  
DepTablePrimToCons(1 ,:)=(/0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0/) ! Density
DepTablePrimToCons(2 ,:)=(/0,1,0,0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0/) ! MomentumX
DepTablePrimToCons(3 ,:)=(/0,1,0,0,0,0,0,1,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0/) ! MomentumY
DepTablePrimToCons(4 ,:)=(/0,1,0,0,0,0,0,0,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0/) ! MomentumZ
DepTablePrimToCons(5 ,:)=(/0,1,0,0,0,0,1,1,1,1,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0/) ! EnergyStagnationDensity
#endif
END SUBROUTINE InitDepTable  

!==================================================================================================================================
!> Build an array containing all available variable names.
!==================================================================================================================================
SUBROUTINE FillDepNames() 
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
!===================================================================================================================================
DepNames(1 ) = "Density"
DepNames(2 ) = "MomentumX"
DepNames(3 ) = "MomentumY"
DepNames(4 ) = "MomentumZ"
DepNames(5 ) = "EnergyStagnationDensity"
DepNames(6 ) = "VelocityX"
DepNames(7 ) = "VelocityY"
DepNames(8 ) = "VelocityZ"
DepNames(9 ) = "Pressure"
DepNames(10) = "Temperature"
DepNames(11) = "VelocityMagnitude" 
DepNames(12) = "VelocitySound"     
DepNames(13) = "Mach"              
DepNames(14) = "EnergyStagnation"  
DepNames(15) = "EnthalpyStagnation"
DepNames(16) = "Entropy"           
DepNames(17) = "TotalTemperature"  
DepNames(18) = "TotalPressure"     
#if PARABOLIC
DepNames(19) = "VorticityX"
DepNames(20) = "VorticityY"
DepNames(21) = "VorticityZ"
DepNames(22) = "VorticityMagnitude"
DepNames(23) = "Helicity"
DepNames(24) = "Lambda2"
DepNames(25) = "Dilatation"
DepNames(26) = "QCriterion"
DepNames(27) = "Schlieren"
#endif
END SUBROUTINE FillDepNames

!==================================================================================================================================
!> Build the actual dependancies. If we call the DGTimeDerivative_weakForm, the primitive variables will be calculated
!> in there and don't need to be computed from the conservative variables.
!> Also, recursively add the dependancies of variabels that a variable depends on to the dependancy of the first variable.
!==================================================================================================================================
SUBROUTINE FillDepTable(withGradients) 
USE MOD_EOS_Posti_Vars
USE MOD_Equation_Vars ,ONLY: StrVarNames,StrVarNamesPrim
USE MOD_StringTools   ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: withGradients
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVar,ivar2,iVarPrim,nVar
!===================================================================================================================================
! if withGradients==True then the DGTimeDerivative_weakForm routine is called and the primitive
! quantities are computed in this operator, which can be used for calculation of other quantities
! directly. Therefore the dependecies of the primitive variables on the conservative variables 
! can be removed.
IF (withGradients) THEN ! remove conservative from all prim variables
  nVar = nVarTotal
  DO iVar=1,nVar ; DO iVarPrim=1,PP_nVarPrim ! search prim
    IF (STRICMP(StrVarNamesPrim(iVarPrim),DepNames(iVar))) THEN
      DepTable(iVar,:)    = 0
      DepTable(iVar,iVar) = 1
    END IF
  END DO ; END DO
END IF

! for each quantity copy from all quantities that this quantity depends on the dependencies.
DO iVar=1,nVarTotal
  DepTable(iVar,iVar) = 1
  DO iVar2=1,iVar-1
    IF (DepTable(iVar,iVar2).EQ.1) THEN
      DepTable(iVar,:) = MAX(DepTable(iVar,:), DepTable(iVar2,:))
    END IF
  END DO
END DO
END SUBROUTINE FillDepTable

#if FV_ENABLED && FV_RECONSTRUCT
SUBROUTINE AppendNeededPrims(mapCalc,mapCalc_FV,nVarCalc) 
USE MOD_EOS_Posti_Vars
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)      :: mapCalc(nVarTotal)
INTEGER,INTENT(OUT)     :: mapCalc_FV(nVarTotal)
INTEGER,INTENT(OUT)     :: nVarCalc
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar
!===================================================================================================================================
mapCalc_FV = mapCalc
DO iVar=1,PP_nVar
  IF (mapCalc(iVar).GT.0) THEN
    mapCalc_FV = MAX(mapCalc_FV, DepTablePrimToCons(iVar,1:nVarTotal))
  END IF
END DO

! renumerate mapCalc_FV
nVarCalc   = 0
DO iVar=1,nVarTotal
  IF (mapCalc_FV(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapCalc_FV(iVar) = nVarCalc
  END IF
END DO
END SUBROUTINE AppendNeededPrims
#endif

SUBROUTINE GetVarnames(varnames,nVarAdd) 
USE MOD_EOS_Posti_Vars
! MODULES                                                                                                                          !
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),POINTER,INTENT(INOUT) :: varnames(:)
INTEGER,INTENT(IN)                       :: nVarAdd
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar
!===================================================================================================================================
CALL FillDepNames()
ALLOCATE(varnames(1:nVarTotal+nVarAdd))
DO iVar=1,nVarTotal
  varnames(iVar) = DepNames(iVar) 
END DO
END SUBROUTINE GetVarnames

FUNCTION GetMaskCons() 
USE MOD_EOS_Posti_Vars
USE MOD_Equation_Vars ,ONLY: StrVarNames
USE MOD_StringTools   ,ONLY: STRICMP
! INPUT / OUTPUT VARIABLES 
INTEGER :: GetMaskCons(nVarTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar,iVar2
!===================================================================================================================================
GetMaskCons = 0
DO iVar=1,nVarTotal
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
INTEGER :: GetMaskPrim(nVarTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar,iVar2
!===================================================================================================================================
GetMaskPrim = 0
DO iVar=1,nVarTotal
  DO iVar2=1,PP_nVarPrim
    IF (STRICMP(StrVarNamesPrim(iVar2),DepNames(iVar))) THEN
      GetMaskPrim(iVar) = 1
    END IF
  END DO
END DO
END FUNCTION GetMaskPrim

FUNCTION GetMaskGrad() 
USE MOD_EOS_Posti_Vars
INTEGER :: GetMaskGrad(nVarTotal)
!===================================================================================================================================
GetMaskGrad = DepTable(:,0)
END FUNCTION GetMaskGrad


PURE SUBROUTINE FillCopy(nVar,Nloc,nElems,U,nElems_calc,indices,Cons,iCons)
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nElems_calc
INTEGER,INTENT(IN) :: indices(nElems_calc)
INTEGER,INTENT(IN) :: Nloc
INTEGER,INTENT(IN) :: nVar
INTEGER,INTENT(IN) :: nElems
INTEGER,INTENT(IN) :: iCons
REAL,INTENT(IN)    :: U(nVar,0:Nloc,0:Nloc,0:Nloc,nElems)
REAL,INTENT(OUT)   :: Cons( 0:Nloc,0:Nloc,0:Nloc,nElems_calc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER         :: iElem,iElem_calc
!==================================================================================================================================
DO iElem_calc=1,nElems_calc
  iElem = indices(iElem_calc)
  Cons(:,:,:,iElem_calc) = U(iCons,:,:,:,iElem)
END DO ! iElem
END SUBROUTINE FillCopy

SUBROUTINE FillPressure(nElems_calc,indices,Nloc,nElems,U,Pressure)
!==================================================================================================================================
! MODULES
USE MOD_Eos_Vars  ,ONLY: KappaM1
! IMPLICIT VARIABLE HANDLING
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

SUBROUTINE CalcQuantities(nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc,withGradients,maskCalc) 
! MODULES
USE MOD_Globals
USE MOD_PreProc         ,ONLY: PP_N
USE MOD_EOS_Posti_Vars  ,ONLY: nVarTotal,DepNames
USE MOD_DG_Vars         ,ONLY: U,UPrim
USE MOD_Mesh_Vars       ,ONLY: nElems
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN) :: nVarCalc
INTEGER,INTENT(IN) :: Nloc
INTEGER,INTENT(IN) :: nElems_loc
INTEGER,INTENT(IN) :: iElems(nElems_loc)
INTEGER,INTENT(IN) :: mapCalc(nVarTotal)
REAL,INTENT(OUT)   :: UCalc(0:Nloc,0:Nloc,0:Nloc,1:nElems_loc,1:nVarCalc)
LOGICAL,INTENT(IN) :: withGradients
INTEGER,INTENT(IN),OPTIONAL :: maskCalc(nVarTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVar,iVarCalc
INTEGER            :: maskCons(nVarTotal)
INTEGER            :: maskPrim(nVarTotal)
!===================================================================================================================================
maskCons = GetMaskCons()
maskPrim = GetMaskPrim()
DO iVar=1,nVarTotal
  iVarCalc = mapCalc(iVar)
  iVarCalc = MERGE(maskCalc(iVar)*iVarCalc, iVarCalc, PRESENT(maskCalc))
  IF (iVarCalc.GT.0) THEN
    SWRITE(*,*) "  ",TRIM(DepNames(iVar))
    IF (maskCons(iVar).GT.0) THEN
      ! ATTENTION: The first 5 variables in DepTable must be the conservative ones
      CALL FillCopy(PP_nVar,PP_N,nElems,U,nElems_loc,iElems,UCalc(:,:,:,:,iVarCalc),iVar)
    ELSE IF(maskPrim(iVar).GT.0) THEN
      IF (withGradients) THEN
        ! ATTENTION: The following 5 variables in DepTable must be the primitive ones
        CALL FillCopy(PP_nVarPrim,PP_N,nElems,UPrim,nElems_loc,iElems,UCalc(:,:,:,:,iVarCalc),iVar-5+1)
      ELSE
        CALL CalcPrimitiveQuantity(iVarCalc,DepNames(iVar),nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc)
      END IF
    ELSE
      CALL CalcDerivedQuantity(iVarCalc,DepNames(iVar),nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc)
    END IF
  END IF
END DO
END SUBROUTINE CalcQuantities

SUBROUTINE CalcPrimitiveQuantity(iVarCalc,DepName,nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc)
! MODULES
USE MOD_PreProc         ,ONLY: PP_N
USE MOD_EOS_Posti_Vars  ,ONLY: nVarTotal
USE MOD_StringTools     ,ONLY: LowCase
USE MOD_DG_Vars         ,ONLY: U
USE MOD_Mesh_Vars       ,ONLY: nElems
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)            :: iVarCalc
CHARACTER(LEN=255),INTENT(IN) :: DepName
INTEGER,INTENT(IN)            :: nVarCalc
INTEGER,INTENT(IN)            :: Nloc
INTEGER,INTENT(IN)            :: nElems_loc
INTEGER,INTENT(IN)            :: iElems(nElems_loc)
INTEGER,INTENT(IN)            :: mapCalc(nVarTotal)
REAL,INTENT(OUT)              :: UCalc(0:Nloc,0:Nloc,0:Nloc,1:nElems_loc,1:nVarCalc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iDens,iMom,iPres
CHARACTER(LEN=255) :: DepName_low
!===================================================================================================================================
CALL LowCase(DepName,DepName_low)
SELECT CASE(DepName_low)
CASE("velocityx")
  iDens = GETiCalc('density',mapCalc)
  iMom  = GETiCalc('momentumx',mapCalc)
  UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iMom) / UCalc(:,:,:,:,iDens)
CASE("velocityy")
  iDens = GETiCalc('density',mapCalc)
  iMom  = GETiCalc('momentumy',mapCalc)
  UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iMom) / UCalc(:,:,:,:,iDens)
CASE("velocityz")
  iDens = GETiCalc('density',mapCalc)
  iMom  = GETiCalc('momentumz',mapCalc)
  UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iMom) / UCalc(:,:,:,:,iDens)
CASE("pressure")
  CALL FillPressure(nElems_loc,iElems,PP_N,nElems,U,UCalc(:,:,:,:,iVarCalc))
CASE("temperature")
  iDens = GETiCalc("density",mapCalc)
  iPres = GETiCalc("pressure",mapCalc)
  CALL FillTemperature(nElems_loc,PP_N,UCalc(:,:,:,:,iDens),UCalc(:,:,:,:,iPres),UCalc(:,:,:,:,iVarCalc))
END SELECT
END SUBROUTINE CalcPrimitiveQuantity

SUBROUTINE CalcDerivedQuantity(iVarCalc,DepName,nVarCalc,Nloc,nElems_loc,iElems,mapCalc,UCalc)
! MODULES
#if PARABOLIC
USE MOD_PreProc         ,ONLY: PP_N
#endif
USE MOD_EOS_Vars        ,ONLY: cp,kappa,R,sKappaM1
USE MOD_EOS_Posti_Vars  ,ONLY: nVarTotal
USE MOD_StringTools     ,ONLY: LowCase
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)            :: iVarCalc
CHARACTER(LEN=255),INTENT(IN) :: DepName
INTEGER,INTENT(IN)            :: nVarCalc
INTEGER,INTENT(IN)            :: Nloc
INTEGER,INTENT(IN)            :: nElems_loc
INTEGER,INTENT(IN)            :: iElems(nElems_loc)
INTEGER,INTENT(IN)            :: mapCalc(nVarTotal)
REAL,INTENT(OUT)              :: UCalc(0:Nloc,0:Nloc,0:Nloc,1:nElems_loc,1:nVarCalc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iDens,iPres,iVel1,iVel2,iVel3,iVelM,iVelS,iEner,iTemp
CHARACTER(LEN=255) :: DepName_low
#if PARABOLIC
INTEGER            :: iVor1,iVor2,iVor3,iVorM
INTEGER            :: i,j,k,iElem
REAL               :: Denominator
#endif
!===================================================================================================================================
CALL LowCase(DepName,DepName_low)
SELECT CASE(DepName_low)
  CASE("velocitymagnitude")
    iVel1 = GETiCalc("velocityx",mapCalc)
    iVel2 = GETiCalc("velocityy",mapCalc)
    iVel3 = GETiCalc("velocityz",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = SQRT(UCalc(:,:,:,:,iVel1)**2 + UCalc(:,:,:,:,iVel2)**2 + UCalc(:,:,:,:,iVel3)**2)
  CASE("velocitysound")
    iDens = GETiCalc("density",mapCalc)
    iPres = GETiCalc("pressure",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = SQRT(Kappa*UCalc(:,:,:,:,iPres)/UCalc(:,:,:,:,iDens))
  CASE("mach")
    iVelM = GETiCalc("velocitymagnitude",mapCalc)
    iVelS = GETiCalc("velocitysound",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iVelM)/UCalc(:,:,:,:,iVelS)
  CASE("energystagnation")
    iDens = GETiCalc("density",mapCalc)
    iEner = GETiCalc("energystagnationdensity",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iEner)/UCalc(:,:,:,:,iDens)
  CASE("enthalpystagnation")
    iDens = GETiCalc("density",mapCalc)
    iPres = GETiCalc("pressure",mapCalc)
    iEner = GETiCalc("energystagnation",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iEner) + UCalc(:,:,:,:,iPres)/UCalc(:,:,:,:,iDens)
  CASE("entropy")
    iDens = GETiCalc("density",mapCalc)
    iTemp = GETiCalc("temperature",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = R*(sKappaM1*LOG(UCalc(:,:,:,:,iTemp))) - LOG(UCalc(:,:,:,:,iDens))
  CASE("totaltemperature")
    iTemp = GETiCalc("temperature",mapCalc)
    iVelM = GETiCalc("velocitymagnitude",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iTemp)+UCalc(:,:,:,:,iVelM)**2/(2*cp)
  CASE("totalpressure")
    iDens = GETiCalc("density",mapCalc)
    iPres = GETiCalc("pressure",mapCalc)
    iVelM = GETiCalc("velocitymagnitude",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iPres)+0.5*UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVelM)**2
#if PARABOLIC      
  CASE("vorticityx")
    CALL FillVorticity(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc),1)
  CASE("vorticityy")
    CALL FillVorticity(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc),2)
  CASE("vorticityz")
    CALL FillVorticity(nElems_loc,iElems,PP_N,UCalc(:,:,:,:,iVarCalc),3)
  CASE("vorticitymagnitude")
    iVor1 = GETiCalc("vorticityx",mapCalc)
    iVor2 = GETiCalc("vorticityy",mapCalc)
    iVor3 = GETiCalc("vorticityy",mapCalc)
    UCalc(:,:,:,:,iVarCalc) = SQRT(UCalc(:,:,:,:,iVor1)**2 + UCalc(:,:,:,:,iVor2)**2 + UCalc(:,:,:,:,iVor3)**2)
  CASE("helicity")
    iVel1 = GETiCalc("velocityx",mapCalc)
    iVel2 = GETiCalc("velocityy",mapCalc)
    iVel3 = GETiCalc("velocityz",mapCalc)
    iVor1 = GETiCalc("vorticityx",mapCalc)
    iVor2 = GETiCalc("vorticityy",mapCalc)
    iVor3 = GETiCalc("vorticityy",mapCalc)
    iVelM = GETiCalc("velocitymagnitude",mapCalc)
    iVorM = GETiCalc("vorticitymagnitude",mapCalc)
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



#if FV_ENABLED && FV_RECONSTRUCT
SUBROUTINE CalcConsFromPrim(mapCalc,nVarCalc,Nloc,nElems_loc,UCalc) 
USE MOD_Globals    
USE MOD_PreProc     
USE MOD_EOS_Posti_Vars ,ONLY: nVarTotal,DepNames
USE MOD_EOS_Vars       ,ONLY: sKappaM1
USE MOD_StringTools    ,ONLY: LowCase
USE MOD_DG_Vars        ,ONLY: U,UPrim
USE MOD_Posti_Vars     ,ONLY: nElems_FV,mapElems_FV
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)          :: mapCalc(nVarTotal)
INTEGER,INTENT(IN)          :: nVarCalc
INTEGER,INTENT(IN)          :: Nloc
INTEGER,INTENT(IN)          :: nElems_loc
REAL,INTENT(INOUT),OPTIONAL :: UCalc(0:Nloc,0:Nloc,0:Nloc,1:nElems_loc,1:nVarCalc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iDens, iVel1, iVel2, iVel3, iPres, iTemp
INTEGER            :: iVar,iVarCalc,iElem,iElem_FV
CHARACTER(LEN=255) :: DepName_low
INTEGER            :: mapCons(nVarTotal)
INTEGER            :: maskCons(nVarTotal)
!===================================================================================================================================
! check wether any cons-quantity is needed
maskCons = GetMaskCons()
mapCons = maskCons * mapCalc
IF (SUM(mapCons).GT.0) THEN
  IF (PRESENT(UCalc)) THEN! fill conservative in UCalc from UCalc
    ! calculate all needed conservative variables
    iDens = GETiCalc("density",mapCalc)
    iVel1 = GETiCalc("velocityx",mapCalc)
    iVel2 = GETiCalc("velocityy",mapCalc)
    iVel3 = GETiCalc("velocityz",mapCalc)
    iPres = GETiCalc("pressure",mapCalc)
    iTemp = GETiCalc("temperature",mapCalc)

    DO iVar=1,nVarTotal
      IF (mapCons(iVar).GT.0) THEN
        iVarCalc = mapCalc(iVar)
        SWRITE(*,*) "  ", TRIM(DepNames(iVar)), " (NVisu_FV)"
        CALL LowCase(DepNames(iVar),DepName_low)
        SELECT CASE(DepName_low)
        CASE("momentumx")
          UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVel1)
        CASE("momentumy")
          UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVel2)
        CASE("momentumz")
          UCalc(:,:,:,:,iVarCalc) = UCalc(:,:,:,:,iDens)*UCalc(:,:,:,:,iVel3)
        CASE("energystagnationdensity")
          ! TODO: use a function from eos.f90
          UCalc(:,:,:,:,iVarCalc) = sKappaM1*UCalc(:,:,:,:,iPres) + 0.5*UCalc(:,:,:,:,iDens)* &
              (UCalc(:,:,:,:,iVel1)**2 + UCalc(:,:,:,:,iVel2)**2 + UCalc(:,:,:,:,iVel3)**2)
        END SELECT
      END IF
    END DO
  ELSE ! fill U from UPrim 
    DO iVar=1,nVarTotal
      IF (mapCons(iVar).GT.0) THEN
        SWRITE(*,*) "  ", TRIM(DepNames(iVar)), " (PP_N)"
        CALL LowCase(DepNames(iVar),DepName_low)
        SELECT CASE(DepName_low)
        CASE("density")
          DO iElem_FV=1,nElems_FV
            iElem = mapElems_FV(iElem_FV)
            U(1,:,:,:,iElem) = UPrim(1,:,:,:,iElem)
          END DO
        CASE("momentumx")
          DO iElem_FV=1,nElems_FV
            iElem = mapElems_FV(iElem_FV)
            U(2,:,:,:,iElem) = UPrim(1,:,:,:,iElem)*UPrim(2,:,:,:,iElem)
          END DO
        CASE("momentumy")
          DO iElem_FV=1,nElems_FV
            iElem = mapElems_FV(iElem_FV)
            U(3,:,:,:,iElem) = UPrim(1,:,:,:,iElem)*UPrim(3,:,:,:,iElem)
          END DO
        CASE("momentumz")
          DO iElem_FV=1,nElems_FV
            iElem = mapElems_FV(iElem_FV)
            U(4,:,:,:,iElem) = UPrim(1,:,:,:,iElem)*UPrim(4,:,:,:,iElem)
          END DO
        CASE("energystagnationdensity")
          ! TODO: use a function from eos.f90
          DO iElem_FV=1,nElems_FV
            iElem = mapElems_FV(iElem_FV)
            U(5,:,:,:,iElem) = sKappaM1*UPrim(5,:,:,:,iElem) + 0.5*UPrim(1,:,:,:,iElem)* &
                (UPrim(2,:,:,:,iElem)**2 + UPrim(3,:,:,:,iElem)**2 + UPrim(4,:,:,:,iElem)**2)
          END DO
        END SELECT
      END IF
    END DO
  END IF
END IF
END SUBROUTINE CalcConsFromPrim
#endif

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


FUNCTION GETiCalc(varname,mapCalc) 
USE MOD_EOS_Posti_Vars ,ONLY: nVarTotal,DepNames
USE MOD_StringTools    ,ONLY: STRICMP
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: varname
INTEGER,INTENT(IN)          :: mapCalc(nVarTotal)
INTEGER                     :: GETiCalc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar
!===================================================================================================================================
DO iVar=1,nVarTotal
  IF (STRICMP(DepNames(iVar),varname)) THEN
    GETiCalc = mapCalc(iVar)
    RETURN
  END IF
END DO
END FUNCTION GETiCalc

END MODULE MOD_EOS_Posti
