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
!> Subroutines  that provide gas properties and conversion between primitive and conservative variables for a ideal gas.
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

INTERFACE InitEOS
  MODULE PROCEDURE InitEOS
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
  MODULE PROCEDURE ConsToPrim_Side
  MODULE PROCEDURE ConsToPrim_Volume
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
  MODULE PROCEDURE PrimToCons_Side
END INTERFACE

INTERFACE PRESSURE_RIEMANN
  MODULE PROCEDURE PRESSURE_RIEMANN
END INTERFACE

PUBLIC::InitEos
PUBLIC::ConsToPrim
PUBLIC::PrimToCons
PUBLIC::PRESSURE_RIEMANN
PUBLIC::DefineParametersEos
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Define parameters for the used eos
!==================================================================================================================================
SUBROUTINE DefineParametersEos()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Equation of State")
CALL prms%CreateLogicalOption('UseNonDimensionalEqn',"Set true to compute R and mu from bulk Mach Reynolds (nondimensional form.",&
                                                '.FALSE.')
CALL prms%CreateRealOption(     'BulkMach',     "Bulk Mach     (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'BulkReynolds', "Bulk Reynolds (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent", '1.4')
CALL prms%CreateRealOption(     'R',            "Specific gas constant", '287.058')
CALL prms%CreateRealOption(     'Pr',           "Prandtl number", '0.72')
CALL prms%CreateRealOption(     'mu0',          "Dynamic Viscosity", '0.')
CALL prms%CreateRealOption(     'Ts',           "Sutherland's law for variable viscosity: Ts", '110.4')
CALL prms%CreateRealOption(     'Tref',         "Sutherland's law for variable viscosity: Tref ", '280.0')
CALL prms%CreateRealOption(     'ExpoSuth',     "Sutherland's law for variable viscosity: Exponent", '1.5')

END SUBROUTINE DefineParametersEos

!==================================================================================================================================
!> Initialize variables needed by the ideal gas equation of state.
!==================================================================================================================================
SUBROUTINE InitEos()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1,KappaP1,cp,cv
USE MOD_EOS_Vars      ,ONLY: R,sKappaM1,sKappaP1
#if PARABOLIC
USE MOD_EOS_Vars      ,ONLY: mu0,Pr,KappaSpr
#if PP_VISC == 1
USE MOD_EOS_Vars      ,ONLY: Ts,cSuth
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
USE MOD_EOS_Vars      ,ONLY: Tref,ExpoSuth
#endif
#endif /*PARABOLIC*/

! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: BulkMach,BulkReynolds
LOGICAL :: UseNonDimensionalEqn=.FALSE.
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT IDEAL GAS...'


UseNonDimensionalEqn=GETLOGICAL('UseNonDimensionalEqn','.FALSE.')
IF(UseNonDimensionalEqn)THEN
  BulkMach    =GETREAL('BulkMach')
  BulkReynolds=GETREAL('BulkReynolds')
END IF

! Gas constants
Kappa    =GETREAL('kappa','1.4')
KappaM1  =Kappa-1.
sKappaM1 =1./KappaM1
KappaP1  =Kappa+1.
sKappaP1 =1./(KappaP1)
IF(.NOT. UseNonDimensionalEqn)THEN
  R=GETREAL('R','287.058')
ELSE
  R=1./(Kappa*BulkMach*BulkMach)
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                              R | Set to 1/(Kappa*BulkMach**2)=',R
END IF

cp=R*kappa/(kappa-1.)
cv=R/(kappa-1.)

#if PARABOLIC
Pr       =GETREAL('Pr','0.72')
KappasPr =Kappa/Pr

! Viscosity
#if   PP_VISC == 0
! Constant viscosity
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0','0.0')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
END IF
#elif PP_VISC == 1
! mu-Sutherland
mu0     =GETREAL('mu0','1.735E-5')
Ts      =GETREAL('Ts','110.4')
Tref    =1./GETREAL('Tref','280.')
ExpoSuth=GETREAL('ExpoSuth','1.5')
Ts      =Ts*Tref
cSuth   =Ts**ExpoSuth*(1+Ts)/(2*Ts*Ts)
#elif PP_VISC == 2
! mu power-law
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0','0.')
  Tref    =GETREAL('Tref')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
  Tref=1.
  SWRITE(UNIT_stdOut,*)'|                           Tref | Set to 1.'
END IF
ExpoSuth=GETREAL('ExpoSuth')
mu0     =mu0/Tref**ExpoSuth
#endif
#endif /*PARABOLIC*/

SWRITE(UNIT_stdOut,'(A)')' INIT IDEAL-GAS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEos


!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single state
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim(prim,cons)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(PP_nVar)     !< vector of conservative variables
REAL,INTENT(OUT) :: prim(PP_nVarPrim) !< vector of primitive variables (density,velocities,temperature,pressure)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! rho
prim(1)=cons(1)
! vel/rho
prim(2:4)=cons(2:4)*sRho
! pressure
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
! temperature
prim(6) = prim(5)*sRho / R
END SUBROUTINE ConsToPrim

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables on a single side
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim_Side(Nloc,prim,cons)
! MODULES
USE MOD_Eos_Vars, ONLY:KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: cons(PP_nVar    ,0:Nloc,0:Nloc) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc) !< vector of primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
INTEGER          :: p,q
!==================================================================================================================================
DO q=0,Nloc; DO p=0,Nloc
  sRho=1./cons(1,p,q)
  ! rho
  prim(1,p,q)=cons(1,p,q)
  ! vel/rho
  prim(2:4,p,q)=cons(2:4,p,q)*sRho
  ! pressure
  prim(5,p,q)=KappaM1*(cons(5,p,q)-0.5*SUM(cons(2:4,p,q)*prim(2:4,p,q)))
  ! temperature
  prim(6,p,q) = prim(5,p,q)*sRho / R
END DO; END DO
END SUBROUTINE ConsToPrim_Side

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables in the whole volume
!==================================================================================================================================
PURE SUBROUTINE ConsToPrim_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Eos_Vars, ONLY:KappaM1,R
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:Nloc,1:nElems) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:Nloc,1:nElems) !< vector of primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
INTEGER          :: i,j,k,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    sRho=1./cons(1,i,j,k,iElem)
    ! rho
    prim(1,i,j,k,iElem)=cons(1,i,j,k,iElem)
    ! vel/rho
    prim(2:4,i,j,k,iElem)=cons(2:4,i,j,k,iElem)*sRho
    ! pressure
    prim(5,i,j,k,iElem)=KappaM1*(cons(5,i,j,k,iElem)-0.5*SUM(cons(2:4,i,j,k,iElem)*prim(2:4,i,j,k,iElem)))
    ! temperature
    prim(6,i,j,k,iElem) = prim(5,i,j,k,iElem)*sRho / R
  END DO; END DO; END DO! i,j,k=0,Nloc
END DO ! iElem
END SUBROUTINE ConsToPrim_Volume

!==================================================================================================================================
!> Transformation from primitive to conservative variables
!==================================================================================================================================
PURE SUBROUTINE PrimToCons(prim,cons)
! MODULES
USE MOD_Eos_Vars,ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: prim(PP_nVarPrim) !< vector of primitive variables
REAL,INTENT(OUT) :: cons(PP_nVar)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! conversion
cons(1)=prim(1)
! rho
cons(2:4)=prim(2:4)*prim(1)
! vel/rho
cons(5)=sKappaM1*prim(5)+0.5*SUM(cons(2:4)*prim(2:4))
! inner energy
END SUBROUTINE PrimToCons

!==================================================================================================================================
!> Transformation from primitive to conservative variables on a single side
!==================================================================================================================================
PURE SUBROUTINE PrimToCons_Side(Nloc,prim,cons)
! MODULES
USE MOD_Eos_Vars,ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN):: Nloc
REAL,INTENT(IN)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc) !< vector of primitive variables
REAL,INTENT(OUT)  :: cons(PP_nVar    ,0:Nloc,0:Nloc)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q
!==================================================================================================================================
DO q=0,Nloc; DO p=0,Nloc
  ! conversion
  cons(1,p,q)=prim(1,p,q)
  ! rho
  cons(2:4,p,q)=prim(2:4,p,q)*prim(1,p,q)
  ! vel/rho
  cons(5,p,q)=sKappaM1*prim(5,p,q)+0.5*SUM(cons(2:4,p,q)*prim(2:4,p,q))
END DO; END DO ! p,q=0,Nloc
! inner energy
END SUBROUTINE PrimToCons_Side


PURE FUNCTION PRESSURE_RIEMANN(U_Prim)
!==================================================================================================================================
! MODULES
USE MOD_Eos_Vars      ,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: U_Prim(PP_nVarPrim)
REAL            :: PRESSURE_RIEMANN
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL     :: kappaFac,ar,br,P_RP
!==================================================================================================================================
kappaFac=2.*Kappa*sKappaM1
IF(U_Prim(2) .LE. 0.)THEN ! rarefaction
  P_RP=U_Prim(5) * MAX(0.0001,(1.+0.5*KappaM1*U_Prim(2)/SQRT(Kappa*U_Prim(5)/U_Prim(1))))**kappaFac
ELSE ! shock
  ar=2.*sKappaP1/U_Prim(1)
  br=KappaM1*sKappaP1*U_Prim(5)
  P_RP=U_Prim(5)+U_Prim(2)/ar*0.5*(U_Prim(2)+SQRT(U_Prim(2)*U_Prim(2)+4.*ar*(U_Prim(5)+br)))
END IF
PRESSURE_RIEMANN=P_RP
END FUNCTION PRESSURE_RIEMANN


END MODULE MOD_EOS
