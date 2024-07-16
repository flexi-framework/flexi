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

INTERFACE PRESSURE_RIEMANN
  MODULE PROCEDURE PRESSURE_RIEMANN
END INTERFACE

PUBLIC::InitEos
PUBLIC::ConsToPrim
PUBLIC::PrimToCons
PUBLIC::ConsToEntropy
PUBLIC::EntropyToCons
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
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent"        , '1.4')
CALL prms%CreateRealOption(     'R',            "Specific gas constant"                            , '287.058')
CALL prms%CreateRealOption(     'Pr',           "Prandtl number"                                   , '0.72')
CALL prms%CreateRealOption(     'mu0',          "Dynamic Viscosity"                                , '0.')
CALL prms%CreateRealOption(     'Ts',           "Sutherland's law for variable viscosity: Ts"      , '110.4')
CALL prms%CreateRealOption(     'Tref',         "Sutherland's law for variable viscosity: Tref"    , '273.15')
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
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT IDEAL GAS...'


UseNonDimensionalEqn=GETLOGICAL('UseNonDimensionalEqn')
IF(UseNonDimensionalEqn)THEN
  BulkMach    =GETREAL('BulkMach')
  BulkReynolds=GETREAL('BulkReynolds')
END IF

! Gas constants
Kappa    = GETREAL('kappa')
KappaM1  = Kappa-1.
sKappaM1 = 1./KappaM1
KappaP1  = Kappa+1.
sKappaP1 = 1./(KappaP1)
IF(.NOT. UseNonDimensionalEqn)THEN
  R = GETREAL('R')
ELSE
  R = 1./(Kappa*BulkMach*BulkMach)
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                              R | Set to 1/(Kappa*BulkMach**2)=',R
END IF

cp=R*kappa/(kappa-1.)
cv=R/(kappa-1.)

#if PARABOLIC
Pr       =GETREAL('Pr')
KappasPr =Kappa/Pr

! Viscosity
#if   PP_VISC == 0
! Constant viscosity
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
END IF
#elif PP_VISC == 1
! mu-Sutherland
! Coefficients from White, F. M., Viscous fluid flow, McGraw-Hill, 2006
mu0     =GETREAL('mu0')
Ts      =GETREAL('Ts')
Tref    =1./GETREAL('Tref')
ExpoSuth=GETREAL('ExpoSuth')
Ts      =Ts*Tref
cSuth   =Ts**ExpoSuth*(1+Ts)/(2*Ts*Ts)
#elif PP_VISC == 2
! mu power-law
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0')
  Tref    =GETREAL('Tref')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
  Tref=1.
  SWRITE(UNIT_stdOut,'(A)')       ' |                           Tref | Set to 1.'
END IF
ExpoSuth=GETREAL('ExpoSuth')
mu0     =mu0/Tref**ExpoSuth
#endif
#endif /*PARABOLIC*/

SWRITE(UNIT_stdOut,'(A)')' INIT IDEAL-GAS DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEos

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single state
!==================================================================================================================================
PPURE SUBROUTINE ConsToPrim(prim,cons)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(CONS)     !< vector of conservative variables
REAL,INTENT(OUT) :: prim(PRIM)     !< vector of primitive variables (density,velocities,temperature,pressure)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(DENS)
! density
prim(DENS)=cons(DENS)
! velocity
prim(VEL1:VEL2)=cons(MOM1:MOM2)*sRho
#if (PP_dim==3)
prim(VEL3)=cons(MOM3)*sRho
#else
prim(VEL3)=0.
#endif
! pressure
prim(PRES)=KappaM1*(cons(ENER)-0.5*SUM(cons(MOMV)*prim(VELV)))
! temperature
prim(TEMP) = prim(PRES)*sRho / R
! kinematic SA viscosity
prim(NUSA) = cons(MUSA)*sRho
END SUBROUTINE ConsToPrim

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables on a single side
!==================================================================================================================================
PPURE SUBROUTINE ConsToPrim_Side(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                           !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: cons(CONS,0:Nloc,0:ZDIM(Nloc)) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PRIM,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  CALL ConsToPrim(prim(:,p,q),cons(:,p,q))
END DO; END DO
END SUBROUTINE ConsToPrim_Side

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables in the whole volume
!==================================================================================================================================
PPURE SUBROUTINE ConsToPrim_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                           !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: cons(CONS,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PRIM,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems) !< vector of primitive variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    CALL ConsToPrim(prim(:,i,j,k,iElem),cons(:,i,j,k,iElem))
  END DO; END DO; END DO! i,j,k=0,Nloc
END DO ! iElem
END SUBROUTINE ConsToPrim_Volume

!==================================================================================================================================
!> Transformation from primitive to conservative variables for a single state
!==================================================================================================================================
PPURE SUBROUTINE PrimToCons(prim,cons)
! MODULES
USE MOD_EOS_Vars,ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: prim(PRIM)     !< vector of primitive variables
REAL,INTENT(OUT) :: cons(CONS)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! density
cons(DENS)=prim(DENS)
! momentum
cons(MOM1:MOM2)=prim(VEL1:VEL2)*prim(DENS)
#if (PP_dim==3)
cons(MOM3)=prim(VEL3)*prim(DENS)
#else
cons(MOM3)=0.
#endif
! energy
cons(ENER)=sKappaM1*prim(PRES)+0.5*SUM(cons(MOMV)*prim(VELV))
! dynamic SA viscosity
cons(MUSA) = prim(NUSA)*prim(DENS)
END SUBROUTINE PrimToCons

!==================================================================================================================================
!> Transformation from primitive to conservative variables on a single side
!==================================================================================================================================
PPURE SUBROUTINE PrimToCons_Side(Nloc,prim,cons)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                           !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: prim(PRIM,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(CONS,0:Nloc,0:ZDIM(Nloc)) !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  CALL PrimToCons(prim(:,p,q),cons(:,p,q))
END DO; END DO ! p,q=0,Nloc
END SUBROUTINE PrimToCons_Side

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PPURE SUBROUTINE PrimToCons_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                               !< local polynomial degree of solution representation
REAL,INTENT(IN)    :: prim(PRIM,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of primitive variables
REAL,INTENT(OUT)   :: cons(CONS,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    CALL PrimToCons(prim(:,i,j,k,iElem),cons(:,i,j,k,iElem))
  END DO; END DO; END DO
END DO
END SUBROUTINE PrimToCons_Volume

!==================================================================================================================================
!> Transformation from conservative variables U to entropy vector, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PPURE SUBROUTINE ConsToEntropy(entropy,cons)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,kappa,sKappaM1
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
REAL                                :: vel(3),s,p,rho_p
!==================================================================================================================================
vel(:) = cons(MOMV)/cons(DENS)
p      = KappaM1*(cons(ENER)-0.5*SUM(cons(MOMV)*vel(:)))
! entropy: log(p) - eq.γ * log(ρ)
s      = LOG(p) - kappa*LOG(cons(DENS))
rho_p  = cons(DENS)/p

! Convert to entropy variables
entropy(DENS)      = (kappa-s)*skappaM1 - rho_p * 0.5*SUM(vel**2)  ! (γ - s) / (γ - 1) - (ρu^2 + ρv^2 + ρw^2) / ρ / 2p,
entropy(MOM1:MOM2) = rho_p*vel(1:2)        ! ρu / p
#if (PP_dim==3)
entropy(MOM3)      = rho_p*vel(3)
#else
entropy(MOM3)      = 0.
#endif
entropy(ENER)      = - rho_p          ! -ρ / p

entropy(MUSA)      = cons(MUSA)

END SUBROUTINE ConsToEntropy

!==================================================================================================================================
!> Transformation from entropy to conservative variables U, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PPURE SUBROUTINE EntropyToCons(entropy,cons)
! MODULES
USE MOD_EOS_Vars,ONLY:KappaM1,kappa,sKappaM1
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
entropy2 = entropy*kappaM1
s        = kappa - entropy2(DENS) + 0.5 * SUM(entropy2(MOMV)**2) / entropy2(ENER)
rhoe     = (kappaM1 / ((-entropy2(ENER))**kappa))**(skappaM1) * EXP(-s*skappaM1)

cons(DENS) = - rhoe * entropy2(ENER) ! ρ = -p * W[5]
cons(MOMV) = rhoe * entropy2(MOMV)
#if (PP_dim==2)
cons(MOM3) = 0.
#endif
cons(ENER) = rhoe * (1 - SUM(entropy2(MOMV)**2) * 0.5/ entropy2(ENER)) !sKappaM1*p+0.5*SUM(cons(MOMV)*vel)

cons(MUSA) = entropy(MUSA)

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
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    CALL ConsToEntropy(entropy(:,i,j,k,iElem),cons(:,i,j,k,iElem))
  END DO; END DO; END DO
END DO
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
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  CALL EntropyToCons(entropy(:,p,q),cons(:,p,q))
END DO; END DO ! p,q=0,Nloc
END SUBROUTINE EntropyToCons_Side


!==================================================================================================================================
!> Riemann solver function to get pressure at BCs
!==================================================================================================================================
PPURE FUNCTION PRESSURE_RIEMANN(U_Prim)
!==================================================================================================================================
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: U_Prim(PRIM)        !< vector of primitive variables
REAL            :: PRESSURE_RIEMANN    !< pressure as the return value of the Riemann problem
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: kappaFac,ar,br,P_RP
!==================================================================================================================================
kappaFac=2.*Kappa*sKappaM1
IF(U_Prim(VEL1) .LE. 0.)THEN ! rarefaction
  P_RP=U_Prim(PRES) * MAX(0.0001,(1.+0.5*KappaM1*U_Prim(VEL1)/SQRT(Kappa*U_Prim(PRES)/U_Prim(DENS))))**kappaFac
ELSE ! shock
  ar=2.*sKappaP1/U_Prim(DENS)
  br=KappaM1*sKappaP1*U_Prim(PRES)
  P_RP=U_Prim(PRES)+U_Prim(VEL1)/ar*0.5*(U_Prim(VEL1)+SQRT(U_Prim(VEL1)*U_Prim(VEL1)+4.*ar*(U_Prim(PRES)+br)))
END IF
PRESSURE_RIEMANN=P_RP
END FUNCTION PRESSURE_RIEMANN


END MODULE MOD_EOS
