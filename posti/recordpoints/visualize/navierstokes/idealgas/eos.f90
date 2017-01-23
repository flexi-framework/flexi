#include "flexi.h"

MODULE MOD_EOS
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitEOS
  MODULE PROCEDURE InitEOS
END INTERFACE
INTERFACE CalcPrims
  MODULE PROCEDURE CalcPrims
END INTERFACE

INTERFACE FinalizeEOS
  MODULE PROCEDURE FinalizeEOS
END INTERFACE

PUBLIC::InitEOS,FinalizeEOS
PUBLIC::CalcPrims
!===================================================================================================================================

CONTAINS

SUBROUTINE InitEOS()
!===================================================================================================================================
! Computes the gradient of the conservative variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools           ,ONLY:GETINT,GETREAL
USE MOD_EOS_Vars
USE MOD_VarNameMappingsRP_Vars,ONLY:Prim
USE MOD_VarNameMappingsRP     ,ONLY:CreateVarMappings
! USE MOD_Viscosity    ,ONLY:InitViscosity
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: nPrim
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT EOS (ideal gas) ...'

nPrim=12
ALLOCATE(VarNamePrim(nPrim))
VarNamePrim(1) ='VelocityX'
VarNamePrim(2) ='VelocityY'
VarNamePrim(3) ='VelocityZ'
VarNamePrim(4) ='VelocityMagnitude'
VarNamePrim(5) ='Pressure'
VarNamePrim(6) ='VelocitySound'
VarNamePrim(7) ='Mach'
VarNamePrim(8) ='Temperature'
VarNamePrim(9) ='EnergyStagnation'
VarNamePrim(10)='EnthalpyStagnation'
VarNamePrim(11)='Entropy'
VarNamePrim(12)='PotTemp'

CALL CreateVarMappings(nPrim,VarNamePrim,Prim)

R      =GETREAL('R','287.058')
Kappa  =GETREAL('kappa','1.4')
KappaM1=Kappa-1
sKappaM1=1./(KappaM1)
sKappa=1./Kappa

P0  =GETREAL('P0','10e5')

WRITE(UNIT_stdOut,'(A)')' INIT EOS DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEOS


SUBROUTINE CalcPrims(nVar,nRP,U_in,U_out)
!===================================================================================================================================
! Computes the primitive variables on output grid
!===================================================================================================================================
! MODULES
USE MOD_Globals             ,ONLY:Abort
USE MOD_VarNameMappingsRP_Vars,ONLY: Prim
USE MOD_EOS_Vars            ,ONLY: kappa,kappaM1,sKappaM1,R
USE MOD_EOS_Vars            ,ONLY: sKappa,P0
USE MOD_Equation_Vars       ,ONLY: PrimMap
USE MOD_Parameters      ,ONLY: usePrims
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: nVar,nRP
REAL,INTENT(IN)              :: U_in(nVar,nRP)
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(OUT)             :: U_out(Prim%nVar_visu,nRP)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                         :: PrimVar(Prim%nVar),Cons(nVar),sRho
INTEGER                      :: iRP,iVar
!===================================================================================================================================
IF(usePrims) THEN
  DO iRP=1,nRP
    PrimVar=0.
    Cons(:)=U_in(:,iRP)
    sRho=1./Cons(1)
    ! Velocities
    PrimVar(1:3)=Cons(PrimMap(2:4))
    ! VelocityMagnitude
    PrimVar(4)=SQRT(SUM(PrimVar(1:3)*PrimVar(1:3)))
    ! Pressure is directly used from state file
    PrimVar(5)=Cons(PrimMap(5))
    ! VelocitySound
    PrimVar(6)=SQRT(Kappa*PrimVar(5)*sRho)
    ! Mach
    PrimVar(7)=PrimVar(4)/PrimVar(6)
    ! Temperature
    PrimVar(8)=PrimVar(5)*sRho/R
    ! EnergyStagnation
    PrimVar(9)=sKappaM1*PrimVar(5)*sRho+PrimVar(4)*PrimVar(4)
    ! EnthalpyStagnation
    PrimVar(10)=PrimVar(9)+PrimVar(5)*sRho
    ! Entropy
    PrimVar(11)= R*(sKappaM1*LOG(PrimVar(8))-LOG(Cons(PrimMap(1)))) 
    ! Potential Temperature 
    PrimVar(12)=PrimVar(8)/(PrimVar(5)/P0)**(1.-sKappa)
    ! write desired PrimVars to output
    U_out(:,iRP)=PrimVar(Prim%Ind)
  END DO  ! iRP
ELSE
DO iRP=1,nRP
  PrimVar=0.
  Cons(:)=U_in(:,iRP)
  sRho=1./Cons(1)
  ! Velocities
  DO iVar=1,3
    PrimVar(iVar)=Cons(iVar+1)*sRho
  END DO
  ! VelocityMagnitude
  PrimVar(4)=SQRT(SUM(PrimVar(1:3)*PrimVar(1:3)))
  ! Pressure
  PrimVar(5)=KappaM1*(Cons(5)-Cons(1)*PrimVar(4)*PrimVar(4)*0.5)
  ! VelocitySound
  PrimVar(6)=SQRT(Kappa*PrimVar(5)*sRho)
  ! Mach
  PrimVar(7)=PrimVar(4)/PrimVar(6)
  ! Temperature
  PrimVar(8)=PrimVar(5)*sRho/R
  ! EnergyStagnation
  PrimVar(9)=Cons(5)*sRho
  ! EnthalpyStagnation
  PrimVar(10)=PrimVar(9)+PrimVar(5)*sRho
  ! Entropy
  PrimVar(11)= R*(sKappaM1*LOG(PrimVar(8))-LOG(Cons(1))) 
  ! Potential Temperature 
  PrimVar(12)=PrimVar(8)/(PrimVar(5)/P0)**(1.-sKappa)
  ! write desired PrimVars to output
  U_out(:,iRP)=PrimVar(Prim%Ind)
END DO  ! iRP
END IF
END SUBROUTINE CalcPrims





SUBROUTINE FinalizeEOS()
!===================================================================================================================================
! Computes the gradient of the conservative variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_EOS_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)') '  EOS FINALIZED'
END SUBROUTINE FinalizeEOS

END MODULE MOD_EOS
