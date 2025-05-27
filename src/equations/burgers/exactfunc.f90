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
MODULE MOD_Exactfunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::DefineParametersExactFunc
PUBLIC::InitExactFunc
PUBLIC::ExactFunc
PUBLIC::CalcSource
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersExactFunc()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Exactfunc")
CALL prms%CreateIntFromStringOption(      'IniExactFunc', "Number of exactfunction to be used, to initialize the solution. ")
CALL addStrListEntry('IniExactFunc','constant' ,1)
CALL addStrListEntry('IniExactFunc','sinex'    ,11)
CALL addStrListEntry('IniExactFunc','jumpx'    ,12)
CALL addStrListEntry('IniExactFunc','rarex'    ,13)
CALL addStrListEntry('IniExactFunc','siney'    ,14)
CALL addStrListEntry('IniExactFunc','jumpy'    ,15)
CALL addStrListEntry('IniExactFunc','rarey'    ,16)
CALL addStrListEntry('IniExactFunc','sinexy'   ,21)
CALL addStrListEntry('IniExactFunc','Uwe'      ,99)
END SUBROUTINE DefineParametersExactFunc

!==================================================================================================================================
!> Get some parameters needed for exact function
!==================================================================================================================================
SUBROUTINE InitExactFunc()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools,   ONLY: GETINTFROMSTR,GETREAL,GETINT
USE MOD_Equation_Vars, ONLY: IniExactFunc,IniRefState
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT FUNCTION...'

! Read in boundary parameters
IniExactFunc = GETINTFROMSTR('IniExactFunc')
IniRefState  = -1 ! only dummy for linadv

! Read in parameters specific to certain init functions
SELECT CASE (IniExactFunc)
  CASE DEFAULT
    ! Nothing to do
END SELECT

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT FUNCTION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitExactFunc

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu,RefStateOpt)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Timedisc_Vars, ONLY: t
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: tIn                    !< input time (either time at RK stage or time at the beginning of
                                                          !< timestep if full boundary order is used (only with RK3)
REAL,INTENT(IN)                 :: x(3)                   !< coordinates to evaluate exact function
INTEGER,INTENT(IN)              :: ExactFunction          !< specifies the exact function to be used
REAL,INTENT(OUT)                :: Resu(PP_nVar)          !< output state in conservative variables
INTEGER,INTENT(IN),OPTIONAL     :: RefStateOpt            !< refstate to be used for exact func (dummy)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: x0(3)
REAL                            :: Pi
!==================================================================================================================================
Pi   = ACOS(-1.)

Resu = 0.

SELECT CASE (ExactFunction)
CASE(1) ! constant
  Resu(1) = 1.0
  Resu(2) = 0.3
CASE(11)  ! sine in x direction
  Resu(1) = SIN(2.*Pi*x(1))
CASE(12)  ! shockx
  IF(x(1).LT.0.50)THEN
    Resu(1) = 1.
  END IF
CASE(13)  ! rarefactionx
  IF(x(1).LT.0.5)THEN
    Resu(1) = -1.
  ELSE
    Resu(1) =  1.
  END IF
CASE(14)  ! sine in y direction
  Resu(2) = SIN(2.*Pi*x(2))
CASE(15)  ! shocky
  IF(x(2).LT.0.5)THEN
    Resu(2) = 1.
  END IF
CASE(16)  ! rarefactionx
  IF(x(2).LT.0.5)THEN
    Resu(2) = -1.
  ELSE
    Resu(2) =  1.
  END IF
CASE(21)  ! sine in xy direction
  Resu(1) = 0.1*SIN(2.*Pi*SUM(x))
  Resu(2) = Resu(1)
CASE(99)  ! Uwe
  Resu(1) = 5*x(2)*(1.-x(2))*(1.-EXP(-t))
CASE DEFAULT
  CALL abort(__STAMP__,&
             'Specified exactfuntion not implemented!')
END SELECT ! ExactFunction

END SUBROUTINE ExactFunc



!==================================================================================================================================
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,t)
! MODULES
USE MOD_Globals,       ONLY:Abort
USE MOD_Equation_Vars, ONLY:IniExactFunc,doCalcSource
USE MOD_Mesh_Vars,     ONLY:nElems,Elem_xGP,sJ
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                        !< solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k
REAL                :: Ut_src(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE DEFAULT
  doCalcSource=.FALSE.
END SELECT ! ExactFunction
END SUBROUTINE CalcSource

END MODULE MOD_Exactfunc
