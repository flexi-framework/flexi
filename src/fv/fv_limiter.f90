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
#if FV_ENABLED
#include "flexi.h"

!==================================================================================================================================
!> Module for slope limiters of FV sub-cells.
!==================================================================================================================================
MODULE MOD_FV_Limiter
! MODULES
IMPLICIT NONE
PRIVATE

#if FV_RECONSTRUCT
INTEGER,PARAMETER :: FV_LIMITERTYPE_NULL      = 0
INTEGER,PARAMETER :: FV_LIMITERTYPE_MINMOD    = 1
INTEGER,PARAMETER :: FV_LIMITERTYPE_SWEBY     = 2
INTEGER,PARAMETER :: FV_LIMITERTYPE_VANALBADA = 3
INTEGER,PARAMETER :: FV_LIMITERTYPE_GMINMOD   = 4
INTEGER,PARAMETER :: FV_LIMITERTYPE_OSPRE     = 5
INTEGER,PARAMETER :: FV_LIMITERTYPE_CENTRAL   = 9

INTERFACE DefineParametersFV_Limiter
  MODULE PROCEDURE DefineParametersFV_Limiter
END INTERFACE

INTERFACE InitFV_Limiter
  MODULE PROCEDURE InitFV_Limiter
END INTERFACE

ABSTRACT INTERFACE
  PPURE SUBROUTINE LimiterInt(sL, sR, s)
    REAL,INTENT(IN)  :: sL(PP_nVarPrim),sR(PP_nVarPrim)
    REAL,INTENT(OUT) :: s(PP_nVarPrim)
  END SUBROUTINE
END INTERFACE

PROCEDURE(LimiterInt),POINTER :: FV_Limiter !< limiting function (see: fv_limiter.f90)

PUBLIC::DefineParametersFV_Limiter
PUBLIC::InitFV_Limiter
PUBLIC::FV_Limiter
!==================================================================================================================================

#endif /* FV_RECONSTRUCT */

CONTAINS

#if FV_RECONSTRUCT

!==================================================================================================================================
!> Define parameters for FV Limiter
!==================================================================================================================================
SUBROUTINE DefineParametersFV_Limiter()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('FV')
CALL prms%CreateIntFromStringOption('FV_LimiterType',"Type of slope limiter of second order reconstruction", '1')
CALL addStrListEntry('FV_LimiterType','none',     FV_LIMITERTYPE_NULL)
CALL addStrListEntry('FV_LimiterType','minmod',   FV_LIMITERTYPE_MINMOD)
CALL addStrListEntry('FV_LimiterType','sweby',    FV_LIMITERTYPE_SWEBY)
CALL addStrListEntry('FV_LimiterType','vanalbada',FV_LIMITERTYPE_VANALBADA)
CALL addStrListEntry('FV_LimiterType','gminmod'  ,FV_LIMITERTYPE_GMINMOD)
CALL addStrListEntry('FV_LimiterType','ospre',    FV_LIMITERTYPE_OSPRE)
CALL addStrListEntry('FV_LimiterType','central',  FV_LIMITERTYPE_CENTRAL)
CALL prms%CreateRealOption('swebyb', "beta parameter for Sweby limiter")
END SUBROUTINE DefineParametersFV_Limiter

!==================================================================================================================================
!> Initialize pointer to chosen limiter type and readin of required parameter
!==================================================================================================================================
SUBROUTINE InitFV_Limiter()
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_FV_Vars     ,ONLY: LimiterType,FV_sweby_beta
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Read in LimiterType parameters
LimiterType = GETINTFROMSTR('FV_LimiterType')
SELECT CASE(LimiterType)
CASE (FV_LIMITERTYPE_NULL) ! NullLimiter
  FV_Limiter => NullLimiter
  SWRITE(UNIT_stdOut,'(A)') '  Using "NullLimiter" limiter.'
CASE (FV_LIMITERTYPE_MINMOD) ! MinMod
  FV_Limiter => MinMod
  SWRITE(UNIT_stdOut,'(A)') '  Using "Minmod" limiter.'
CASE (FV_LIMITERTYPE_SWEBY) ! Sweby
  FV_Limiter => Sweby
  FV_sweby_beta = GETREAL('swebyb')
  SWRITE(UNIT_stdOut,'(A,F8.6)') '  Using "Sweby" limiter with beta = ', FV_sweby_beta
CASE (FV_LIMITERTYPE_VANALBADA) ! van Albada
  FV_Limiter => VanAlbada
  SWRITE(UNIT_stdOut,'(A)') '  Using "van Albada" limiter.'
CASE (FV_LIMITERTYPE_GMINMOD) ! GMinMod
  FV_Limiter => GMinMod
  SWRITE(UNIT_stdOut,'(A)') '  Using "GMinmod" limiter.'
CASE (FV_LIMITERTYPE_OSPRE) ! Ospre
  FV_Limiter => Ospre
  SWRITE(UNIT_stdOut,'(A)') '  Using "Ospre" limiter.'
CASE (FV_LIMITERTYPE_CENTRAL) ! Central
  FV_Limiter => CentralLimiter
  SWRITE(UNIT_stdOut,'(A,F8.6)') '  Using "Central" limiter.'
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'FV Limiter-Type unknown.')
END SELECT
END SUBROUTINE InitFV_Limiter

!==================================================================================================================================
!> Limiter sets slope to zero.
!==================================================================================================================================
PPURE SUBROUTINE NullLimiter(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(PP_nVarPrim) !< left slope
REAL,INTENT(IN)  :: sR(PP_nVarPrim) !< right slope
REAL,INTENT(OUT) :: s(PP_nVarPrim)  !< limited slope
!==================================================================================================================================
! NullLimiter
s = 0.
END SUBROUTINE NullLimiter

!==================================================================================================================================
!> MinMod slope limiter.
!==================================================================================================================================
PPURE SUBROUTINE MinMod(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(PP_nVarPrim) !< left slope
REAL,INTENT(IN)  :: sR(PP_nVarPrim) !< right slope
REAL,INTENT(OUT) :: s(PP_nVarPrim)  !< limited slope
!==================================================================================================================================
! MinMod
s = MERGE(sL,sR, ABS(sL) .LT. ABS(sR))
s = MERGE(s,0., sL*sR .GT. 0.)
!if (sL(1)*sR(1) .LT.0) then
  !if (abs(sl(1)*sR(1)).gt.0.00001) then
     !WRITE (*,*) sl(1), sr(1), s(1)
  !end if
!end if
END SUBROUTINE MinMod

!==================================================================================================================================
!> Sweby slope limiter.
!==================================================================================================================================
PPURE SUBROUTINE Sweby(sL, sR, s)
! MODULES
USE MOD_FV_Vars ,ONLY: FV_sweby_beta
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(PP_nVarPrim) !< left slope
REAL,INTENT(IN)  :: sR(PP_nVarPrim) !< right slope
REAL,INTENT(OUT) :: s(PP_nVarPrim)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: sa(PP_nVarPrim),sb(PP_nVarPrim)
!==================================================================================================================================
CALL MinMod(sL*FV_sweby_beta,sR,sa)
CALL MinMod(sL,sR*FV_sweby_beta,sb)
s = SIGN(MAX(ABS(sa),ABS(sb)),sL)
END SUBROUTINE Sweby

!==================================================================================================================================
!> van Albada slope limiter.
!==================================================================================================================================
PPURE SUBROUTINE VanAlbada(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(PP_nVarPrim) !< left slope
REAL,INTENT(IN)  :: sR(PP_nVarPrim) !< right slope
REAL,INTENT(OUT) :: s(PP_nVarPrim)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! van Albada
s = (sL*sR*(sL+sR))/MAX(sL**2+sR**2,1e-13)
s = MERGE(s,0., sL*sR .GT. 0.)
END SUBROUTINE VanAlbada

!==================================================================================================================================
!> GMinmod slope limiter.
!==================================================================================================================================
PURE SUBROUTINE GMinMod(sL, sR,s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(PP_nVarPrim) !< left slope
REAL,INTENT(IN)  :: sR(PP_nVarPrim) !< right slope
REAL,INTENT(OUT) :: s(PP_nVarPrim)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER :: beta = 2.0d0
REAL :: s1(PP_nVarPrim),s2(PP_nVarPrim)
!==================================================================================================================================
  s1 = sign(1.0d0, sR)
  s2 = sign(1.0d0, sL)
  s = 0.5d0*(s1*max(0.0d0, min(beta*sR*s1, sL*s1)) + &
             s2*max(0.0d0, min(beta*sL*s2, sR*s2)))
END SUBROUTINE GMinMod

!==================================================================================================================================
!> OSPRE slope limiter.
!==================================================================================================================================
PURE SUBROUTINE Ospre(sL, sR,s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(PP_nVarPrim) !< left slope
REAL,INTENT(IN)  :: sR(PP_nVarPrim) !< right slope
REAL,INTENT(OUT) :: s(PP_nVarPrim)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: d(PP_nVarPrim)
!==================================================================================================================================
  d = sL**2 + sL*sR + sR**2
  s = (d/(d**2 + 1.0d-28))*1.5d0*sL*sR*(sL + sR)
END SUBROUTINE Ospre

!==================================================================================================================================
!> central limiter s = (sL + sR)/2  (ATTENTION: unstable and not TVD)
!==================================================================================================================================
PPURE SUBROUTINE CentralLimiter(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(PP_nVarPrim) !< left slope
REAL,INTENT(IN)  :: sR(PP_nVarPrim) !< right slope
REAL,INTENT(OUT) :: s(PP_nVarPrim)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
s = (sL+sR)/2.0
END SUBROUTINE CentralLimiter

#endif /* FV_RECONSTRUCT */

END MODULE MOD_FV_Limiter
#endif /* FV_ENABLED */
