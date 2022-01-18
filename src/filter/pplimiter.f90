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
#if EQNSYSNR == 2
#include "eos.h"
#endif
!==================================================================================================================================
! Positivity-preserving limiter from:
! Xiangxiong Zhang, Chi-Wang Shu, Journal of Computational Physics (2010):
! "On positivity-preserving high order discontinuous Galerkin schemes for compressible Euler equations on rectangular meshes",
! Volume 229, Issue 23, Pages 8918-8934,
! https://doi.org/10.1016/j.jcp.2010.08.016.
!==================================================================================================================================
MODULE MOD_PPLimiter
#if PP_LIMITER
  !----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitPPLimiter
  MODULE PROCEDURE InitPPLimiter
END INTERFACE

INTERFACE PPLimiter
  MODULE PROCEDURE PPLimiter
END INTERFACE

INTERFACE FinalizePPLimiter
  MODULE PROCEDURE FinalizePPLimiter
END INTERFACE

INTERFACE PPLimiter_Info
  MODULE PROCEDURE PPLimiter_Info
END INTERFACE

PUBLIC:: DefineParametersPPLimiter
PUBLIC:: InitPPLimiter
PUBLIC:: PPLimiter
PUBLIC:: FinalizePPLimiter
PUBLIC:: PPLimiter_Info
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters needed for filtering
!==================================================================================================================================
SUBROUTINE DefineParametersPPLimiter()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Positivity-Preserving Limiter")
CALL prms%CreateLogicalOption('DoPPLimiter',       "Toggles if positivity preserving limiter is used.",       ".FALSE.")
CALL prms%CreateRealOption(   'PPThresholdFactor', "Factor for PP limiter thresholds relative to PPRefState", "1.E-10")
CALL prms%CreateIntOption(    'PPRefState',        "Index of reference state for PP limiter threshold",       "1")
END SUBROUTINE DefineParametersPPLimiter

!==================================================================================================================================
!> Initialize  information and  operators
!==================================================================================================================================
SUBROUTINE InitPPLimiter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars
USE MOD_ReadInTools        ,ONLY: GETREAL,GETLOGICAL,GETINT
USE MOD_IO_HDF5            ,ONLY: AddToFieldData
USE MOD_IO_HDF5            ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars          ,ONLY: nElems,sJ
USE MOD_Interpolation_Vars ,ONLY: wGP
#if FV_ENABLED
USE MOD_FV_Vars            ,ONLY: FV_w
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,i,j,k
REAL                         :: Vol
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT POSITIVITY-PRESERVING LIMITER...'

! Read in variables
DoPPLimiter = GETLOGICAL('DoPPLimiter','.FALSE.')
PPEpsFac    = GETREAL('PPThresholdFactor','1.E-10')
iPPRefState = GETINT('PPRefState','1')
!PPepsDens and PPepsPres are calculated in limiter, since RefStatePrim is needed from InitEquation

! Prepare PP Limiter
ALLOCATE(PP_Elems(nElems))
PP_Elems=0
! Add information about PP Limiter to Element Data output
CALL AddToElemData(ElementOut,'PP_Elems',IntArray=PP_Elems)

! arrays for PP Limiter statistics
ALLOCATE(PP_Elems_counter(nElems))
ALLOCATE(PP_Elems_Amount(nElems))
PP_Elems_counter  = 0
PP_Switch_counter = 0
PP_Elems_Amount = 0
CALL AddToElemData(ElementOut,'PP_Elems_Amount',RealArray=PP_Elems_Amount)

IF (DoPPLimiter) THEN
  ! Allocate volume and weights if not already done by LAF Filter
#if FV_ENABLED
  SDEALLOCATE(IntegrationWeight) ! different array size in last dim
#endif
  IF( .NOT. ALLOCATED(IntegrationWeight)) THEN
    ALLOCATE(IntegrationWeight(0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_ENABLED))
    DO iElem=1,nElems
      Vol = 0.
      DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
        IntegrationWeight(i,j,k,iElem,0) = wGP(i) * wGP(j) &
#if PP_dim == 3
                                                         * wGP(k) &
#endif
                                                                  / sJ(i,j,k,iElem,0)
        Vol = Vol + IntegrationWeight(i,j,k,iElem,0)
      END DO; END DO; END DO
      IntegrationWeight(:,:,:,iElem,0) = IntegrationWeight(:,:,:,iElem,0) / Vol
    END DO !iElem

#if FV_ENABLED
    DO iElem=1,nElems
      Vol = 0.
      DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
        IntegrationWeight(i,j,k,iElem,1) = FV_w**PP_dim / sJ(i,j,k,iElem,1)
        Vol = Vol + IntegrationWeight(i,j,k,iElem,1)
      END DO; END DO; END DO
      IntegrationWeight(:,:,:,iElem,1) = IntegrationWeight(:,:,:,iElem,1) / Vol
    END DO !iElem
#endif
  END IF
END IF
SWRITE(UNIT_stdOut,'(A)')' INIT POSITIVITY-PRESERVING LIMITER DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitPPLimiter

!==================================================================================================================================
!> Hyperbolicity Preserving Limiter, limits polynomial towards admissible cellmean
!==================================================================================================================================
SUBROUTINE PPLimiter()
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars             ,ONLY: U,UPrim
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Filter_Vars         ,ONLY: iPPRefState,PPEpsFac
USE MOD_Filter_Vars         ,ONLY: PPepsDens,PPepsPres
USE MOD_Filter_Vars         ,ONLY: PP_Elems,PP_Switch_counter,PP_Elems_counter,PP_Elems_Amount
USE MOD_EOS                 ,ONLY: ConsToPrim,PrimtoCons
USE MOD_Equation_Vars       ,ONLY: RefStatePrim
USE MOD_Filter_Vars         ,ONLY: IntegrationWeight
#if FV_ENABLED
USE MOD_FV_Vars             ,ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: iElem,i,j,k
REAL                   :: UMean(PP_nVar),UMeanPrim(PP_nVarPrim),rhoMin,pMin
REAL                   :: Theta1,Theta2,tLoc
!==================================================================================================================================
!InitFilter (InitPPLimiter) is called before InitEquation (where RefStatePrim is filled), so we have to do this here:
PPepsDens=PPEpsFac*RefStatePrim(DENS,iPPRefState)
PPepsPres=PPEpsFac*RefStatePrim(PRES,iPPRefState) !Isothermal expansion
!PPepsPres=PPEpsFac**Kappa*RefStatePrim(PRES,iPPRefState) !Adiabatic expansion

CALL ConsToPrim(PP_N,UPrim,U)

! reset PP element array
PP_Elems=0

DO iElem=1,nElems
  ! Set initial values
  rhoMin = MINVAL(UPrim(DENS,:,:,:,iElem))
  pMin   = MINVAL(UPrim(PRES,:,:,:,iElem))

  ! Do nothing if density and pressure are above limit
  IF ((rhoMin.GE.PPepsDens).AND.(pMin.GE.PPepsPres)) CYCLE

  ! Set element as active (for analyze/visu purposes only)
  PP_Elems(iElem)=1

  ! Calculate mean values in element
  UMean = 0.
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    UMean = UMean + U(:,i,j,k,iElem)*IntegrationWeight(i,j,k,iElem,FV_Elems(iElem))
  END DO; END DO; END DO

  ! check if mean is admissible
  CALL ConsToPrim(UMeanPrim,UMean)
  IF ((UMean(DENS).LE.PPepsDens).OR.(UMeanPrim(PRES).LE.PPepsPres)) THEN
    ! mean is not admissible: Make mean admissible and set to a constant value
    UMeanPrim(DENS) = MAX(UMeanPrim(DENS),PPepsDens)
    UMeanPrim(PRES) = MAX(UMeanPrim(PRES),PPepsPres)
    CALL PrimToCons(UMeanPrim,UMean)
    DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
      U(:,i,j,k,iElem) = UMean
    END DO;END DO;END DO
    CYCLE
  END IF

  ! Here the actual limiting starts
  ! Step 1
  IF (rhoMin.LT.PPepsDens) THEN
    Theta1 = (UMean(DENS)-PPepsDens) / (Umean(DENS)-rhoMin)
    DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
      U(DENS,i,j,k,iElem) = Theta1*(U(DENS,i,j,k,iElem)-UMean(DENS)) + UMean(DENS)
    END DO;END DO;END DO
  END IF

  ! Step 2
  Theta2 = 1.
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    tLoc   = CalcT(U(:,i,j,k,iElem),UMean)
    Theta2 = MIN(Theta2,tLoc)
  END DO;END DO;END DO
  ! Limit Solution U
  IF(Theta2.LT.1.) THEN
    DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
      U(:,i,j,k,iElem) = Theta2*(U(:,i,j,k,iElem)-UMean) + UMean
    END DO;END DO;END DO
  END IF
END DO !iElem
! collect statistics
PP_Elems_counter  = PP_Elems_counter  + PP_Elems
PP_Switch_counter = PP_Switch_counter + 1
PP_Elems_Amount   = REAL(PP_Elems_Counter)/PP_Switch_counter
END SUBROUTINE PPLimiter

!==================================================================================================================================
!> Computes t, such that t*U+(1-t)*Umean is admissible (i.e. has p>eps)
!==================================================================================================================================
PPURE FUNCTION CalcT(ULoc,UMean) RESULT (tLoc)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars      ,ONLY: KappaM1
USE MOD_Filter_Vars   ,ONLY: PPepsPres
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)         :: ULoc( PP_nVar)                       ! Unfiltered solution
REAL,INTENT(IN)         :: UMean(PP_nVar)                       ! Mean of the solution
REAL                    :: tLoc                                 ! Result
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: p
REAL                    :: a,b,c
REAL                    :: UDiff(PP_nVar)
!==================================================================================================================================
! Calculate new pressure
p = KappaM1*(ULoc(ENER)-0.5*DOT_PRODUCT(ULoc(MMV2),ULoc(MMV2))/ULoc(DENS))

! If pressure is greater than limit, return and do nothing
IF (p .GE. PPepsPres) THEN
  tLoc = 1.
  RETURN
END IF

! Otherwise calculate the necessary value for Theta2
UDiff=ULoc-UMean

a = UDiff(ENER)*UDiff(DENS)                                                             - 0.5*DOT_PRODUCT(UDiff(MMV2),UDiff(MMV2))
b = UMean(ENER)*UDiff(DENS) + UMean(DENS)*UDiff(ENER) - (PPepsPres/KappaM1)*UDiff(DENS)     - DOT_PRODUCT(UMean(MMV2),UDiff(MMV2))
c = UMean(ENER)*UMean(DENS)                           - (PPepsPres/KappaM1)*UMean(DENS) - 0.5*DOT_PRODUCT(UMean(MMV2),UMean(MMV2))

! f = at^2+bt+c; f(0)>0; f(1)<0
! => exactly one root in [0,1], with df/dt<0
! Check which solution for t has to be returned
IF(ABS(a).LT.1.E-12) THEN
  ! linear case
  tLoc = -c/b
ELSE
  ! non-linear cases:
  ! a<0: parabola has a MAXIMUM, root with df/dt<0 is the larger of the two
  ! a>0: parabola has a MINIMUM, root with df/dt<0 is the smaller of the two
  ! formula is the same for both cases due to sign of denominator
  tLoc = -0.5*(b + SQRT(b*b-4.*a*c))/a
END IF

! Sanity check
IF(ISNAN(tLoc).OR.(tLoc.GT.1.).OR.(tLoc.LT.0.)) tLoc=0.

END FUNCTION CalcT

!==================================================================================================================================
!> Print information on the amount of PP subcells
!==================================================================================================================================
SUBROUTINE PPLimiter_Info(iter)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars    ,ONLY: nGlobalElems
USE MOD_Analyze_Vars ,ONLY: totalPP_nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN) :: iter !< number of iterations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,totalPP_nElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  ! totalPP_nElems is counted in PrintStatusLine
ELSE
  CALL MPI_REDUCE(totalPP_nElems,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(A,F8.3,A,I0,A)')' PP Elems amount %: ',totalPP_nElems/REAL(nGlobalElems)/iter*100,', ',totalPP_nElems,' elems'
totalPP_nElems = 0
END SUBROUTINE PPLimiter_Info

!==================================================================================================================================
!> Deallocate filter arrays
!==================================================================================================================================
SUBROUTINE FinalizePPLimiter()
! MODULES
USE MOD_Filter_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(PP_Elems)
SDEALLOCATE(PP_Elems_counter)
SDEALLOCATE(PP_Elems_Amount)
END SUBROUTINE FinalizePPLimiter

#endif
END MODULE MOD_PPLimiter
