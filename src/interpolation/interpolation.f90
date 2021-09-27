!==================================================================================================================================
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
!==================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
!> Contains routines to prepare for interpolation procedures:
!> - Initialize interpolation variables.
!> - Calculate node positions and weights.
!> - Build Vandermonde matrices.
!> - Build derivative matrices.
!> Also contains routines to map the solution between physical and reference space.
!==================================================================================================================================
MODULE MOD_Interpolation
! MODULES
USE MOD_Basis
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part
! ----------------------------------------------------------------------------------------------------------------------------------
! Public Part
! ----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitInterpolation
   MODULE PROCEDURE InitInterpolation
END INTERFACE

INTERFACE InitInterpolationBasis
   MODULE PROCEDURE InitInterpolationBasis
END INTERFACE

INTERFACE GetNodesAndWeights
   MODULE PROCEDURE GetNodesAndWeights
END INTERFACE

INTERFACE GetVandermonde
   MODULE PROCEDURE GetVandermonde
END INTERFACE

INTERFACE GetDerivativeMatrix
   MODULE PROCEDURE GetDerivativeMatrix
END INTERFACE

INTERFACE FinalizeInterpolation
   MODULE PROCEDURE FinalizeInterpolation
END INTERFACE

PUBLIC::InitInterpolation
PUBLIC::InitInterpolationBasis
PUBLIC::GetNodesAndWeights
PUBLIC::GetVandermonde
PUBLIC::GetDerivativeMatrix
PUBLIC::FinalizeInterpolation
!==================================================================================================================================


PUBLIC::DefineParametersInterpolation

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersInterpolation()
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
CALL prms%SetSection("Interpolation")
CALL prms%CreateIntOption('N'    , "Polynomial degree of computation to represent to solution")
END SUBROUTINE DefineParametersInterpolation


!=================================================================================================================================
!> Initialize interpolation. Call basis initialization and calculate Vandermonde matrices between NUnder/N/Over.
!=================================================================================================================================
SUBROUTINE InitInterpolation(NIn)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars
USE MOD_ReadInTools,        ONLY:GETINT,CountOption
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: NIn  !< optional polynomial degree
!----------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
#if PP_N != N
INTEGER :: Ntmp
#endif
!==================================================================================================================================
IF (InterpolationInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'InitInterpolation already called.')
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT INTERPOLATION...'

! Access ini-file
#if PP_N == N
IF(PRESENT(Nin))THEN
  PP_N = NIn
ELSE
  PP_N = GETINT('N')
END IF
#else
IF(PRESENT(Nin))THEN
  Ntmp = NIn
ELSE
  Ntmp=PP_N
  IF(CountOption('N').EQ.1) Ntmp=GETINT('N')
END IF
IF(PP_N.NE.Ntmp) THEN
  CALL CollectiveStop(__STAMP__,&
  'N in ini-file is different from hard-compiled N in Flexi. Ini/Compiled:',Ntmp,REAL(PP_N))
END IF
#endif

! Compute Nodes and weights for Gauss or GaussLobatto-Nodes
SWRITE(UNIT_stdOut,'(A)') ' NodeType: '//NodeType
CALL InitInterpolationBasis(PP_N, xGP ,wGP ,wBary ,L_Minus ,L_Plus ,Vdm_Leg ,sVdm_Leg)

InterpolationInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT INTERPOLATION DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitInterpolation



!==================================================================================================================================
!> Initialize basis for Gauss-points of order N_in.
!> Calculate positions of Gauss-points, integration weights and barycentric weights. Prepare basis evaluation at -1 and +1.
!> Calculate Vandermonde Nodal <--> Modal.
!==================================================================================================================================
SUBROUTINE InitInterpolationBasis(N_in,xGP,wGP,wBary,L_Minus,L_Plus,Vdm_Leg,sVdm_Leg)
! MODULES
USE MOD_Interpolation_Vars,      ONLY:NodeType
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys,buildLegendreVdm
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_in               !< Polynomial degree
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT):: xGP                !< Gauss points in [-1,1]
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT):: wGP                !< Integration weights
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT):: wBary              !< Barycentric weights
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT):: L_Minus            !< Lagrange polynomials at -1
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT):: L_Plus             !< Lagrange polynomials at +1
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT):: Vdm_Leg            !< Vandermonde Modal->Nodal
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT):: sVdm_Leg           !< Vandermonde Nodal->Modal
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Allocate global variables, needs to go somewhere else later
ALLOCATE(xGP(0:N_in), wGP(0:N_in), wBary(0:N_in))
ALLOCATE(L_Minus(0:N_in), L_Plus(0:N_in))

CALL GetNodesAndWeights(N_in,NodeType,xGP,wGP,wBary)

!! interpolate to left and right face (1 and -1)
CALL LagrangeInterpolationPolys(1.,N_in,xGP,wBary,L_Plus)
CALL LagrangeInterpolationPolys(-1.,N_in,xGP,wBary,L_Minus)

! Vandermonde from NODAL <--> MODAL
ALLOCATE(Vdm_Leg(0:N_in,0:N_in),sVdm_Leg(0:N_in,0:N_in))
CALL buildLegendreVdm(N_in,xGP,Vdm_Leg,sVdm_Leg)

END SUBROUTINE InitInterpolationBasis


!==================================================================================================================================
!> Compute 1D nodes and weights for several node types in interval [-1,1]
!==================================================================================================================================
SUBROUTINE GetNodesAndWeights(N_in,NodeType_in,xIP,wIP,wIPBary)
! MODULES
USE MOD_Globals
USE MOD_StringTools, ONLY: LowCase
USE MOD_Basis,       ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,ChebyGaussLobNodesAndWeights,BarycentricWeights
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in            !< Polynomial degree
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_in     !< Type of 1D points
REAL,INTENT(OUT)                   :: xIP(0:N_in)     !< Position of nodes
REAL,INTENT(OUT),OPTIONAL          :: wIP(0:N_in)     !< Integration weights
REAL,INTENT(OUT),OPTIONAL          :: wIPBary(0:N_in) !< Barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i
REAL                               :: x(0:(N_in+1)/2-1)
REAL                               :: w(0:(N_in+1)/2-1)
INTEGER                            :: NN
CHARACTER(LEN=255)                 :: NodeType_loc    !< Type of 1D points
!==================================================================================================================================
CALL LowCase(NodeType_in,NodeType_loc)
IF(PRESENT(wIP))THEN
  SELECT CASE(TRIM(NodeType_loc))
  CASE('gauss')
    CALL LegendreGaussNodesAndWeights(N_in,xIP,wIP)
  CASE('gauss-lobatto')
    CALL LegGaussLobNodesAndWeights(N_in,xIP,wIP)
  CASE('chebyshev-gauss-lobatto')
    CALL ChebyGaussLobNodesAndWeights(N_in,xIP,wIP)
  CASE('visu')
    DO i=0,N_in
      xIP(i) = 2.*REAL(i)/REAL(N_in) - 1.
    END DO
    ! Trapez rule for integration !!!
    wIP(:) = 2./REAL(N_in)
    wIP(0) = 0.5*wIP(0)
    wIP(N_in) = 0.5*wIP(N_in)
  CASE('visu_inner')
    DO i=0,N_in
      xIP(i) = 1./REAL(N_in+1)+2.*REAL(i)/REAL(N_in+1) - 1.
    END DO
    ! first order intergration !!!
    wIP=2./REAL(N_in+1)
  CASE('fv_gauss')
    CALL Abort(__STAMP__,&
        'Not implemented!')
  CASE('fv_gausslob')
    CALL Abort(__STAMP__,&
        'Not implemented!')
  CASE('visu_fvequi')
    CALL Abort(__STAMP__,&
        'Not implemented!')
  CASE DEFAULT
    CALL Abort(__STAMP__,&
      'NodeType "'//TRIM(NodeType_loc)//'" in GetNodesAndWeights not found!')
  END SELECT
ELSE
  SELECT CASE(TRIM(NodeType_loc))
  CASE('gauss')
    CALL LegendreGaussNodesAndWeights(N_in,xIP)
  CASE('gauss-lobatto')
    CALL LegGaussLobNodesAndWeights(N_in,xIP)
  CASE('chebyshev-gauss-lobatto')
    CALL ChebyGaussLobNodesAndWeights(N_in,xIP)
  CASE('visu')
    DO i=0,N_in
      xIP(i) = 2.*REAL(i)/REAL(N_in) - 1.
    END DO
  CASE('visu_inner')
    DO i=0,N_in
      xIP(i) = 1./REAL(N_in+1)+2.*REAL(i)/REAL(N_in+1) - 1.
    END DO
  CASE('fv_gauss')
    IF (MOD(N_in,2).EQ.0) THEN
      CALL Abort(__STAMP__,&
        'NodeType "'//TRIM(NodeType_in)//'" needs N to be odd!')
    END IF
    NN=(N_in+1)/2 -1
    CALL LegendreGaussNodesAndWeights(NN,x,w)

    DO i=0,NN
      xIP(i*2  ) = -1.+SUM(w(0:i-1))
      xIP(i*2+1) = -1.+SUM(w(0:i  ))
    END DO
  CASE('fv_gausslob')
    IF (MOD(N_in,2).EQ.0) THEN
      CALL Abort(__STAMP__,&
        'NodeType "'//TRIM(NodeType_in)//'" needs N to be odd!')
    END IF
    NN=(N_in+1)/2 -1
    CALL LegGaussLobNodesAndWeights(NN,x,w)

    DO i=0,NN
      xIP(i*2  ) = -1.+SUM(w(0:i-1))
      xIP(i*2+1) = -1.+SUM(w(0:i  ))
    END DO
  CASE('visu_fvequi')
    IF (MOD(N_in,2).EQ.0) THEN
      CALL Abort(__STAMP__,&
        'NodeType "'//TRIM(NodeType_in)//'" needs N to be odd!')
    END IF
    NN=(N_in-1)/2
    DO i=0,NN
      xIP(i*2  ) = -1.+ i   *(2./(NN+1))
      xIP(i*2+1) = -1.+(i+1)*(2./(NN+1))
    END DO

  CASE DEFAULT
    CALL Abort(__STAMP__,&
      'NodeType "'//TRIM(NodeType_in)//'" in GetNodesAndWeights not found!')
  END SELECT
END IF !present wIP
IF(PRESENT(wIPBary)) CALL BarycentricWeights(N_in,xIP,wIPBary)
END SUBROUTINE GetNodesAndWeights


!==================================================================================================================================
!> Build a Vandermonde-Matrix from/to different node types and polynomial degrees.
!==================================================================================================================================
SUBROUTINE GetVandermonde(N_in,NodeType_in,N_out,NodeType_out,Vdm_In_Out,Vdm_Out_In,modal)
! MODULES
USE MOD_Preproc
USE MOD_StringTools,       ONLY:STRICMP
USE MOD_Basis,             ONLY:BarycentricWeights,InitializeVandermonde
USE MOD_Interpolation_Vars,ONLY:NodeType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in                       !< Input polynomial degree
INTEGER,INTENT(IN)                 :: N_out                      !< Output polynomial degree
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_in                !< Type of 1D input points
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_out               !< Type of 1D output points
LOGICAL,INTENT(IN),OPTIONAL        :: modal                      !< Switch if a modal Vandermonde should be built
REAL,INTENT(OUT)                   :: Vdm_In_out(0:N_out,0:N_in) !< Vandermonde In->Out
REAL,INTENT(OUT),OPTIONAL          :: Vdm_Out_In(0:N_in,0:N_out) !< Vandermonde Out->in
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i
REAL                               :: xIP_in(0:N_in)
REAL                               :: xIP_out(0:N_out)
REAL                               :: wBary_in(0:N_in)
REAL                               :: wBary_out(0:N_out)
REAL                               ::  Vdm_Leg_in( 0:N_in,0:N_in)
REAL                               :: sVdm_Leg_in( 0:N_in,0:N_in)
REAL                               ::  Vdm_Leg_out(0:N_out,0:N_out)
REAL                               :: sVdm_Leg_out(0:N_out,0:N_out)
LOGICAL                            :: modalLoc
!==================================================================================================================================
modalLoc=.FALSE.
IF(PRESENT(modal)) modalLoc=modal

! Check if change Basis is needed

IF(STRICMP(NodeType_out,NodeType_in).AND.(N_in.EQ.N_out))THEN
  Vdm_In_Out=0.
  DO i=0,N_in
    Vdm_In_out(i,i)=1.
  END DO
  IF(PRESENT(Vdm_Out_In))THEN
    Vdm_Out_In=0.
    DO i=0,N_Out
      Vdm_Out_In(i,i)=1.
    END DO
  END IF
ELSE
  ! Input points
  CALL GetNodesAndWeights(N_in,NodeType_in,xIP_in)
  CALL BarycentricWeights(N_in,xIP_in,wBary_in)
  ! Output points
  CALL GetNodesAndWeights(N_out,NodeType_out,xIP_out)

  IF(modalLoc)THEN
    CALL buildLegendreVdm(N_In, xIP_in, Vdm_Leg_in, sVdm_Leg_in)
    CALL buildLegendreVdm(N_Out,xIP_out,Vdm_Leg_out,sVdm_Leg_out)
  END IF

  IF((N_Out.LT.N_In).AND.modalLoc)THEN
    Vdm_In_Out=MATMUL(Vdm_Leg_Out(0:N_Out,0:N_Out),sVdm_Leg_In(0:N_Out,0:N_In))
  ELSE
    CALL InitializeVandermonde(N_in,N_out,wBary_in,xIP_in,xIP_out,Vdm_In_Out)
  END IF
  IF(PRESENT(Vdm_Out_In))THEN
    IF((N_In.LT.N_Out).AND.modalLoc)THEN
      Vdm_Out_In=MATMUL(Vdm_Leg_In(0:N_In,0:N_In),sVdm_Leg_Out(0:N_In,0:N_Out))
    ELSE
      CALL BarycentricWeights(N_out,xIP_out,wBary_out)
      CALL InitializeVandermonde(N_out,N_in,wBary_out,xIP_out,xIP_in,Vdm_Out_In)
    END IF
  END IF
END IF
END SUBROUTINE GetVandermonde


!==================================================================================================================================
!> Compute polynomial derivative matrix. D(i,j) = Derivative of basis function j evaluated at node point i.
!==================================================================================================================================
SUBROUTINE GetDerivativeMatrix(N_in,NodeType_in,D)
! MODULES
USE MOD_Basis,             ONLY:PolynomialDerivativeMatrix
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in                       !< Polynomial degree
CHARACTER(LEN=255),INTENT(IN)      :: NodeType_in                !< Type of 1D input points
REAL,INTENT(OUT)                   :: D(0:N_in,0:N_in)           !< Derivative matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: xIP(0:N_in)
!==================================================================================================================================
CALL GetNodesAndWeights(N_in,NodeType_in,xIP)
CALL PolynomialDerivativeMatrix(N_in,xIP,D)
END SUBROUTINE GetDerivativeMatrix


!============================================================================================================================
!> Deallocate all global interpolation variables
!============================================================================================================================
SUBROUTINE FinalizeInterpolation()
! MODULES
USE MOD_Interpolation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
! Deallocate global variables
SDEALLOCATE(xGP)
SDEALLOCATE(wGP)
SDEALLOCATE(wBary)
SDEALLOCATE(L_Minus)
SDEALLOCATE(L_Plus)
SDEALLOCATE(Vdm_Leg)
SDEALLOCATE(sVdm_Leg)

InterpolationInitIsDone = .FALSE.
END SUBROUTINE FinalizeInterpolation

END MODULE MOD_Interpolation
