!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
#if FV_ENABLED

!==================================================================================================================================
!> Initialize a lot of basic variables for the Finite Volume sub-cells. Especially variables in reference space
!==================================================================================================================================
MODULE MOD_FV_Basis
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER      :: FV_NODETYPE_SAME              =-1
INTEGER,PARAMETER      :: FV_NODETYPE_EQUIDISTANT       = 0
INTEGER,PARAMETER      :: FV_NODETYPE_LEGENDRE_GAUSS    = 1
INTEGER,PARAMETER      :: FV_NODETYPE_LEGENDRE_LOBATTO  = 2
INTEGER,PARAMETER      :: FV_NODETYPE_CHEBYSHEV_LOBATTO = 3

INTERFACE DefineParametersFV_Basis
  MODULE PROCEDURE DefineParametersFV_Basis
END INTERFACE

INTERFACE InitFV_Basis
  MODULE PROCEDURE InitFV_Basis
END INTERFACE

INTERFACE FV_Build_X_w_BdryX
  MODULE PROCEDURE FV_Build_X_w_BdryX
END INTERFACE

INTERFACE FV_Build_VisuVdm
  MODULE PROCEDURE FV_Build_VisuVdm
END INTERFACE

INTERFACE FV_Build_Vdm_Gauss_FVboundary
  MODULE PROCEDURE FV_Build_Vdm_Gauss_FVboundary
END INTERFACE

INTERFACE FV_GetVandermonde
  MODULE PROCEDURE FV_GetVandermonde
END INTERFACE

INTERFACE FinalizeFV_Basis
  MODULE PROCEDURE FinalizeFV_Basis
END INTERFACE


PUBLIC::DefineParametersFV_Basis
PUBLIC::InitFV_Basis
PUBLIC::FV_Build_X_w_BdryX
PUBLIC::FV_Build_VisuVdm
PUBLIC::FV_Build_Vdm_Gauss_FVboundary
PUBLIC::FV_GetVandermonde
PUBLIC::FinalizeFV_Basis
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FV basis
!==================================================================================================================================
SUBROUTINE DefineParametersFV_Basis()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('FV_Basis')
CALL prms%CreateIntFromStringOption('FV_CellType',"Specify type of FV cell distribution to be used:&
                                                                            & equidistant\n&
                                                                            & same\n&
                                                                            & legendre_gauss\n&
                                                                            & legendre_lobatto\n&
                                                                            & chebychev_lobatto"&
                                                                            , value='equidistant')
CALL addStrListEntry('FV_CellType','same',              FV_NODETYPE_SAME             )
CALL addStrListEntry('FV_CellType','equidistant',       FV_NODETYPE_EQUIDISTANT      )
CALL addStrListEntry('FV_CellType','legendre_gauss',    FV_NODETYPE_LEGENDRE_GAUSS   )
CALL addStrListEntry('FV_CellType','legendre_lobatto',  FV_NODETYPE_LEGENDRE_LOBATTO )
CALL addStrListEntry('FV_CellType','chebyshev_lobatto', FV_NODETYPE_CHEBYSHEV_LOBATTO)
END SUBROUTINE DefineParametersFV_Basis

!==================================================================================================================================
!> Initialize sub-cells width/points/distances/... Vandermondes to switch between DG and FV ...
!==================================================================================================================================
SUBROUTINE InitFV_Basis()
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone,NodeType
USE MOD_ReadInTools        ,ONLY: GETINTFROMSTR
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FV Basis...'
IF(InterpolationInitIsDone.AND.FVInitBasisIsDone)THEN
  CALL CollectiveStop(__STAMP__, &
    'InitFV_Basis not ready to be called or already called.')
END IF

#if PARABOLIC
#if !(FV_RECONSTRUCT)
CALL CollectiveStop(__STAMP__, &
  'FV_RECONSTRUCT=0 and PARABOLIC=T is not allowed. Switch off PARABOLIC or switch on FV_RECONSTRUCT!')
#endif
#endif

! The indicator value is used to decide where FV sub-cells are needed
! TODO: Force overwrite not so nice....
#if FV_ENABLED == 2 || FV_ENABLED == 3
  FV_CellType = FV_NODETYPE_SAME
#else
  FV_CellType = GETINTFROMSTR('FV_CellType')
#endif /*FV_BLENDING*/

! DG reference element [-1,1] is subdivided into (N+1) sub-cells
!
! -1                                       1
! |---------|---------|---------|----------|
! ^         ^         ^         ^          ^   FV_BdryX: boundary positions of subcells
!      ^         ^         ^         ^         FV_X    : positions of support points of subcells
!  <------->                                   FV_w    : width of subcell

! allocate arrays for precomputed stuff on reference-element
ALLOCATE(FV_BdryX(0:PP_N+1))  ! 1D boundary positions of FV-Subcells
ALLOCATE(FV_X(0:PP_N))        ! 1D positions of support points of FV-Subcells
ALLOCATE(FV_w(0:PP_N))        ! 1D width of FV-Subcells
ALLOCATE(FV_w_inv(0:PP_N))    ! Inverse of 1D width of FV-Subcells

! precompute stuff on reference-element
CALL FV_Build_X_w_BdryX(PP_N, FV_X, FV_w, FV_BdryX, FV_CellType)

! equidistant FV needs Vdm-Matrix to interpolate/reconstruct between FV- and DG-Solution
ALLOCATE(FV_Vdm(0:PP_N,0:PP_N))    ! used for DG -> FV
ALLOCATE(FV_sVdm(0:PP_N,0:PP_N))   ! used for FV -> DG
! Build Vandermondes for switch between DG and FV
CALL FV_GetVandermonde(PP_N,NodeType,FV_Vdm,FV_sVdm)

! calculate inverse FV-widths (little speedup)
FV_w_inv      = 1.0 / FV_w

FVInitBasisIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT FV DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitFV_Basis


!==================================================================================================================================
!> Build Vandermondes to switch solution between DG and FV sub-cells.
!==================================================================================================================================
SUBROUTINE FV_GetVandermonde(N_in,NodeType_in,FV_Vdm,FV_sVdm)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Interpolation ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_Basis         ,ONLY: InitializeVandermonde
USE MOD_Mathtools     ,ONLY: INVERSE
USE MOD_FV_Vars       ,ONLY: FV_CellType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: N_in                    !< Number of 1D input points / output points
CHARACTER(LEN=255),INTENT(IN) :: NodeType_in             !< Type of 1D input points
REAL,INTENT(OUT)              :: FV_Vdm(0:N_in,0:N_in)   !< Vandermonde matrix for converstion from DG to FV
REAL,INTENT(OUT),OPTIONAL     :: FV_sVdm(0:N_in,0:N_in)  !< Vandermonde matrix for converstion from FV to DG
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                   :: FV_X(0:N_in),FV_w(0:N_in),FV_BdryX(0:N_In+1)
REAL,DIMENSION(0:N_In) :: xGP,wGP,wBary
REAL                   :: SubxGP(0:N_In)
REAL                   :: VDM(0:N_In,0:N_In)
INTEGER                :: i,j,k
!==================================================================================================================================
CALL GetNodesAndWeights(N_in,NodeType_in,xGP,wGP,wBary)

! one DG cell [-1,1] is divided in FV-Subcells
!
! -1                  i-th Subcell                              1
!  |------- ... ----|---x---x---x---|---- ... -----|------------|
!                     x = xGP in Subcell
!
!  xFV_i(k) = k-th Gauss point in i-th Subcell
!
! We have to integrate the solution U in each Subcell (loop over i) and compute with this the integral mean value, which
! then is the Finite-Volume solution U_FV of this Subcell.
! Therefore we must compute the solution in each Gauss point xFV_i of the Subcell, which is done by evaluating all Lagrange
! polynomials l_j in all Gauss points xFV_i => stored in the matrix VDM_i = [ l_j(xFV_i(k)) ]_kj
! Multiplying this Vandermonde with the DG-solution U gives us the solution in the Gauss points xFV_i of the i-th Subcell:
!    u_i = VDM_i . U         meaning u_i(k) = U(xFV_i(k))
! The integration over the Subcell is done by Gauss-Quadrature:
!    \int_{Subcell(i)} U dx = [ \sum_{k=0..N} u_i(k) . wGP(k) ] * a       with a = (width of Subcell(i)) / 2
! ( The /2 comes from the width of the reference element. )
! We get the integral mean value and therewith the Finite-Volume solution of the i-th Subcell  by dividing this integral by
! the width of the i-th Subcell :
!    u_FV(i) = [ \sum_{k=0..N} wGP(k) * u_i(k) ] / 2
!
! All this can be write in Matrix-Vector-Multiplication:
!    u_FV(i) = 1/2 * wGP^T . VDM_i . U
! If we combine the vectors (wGP^T . VDM_i) for the i-th Subcell for all Subcells in one matrix we get :
!    u_FV = 1/2 * FV_Vdm . U          with FV_Vdm = ( wGP^T.VDM_0 // wGP^T.VDM_1 // ... // wGP^T. VDM_N )
! With the inverse of the matrix FV_Vdm we then can reconstruct the DG-Solution from a FV-Solution
!
! Algorithm :
FV_Vdm = 0.0
CALL FV_Build_X_w_BdryX(N_in, FV_X, FV_w, FV_BdryX, FV_CellType)
!FV_BdryX = (FV_BdryX + 1.)/2. -1.
DO i=0,N_in
  ! 1. Compute the Gauss points xFV_i in the i-th Subcell
  SubxGP = FV_BdryX(i) + (xGP + 1.)/2. * (FV_BdryX(i+1) - FV_BdryX(i))
  ! 2. Evaluate the all Lagrange-Polys in all Gauss points xFV_i of the i-th Subcell  =>  store in VDM
  CALL InitializeVandermonde(N_in,N_in,wBary,xGP,SubxGP,VDM(:,:))
  ! 3. Multiply wGP^T with VDM and store it in the i-th row of the matrix FV_Vdm
  DO j=0,N_in
    DO k=0,N_in
      FV_Vdm(i,j) = FV_Vdm(i,j) + wGP(k) * VDM(k,j)
    END DO
  END DO
END DO
! 4. don't forget the 1/2
FV_Vdm = FV_Vdm * 0.5

! Compute the inverse of FV_Vdm
IF (PRESENT(FV_sVdm)) THEN
  FV_sVdm = INVERSE(FV_Vdm)
END IF
END SUBROUTINE FV_GetVandermonde

!==================================================================================================================================
!> Build positions FV_X, widths, and boundary positions
!==================================================================================================================================
SUBROUTINE FV_Build_X_w_BdryX(N, FV_X, FV_w, FV_BdryX, FV_CellType)
! MODULES
USE MOD_Globals
USE MOD_Basis
USE MOD_Interpolation      ,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars ,ONLY: NodeType, NodeTypeG, NodeTypeGL, NodeTypeCL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N               !< polynomial degree of DG elements / number of sub-cells per direction (N+1)
REAL,INTENT(OUT)   :: FV_X(0:N)       !< cell-centers of the sub-cells in reference space
REAL,INTENT(OUT)   :: FV_w(0:N)       !< width of the sub-cells in reference space
REAL,INTENT(OUT)   :: FV_BdryX(0:N+1) !< positions of the boundaries of the sub-cells in reference space
INTEGER,INTENT(IN) :: FV_CellType
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!==================================================================================================================================
! NOTE: calc different node distributions for fv subcell: default is the same as DG
SELECT CASE (FV_CellType)
CASE (FV_NODETYPE_EQUIDISTANT)
  FV_w(:)  = 2.0 / (N+1)
CASE (FV_NODETYPE_LEGENDRE_GAUSS)
  CALL GetNodesAndWeights(N,NodeTypeG,FV_X,FV_w)
CASE (FV_NODETYPE_LEGENDRE_LOBATTO)
  CALL GetNodesAndWeights(N,NodeTypeGL,FV_X,FV_w)
CASE (FV_NODETYPE_CHEBYSHEV_LOBATTO)
  CALL GetNodesAndWeights(N,NodeTypeCL,FV_X,FV_w)
CASE(FV_NODETYPE_SAME)
  ! Set FV positions FV_X as the GLL/G nodes, and the FV cell width to be the GLL/G weights.
  CALL GetNodesAndWeights(N,NodeType,FV_X,FV_w)
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,"Nodetype for FV subcells unknown!")
END SELECT

! calculate boundaries and nodes (midpoints) of FV-Subcells
FV_BdryX(0) = -1.0
DO i=1,N+1
  FV_BdryX(i) = FV_BdryX(i-1) + FV_w(i-1)
  ! The following is for FVx at midpoint of cell.
  FV_X(i-1)   = (FV_BdryX(i-1) + FV_BdryX(i))/2.0
END DO
END SUBROUTINE FV_Build_X_w_BdryX

!==================================================================================================================================
!> Build Vandermonde to convert solution from Gauss points (poly degree N) to left and right face of each subcells
!> (inner element faces are doubled).
!==================================================================================================================================
SUBROUTINE FV_Build_VisuVdm(N, Vdm)
! MODULES
USE MOD_Basis
USE MOD_FV_Vars       ,ONLY: FV_CellType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N                    !< polynomial degree of DG elements / number of sub-cells per direction (N+1)
REAL,INTENT(OUT)   :: Vdm(0:(N+1)*2-1,0:N) !< Vandermonde matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:PP_N) :: FV_X,xGP,wBary
REAL                   :: X  (0:(N+1)*2-1)
REAL                   :: FV_BdryX(0:N+1),FV_w(0:N)
INTEGER                :: i,l,k
!==================================================================================================================================
#if (PP_NodeType == 1)
  CALL LegendreGaussNodesAndWeights(N, xGP, wBary)
#elif (PP_NodeType == 2)
  CALL LegGaussLobNodesAndWeights(N, xGP, wBary)
#endif
CALL BarycentricWeights(N, xGP, wBary)

CALL FV_Build_X_w_BdryX(N, FV_X, FV_w, FV_BdryX, FV_CellType)
k = 0
DO i=0,N
  DO l=0,1
    X(k) = FV_BdryX(i) + l*(FV_BdryX(i+1) - FV_BdryX(i))
    k = k+1
  END DO
END DO

CALL InitializeVandermonde(N,(N+1)*2-1,wBary,xGP,X,Vdm)
END SUBROUTINE FV_Build_VisuVdm

!==================================================================================================================================
!> Build Vandermonde to convert solution from Gauss points (poly degree N) to FV_BdryX points.
!==================================================================================================================================
SUBROUTINE FV_Build_Vdm_Gauss_FVboundary(N, Vdm)
! MODULES
USE MOD_Basis
USE MOD_FV_Vars       ,ONLY: FV_CellType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N              !< polynomial degree of DG elements / number of sub-cells per direction (N+1)
REAL,INTENT(OUT)   :: Vdm(0:N+1,0:N) !< Vandermonde matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:PP_N) :: FV_X,xGP,wBary
REAL                   :: FV_BdryX(0:N+1),FV_w(0:N)
!==================================================================================================================================
#if (PP_NodeType == 1)
  CALL LegendreGaussNodesAndWeights(N, xGP, wBary)
#elif (PP_NodeType == 2)
  CALL LegGaussLobNodesAndWeights(N, xGP, wBary)
#endif
CALL BarycentricWeights(N, xGP, wBary)
CALL FV_Build_X_w_BdryX(N, FV_X, FV_w, FV_BdryX, FV_CellType)
CALL InitializeVandermonde(N,N+1,wBary,xGP,FV_BdryX,Vdm)
END SUBROUTINE FV_Build_Vdm_Gauss_FVboundary

!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeFV_Basis()
USE MOD_FV_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(FV_BdryX)
SDEALLOCATE(FV_X)
SDEALLOCATE(FV_w)
SDEALLOCATE(FV_w_inv)
SDEALLOCATE(FV_Vdm)
SDEALLOCATE(FV_sVdm)
FVInitBasisIsDone=.FALSE.
END SUBROUTINE FinalizeFV_Basis


END MODULE MOD_FV_Basis
#endif /* FV_ENABLED */
