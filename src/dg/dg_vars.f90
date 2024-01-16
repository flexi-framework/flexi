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

!==================================================================================================================================
!> \brief Contains global variables provided by the DG modules.
!>
!> This module contains the building blocks for the DGSEM operator, i.e. the interpolation and differentiation
!> operators and the variables for solution storage as well as auxilliary variables.
!==================================================================================================================================
MODULE MOD_DG_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! DG basis, contains the differentiation and interpolation operators
REAL,ALLOCATABLE                      :: D(:,:)                 !< Differentiation matrix of size [0..N,0..N], contains the
                                                                !< first derivative of each Lagrange polynomial at each node.

REAL,ALLOCATABLE                      :: D_T(:,:)               !< Transpose of differentiation matrix, size [0..N,0..N].

REAL,ALLOCATABLE                      :: D_Hat(:,:)             !< Differentiation matrix premultiplied by
                                                                !< mass matrix, \f$ \hat{D} = M^{-1} D^T M \f$, size [0..N,0..N].

REAL,ALLOCATABLE                      :: D_Hat_T(:,:)           !< Transpose of differentiation matrix premultiplied by
                                                                !< mass matrix, size [0..N,0..N].

REAL,ALLOCATABLE                      :: L_HatMinus(:)          !< Lagrange polynomials evaluated at \f$\xi=-1\f$
                                                                !< premultiplied by mass matrix

REAL,ALLOCATABLE                      :: L_HatPlus(:)           !< Lagrange polynomials evaluated at \f$\xi=+1\f$
                                                                !< premultiplied by mass matrix

#ifdef SPLIT_DG
REAL,ALLOCATABLE                      :: DVolSurf(:,:)          !< Transpose of differentiation matrix used for calculating the strong form
                                                                !< DVolSurf = D_T
                                                                !< DVolSurf(0,0) = D_T(0,0) + 1/(2 * wGP(0))
                                                                !< DVolSurf(N,N) = D_T(N,N) - 1/(2 * wGP(N))
#endif /*SPLIT_DG*/

!----------------------------------------------------------------------------------------------------------------------------------
! DG solution (JU or U) vectors)
REAL,ALLOCATABLE,TARGET               :: U(:,:,:,:,:)           !< Solution variable for each equation, node and element,
                                                                !< size [1..NVar,0..N,0..N,0..N,nElems].
#ifdef PP_EntropyVars
REAL,ALLOCATABLE,TARGET               :: V   (:,:,:,:,:)        !< Entropy variables for each node and element,
                                                                !< size \([1..PP_nVar,0..N,0..N,0..N,nElems]\).
REAL,ALLOCATABLE                      :: V_master(:,:,:,:)      !< 2D Solution on face nodes for the master sides,
                                                                !< size \([1..nVar,0..N,0..N,all\_master\_sides]\)

REAL,ALLOCATABLE                      :: V_slave(:,:,:,:)       !< 2D Solution on face nodes for the slave sides,
#endif

!----------------------------------------------------------------------------------------------------------------------------------
! DG time derivative or Residual U_t
REAL,ALLOCATABLE                      :: Ut(:,:,:,:,:)          !< Residual/time derivative, size [1..NVar,0..N,0..N,0..NZ,nElems].

!----------------------------------------------------------------------------------------------------------------------------------
! auxilliary counters: number of entries in U, Ut, gradUx, gradUy, gradUz, used of optimization
INTEGER                               :: nTotalU                !< Total number of entries in U / size of U.
INTEGER                               :: nDOFFace               !< Degrees of freedom on single face(per equation)
                                                                !< $ nDOFace=(N+1)^2 $.
INTEGER                               :: nDOFElem               !< Degrees of freedom in single element(per equation)
                                                                !< $ nDOFElem=(N+1)^3 $.

!----------------------------------------------------------------------------------------------------------------------------------
! interior face values for all elements
REAL,ALLOCATABLE                      :: U_master(:,:,:,:)      !< 1D/2D Solution on face nodes for the master sides,
                                                                !< size [1..NVar,0..N,0..NZ,all_master_sides]

REAL,ALLOCATABLE                      :: U_slave(:,:,:,:)       !< 1D/2D Solution on face nodes for the slave sides,
                                                                !< size [1..NVar,0..N,0..NZ,all_slave_sides]

REAL,ALLOCATABLE                      :: Flux_master(:,:,:,:)   !< Fluxes on face, size [1..NVar,0..N,0..NZ,allsides].
REAL,ALLOCATABLE                      :: Flux_slave (:,:,:,:)   !< Fluxes on face, size [1..NVar,0..N,0..NZ,allsides].

!----------------------------------------------------------------------------------------------------------------------------------
! Variables in case of primitive lifting
REAL,ALLOCATABLE                      :: UPrim(:,:,:,:,:)       !< Solution in primitive variables per equation, node and element,
                                                                !< size [1..NVar,0..N,0..N,0..NZ,nElems].
REAL,ALLOCATABLE                      :: UPrim_master(:,:,:,:)  !< 2D Solution in Primitive variables on face, master side,
                                                                !< size [1..NVar,0..N,0..NZ,all_master_sides]
REAL,ALLOCATABLE                      :: UPrim_slave(:,:,:,:)   !< 2D Solution in Primitive variables on face, slave side,
                                                                !<size [1..NVar,0..N,0..NZ,all_slave_sides]
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                               :: SplitDG                !< Shows which split formulation is used
! Variables for boundary flux calculation
REAL,ALLOCATABLE                      :: UPrim_Boundary(:,:,:)

!----------------------------------------------------------------------------------------------------------------------------------
! Auxilliary variables
LOGICAL                               :: DGInitIsDone=.FALSE.   !< Switch to check DGInit status
!==================================================================================================================================
END MODULE MOD_DG_Vars
