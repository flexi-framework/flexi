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
!> \brief Computes the DGSEM spatial operator and updates residual Ut

!> Contains the routines to
!> - initialize and finalize the DG global variables and the DG basis
!> - compute the DG spatial operators/residuals(Ut) using U from the volume, surface and source contribution, incl.
!> lifting for the gradients and parallelization
!==================================================================================================================================
MODULE MOD_DG
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: InitDG
PUBLIC:: DGTimeDerivative_weakForm
PUBLIC:: FinalizeDG
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Allocate all global DG variables like U (solution in volume), U_slave/U_master (solution on faces), Flux, Ut (DG time derivative),
!> also fill the initial solution and call init DG basis. Operator building are also initialized by calling InitDGBasis.
!==================================================================================================================================
SUBROUTINE InitDG()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars
USE MOD_Interpolation_Vars,   ONLY: xGP,wGP,L_minus,L_plus
USE MOD_Interpolation_Vars,   ONLY: InterpolationInitIsDone
USE MOD_Restart_Vars,         ONLY: DoRestart,RestartInitIsDone
USE MOD_Mesh_Vars,            ONLY: nElems,nSides,Elem_xGP,MeshInitIsDone
USE MOD_ChangeBasisByDim,     ONLY: ChangeBasisVolume
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Check if all the necessary initialization is done before
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.DGInitIsDone) &
  CALL CollectiveStop(__STAMP__,'InitDG not ready to be called or already called.')
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

! Pre-compute the dg operator building blocks (differentiation matrices and prolongation operators)
CALL InitDGBasis(PP_N, xGP,wGP,L_minus,L_plus,D ,D_T ,D_Hat ,D_Hat_T ,L_HatMinus ,L_HatPlus)

! Allocate the local DG solution (JU or U): element-based
ALLOCATE(U(        PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
! Allocate the time derivative / solution update /residual vector dU/dt: element-based
ALLOCATE(Ut(       PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
U=0.
Ut=0.

#if PP_EntropyVars==1
ALLOCATE(V   (PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
V=0.
ALLOCATE(V_master(PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(V_slave( PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
V_master=0.
V_slave=0.
#endif /*if PP_EntropyVars==1*/

! Allocate the 2D solution vectors on the sides, one array for the data belonging to the proc (the master)
! and one for the sides which belong to another proc (slaves): side-based
ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
U_master=0.
U_slave=0.

! Repeat the U, U_Minus, U_Plus structure for the primitive quantities
ALLOCATE(UPrim(       PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
UPrim=0.
UPrim_master=0.
UPrim_slave=0.

! Allocate the UPrim_boundary for the boundary fluxes
ALLOCATE(UPrim_boundary(PP_nVarPrim,0:PP_N,0:PP_NZ))

! Allocate two fluxes per side (necessary for coupling of FV and DG)
ALLOCATE(Flux_master(PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(Flux_slave (PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
Flux_master=0.
Flux_slave=0.

! variables for performance tricks
nDOFFace=(PP_N+1)**(PP_dim-1)
nDOFElem=(PP_N+1)**PP_dim
nTotalU=PP_nVar*nDOFElem*nElems

! Fill the solution vector U with the initial solution by interpolation, if not filled through restart already
IF(.NOT.DoRestart)THEN
  CALL FillIni(PP_N,Elem_xGP,U)
END IF

DGInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitDG


!==================================================================================================================================
!> Allocate and initialize the building blocks for the DG operator: Differentiation matrices and prolongation operators
!==================================================================================================================================
SUBROUTINE InitDGBasis(N_in,xGP,wGP,L_Minus,L_Plus,D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatPlus)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Interpolation,    ONLY: GetNodesAndWeights
USE MOD_Basis,            ONLY: PolynomialDerivativeMatrix,LagrangeInterpolationPolys,PolynomialMassMatrix
#ifdef SPLIT_DG
USE MOD_DG_Vars,          ONLY: DVolSurf ! Transpose of differentiation matrix used for calculating the strong form
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                             :: N_in                   !< Polynomial degree
REAL,DIMENSION(0:N_in),INTENT(IN)              :: xGP                    !< Gauss/Gauss-Lobatto Nodes
REAL,DIMENSION(0:N_in),INTENT(IN)              :: wGP                    !< Gauss/Gauss-Lobatto Weights
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Minus                !< Values of lagrange polynomials at \f$ \xi = -1 \f$
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Plus                 !< Values of lagrange polynomials at \f$ \xi = +1 \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D                      !< Differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_T                    !< Transpose of differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat                  !< Differentiation matrix premultiplied by mass matrix,
                                                                         !< \f$ \hat{D} = M^{-1} D^T M \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat_T                !< Transpose of D_Hat matrix \f$ \hat{D}^T \f$
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatMinus             !< Values of lagrange polynomials at \f$ \xi = -1 \f$
                                                                         !< premultiplied with mass matrix
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatPlus              !< Values of lagrange polynomials at \f$ \xi = +1 \f$
                                                                         !< premultiplied with mass matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_in,0:N_in)              :: M,Minv
#if ((PP_NodeType==1) && defined(SPLIT_DG))
REAL,DIMENSION(2,0:N_in)                   :: Vf
REAL,DIMENSION(2,2)                        :: B
#endif /*((PP_NodeType==1) && defined(SPLIT_DG))*/
!==================================================================================================================================

ALLOCATE(L_HatMinus(0:N_in), L_HatPlus(0:N_in))
ALLOCATE(D(    0:N_in,0:N_in), D_T(    0:N_in,0:N_in))
ALLOCATE(D_Hat(0:N_in,0:N_in), D_Hat_T(0:N_in,0:N_in))
! Compute Differentiation matrix D for given Gausspoints
CALL PolynomialDerivativeMatrix(N_in,xGP,D)
D_T=TRANSPOSE(D)

! Build Mass Matrix
CALL PolynomialMassMatrix(N_in,xGP,wGP,M,Minv)

! Build D_Hat matrix. D^ = - (M^(-1) * D^T * M)
D_Hat  = -MATMUL(Minv,MATMUL(TRANSPOSE(D),M))
D_Hat_T= TRANSPOSE(D_Hat)

#ifdef SPLIT_DG
! Use a modified D matrix for the strong form volume integral, that incorporates the inner fluxes that are subtracted from the
! surfaces
ALLOCATE(DVolSurf(0:N_in,0:N_in))
#if (PP_NodeType==1)
! Chan, J.: Efficient Entropy Stable Gauss Collocation Methods, JSC, 2019.
! S = D - 0.5 V^T B V
! Lagrange polynomials evaluated at the cell boundaries
Vf(1,:) = L_Minus
Vf(2,:) = L_Plus
B(1,:)  = (/-1.,0./)
B(2,:)  = (/ 0.,1./)
DVolSurf = D - 1./2.*MATMUL(Minv,MATMUL(MATMUL(TRANSPOSE(Vf),B),Vf))
! Transpose for efficiency
DVolSurf = TRANSPOSE(DVolSurf)
#else
DVolSurf = D_T
! Modify the D matrix here, the integral over the inner fluxes at the boundaries will then be automatically done in the volume
! integral. The factor 1/2 is needed since we incorporate a factor of 2 in the split fluxes themselves!
! For Gauss-Lobatto points, these inner flux contributions cancel exactly with entries in the DVolSurf matrix, resulting in zeros.
DVolSurf(   0,   0) = DVolSurf(   0   ,0) + 1.0/(2.0 * wGP(   0))  ! = 0. (for LGL)
DVolSurf(N_in,N_in) = DVolSurf(N_in,N_in) - 1.0/(2.0 * wGP(N_in))  ! = 0. (for LGL)
#endif
#endif /*SPLIT_DG*/

! interpolate to left and right face (1 and -1 in reference space) and pre-divide by mass matrix
L_HatPlus  = MATMUL(Minv,L_Plus)
L_HatMinus = MATMUL(Minv,L_Minus)

END SUBROUTINE InitDGBasis


!==================================================================================================================================
!> \brief Computes the residual Ut = \f$ \frac {d\vec{U}} {dt} \f$ from the current solution U employing the DG method.
!> Computes the weak DGSEM space operator from surface, volume and source contributions. To do this we need to:
!> - Prolong the solution from the volume to the interface
!> - Invoke the lifting operator to calculate the gradients
!> - Perform the volume integral
!> - Perform the surface integral
!> - If needed, add source terms to the residual
!==================================================================================================================================
SUBROUTINE DGTimeDerivative_weakForm(t)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_DG_Vars             ,ONLY: Ut,U,U_slave,U_master,Flux_master,Flux_slave,L_HatPlus,L_HatMinus
USE MOD_DG_Vars             ,ONLY: UPrim,UPrim_master,UPrim_slave
!USE MOD_DG_Vars,             ONLY: nTotalU
USE MOD_VolInt
USE MOD_SurfIntCons         ,ONLY: SurfIntCons
USE MOD_ProlongToFaceCons   ,ONLY: ProlongToFaceCons
USE MOD_FillFlux            ,ONLY: FillFlux
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
USE MOD_Interpolation_Vars  ,ONLY: L_Minus,L_Plus
USE MOD_Overintegration_Vars,ONLY: OverintegrationType
USE MOD_Overintegration,     ONLY: Overintegration
USE MOD_ChangeBasisByDim    ,ONLY: ChangeBasisVolume
USE MOD_TestCase            ,ONLY: TestcaseSource
USE MOD_TestCase_Vars       ,ONLY: doTCSource
USE MOD_Equation            ,ONLY: GetPrimitiveStateSurface,GetConservativeStateSurface
USE MOD_EOS                 ,ONLY: ConsToPrim
USE MOD_Exactfunc           ,ONLY: CalcSource
USE MOD_Equation_Vars       ,ONLY: doCalcSource
USE MOD_Sponge              ,ONLY: Sponge
USE MOD_Sponge_Vars         ,ONLY: doSponge
USE MOD_Filter              ,ONLY: Filter_Pointer
USE MOD_Filter_Vars         ,ONLY: FilterType,FilterMat
USE MOD_FillMortarCons      ,ONLY: U_MortarCons,Flux_MortarCons
USE MOD_FillMortarPrim      ,ONLY: U_MortarPrim
#if PARABOLIC
USE MOD_Lifting             ,ONLY: Lifting
USE MOD_Lifting_Vars
#endif /*PARABOLIC*/
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI                 ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,           ONLY: nSides
#endif /*USE_MPI*/
#if FV_ENABLED
USE MOD_FV_Vars             ,ONLY: FV_Elems_master,FV_Elems_slave,FV_Elems_Sum
USE MOD_FV_Mortar           ,ONLY: FV_Elems_Mortar
USE MOD_FV                  ,ONLY: FV_DGtoFV,FV_ConsToPrim
USE MOD_FV_VolInt           ,ONLY: FV_VolInt
#if ((FV_ENABLED >= 2) && (PP_NodeType == 1))
USE MOD_FV_Vars             ,ONLY: FV_U_master,FV_U_slave,FV_UPrim_master,FV_UPrim_slave
USE MOD_FV_Vars             ,ONLY: FV_Flux_master,FV_Flux_slave
#endif /*((FV_ENABLED >= 2) && (PP_NodeType == 1))*/
#if USE_MPI
USE MOD_MPI                 ,ONLY: StartExchange_FV_Elems
#endif /*USE_MPI*/
#if FV_RECONSTRUCT
USE MOD_FV_Vars             ,ONLY: gradUxi,gradUeta,gradUzeta
#if VOLINT_VISC
USE MOD_FV_Vars             ,ONLY: gradUxi_central,gradUeta_central,gradUzeta_central
USE MOD_FV_Vars             ,ONLY: FV_surf_gradU_master,FV_surf_gradU_slave
USE MOD_FV_Reconstruction   ,ONLY: FV_SurfCalcGradients_Parabolic
#endif /* VOLINT_VISC */
USE MOD_FV_Vars             ,ONLY: FV_surf_gradU,FV_multi_master,FV_multi_slave
USE MOD_FV_ProlongToFace    ,ONLY: FV_ProlongToDGFace
USE MOD_FV_Mortar           ,ONLY: FV_gradU_mortar
USE MOD_FV_Reconstruction   ,ONLY: FV_PrepareSurfGradient,FV_SurfCalcGradients,FV_SurfCalcGradients_BC,FV_CalcGradients
#endif /* FV_RECONSTRUCT */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars       ,ONLY: ComputeEddyViscosity, muSGS, muSGS_master, muSGS_slave
USE MOD_ProlongToFace       ,ONLY: ProlongToFace
USE MOD_TimeDisc_Vars       ,ONLY: CurrentStage
#endif
#if PP_EntropyVars==1
USE MOD_DG_Vars             ,ONLY: V,V_slave,V_master
USE MOD_EOS                 ,ONLY: ConsToEntropy
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< Current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! -----------------------------------------------------------------------------
! MAIN STEPS        []=FV only
! -----------------------------------------------------------------------------
! 1.   Filter solution vector
! 2.   Convert volume solution to primitive
! 3.   Prolong to face (fill U_master/slave)
! 4.   ConsToPrim of face data (U_master/slave)
![5. ] Second order reconstruction for FV
! 6.   Lifting
! 7.   IF EDDYVISCOSITY: Prolong muSGS to face and send from slave to master
! 8.   Volume integral (DG only)
![9. ] FV volume integral
![10.] Volume integral (viscous contribution) if FV-blending
! 11.  Fill flux (Riemann solver) + surface integral
! 12.  Ut = -Ut
! 13.  Sponge and source terms
! 14.  Perform overintegration and apply Jacobian
! -----------------------------------------------------------------------------

! (0. Nullify arrays)
! NOTE: UT and U are nullified in DGInit, and Ut is set directly in the volume integral, so in this implementation,
!       ARRAYS DO NOT NEED TO BE NULLIFIED, OTHERWISE THEY HAVE TO!
! CALL VNullify(nTotalU,Ut)

! 1. Filter the solution vector if applicable, filter_pointer points to cut-off filter or LAF filter (see filter.f90)
IF(FilterType.GT.0) CALL Filter_Pointer(U,FilterMat)

! 2. Convert Volume solution to primitive
CALL ConsToPrim(PP_N,UPrim,U)

! Compute entropy variables
#if PP_EntropyVars == 1
Call ConsToEntropy(PP_N,V,U)
#endif

! 3. Prolong the solution to the face integration points for flux computation (and do overlapping communication)
! -----------------------------------------------------------------------------------------------------------
! General idea: The slave sends its surface data to the master, where the flux is computed and sent back to the slaves.
! Steps:
! * (these steps are done for all slave MPI sides first and then for all remaining sides):
! 3.1)  Prolong solution to faces and store in U_master/slave. Use them to build mortar data (split into 2/4 smaller sides).
!       Then U_slave can be communicated from the slave to master MPI side.
![3.2)] The information which element is a DG or FV subcells element is stored in FV_Elems per element. To know which of the
!       data inside the face-arrays U_master/slave is DG or FV the FV_Elems is 'prolongated' (copied) to FV_Elems_master/slave.
!       These directly correspond to U_master/slave and must be handled in the same way as U_master/slave. Therefore they are
!       'mortarized' and then FV_Elems_slave is transmitted like U_slave.
![3.3)] The reconstruction of slopes over element interfaces requires, besides U_slave and FV_Elems_slave, some more
!       information that has to be transmitted from the slave to the master MPI side (same direction as U_slave and
!       FV_Elems_slave), which does the whole flux computation.
!       This additional data is different for the two the element types (DG or FV) and is stored in the multipurpose array
!       FV_multi_master/slave.
! 3.4)  Finish all started MPI communications (after step 2. due to latency hiding)

#if USE_MPI
! Step 3 for all slave MPI sides
! 3.1)
CALL StartReceiveMPIData(U_slave   ,DataSizeSide,1,nSides,MPIRequest_U(:,SEND)   ,SendID=2) ! Receive MINE / U_slave: slave -> master
#if (FV_ENABLED == 2) && (PP_NodeType==1)
CALL StartReceiveMPIData(FV_U_slave,DataSizeSide,1,nSides,MPIRequest_FV_U(:,SEND),SendID=2) ! Receive MINE / FV_U_slave: slave -> master
#endif

#if (FV_ENABLED == 2) && (PP_NodeType==1)
#if PP_EntropyVars == 0
CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL ProlongToFaceCons(PP_N,U,FV_U_master,FV_U_slave,L_Minus,L_Plus,doMPISides=.TRUE.,pureFV=.TRUE.)
#endif /*if PP_EntropyVars == 0*/
#else /*FV_ENABLED*/
#if PP_EntropyVars == 1
CALL ProlongToFaceCons(PP_N,V,V_master,V_slave,U_master,U_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
#else
CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
#endif /*if PP_EntropyVars == 1*/
#endif /*FV_ENABLED*/

CALL U_MortarCons(U_master,U_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,RECV)   ,SendID=2) ! SEND YOUR / U_slave: slave -> master
#if (FV_ENABLED == 2) && (PP_NodeType==1)
CALL U_MortarCons(FV_U_master,FV_U_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(FV_U_slave,DataSizeSide,1,nSides,MPIRequest_FV_U(:,RECV),SendID=2) ! SEND YOUR / FV_U_slave: slave -> master
#endif
#if FV_ENABLED
! 3.2)
CALL FV_Elems_Mortar(FV_Elems_master,FV_Elems_slave,doMPISides=.TRUE.)
CALL StartExchange_FV_Elems(FV_Elems_slave,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=2)
                                                                 ! Receive MINE, Send YOUR / FV_Elems_slave: slave -> master
#if FV_RECONSTRUCT
! 3.3)
CALL StartReceiveMPIData(FV_multi_slave,DataSizeSidePrim,1,nSides,MPIRequest_FV_gradU(:,SEND),SendID=2)
                                                                 ! Receive MINE / FV_multi_slave: slave -> master
CALL FV_PrepareSurfGradient(UPrim,FV_multi_master,FV_multi_slave,doMPiSides=.TRUE.)
CALL U_MortarPrim(FV_multi_master,FV_multi_slave,doMPiSides=.TRUE.)
CALL StartSendMPIData(   FV_multi_slave,DataSizeSidePrim,1,nSides,MPIRequest_FV_gradU(:,RECV),SendID=2)
                                                                 ! SEND YOUR / FV_multi_slave: slave -> master
#endif /* FV_RECONSTRUCT */
#endif /* FV_ENABLED */
#endif /*USE_MPI*/

! Step 3 for all remaining sides
! 3.1)
#if (FV_ENABLED == 2) && (PP_NodeType==1)
#if PP_EntropyVars == 0
CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL ProlongToFaceCons(PP_N,U,FV_U_master,FV_U_slave,L_Minus,L_Plus,doMPISides=.FALSE.,pureFV=.TRUE.)
#endif /*if PP_EntropyVars == 0*/
#else /*FV_ENABLED*/
#if PP_EntropyVars == 1
CALL ProlongToFaceCons(PP_N,V,V_master,V_slave,U_master,U_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#else
CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#endif /*if PP_EntropyVars == 1*/
#endif

CALL U_MortarCons(U_master,U_slave,doMPISides=.FALSE.)
#if (FV_ENABLED == 2) && (PP_NodeType==1)
CALL U_MortarCons(FV_U_master,FV_U_slave,doMPISides=.FALSE.)
#endif
#if FV_ENABLED
! 3.2)
CALL FV_Elems_Mortar(FV_Elems_master,FV_Elems_slave,doMPISides=.FALSE.)
#if FV_RECONSTRUCT
! 3.3)
CALL FV_PrepareSurfGradient(UPrim,FV_multi_master,FV_multi_slave,doMPiSides=.FALSE.)
CALL U_MortarPrim(FV_multi_master,FV_multi_slave,doMPiSides=.FALSE.)
#endif
#endif

#if USE_MPI
! 3.4) complete send / receive of side data from step 3.
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! U_slave: slave -> master
#if FV_ENABLED
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems) ! FV_Elems_slave: slave -> master
! communicate FV_Elems_master from master to slave to define FV_Elems_Sum properly also for YOUR sides (needed e.g. by Lifting_FillFlux_NormVec)
CALL StartExchange_FV_Elems(FV_Elems_master,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=1)
                                                            ! Receive YOUR, Send MINE / FV_Elems_master: master -> slave
#if (FV_ENABLED == 2) && (PP_NodeType==1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_U)     ! FV_U_slave: slave -> master
#endif
#if FV_RECONSTRUCT
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_gradU) ! FV_multi_slave: slave -> master
#endif /*FV_RECONSTRUCT*/
#endif /*FV_ENABLED*/
#endif /*USE_MPI*/

! 4. Convert face data from conservative to primitive variables
!    Attention: For FV with 2nd order reconstruction U_master/slave and therewith UPrim_master/slave are still only 1st order
! TODO: Linadv?
CALL GetPrimitiveStateSurface(U_master,U_slave,UPrim_master,UPrim_slave)
#if (FV_ENABLED == 2) && (PP_NodeType==1)
CALL GetPrimitiveStateSurface(FV_U_master,FV_U_slave,FV_UPrim_master,FV_UPrim_slave)
#endif
#if FV_ENABLED
#if USE_MPI
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems) ! FV_Elems_master: master -> slave
#endif
! Build four-states-array for the 4 different combinations DG/DG(0), FV/DG(1), DG/FV(2) and FV/FV(3) a face can be.
FV_Elems_Sum = FV_Elems_master + 2*FV_Elems_slave
#endif /*FV_ENABLED*/

#if FV_ENABLED && FV_RECONSTRUCT
! [ 5. Second order reconstruction (computation of slopes) ]
!-------------------------------------------------------
! General idea at the faces: With the slave data from step 3.) reconstruct the slope over the interface on the master side
!    and send it back to the slave
! Steps:
! * (steps 5.2 and 5.4 are done for all master MPI sides first and then for all remaining sides)
! 5.1) Convert FV_multi_master/slave (only the DG parts of it) from DG nodes to FV nodes (equidistant)
! 5.2) Reconstruct the slope over the interface (and send it from master to slave)
! 5.3) On the slave side combine the slopes from the 2/4 small mortar sides to the big mortar side (when communication finished)
! 5.4) Calculate slopes at boundary conditions
! 5.5) Use the slope to prolongate the solution to UPrim_master/slave (ATTENTION: U_master/slave are only 1st order!)
! 5.6) Calculate the inner (volume) slopes

! 5.1)
CALL FV_DGtoFV(PP_nVarPrim,FV_multi_master,FV_multi_slave)

#if USE_MPI
! 5.2)
CALL StartReceiveMPIData(FV_surf_gradU,DataSizeSidePrim,1,nSides,MPIRequest_Flux(:,SEND),SendID=1)
                                                         ! Receive YOUR / FV_surf_gradU: master -> slave
CALL FV_SurfCalcGradients(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,&
    FV_surf_gradU,doMPISides=.TRUE.)
CALL StartSendMPIData(   FV_surf_gradU,DataSizeSidePrim,1,nSides,MPIRequest_Flux(:,RECV),SendID=1)
                                                         ! Send MINE  /   FV_surf_gradU: master -> slave
! 5.4)
CALL FV_ProlongToDGFace(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,FV_surf_gradU,doMPISides=.TRUE.)
#endif /*USE_MPI*/

! Calculate FV-Gradients over inner Sides
! 5.2)
CALL FV_SurfCalcGradients(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,&
    FV_surf_gradU,doMPISides=.FALSE.)
! 5.3)
CALL FV_gradU_mortar(FV_surf_gradU,doMPISides=.FALSE.)
#if USE_MPI
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux)   ! FV_surf_gradU: master -> slave
CALL FV_gradU_mortar(FV_surf_gradU,doMPISides=.TRUE.)
#endif /*USE_MPI*/
! 5.4)
CALL FV_SurfCalcGradients_BC(UPrim_master,FV_surf_gradU,t)
! 5.5)
CALL FV_ProlongToDGFace(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,FV_surf_gradU,doMPISides=.FALSE.)
! 5.6)
#if VOLINT_VISC && USE_MPI
CALL StartReceiveMPIData(FV_surf_gradU_slave,DataSizeSideGradParabolic,1,nSides,MPIRequest_FV_gradU(:,RECV),SendID=2)
#endif
CALL FV_CalcGradients(UPrim,FV_surf_gradU,gradUxi,gradUeta,gradUzeta &
#if VOLINT_VISC
    ,gradUxi_central,gradUeta_central,gradUzeta_central &
#endif /* VOLINT_VISC */
    )
#if VOLINT_VISC && USE_MPI
CALL StartSendMPIData(FV_surf_gradU_slave,DataSizeSideGradParabolic,1,nSides,MPIRequest_FV_gradU(:,SEND),SendID=2)
#endif
#endif /* FV_ENABLED && FV_RECONSTRUCT */

#if PARABOLIC
! 6. Lifting
! Compute the gradients using Lifting (BR1 scheme,BR2 scheme ...)
! The communication of the gradients is initialized within the lifting routines
CALL Lifting(UPrim,UPrim_master,UPrim_slave,t)

#if EDDYVISCOSITY
! 7. [ After the lifting we can now compute the eddy viscosity, which then has to be evaluated at the boundary. ]
! 7.1) - [ Open receive channel ]
! 7.2) - Compute SGS viscosity
! 7.3).  Prolong muSGS to face and send from slave to master, first MPI sides then inner sides
! 7.4) - After step 6. receive SGS data
IF(CurrentStage.EQ.1) THEN
#if USE_MPI
  CALL StartReceiveMPIData(muSGS_slave,DataSizeSideSGS,1,nSides,MPIRequest_SGS(:,RECV),SendID=2)
#endif /*USE_MPI*/
  CALL ComputeEddyViscosity()
#if USE_MPI
  CALL ProlongToFace(1,PP_N,muSGS,muSGS_master,muSGS_slave,L_Minus,L_Plus,.TRUE.)
  CALL StartSendMPIData   (muSGS_slave,DataSizeSideSGS,1,nSides,MPIRequest_SGS(:,SEND),SendID=2)
#endif /*USE_MPI*/
  CALL ProlongToFace(1,PP_N,muSGS,muSGS_master,muSGS_slave,L_Minus,L_Plus,.FALSE.)
END IF
#endif /* EDDYVISCOSITY */
#endif /*PARABOLIC*/

! 8. Compute volume integral contribution and add to Ut
CALL VolInt(Ut)

#if FV_ENABLED
! [ 9. Volume integral (advective and viscous) for all FV elements ]
CALL FV_VolInt(UPrim,Ut)
#endif /*FV_ENABLED*/


#if (FV_ENABLED >= 2) && PARABOLIC
! [10. Compute viscous volume integral contribution separately and add to Ut (FV-blending only)]
CALL VolInt_Visc(Ut)
#endif

#if PARABOLIC && USE_MPI
#if EDDYVISCOSITY
IF(CurrentStage.EQ.1) THEN
  CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_SGS)  ! muSGS_slave: slave -> master
END IF
#endif /* EDDYVISCOSITY */
! Complete send / receive for gradUx, gradUy, gradUz, started in the lifting routines
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU) ! gradUx,y,z: slave -> master
#endif /*PARABOLIC && USE_MPI*/


! 11. Fill flux and Surface integral
! General idea: U_master/slave and gradUx,y,z_master/slave are filled and can be used to compute the Riemann solver
!               and viscous flux at the faces. This is done for the MPI master sides first, to start communication early
!               and then for all other sides.
!               At mixed DG/FV interfaces the flux is computed in FV points, therefore the DG part of the solution and gradients
!               at mixed interfaces must be converted from DG to FV representation.
!               After communication from master to slave the flux can be integrated over the faces.
! Steps:
! * (step 11.2 is done for all MPI master sides first and then for all remaining sides)
! * (step 11.3 and 10.4 are done for all other sides first and then for the MPI master sides)
![11.1)] Change basis of DG solution and gradients at mixed FV/DG interfaces to the FV grid
![11.2)] Convert primitive face solution to conservative at FV faces
! 11.3)  Fill flux (Riemann solver + viscous flux)
! 11.4)  Combine fluxes from the 2/4 small mortar sides to the flux on the big mortar side (when communication finished)
! 11.5)  Compute surface integral
#if FV_ENABLED
! 11.1)
#if VOLINT_VISC
CALL FV_DGtoFV(PP_nVarLifting,gradUx_master,gradUx_slave)
CALL FV_DGtoFV(PP_nVarLifting,gradUy_master,gradUy_slave)
CALL FV_DGtoFV(PP_nVarLifting,gradUz_master,gradUz_slave)
#if USE_MPI
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_gradU) ! FV_surf_gradU_slave: slave -> master
#endif
CALL FV_SurfCalcGradients_Parabolic()
#endif /* VOLINT_VISC */

CALL FV_DGtoFV(PP_nVar    ,U_master     ,U_slave     )
CALL FV_ConsToPrim(PP_nVarPrim,PP_nVar,UPrim_master,UPrim_slave,U_master,U_slave)
#if FV_RECONSTRUCT
! 10.2)
CALL GetConservativeStateSurface(UPrim_master, UPrim_slave, U_master, U_slave, FV_Elems_master, FV_Elems_slave, 1)
#endif /*FV_RECONSTRUCT*/
#endif /*FV_ENABLED*/

#if USE_MPI
! 11.3)
CALL StartReceiveMPIData(Flux_slave, DataSizeSide, 1,nSides,MPIRequest_Flux( :,SEND),SendID=1)
                                                                              ! Receive YOUR / Flux_slave: master -> slave
CALL FillFlux(t,Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   Flux_slave, DataSizeSide, 1,nSides,MPIRequest_Flux( :,RECV),SendID=1)
                                                                              ! Send MINE  /   Flux_slave: master -> slave

#if ((FV_ENABLED == 2) && (PP_NodeType == 1))
CALL StartReceiveMPIData(FV_Flux_slave, DataSizeSide, 1,nSides,MPIRequest_FV_Flux( :,SEND),SendID=1)
                                                                              ! Receive YOUR / Flux_slave: master -> slave
CALL FillFlux(t,FV_Flux_master,FV_Flux_slave,FV_U_master,FV_U_slave,FV_UPrim_master,FV_UPrim_slave,doMPISides=.TRUE.,pureFV=.TRUE.)
CALL StartSendMPIData(   FV_Flux_slave, DataSizeSide, 1,nSides,MPIRequest_FV_Flux( :,RECV),SendID=1)
                                                                              ! Send MINE  /   Flux_slave: master -> slave
#endif /*((FV_ENABLED == 2) && (PP_NodeType == 1))*/
#endif /*USE_MPI*/

! 11.3)
CALL FillFlux(t,Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave,doMPISides=.FALSE.)
! 11.4)
CALL Flux_MortarCons(Flux_master,Flux_slave,doMPISides=.FALSE.,weak=.TRUE.)

#if ((FV_ENABLED == 2) && (PP_NodeType == 1))
! 11.3)
CALL FillFlux(t,FV_Flux_master,FV_Flux_slave,FV_U_master,FV_U_slave,FV_UPrim_master,FV_UPrim_slave,doMPISides=.FALSE.)
! 11.4)
CALL Flux_MortarCons(FV_Flux_master,FV_Flux_slave,doMPISides=.FALSE.,weak=.TRUE.)
! 11.5)
CALL SurfIntCons(PP_N,Flux_master,Flux_slave,FV_Flux_master,FV_Flux_slave,Ut,.FALSE.,L_HatMinus,L_hatPlus)
#else
! 11.5)
CALL SurfIntCons(PP_N,Flux_master,Flux_slave,Ut,.FALSE.,L_HatMinus,L_hatPlus)
#endif

#if USE_MPI
! 11.4)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux )                       ! Flux_slave: master -> slave
CALL Flux_MortarCons(Flux_master,Flux_slave,doMPISides=.TRUE.,weak=.TRUE.)
#if ((FV_ENABLED == 2) && (PP_NodeType == 1))
! 11.4)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Flux )                       ! Flux_slave: master -> slave
CALL Flux_MortarCons(FV_Flux_master,FV_Flux_slave,doMPISides=.TRUE.,weak=.TRUE.)
! 11.5)
CALL SurfIntCons(PP_N,Flux_master,Flux_slave,FV_Flux_master,FV_Flux_slave,Ut,.TRUE.,L_HatMinus,L_HatPlus)
#else
! 11.5)
CALL SurfIntCons(PP_N,Flux_master,Flux_slave,Ut,.TRUE.,L_HatMinus,L_HatPlus)
#endif
#endif /*USE_MPI*/

! 12. Swap to right sign :)
Ut=-Ut

! 13. Compute source terms and sponge (in physical space, conversion to reference space inside routines)
IF(doCalcSource) CALL CalcSource(Ut,t)
IF(doSponge)     CALL Sponge(Ut)
IF(doTCSource)   CALL TestcaseSource(Ut)

! 14. Perform overintegration and apply Jacobian
! Perform overintegration (projection filtering type overintegration)
SELECT CASE (OverintegrationType)
  CASE (OVERINTEGRATIONTYPE_CONSCUTOFF )
    CALL Overintegration(Ut)
#if FV_ENABLED
    CALL ApplyJacobianCons(Ut,toPhysical=.TRUE.,FVE=1)
#endif /*FV_ENABLED*/
  CASE (OVERINTEGRATIONTYPE_CUTOFF)
    CALL Overintegration(Ut)
    CALL ApplyJacobianCons(Ut,toPhysical=.TRUE.)
  CASE DEFAULT
    CALL ApplyJacobianCons(Ut,toPhysical=.TRUE.)
END SELECT

END SUBROUTINE DGTimeDerivative_weakForm


!==================================================================================================================================
!> Fills the solution array U with a initial solution provided by the ExactFunc subroutine though interpolation
!==================================================================================================================================
SUBROUTINE FillIni(Nloc,xGP,U)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Exactfunc     ,ONLY: ExactFunc
USE MOD_Mesh_Vars     ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: Nloc                                    !< Polynomial degree of solution
REAL,INTENT(IN)                 :: xGP(3,    0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)  !< Coordinates of Gauss-points
REAL,INTENT(OUT)                :: U(PP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)  !< Solution array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================

! Evaluate the initial solution at the nodes and fill the solution vector U.
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    CALL ExactFunc(IniExactFunc,0.,xGP(1:3,i,j,k,iElem),U(:,i,j,k,iElem))
  END DO; END DO; END DO
END DO

END SUBROUTINE FillIni


!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeDG()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_DG_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(D)
SDEALLOCATE(D_T)
SDEALLOCATE(D_Hat)
SDEALLOCATE(D_Hat_T)
#if SPLIT_DG
SDEALLOCATE(DVolSurf)
#endif
SDEALLOCATE(L_HatMinus)
SDEALLOCATE(L_HatPlus)
SDEALLOCATE(U)
SDEALLOCATE(Ut)
SDEALLOCATE(U_master)
SDEALLOCATE(U_slave)
#if PP_EntropyVars == 1
SDEALLOCATE(V)
SDEALLOCATE(V_master)
SDEALLOCATE(V_slave)
#endif
SDEALLOCATE(Flux_master)
SDEALLOCATE(Flux_slave)
SDEALLOCATE(UPrim)
SDEALLOCATE(UPrim_master)
SDEALLOCATE(UPrim_slave)
SDEALLOCATE(UPrim_boundary)

DGInitIsDone = .FALSE.

END SUBROUTINE FinalizeDG

END MODULE MOD_DG
