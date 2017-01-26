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
!> \brief Computes the DGSEM spatial operator and updates residual Ut

!> Contains the routines to 
!> - initialize and finalize the DG global variables and the DG basis
!> - compute the DG spatial operators/residuals(Ut) using U from the volume, surface and source contribution, incl. 
!> lifting for the gradients and parallelization
!==================================================================================================================================
MODULE MOD_DG
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
INTERFACE FillIni
  MODULE PROCEDURE FillIni
END INTERFACE


! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitDG
  MODULE PROCEDURE InitDG
END INTERFACE


INTERFACE DGTimeDerivative_weakForm
  MODULE PROCEDURE DGTimeDerivative_weakForm
END INTERFACE


INTERFACE FinalizeDG
  MODULE PROCEDURE FinalizeDG
END INTERFACE


PUBLIC::InitDG,DGTimeDerivative_weakForm,FinalizeDG
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
USE MOD_Overintegration_Vars, ONLY: NOver,VdmNOverToN,xGPO,wGPO,OverintegrationType
USE MOD_Interpolation_Vars,   ONLY: xGP,wGP,L_minus,L_plus
USE MOD_Interpolation_Vars,   ONLY: InterpolationInitIsDone
USE MOD_Restart_Vars,         ONLY: DoRestart,RestartInitIsDone
USE MOD_Mesh_Vars,            ONLY: nElems,nSides,Elem_xGP,Elem_xGPO,MeshInitIsDone
USE MOD_ChangeBasis,          ONLY: ChangeBasis3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE :: L_minusO(:),L_plusO(:),L_HatMinusO(:),L_HatPlusO(:),D_O(:,:),D_TO(:,:),D_HatO(:,:) !< dummy variables for 
                                                                                                       !< selective overintegration
!==================================================================================================================================

! Check if all the necessary initialization is done before
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.DGInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    'InitDG not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

! Pre-compute the dg operator building blocks (differentiation matrices and prolongation operators) 
CALL InitDGBasis(PP_N, xGP,wGP,L_minus,L_plus,D ,D_T ,D_Hat ,D_Hat_T ,L_HatMinus ,L_HatPlus)

! Allocate the local DG solution (JU or U): element-based 
ALLOCATE(U(        PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
! Allocate the time derivative / solution update /residual vector dU/dt: element-based
ALLOCATE(Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
U=0.
Ut=0.

! Allocate the 2D solution vectors on the sides, one array for the data belonging to the proc (the master) 
! and one for the sides which belong to another proc (slaves): side-based
ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_N,1:nSides))
U_master=0.
U_slave=0.

! Repeat the U, U_Minus, U_Plus structure for the primitive quantities
ALLOCATE(UPrim(       PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(UPrim_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_N,1:nSides))
UPrim=0.
UPrim_master=0.
UPrim_slave=0.

! Allocate two fluxes per side (necessary for coupling of FV and DG)
ALLOCATE(Flux_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Flux_slave (PP_nVar,0:PP_N,0:PP_N,1:nSides))
Flux_master=0.
Flux_slave=0.

! variables for performance tricks
nDOFElem=(PP_N+1)**3
nTotalU=PP_nVar*nDOFElem*nElems

! Allocate variables for selective Overintegration of advective fluxes
IF(OverintegrationType.EQ.SELECTIVE)THEN 
  ALLOCATE(L_minusO(0:NOver)) ! dummy variable
  ALLOCATE(L_plusO(0:NOver))  ! dummy variable
  ! The following call of InitDGBasis just creates D_Hat_TO. All other outputs of this
  ! routine are thrown away.
  CALL InitDGBasis(NOver,xGPO,wGPO,L_minusO,L_plusO,D_O,D_TO,D_HatO,D_Hat_TO,L_HatMinusO,L_HatPlusO)
  ! Deallocate dummy variables
  SDEALLOCATE(D_O)
  SDEALLOCATE(D_TO)
  SDEALLOCATE(D_HatO)
  SDEALLOCATE(L_MinusO)
  SDEALLOCATE(L_PlusO)
  SDEALLOCATE(L_HatMinusO)
  SDEALLOCATE(L_HatPlusO)

  ALLOCATE(UO    (PP_nVar    ,0:NOver,0:NOver,0:NOver,nElems))
  ALLOCATE(UPrimO(PP_nVarPrim,0:NOver,0:NOver,0:NOver,nElems))
  ALLOCATE(UtO   (PP_nVar    ,0:NOver,0:NOver,0:NOver,nElems))
  UO=0.
  UPrimO=0.
  UtO=0.
  nDOFElemO=(NOver+1)**3 ! variable for performance tricks

  ALLOCATE(U_masterO(PP_nVar,0:NOver,0:NOver,1:nSides))
  ALLOCATE(U_slaveO( PP_nVar,0:NOver,0:NOver,1:nSides))
  ALLOCATE(UPrim_masterO(PP_nVarPrim,0:NOver,0:NOver,1:nSides))
  ALLOCATE(UPrim_slaveO( PP_nVarPrim,0:NOver,0:NOver,1:nSides))
  ALLOCATE(FluxO(    PP_nVar,0:NOver,0:NOver,1:nSides))
  U_masterO=0.
  U_slaveO=0.
  UPrim_masterO=0.
  UPrim_slaveO=0.
  FluxO=0.
END IF

! Fill the solution vector U with the initial solution by interpolation, if not filled through restart already
IF(.NOT.DoRestart)THEN
  IF(OverintegrationType.EQ.SELECTIVE) THEN 
    ! Interpolate onto the NOver grid, then project onto N
    CALL FillIni(NOver,Elem_xGPO,UO)
    CALL ChangeBasis3D(PP_nVar,nElems,NOver,PP_N,VdmNOverToN,UO,U,.FALSE.)
  ELSE
    CALL FillIni(PP_N,Elem_xGP,U)
  END IF
END IF

DGInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDG



!==================================================================================================================================
!> Allocate and initialize the building blocks for the DG operator: Differentiation matrices and prolongation operators
!==================================================================================================================================
SUBROUTINE InitDGbasis(N_in,xGP,wGP,L_Minus,L_Plus,D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatPlus)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Interpolation,    ONLY: GetNodesAndWeights
USE MOD_Basis,            ONLY: PolynomialDerivativeMatrix,LagrangeInterpolationPolys
#ifdef SPLIT_DG
USE MOD_DG_Vars,          ONLY: DVolSurf ! Transpose of differentiation matrix used for calculating the strong form
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                             :: N_in                   !< Polynomial degree
REAL,DIMENSION(0:N_in),INTENT(IN)              :: xGP                    !< Gauss/Gauss-Lobatto Nodes 
REAL,DIMENSION(0:N_in),INTENT(IN)              :: wGP                    !< Gauss/Gauss-Lobatto Weights
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Minus                !< Values of lagrange polynomials at \f$ xi = -1 \f$  
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Plus                 !< Values of lagrange polynomials at \f$ xi = +1 \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D                      !< Differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_T                    !< Transpose of differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat                  !< Differentiation matrix premultiplied by mass matrix,
                                                                         !< \f$ \hat{D} = M^{-1} D^T M \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat_T                !< Transpose of D_Hat matrix
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatMinus             !< Values of lagrange polynomials at \f$ xi = -1 \f$
                                                                         !< premultiplied with mass matrix
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatPlus              !< Values of lagrange polynomials at \f$ xi = +1 \f$
                                                                         !< premultiplied with mass matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_in,0:N_in)              :: M,Minv
INTEGER                                    :: iMass
!==================================================================================================================================

ALLOCATE(L_HatMinus(0:N_in), L_HatPlus(0:N_in))
ALLOCATE(D(    0:N_in,0:N_in), D_T(    0:N_in,0:N_in))
ALLOCATE(D_Hat(0:N_in,0:N_in), D_Hat_T(0:N_in,0:N_in))
! Compute Differentiation matrix D for given Gausspoints
CALL PolynomialDerivativeMatrix(N_in,xGP,D)
D_T=TRANSPOSE(D)

! Build D_Hat matrix. (D^ = M^(-1) * D^T * M
M=0.
Minv=0.
DO iMass=0,N_in
  M(iMass,iMass)=wGP(iMass)
  Minv(iMass,iMass)=1./wGP(iMass)
END DO
D_Hat  = -MATMUL(Minv,MATMUL(TRANSPOSE(D),M))
D_Hat_T= TRANSPOSE(D_hat)

#ifdef SPLIT_DG
ALLOCATE(DVolSurf(0:N_in,0:N_in))
DVolSurf = D_T
DVolSurf(0,0) = DVolSurf(0,0) + 1.0/(2.0 * wGP(0))
DVolSurf(N_in,N_in) = DVolSurf(N_in,N_in) - 1.0/(2.0 * wGP(N_in))
#endif /*SPLIT_DG*/

! interpolate to left and right face (1 and -1 in reference space) and pre-divide by mass matrix
L_HatPlus  = MATMUL(Minv,L_Plus)
L_HatMinus = MATMUL(Minv,L_Minus)
END SUBROUTINE InitDGbasis



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
USE MOD_DG_Vars             ,ONLY: Ut,UtO,U,U_slave,U_master,UO,UPrimO,Flux_master,Flux_slave,L_HatPlus,L_HatMinus
USE MOD_DG_Vars             ,ONLY: D_Hat_TO,nDOFElemO
USE MOD_DG_Vars             ,ONLY: UPrim,UPrim_master,UPrim_slave
USE MOD_DG_Vars,             ONLY: nTotalU
USE MOD_VolInt
USE MOD_SurfIntCons         ,ONLY: SurfIntCons
USE MOD_ProlongToFaceCons   ,ONLY: ProlongToFaceCons
USE MOD_FillFlux            ,ONLY: FillFlux
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
USE MOD_Interpolation_Vars  ,ONLY: L_Minus,L_Plus
USE MOD_Overintegration_Vars,ONLY: NOver,VdmNOverToN,VdmNToNOver,OverintegrationType
USE MOD_Overintegration,     ONLY: Overintegration
USE MOD_ChangeBasis         ,ONLY: ChangeBasis3D
USE MOD_Testcase            ,ONLY: TestcaseSource
USE MOD_Testcase_Vars       ,ONLY: doTCSource
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Mesh_Vars,           ONLY: Metrics_fTildeO,Metrics_gTildeO,Metrics_hTildeO
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
#if MPI
USE MOD_MPI_Vars
USE MOD_MPI                 ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,           ONLY: nSides
#endif /*MPI*/
#if FV_ENABLED
USE MOD_FV_Vars             ,ONLY: FV_Elems_master,FV_Elems_slave,FV_Elems_Sum
USE MOD_FV_Mortar           ,ONLY: FV_Elems_Mortar
USE MOD_FV                  ,ONLY: FV_DGtoFV
USE MOD_FV_VolInt           ,ONLY: FV_VolInt
#if MPI
USE MOD_MPI                 ,ONLY: StartExchange_FV_Elems
#endif /*MPI*/
#if FV_RECONSTRUCT
USE MOD_FV_Vars             ,ONLY: gradUxi,gradUeta,gradUzeta,gradUxi_central,gradUeta_central,gradUzeta_central
USE MOD_FV_Vars             ,ONLY: FV_surf_gradU_master,FV_surf_gradU_slave,FV_multi_master,FV_multi_slave
USE MOD_FV_ProlongToFace    ,ONLY: FV_ProlongToDGFace
USE MOD_FV_Mortar           ,ONLY: FV_gradU_mortar
USE MOD_FV_Reconstruction   ,ONLY: FV_PrepareSurfGradient,FV_SurfCalcGradients,FV_SurfCalcGradients_BC,FV_CalcGradients
#endif /* FV_RECONSTRUCT */
#endif /* FV_ENABLED */
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars       ,ONLY:SGS_Ind,SGS_Ind_Master,SGS_Ind_Slave 
USE MOD_EddyVisc_Vars       ,ONLY:EddyViscType
USE MOD_ProlongToFace1      ,ONLY: ProlongToFace1
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< Current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! -----------------------------------------------------------------------------
! MAIN STEPS        []=FV only,  {}=selective overintegration only
! -----------------------------------------------------------------------------
! 0.  Convert volume solution to primitive 
! 1.  Pronlong to face (fill U_master/slave)
!{2.} Prepare terms for advective volume integral (DG selective only)
! 3.  ConsToPrim of face data (U_master/slave)
![4.] Second order reconstruction for FV 
! 5.  Lifting
! 6.  Volume integral (DG only)
![7.] FV volume integral
! 8.  Fill flux (Riemann solver) + surface integral
!{9.} Add advective volume (Ut0) integral to Ut for selective overintegration
! 10. Ut = -Ut
! 11. Sponge and source terms
! 12. Perform overintegration and apply Jacobian
! -----------------------------------------------------------------------------

! Nullify arrays
! TODO fix!!!
! NOTE: UT and U are nullified in DGInit, and Ut is set directly in the volume integral, so in this implementation, 
!       ARRAYS DO NOT NEED TO BE NULLIFIED, OTHERWISE THEY HAVE TO!
CALL VNullify(nTotalU,Ut)

! Filter the solution vector if applicable, filter_pointer points to cut-off filter or LAF filter (see filter.f90)
IF(FilterType.GT.0) CALL Filter_Pointer(U,FilterMat) 

! 0. Convert Volume solution to primitive 
CALL ConsToPrim(PP_N,UPrim,U)

! 1. Prolong the solution to the face integration points for flux computation (and do overlapping communication)
! -----------------------------------------------------------------------------------------------------------
! General idea: The slave sends its surface data to the master, where the flux is computed and sent back to the slaves.
! Steps:
! * (these steps are done for all slave MPI sides first and then for all remaining sides): 
! 1.1)  Prolong solution to faces and store in U_master/slave. Use them to build mortar data (split into 2/4 smaller sides).
!       Then U_slave can be communicated from the slave to master MPI side.
![1.2)] The information which element is a DG or FV subcells element is stored in FV_Elems per element. To know which of the
!       data inside the face-arrays U_master/slave is DG or FV the FV_Elems is 'prolongated' (copied) to FV_Elems_master/slave.
!       These directly correspond to U_master/slave and must be handled in the same way as U_master/slave. Therefore they are 
!       'mortarized' and then FV_Elems_slave is transmitted like U_slave.
![1.3)] The reconstruction of slopes over element interfaces requires, besides U_slave and FV_Elems_slave, some more 
!       information that has to be transmitted from the slave to the master MPI side (same direction as U_slave and 
!       FV_Elems_slave), which does the whole flux computation.
!       This additional data is different for the two the element types (DG or FV) and is stored in the multipurpose array
!       FV_multi_master/slave. 
! 1.4)  Finish all started MPI communications (after step 2. due to latency hiding) 



#if MPI
! Step 1 for all slave MPI sides
! 1.1)
CALL StartReceiveMPIData(U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE / U_slave: slave -> master
CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_MortarCons(U_master,U_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR / U_slave: slave -> master
#if FV_ENABLED
! 1.2)
CALL FV_Elems_Mortar(FV_Elems_master,FV_Elems_slave,doMPISides=.TRUE.) 
CALL StartExchange_FV_Elems(FV_Elems_slave,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=2) 
                                                                 ! Receive MINE, Send YOUR / FV_Elems_slave: slave -> master 
#if FV_RECONSTRUCT
! 1.3)
CALL StartReceiveMPIData(FV_multi_slave,DataSizeSidePrim,1,nSides,MPIRequest_FV_gradU(:,SEND),SendID=2) 
                                                                 ! Receive MINE / FV_multi_slave: slave -> master
CALL FV_PrepareSurfGradient(UPrim,FV_multi_master,FV_multi_slave,doMPiSides=.TRUE.)
CALL U_MortarPrim(FV_multi_master,FV_multi_slave,doMPiSides=.TRUE.)
CALL StartSendMPIData(   FV_multi_slave,DataSizeSidePrim,1,nSides,MPIRequest_FV_gradU(:,RECV),SendID=2)
                                                                 ! SEND YOUR / FV_multi_slave: slave -> master
#endif /* FV_RECONSTRUCT */
#endif /* FV_ENABLED */
#endif /* MPI */

! Step 1 for all remaining sides
! 1.1)
CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL U_MortarCons(U_master,U_slave,doMPISides=.FALSE.)
#if FV_ENABLED
! 1.2)
CALL FV_Elems_Mortar(FV_Elems_master,FV_Elems_slave,doMPISides=.FALSE.)
#if FV_RECONSTRUCT
! 1.3)
CALL FV_PrepareSurfGradient(UPrim,FV_multi_master,FV_multi_slave,doMPiSides=.FALSE.)
CALL U_MortarPrim(FV_multi_master,FV_multi_slave,doMPiSides=.FALSE.)
#endif
#endif

!{2.}Prepare data for selective volume integral (for DG elements only)
!    (latency hiding before finishing communication of side data in 1.4) )
IF(OverintegrationType.EQ.SELECTIVE)THEN
  CALL ChangeBasis3D(PP_nVar,nElems,PP_N,NOver,VdmNToNOver,U,UO,.FALSE.)
  CALL ConsToPrim(NOver,UPrimO,UO)
END IF

#if MPI
! 1.4) complete send / receive of side data from step 1.
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! U_slave: slave -> master 
#if FV_ENABLED
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems) ! FV_Elems_slave: slave -> master 
#if FV_RECONSTRUCT
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_gradU) ! FV_multi_slave: slave -> master
#endif
#endif
#endif

! 3. Convert face data from conservative to primitive variables 
!    Attention: For FV with 2nd order reconstruction U_master/slave and therewith UPrim_master/slave are still only 1st order
! TODO: Linadv?
CALL GetPrimitiveStateSurface(U_master,U_slave,UPrim_master,UPrim_slave)
#if FV_ENABLED
! Build four-states-array for the 4 different combinations DG/DG(0), FV/DG(1), DG/FV(2) and FV/FV(3) a face can be.
FV_Elems_Sum = FV_Elems_master + 2*FV_Elems_slave
#endif

#if FV_ENABLED && FV_RECONSTRUCT
! [ 4. Second order reconstruction (computation of slopes) ]
!-------------------------------------------------------
! General idea at the faces: With the slave data from step 1.) reconstruct the slope over the interface on the master side
!    and send it back to the slave
! Steps:
! * (steps 4.2 and 4.4 are done for all master MPI sides first and then for all remaining sides)
! 4.1) Convert FV_multi_master/slave (only the DG parts of it) from DG nodes to FV nodes (equidistant)
! 4.2) Reconstruct the slope over the interface (and send it from master to slave)
! 4.3) On the slave side combine the slopes from the 2/4 small mortar sides to the big mortar side (when communication finished)
! 4.4) Use the slope to prolongate the solution to UPrim_master/slave (ATTENTION: U_master/slave are only 1st order!)
! 4.5) Calculate slopes at boundary conditions 
! 4.6) Calculate the inner (volume) slopes

! 4.1)
CALL FV_DGtoFV(PP_nVarPrim,FV_multi_master,FV_multi_slave)

#if MPI
! 4.2)
CALL StartReceiveMPIData(FV_surf_gradU_slave,DataSizeSidePrim,1,nSides,MPIRequest_Flux(:,SEND),SendID=1)
                                                         ! Receive YOUR / FV_surf_gradU_slave: master -> slave
CALL FV_SurfCalcGradients(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,&
    FV_surf_gradU_master,FV_surf_gradU_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   FV_surf_gradU_slave,DataSizeSidePrim,1,nSides,MPIRequest_Flux(:,RECV),SendID=1)
                                                         ! Send MINE  /   FV_surf_gradU_slave: master -> slave
! 4.4)
CALL FV_ProlongToDGFace(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,FV_surf_gradU_master,FV_surf_gradU_slave,&
    doMPISides=.TRUE.) 
#endif /*MPI*/

! Calculate FD-Gradients over inner Sides
! 4.2)
CALL FV_SurfCalcGradients(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,&
    FV_surf_gradU_master,FV_surf_gradU_slave,doMPISides=.FALSE.)
! 4.3) 
CALL FV_gradU_mortar(FV_surf_gradU_master,FV_surf_gradU_slave,doMPISides=.FALSE.)
#if MPI
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux)   ! FV_surf_gradU_slave: master -> slave
CALL FV_gradU_mortar(FV_surf_gradU_master,FV_surf_gradU_slave,doMPISides=.TRUE.)
#endif
! 4.4)
CALL FV_ProlongToDGFace(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,FV_surf_gradU_master,FV_surf_gradU_slave,&
    doMPISides=.FALSE.)
! 4.5) 
CALL FV_SurfCalcGradients_BC(UPrim_master,FV_surf_gradU_master,t)
! 4.6) 
CALL FV_CalcGradients(UPrim,FV_surf_gradU_master,FV_surf_gradU_slave,gradUxi,gradUeta,gradUzeta,&
    gradUxi_central,gradUeta_central,gradUzeta_central)
#endif /* FV_ENABLED && FV_RECONSTRUCT */

#if PARABOLIC
! 5. Lifting 
! Compute the gradients using Lifting (BR1 scheme,BR2 scheme ...)
! The communication of the gradients is initialized within the lifting routines
CALL Lifting(UPrim,UPrim_master,UPrim_slave,t)
#endif /*PARABOLIC*/

! 6. Compute volume integral contribution and add to Ut
!    Compute separately in case of selective overintegration
IF(OverintegrationType.NE.SELECTIVE)THEN
  CALL VolInt(Ut)
ELSE
  CALL VolIntAdv(NOver,nDOFElemO,D_Hat_TO,Metrics_fTildeO,Metrics_gTildeO,Metrics_hTildeO,UO,UPrimO,UtO)
#if PARABOLIC
  CALL VolIntVisc(Ut)
#endif
END IF

#if FV_ENABLED
! [ 7. Volume integral (advective and viscous) for all FV elements ]
CALL FV_VolInt(UPrim,Ut)
#endif

#if EDDYVISCOSITY
IF(EddyViscType.EQ.2) THEN
#if MPI
! 4.2)
CALL StartReceiveMPIData(SGS_Ind_master,DataSizeSideScalar,1,nSides,MPIRequest_SGS_Ind(:,SEND),SendID=1)
                                                         ! Receive YOUR / FV_surf_gradU_slave: master -> slave
CALL ProlongToFace1(PP_N,SGS_Ind(1:2,:,:,:,:),SGS_Ind_master(:,:,:,:),SGS_Ind_Slave(:,:,:,:),L_Minus,L_Plus,.TRUE.)
CALL StartSendMPIData(SGS_Ind_master,DataSizeSideScalar,1,nSides,MPIRequest_SGS_Ind(:,RECV),SendID=1)
#endif
! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace1(PP_N,SGS_Ind(1:2,:,:,:,:),SGS_Ind_master(:,:,:,:),SGS_Ind_Slave(:,:,:,:),L_Minus,L_Plus,.FALSE.)
#if MPI  
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_SGS_Ind)  ! U_slave: slave -> master 
#endif

END IF
#endif



#if PARABOLIC && MPI
! Complete send / receive for gradUx, gradUy, gradUz, started in the lifting routines
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU) ! gradUx,y,z: slave -> master
#endif /*PARABOLIC && MPI*/


! 8. Fill flux and Surface integral
! General idea: U_master/slave and gradUx,y,z_master/slave are filled and can be used to compute the Riemann solver
!               and viscous flux at the faces. This is done for the MPI master sides first, to start communication early
!               and then for all other sides.
!               At mixed DG/FV interfaces the flux is computed in FV points, therefore the DG part of the solution and gradients
!               at mixed interfaces must be converted from DG to FV representation.
!               After communication from master to slave the flux can be integrated over the faces.
! Steps:
! * (step 8.2 is done for all MPI master sides first and then for all remaining sides)
! * (step 8.3 and 8.4 are done for all other sides first and then for the MPI master sides) 
![8.1)] Change basis of DG solution and gradients at mixed FV/DG interfaces to the FV grid
![8.2)] Convert primitive face solution to conservative at FV faces
! 8.3)  Fill flux (Riemann solver + viscous flux)
! 8.4)  Combine fluxes from the 2/4 small mortar sides to the flux on the big mortar side (when communication finished)
! 8.5)  Compute surface integral 
#if FV_ENABLED
! 8.1) 
#if PARABOLIC
CALL FV_DGtoFV(PP_nVarPrim,gradUx_master,gradUx_slave)
CALL FV_DGtoFV(PP_nVarPrim,gradUy_master,gradUy_slave)
CALL FV_DGtoFV(PP_nVarPrim,gradUz_master,gradUz_slave)
#endif
CALL FV_DGtoFV(PP_nVar    ,U_master     ,U_slave     )
CALL FV_DGtoFV(PP_nVarPrim,UPrim_master ,UPrim_slave )
! 8.2) 
CALL GetConservativeStateSurface(UPrim_master, UPrim_slave, U_master, U_slave, FV_Elems_master, FV_Elems_slave, 1)
#endif

#if MPI
! 8.3)
CALL StartReceiveMPIData(Flux_slave, DataSizeSide, 1,nSides,MPIRequest_Flux( :,SEND),SendID=1)
                                                                              ! Receive YOUR / Flux_slave: master -> slave
CALL FillFlux(t,Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   Flux_slave, DataSizeSide, 1,nSides,MPIRequest_Flux( :,RECV),SendID=1)
                                                                              ! Send MINE  /   Flux_slave: master -> slave
#endif /* MPI*/

! 8.3)
CALL FillFlux(t,Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave,doMPISides=.FALSE.)
! 8.4)
CALL Flux_MortarCons(Flux_master,Flux_slave,doMPISides=.FALSE.,weak=.TRUE.)
! 8.5)
CALL SurfIntCons(PP_N,Flux_master,Flux_slave,Ut,.FALSE.,L_HatMinus,L_hatPlus)

#if MPI
! 8.4)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux )                                        ! Flux_slave: master -> slave 
CALL Flux_MortarCons(Flux_master,Flux_slave,doMPISides=.TRUE.,weak=.TRUE.)
! 8.5)
CALL SurfIntCons(PP_N,Flux_master,Flux_slave,Ut,.TRUE.,L_HatMinus,L_HatPlus)
#endif /*MPI*/

! 9. Add advection volume integral to residual for selective overintegration
IF(OverintegrationType.EQ.SELECTIVE)THEN
  CALL ChangeBasis3D(PP_nVar,nElems,NOver,PP_N,VdmNOverToN,UtO,Ut,.TRUE.)
END IF

! 10. Swap to right sign :) 
Ut=-Ut

! 11. Compute source terms and sponge (in physical space, conversion to reference space inside routines)
IF(doCalcSource) CALL CalcSource(Ut,t)
IF(doSponge)     CALL Sponge(Ut)
IF(doTCSource)   CALL TestcaseSource(Ut)

! 12. Perform overintegration and apply Jacobian 
! Perform overintegration (projection filtering type overintegration)
IF(OverintegrationType.GT.0) THEN
  CALL Overintegration(Ut)
END IF
! Apply Jacobian (for OverintegrationType==CUTOFFCONS this is already done within the Overintegration, but for DG only)
IF (OverintegrationType.EQ.CUTOFFCONS) THEN 
#if FV_ENABLED
  CALL ApplyJacobianCons(Ut,toPhysical=.TRUE.,FVE=1)
#endif
ELSE
  CALL ApplyJacobianCons(Ut,toPhysical=.TRUE.)
END IF

END SUBROUTINE DGTimeDerivative_weakForm



!==================================================================================================================================
!> Fills the solution array U with a initial solution provided by the ExactFunc subroutine though interpolation
!==================================================================================================================================
SUBROUTINE FillIni(NLoc,xGP,U)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Exactfunc     ,ONLY: ExactFunc
USE MOD_Mesh_Vars     ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: NLoc                                    !< Polynomial degree of solution 
REAL,INTENT(IN)                 :: xGP(3,    0:NLoc,0:NLoc,0:NLoc,nElems)  !< Coordinates of Gauss-points
REAL,INTENT(OUT)                :: U(PP_nVar,0:NLoc,0:NLoc,0:NLoc,nElems)  !< Solution array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================

! Evaluate the initial solution at the nodes and fill the solutin vector U. 
DO iElem=1,nElems
  DO k=0,NLoc; DO j=0,NLoc; DO i=0,NLoc
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
SDEALLOCATE(L_HatMinus)
SDEALLOCATE(L_HatPlus)
SDEALLOCATE(U)
SDEALLOCATE(Ut)
SDEALLOCATE(UO)
SDEALLOCATE(UtO)
SDEALLOCATE(U_master)
SDEALLOCATE(U_slave)
SDEALLOCATE(Flux_master)
SDEALLOCATE(Flux_slave)
SDEALLOCATE(U_masterO)
SDEALLOCATE(U_slaveO)
SDEALLOCATE(FluxO)
SDEALLOCATE(UPrim)
SDEALLOCATE(UPrim_master)
SDEALLOCATE(UPrim_slave)
DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG


END MODULE MOD_DG
