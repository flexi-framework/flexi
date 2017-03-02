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
#if PARABOLIC
#include "flexi.h"

!==================================================================================================================================
!> \brief Contains the BR2 lifting procedure (initialization and lifting operator) for computing the lifted solution gradients
!> according to Bassi, Rebay et al., "A high-order accurate discontinuous Finite Element method for inviscid an viscous 
!> turbomachinery flows", 1997. The lifted gradients are required for the viscous fluxes.
!>
!> The BR1 scheme has been found to be unstable for purely elliptic equations. Consequently, the BR2 scheme provides a stable method
!> with local lifting operators. In contrast to the BR1 scheme, the BR2 scheme requires a strong form of the lifting operator. The
!> surface gradients are lifted with an additional penalty term \f$ \eta_{BR2} \f$. Stability was shown for \f$ \eta_{BR2} > \f$
!> number of element faces.
!> Similar to the BR1 scheme, both, conservative and non-conservative volume integrals are available.
!> 
!> NB: The lifting procedure is only necessary for a DG solution. If the element contains a FV solution, the lifting procedure is
!> skipped and instead central gradients are calculated for the comuptation of viscous fluxes
!==================================================================================================================================
MODULE MOD_Lifting
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitLifting
  MODULE PROCEDURE InitLifting
END INTERFACE

INTERFACE Lifting
  MODULE PROCEDURE Lifting
END INTERFACE

INTERFACE FinalizeLifting
  MODULE PROCEDURE FinalizeLifting
END INTERFACE

PUBLIC::DefineParametersLifting,InitLifting,Lifting,FinalizeLifting
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersLifting()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Lifting")
CALL prms%CreateLogicalOption('doConservativeLifting', "Set true to compute the volume contribution to the gradients in "//&
                                                       "conservative form, i.e. deriving the solution multiplied by the metric "//&
                                                       "terms instead of deriving the solution and multiplying by the metrics.",&
                                                       '.FALSE.')
CALL prms%CreateRealOption(   'etaBR2',                "Lifting penalty for BR2. Increase improves stability at the cost of "//&
                                                       "performance and reduces jumps between two cells.", '2.')
END SUBROUTINE DefineParametersLifting


!==================================================================================================================================
!> \brief Initialize the BR2 lifting: get parameters and allocate the arrays required for the BR2 lifting procedure.
!>
!> Important parameters:
!> - doConservativeLifting: If true, the volume contribution to the gradients is in conservative form, i.e. the solution is derived
!>   multiplied by the metrics terms
!> - etaBR2: Penalty term for the surface contribution of the BR2 lifting. Note, stability is shown only for \f$ \eta_{BR2} > \f$
!>   number of element faces
!>
!> Default are non-conservative lifting with etaBR2 = 2
!>
!> Note that the gradient arrays in x/y/z directions in the volume and on the surfaces contain the gradients of the primitive
!> variables, i.e. they must be allocated for PP_nVarPrim variables, i.e. \f$ (\rho, u_1,u_2,u_3,p,T) \f$.
!==================================================================================================================================
SUBROUTINE InitLifting()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars
USE MOD_Overintegration_Vars, ONLY: OverintegrationType,NOver
USE MOD_DG_Vars,              ONLY: DGInitIsDone
USE MOD_Mesh_Vars,            ONLY: nSides,nBCSides,nElems
USE MOD_ReadInTools,          ONLY: GETREAL,GETLOGICAL
#if USE_MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF((.NOT.DGInitIsDone).OR.LiftingInitIsDone)THEN
   SWRITE(*,*) "InitDG not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR2 ...'

doWeakLifting=.FALSE.
doConservativeLifting=GETLOGICAL('doConservativeLifting','.FALSE.')
etaBR2=GETREAL('etaBR2','2.')

! We store the interior gradients at the each element face
ALLOCATE(gradUx_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUx_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)) ! needed in 2D for interface compatibility in e.g. calcbodyforces
ALLOCATE(FluxX        (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(FluxY        (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(FluxZ        (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
gradUx_slave=0.
gradUy_slave=0.
gradUz_slave=0.
gradUx_master=0.
gradUy_master=0.
gradUz_master=0.
FluxX=0.
FluxY=0.
FluxZ=0.

IF(OverintegrationType.EQ.SELECTIVE)THEN
  ALLOCATE(gradUx_masterO(PP_nVarPrim,0:NOver,0:PP_NOverZ,1:nBCSides))
  ALLOCATE(gradUy_masterO(PP_nVarPrim,0:NOver,0:PP_NOverZ,1:nBCSides))
  ALLOCATE(gradUz_masterO(PP_nVarPrim,0:NOver,0:PP_NOverZ,1:nBCSides))
  gradUx_masterO=0.
  gradUy_masterO=0.
  gradUz_masterO=0.
ENDIF

! The gradients of the conservative variables are stored at each volume integration point
ALLOCATE(gradUx(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUy(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUz(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
gradUx=0.
gradUy=0.
gradUz=0. ! gradUz must be kept also in 2D case for interface compatibility with Lifting_Volint

LiftingInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LIFTING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLifting



!==================================================================================================================================
!> \brief Computes the DG gradients using the BR2 scheme in x/y/z direction.
!>
!> To calculate the lifted gradients we need to:
!> - For reducing the overhead of MPI communication, we calculate the surface contribution, i.e. the numerical fluxes, across the
!>   element faces at MPI boundaries first and communicate them right away. The fluxes are temporarilly stored in the gradient
!>   arrays.
!> - Then, the remainind surface fluxes across the inner element faces are computed
!> - Calculate the volume integral in conservative or non-conservative form
!> - The Jacobian is applied to the gradients, but only for the DG elements!
!> - Prolong the volume contribution to the interface: Note, this is different from the BR1 scheme. Here we prolong the volume
!>   contribution to the surface before any surface gradient has bee n applied to the colume terms
!> - The surface integral of the fluxes is calculated, here the strong form is used and the etaBR2 factor is applied.
!>   In the same step, the surface contribution is added to the prolonged volume contribution
!==================================================================================================================================
SUBROUTINE Lifting(UPrim,UPrim_master,UPrim_slave,t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars
USE MOD_Lifting_SurfInt,   ONLY: Lifting_SurfInt
USE MOD_Lifting_VolInt,    ONLY: Lifting_VolInt
USE MOD_ProlongToFacePrim, ONLY: ProlongToFacePrim
USE MOD_Lifting_FillFlux,  ONLY: Lifting_FillFlux,Lifting_FillFlux_BC
USE MOD_ApplyJacobianPrim, ONLY: ApplyJacobianPrim
USE MOD_Interpolation_Vars,ONLY: L_Minus,L_Plus
USE MOD_FillMortarPrim,    ONLY: U_MortarPrim,Flux_MortarPrim
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI,               ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_Mesh_Vars,         ONLY: nSides
USE MOD_Mesh_Vars,         ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: UPrim(       PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution vector for which lifted gradients will be computed
REAL,INTENT(IN) :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the master sides
REAL,INTENT(IN) :: UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the slave sides
REAL,INTENT(IN) :: t                                                      !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! fill the global surface flux list
!fluxX=0. !don't nullify fluxes if not really needed (very expensive)
!fluxY=0. !don't nullify fluxes if not really needed (very expensive)
!fluxZ=0. !don't nullify fluxes if not really needed (very expensive)
#if USE_MPI
! Receive YOUR
CALL StartReceiveMPIData(FluxX,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=1)
CALL StartReceiveMPIData(FluxY,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=1)
#if PP_dim == 3
CALL StartReceiveMPIData(FluxZ,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=1)
#endif
! Compute lifting MPI fluxes
CALL Lifting_FillFlux(1,UPrim_master,UPrim_slave,FluxX,doMPISides=.TRUE.)
CALL Lifting_FillFlux(2,UPrim_master,UPrim_slave,FluxY,doMPISides=.TRUE.)
#if PP_dim == 3
CALL Lifting_FillFlux(3,UPrim_master,UPrim_slave,FluxZ,doMPISides=.TRUE.)
#endif
! Start Send MINE
CALL StartSendMPIData(   FluxX,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=1)
CALL StartSendMPIData(   FluxY,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=1)
#if PP_dim == 3
CALL StartSendMPIData(   FluxZ,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=1)
#endif
#endif /*USE_MPI*/

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC(t,UPrim_master, FluxX, FluxY, FluxZ)
CALL Lifting_FillFlux(1,UPrim_master,UPrim_slave,FluxX,doMPISides=.FALSE.)
CALL Lifting_FillFlux(2,UPrim_master,UPrim_slave,FluxY,doMPISides=.FALSE.)
#if PP_dim == 3
CALL Lifting_FillFlux(3,UPrim_master,UPrim_slave,FluxZ,doMPISides=.FALSE.)
#endif

CALL Flux_MortarPrim(FluxX,FluxX,doMPISides=.FALSE.,weak=.FALSE.)
CALL Flux_MortarPrim(FluxY,FluxY,doMPISides=.FALSE.,weak=.FALSE.)
#if PP_dim == 3
CALL Flux_MortarPrim(FluxZ,FluxZ,doMPISides=.FALSE.,weak=.FALSE.)
#endif

! compute volume integral contribution and add to gradU
IF(doConservativeLifting)THEN
  CALL Lifting_VolInt(1,UPrim,GradUx)
  CALL Lifting_VolInt(2,UPrim,GradUy)
#if PP_dim == 3
  CALL Lifting_VolInt(3,UPrim,GradUz)
#endif
ELSE
  CALL Lifting_VolInt(UPrim,GradUx,GradUy,GradUz)
END IF

! Account for the jacobian
! The Lifting already has the right sign
! For FV elements no Jacobian is needed here, since already applied in Lifting_Volint (hidden in FV_Metrics_f/g/hTilde_sJ)
CALL ApplyJacobianPrim(gradUx,toPhysical=.TRUE.,FVE=0)
CALL ApplyJacobianPrim(gradUy,toPhysical=.TRUE.,FVE=0)
#if PP_dim == 3
CALL ApplyJacobianPrim(gradUz,toPhysical=.TRUE.,FVE=0)
#endif

! The volume contribution of the gradients must be interpolated to the face of the grid cells
#if USE_MPI
! Prolong to face for MPI sides - send direction
CALL ProlongToFacePrim(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL ProlongToFacePrim(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
#if PP_dim == 3
CALL ProlongToFacePrim(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
#endif
#endif /*USE_MPI*/
! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFacePrim(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL ProlongToFacePrim(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#if PP_dim == 3
CALL ProlongToFacePrim(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#endif

#if USE_MPI
! Complete send / receive
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU)
CALL Flux_MortarPrim(FluxX,FluxX,doMPISides=.TRUE.,weak=.FALSE.)
CALL Flux_MortarPrim(FluxY,FluxY,doMPISides=.TRUE.,weak=.FALSE.)
#if PP_dim == 3
CALL Flux_MortarPrim(FluxZ,FluxZ,doMPISides=.TRUE.,weak=.FALSE.)
#endif

CALL StartReceiveMPIData(gradUx_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=2)
CALL StartReceiveMPIData(gradUy_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=2)
#if PP_dim == 3
CALL StartReceiveMPIData(gradUz_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=2)
#endif

CALL Lifting_SurfInt(FluxX,gradUx,gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL U_MortarPrim(gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUx_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=2)
CALL Lifting_SurfInt(FluxY,gradUy,gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL U_MortarPrim(gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUy_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=2)
#if PP_dim == 3
CALL Lifting_SurfInt(FluxZ,gradUz,gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL U_MortarPrim(gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(gradUz_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=2)
#endif
#endif /*USE_MPI*/

! Add the surface lifting flux to the prolonged volume contributions of the gradients and computes the surface integral
CALL Lifting_SurfInt(FluxX,gradUx,gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL Lifting_SurfInt(FluxY,gradUy,gradUy_master,gradUy_slave,doMPISides=.FALSE.)
#if PP_dim == 3
CALL Lifting_SurfInt(FluxZ,gradUz,gradUz_master,gradUz_slave,doMPISides=.FALSE.)
#endif
CALL U_MortarPrim(gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL U_MortarPrim(gradUy_master,gradUy_slave,doMPISides=.FALSE.)
#if PP_dim == 3
CALL U_MortarPrim(gradUz_master,gradUz_slave,doMPISides=.FALSE.)
#endif

END SUBROUTINE Lifting



!==================================================================================================================================
!> Deallocate BR2 arrays (volume and surface gradients and gradient fluxes)
!==================================================================================================================================
SUBROUTINE FinalizeLifting()
! MODULES
USE MOD_Lifting_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(gradUx_slave)
SDEALLOCATE(gradUy_slave)
SDEALLOCATE(gradUz_slave)
SDEALLOCATE(gradUx_master)
SDEALLOCATE(gradUy_master)
SDEALLOCATE(gradUz_master)
SDEALLOCATE(gradUx)
SDEALLOCATE(gradUy)
SDEALLOCATE(gradUz)
SDEALLOCATE(FluxX)
SDEALLOCATE(FluxY)
SDEALLOCATE(FluxZ)
SDEALLOCATE(gradUx_masterO)
SDEALLOCATE(gradUy_masterO)
SDEALLOCATE(gradUz_masterO)
LiftingInitIsDone = .FALSE.
END SUBROUTINE FinalizeLifting

END MODULE MOD_Lifting
#endif /* PARABOLIC */
