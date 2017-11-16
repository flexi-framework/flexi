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
!> \brief Contains the BR1 lifting procedures (initialization and lifting operator) for computing the lifted solution gradients
!> according to Bassi & Rebay 1997. The lifted gradients are required for the viscous fluxes.
!>
!> Local gradients of the DG polynomial are not well suited to compute the viscous fluxes due to their discontinuous nature. 
!> Instead, modified gradients \f$ Q \approx \nabla_x U \f$ are introduced and this equation is discretized using the DG method.
!> The equation for the gradients can be discretized in the weak or the strong form, which implies integration by parts
!> once (weak) or twice (strong), see "On the Quadrature and Weak Form Choices in Collocation Type Discontinuous Galerkin Spectral
!> Element Methods" (Gassner & Kopriva 2010) for details. If the strong form is chosen, the volume integral can be computed in a
!> conservative and a non conservative way, see "Implementing Spectral Methods for Partial Differential Equations" (Kopriva 2009)
!> for details. In the non conservative form we derive the solution and multiply by the metric terms while the conservative
!> formulation derives the solution multiplied by the metric terms.
!> The numerical flux for the BR1 procedure is simply choosen as the arithmetic mean of the solution.
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
!==================================================================================================================================
CALL prms%SetSection("Lifting")
CALL prms%CreateLogicalOption('doWeakLifting',         "Set true to perform lifting in weak form.", '.FALSE.')
CALL prms%CreateLogicalOption('doConservativeLifting', "Set true to compute the volume contribution to the gradients in "//&
                                                       "conservative form, i.e. deriving the solution multiplied by the metric "//&
                                                       "terms instead of deriving the solution and multiplying by the metrics. "//&
                                                       "Only available for doWeakLifting=.FALSE.",&
                                                       '.FALSE.')
END SUBROUTINE DefineParametersLifting


!==================================================================================================================================
!> \brief Initialize BR1 lifting: get parameters and allocate arrays required for the BR1 lifting procedure.
!>
!> Important parameters:
!> - doWeakLifting will set the lifting procedure to be performed in weak or strong form
!> - In the strong form, the lifting can be performed in a conservative or non conservative  version
!>   using the doConservativeLifting parameter
!>
!> Default ist the non conservative form since this version has the fastest implementation.
!>
!> The arrays containing the lifted gradients in x/y/z direction in the volume as well as on the element faces
!> will be allocated and nullified if necessary.
!>
!> Note that the gradient arrays in x/y/z directions in the volume and on the surfaces contain the gradients of the primitive
!> variables, i.e. they must be allocated for PP_nVarPrim variables, i.e. \f$ (\rho, u_1,u_2,u_3,p,T) \f$.
!==================================================================================================================================
SUBROUTINE InitLifting()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars
USE MOD_DG_Vars,             ONLY: DGInitIsDone
USE MOD_Mesh_Vars,           ONLY: nSides
USE MOD_Mesh_Vars,           ONLY: nElems
USE MOD_ReadinTools,         ONLY: GETLOGICAL
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
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR1...'

doWeakLifting=GETLOGICAL('doWeakLifting','.FALSE.')
IF(.NOT.doWeakLifting)&
  doConservativeLifting=GETLOGICAL('doConservativeLifting','.FALSE.')

! We store the interior gradients at the each element face
ALLOCATE(gradUx_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUx_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
gradUx_slave=0.
gradUy_slave=0.
gradUz_slave=0.
gradUx_master=0.
gradUy_master=0.
gradUz_master=0.

! The gradients of the conservative variables are stored at each volume integration point
ALLOCATE(gradUx(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUy(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUz(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
gradUx=0.
gradUy=0.
gradUz=0.


LiftingInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LIFTING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLifting



!==================================================================================================================================
!> \brief Computes the DG gradients using the BR1 scheme in x/y/z direction.
!>
!> This routine will be called after the current solution U has been prolonged to the element faces.
!> To compute the lifted gradients of the solution the following steps are taken:
!> - Compute and communicate the surface fluxes on the mpi interface, do this first to use latency hiding. The fluxes will be
!>   temporarily stored in the gradUxyz_master arrays. Surface fluxes are different if strong or weak form is used.   
!> - Compute the volume integral. The gradients will also get nullified in this routine.
!>   There are different versions of the VolInt routine depending on the usage of the conservative (weak or strong) 
!>   or non conservative (strong only) form. 
!> - The surface fluxes for all remaining sides (boundaries and inner) will be computed. 
!> - The surface integral is performed, first for all inner and boundary sides. Then the communication of the fluxes is finished
!>   and the surface integral of the remaining sides is performed.
!>   In the surface integral there is a distinction between weak and strong formulation. The fluxes for the strong form are
!>   different on the slave or master sides since we substract the inner solution, this is accounted for in the SurfInt routine. 
!> - The gradients are transformed back to physical space to be used by the DG routines.
!> - The computed volume gradients are prolonged to the surfaces at the end of the routine.
!==================================================================================================================================
SUBROUTINE Lifting(UPrim,UPrim_master,UPrim_slave,t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars
USE MOD_Lifting_VolInt,     ONLY: Lifting_VolInt
USE MOD_Lifting_FillFlux,   ONLY: Lifting_FillFlux,Lifting_FillFlux_BC,Lifting_FillFlux_NormVec
USE MOD_DG_Vars,            ONLY: L_hatMinus,L_hatPlus
USE MOD_Lifting_SurfInt,    ONLY: Lifting_SurfInt
USE MOD_ProlongToFacePrim,  ONLY: ProlongToFacePrim
USE MOD_ApplyJacobianPrim,  ONLY: ApplyJacobianPrim
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_FillMortarPrim,     ONLY: U_MortarPrim,Flux_MortarPrim
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI,                ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_Mesh_Vars,          ONLY: nSides,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim(  PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution vector for which lifted gradients will be computed
REAL,INTENT(INOUT) :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< solution on the master sides
REAL,INTENT(INOUT) :: UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< solution on the slave sides
REAL,INTENT(IN)    :: t                                                 !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! fill the global surface flux list
! #### use gradUz_slave for storing the fluxes, NormVec is applied later ####
! ## DO NOT USE SAME STORAGE for transformed/untransformed fluxes, since NormVec can be applied before communication is finished ##
#if USE_MPI
! Receive YOUR
CALL StartReceiveMPIData(gradUz_slave,DataSizeSidePrim,1,nSides,MPIRequest_Flux(:,RECV),SendID=1)
! Compute lifting MPI fluxes
CALL Lifting_FillFlux(   UPrim_master,UPrim_slave,gradUz_slave,doMPISides=.TRUE.)
! Start Send MINE
CALL StartSendMPIData(   gradUz_slave,DataSizeSidePrim,1,nSides,MPIRequest_Flux(:,SEND),SendID=1)
#endif /*USE_MPI*/


! compute volume integral contribution and add to ut
IF(doWeakLifting.OR.doConservativeLifting)THEN
  CALL Lifting_VolInt(1,UPrim,GradUx)
  CALL Lifting_VolInt(2,UPrim,GradUy)
#if (PP_dim==3)
  CALL Lifting_VolInt(3,UPrim,GradUz)
#endif
ELSE
  CALL Lifting_VolInt(UPrim,GradUx,GradUy,GradUz)
END IF

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC(t,UPrim_master,                  gradUz_slave)
CALL Lifting_FillFlux(     UPrim_master,UPrim_slave,      gradUz_slave,doMPISides=.FALSE.)
! at this point BC, inner and MPI MINE are filled
CALL Lifting_FillFlux_NormVec(gradUz_slave,gradUx_master,gradUy_master,gradUz_master,doMPISides=.FALSE.)

! Attention: we only have one Flux (gradUx/y/z_master) for the Lifting 
!            => input it to Flux_Mortar for both fluxes (master/slave)
CALL Flux_MortarPrim(gradUx_master,gradUx_master,doMPISides=.FALSE.,weak=doWeakLifting)
CALL Flux_MortarPrim(gradUy_master,gradUy_master,doMPISides=.FALSE.,weak=doWeakLifting)
#if (PP_dim==3)
CALL Flux_MortarPrim(gradUz_master,gradUz_master,doMPISides=.FALSE.,weak=doWeakLifting)
#endif


! compute surface integral contribution and add to ut
CALL Lifting_SurfInt(PP_N,gradUx_master,gradUx,.FALSE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
CALL Lifting_SurfInt(PP_N,gradUy_master,gradUy,.FALSE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#if (PP_dim==3)
CALL Lifting_SurfInt(PP_N,gradUz_master,gradUz,.FALSE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#endif
#if USE_MPI
! Complete send / receive
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux)
CALL Lifting_FillFlux_NormVec(gradUz_slave,gradUx_master,gradUy_master,gradUz_master,doMPISides=.TRUE.)
! Attention: we only have one Flux (gradUx/y/z_master) for the Lifting 
!            => input it to Flux_Mortar for both fluxes (master/slave)
CALL Flux_MortarPrim(gradUx_master,gradUx_master,doMPISides=.TRUE.,weak=doWeakLifting)
CALL Flux_MortarPrim(gradUy_master,gradUy_master,doMPISides=.TRUE.,weak=doWeakLifting)
#if (PP_dim==3)
CALL Flux_MortarPrim(gradUz_master,gradUz_master,doMPISides=.TRUE.,weak=doWeakLifting)
#endif
! compute surface integral contribution and add to ut
CALL Lifting_SurfInt(PP_N,gradUx_master,gradUx,.TRUE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
CALL Lifting_SurfInt(PP_N,gradUy_master,gradUy,.TRUE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#if (PP_dim==3)
CALL Lifting_SurfInt(PP_N,gradUz_master,gradUz,.TRUE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#endif
#endif /*USE_MPI*/

! Account for the jacobian
! The Lifting already has the right sign
! For FV elements no Jacobian is needed here, since already applied in Lifting_Volint (hidden in FV_Metrics_f/g/hTilde_sJ)
CALL ApplyJacobianPrim(gradUx,toPhysical=.TRUE.,FVE=0)
CALL ApplyJacobianPrim(gradUy,toPhysical=.TRUE.,FVE=0)
#if (PP_dim==3)
CALL ApplyJacobianPrim(gradUz,toPhysical=.TRUE.,FVE=0)
#endif

! We need the gradients at the face of the grid cells
#if USE_MPI
! Prolong to face for MPI sides - send direction
CALL StartReceiveMPIData(gradUx_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=2)
CALL StartReceiveMPIData(gradUy_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=2)
#if (PP_dim==3)
CALL StartReceiveMPIData(gradUz_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=2)
#endif

CALL ProlongToFacePrim(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_MortarPrim(gradUx_master,gradUx_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   gradUx_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=2)
CALL ProlongToFacePrim(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_MortarPrim(gradUy_master,gradUy_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   gradUy_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=2)
#if (PP_dim==3)
CALL ProlongToFacePrim(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_MortarPrim(gradUz_master,gradUz_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(   gradUz_slave,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=2)
#endif
#endif /*USE_MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFacePrim(PP_N,gradUx,gradUx_master,gradUx_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL ProlongToFacePrim(PP_N,gradUy,gradUy_master,gradUy_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#if (PP_dim==3)
CALL ProlongToFacePrim(PP_N,gradUz,gradUz_master,gradUz_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#endif
CALL U_MortarPrim(gradUx_master,gradUx_slave,doMPISides=.FALSE.)
CALL U_MortarPrim(gradUy_master,gradUy_slave,doMPISides=.FALSE.)
#if (PP_dim==3)
CALL U_MortarPrim(gradUz_master,gradUz_slave,doMPISides=.FALSE.)
#endif

END SUBROUTINE Lifting

!==================================================================================================================================
!> Deallocate BR1 arrays (volume and surface gradients and gradient fluxes)
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
LiftingInitIsDone = .FALSE.
END SUBROUTINE FinalizeLifting

END MODULE MOD_Lifting
#endif /*PARABOLIC*/
