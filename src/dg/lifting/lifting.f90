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
!> Lifting
!==================================================================================================================================
MODULE MOD_Lifting
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!==================================================================================================================================
INTERFACE Lifting
#if PP_Lifting==1
  MODULE PROCEDURE Lifting_BR1
#elif PP_Lifting==2
  MODULE PROCEDURE Lifting_BR2
#endif
END INTERFACE

INTERFACE Lifting_VolInt
  MODULE PROCEDURE Lifting_VolInt_Conservative
  MODULE PROCEDURE Lifting_VolInt_Nonconservative
END INTERFACE
!==================================================================================================================================

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitLifting
  MODULE PROCEDURE InitLifting
END INTERFACE

INTERFACE FinalizeLifting
  MODULE PROCEDURE FinalizeLifting
END INTERFACE

PUBLIC::DefineParametersLifting,InitLifting,FinalizeLifting
PUBLIC::Lifting
!==================================================================================================================================

CONTAINS
#if PP_Lifting==1
#include "lifting_br1.t90"
#elif PP_Lifting==2
#include "lifting_br2.t90"
#endif

#include "lifting_fillflux.t90"
#include "lifting_volint.t90"

!==================================================================================================================================
!> DefineLifting
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
CALL prms%CreateLogicalOption('doWeakLifting',         "Set true to perform lifting in weak form.", '.FALSE.')
CALL prms%CreateLogicalOption('doConservativeLifting', "Set true to compute the volume contribution to the gradients in "//&
                                                       "conservative form, i.e. deriving the solution multiplied by the metric "//&
                                                       "terms instead of deriving the solution and multiplying by the metrics.",&
                                                       '.FALSE.')
#if PP_Lifting==2
CALL prms%CreateRealOption(   'etaBR2',                "Lifting penalty for BR2. Increase improves stability at the cost of "//&
                                                       "performance and reduces jumps between two cells.", '2.')
CALL prms%CreateRealOption(   'etaBR2_wall',           "Lifting penalty for BR2 at wall boundaries. Can be choosen different from"//&
                                                       "to decrease wall velocities.", '-1')
#endif
END SUBROUTINE DefineParametersLifting

!==================================================================================================================================
!> \brief Initialize the BR1 and BR2 lifting: get parameters and allocate the arrays required for the BR1/BR2 lifting procedure.
!>
!> Important parameters:
!> - doConservativeLifting: If true, the volume contribution to the gradients is in conservative form, i.e. the solution is derived
!>   multiplied by the metrics terms
!> - etaBR2: Penalty term for the surface contribution of the BR2 lifting. Note, stability is shown only for \f$ \eta_{BR2} > \f$
!>   number of element faces
!>
!> Important parameters for BR1 (the BR2 is always strong):
!> - doWeakLifting will set the lifting procedure to be performed in weak or strong form
!> - In the strong form, the lifting can be performed in a conservative or non conservative  version
!>   using the doConservativeLifting parameter
!>
!> Default ist the non conservative form since this version has the fastest implementation.
!> Default is etaBR2 = 2
!>
!> Note that the gradient arrays in x/y/z directions in the volume and on the surfaces contain the gradients of the primitive
!> variables, i.e. they must be allocated for PP_nVarPrim variables, i.e. \f$ (\rho, u_1,u_2,u_3,p,T) \f$.
!==================================================================================================================================
SUBROUTINE InitLifting()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars
USE MOD_DG_Vars,              ONLY: DGInitIsDone
USE MOD_Mesh_Vars,            ONLY: nSides,nElems
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
#if PP_Lifting==1
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR1...'
doWeakLifting=GETLOGICAL('doWeakLifting','.FALSE.')
IF(.NOT.doWeakLifting)&
  doConservativeLifting=GETLOGICAL('doConservativeLifting','.FALSE.')

#elif PP_Lifting==2
SWRITE(UNIT_stdOut,'(A)') ' INIT LIFTING WITH BR2 ...'
doWeakLifting=.FALSE.
doConservativeLifting=GETLOGICAL('doConservativeLifting','.FALSE.')
etaBR2=GETREAL('etaBR2','2.')
etaBR2_wall=GETREAL('etaBR2_wall','-1.')
IF(etaBR2_wall .EQ. -1.) etaBR2_wall = etaBR2  !default etaBR2_wall == etaBR2

ALLOCATE(FluxX(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(FluxY(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(FluxZ(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
FluxX=0.
FluxY=0.
FluxZ=0.
#endif

! We store the interior gradients at the each element face
ALLOCATE(gradUx_slave (PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_slave (PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_slave (PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUx_master(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_master(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_master(PP_nVarLifting,0:PP_N,0:PP_NZ,1:nSides))
gradUx_slave=0.
gradUy_slave=0.
gradUz_slave=0.
gradUx_master=0.
gradUy_master=0.
gradUz_master=0.

! The gradients of the conservative variables are stored at each volume integration point
ALLOCATE(gradUx(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUy(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUz(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems))
gradUx=0.
gradUy=0.
gradUz=0.

ALLOCATE(diffFluxX_L(PP_nVar,0:PP_N,0:PP_NZ))
ALLOCATE(diffFluxX_R(PP_nVar,0:PP_N,0:PP_NZ))
ALLOCATE(diffFluxY_L(PP_nVar,0:PP_N,0:PP_NZ))
ALLOCATE(diffFluxY_R(PP_nVar,0:PP_N,0:PP_NZ))
ALLOCATE(diffFluxZ_L(PP_nVar,0:PP_N,0:PP_NZ))
ALLOCATE(diffFluxZ_R(PP_nVar,0:PP_N,0:PP_NZ))
diffFluxX_L=0.
diffFluxX_R=0.
diffFluxY_L=0.
diffFluxY_R=0.
diffFluxZ_L=0.
diffFluxZ_R=0.

LiftingInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LIFTING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLifting

!==================================================================================================================================
!> FinalizeLifting
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
#if PP_Lifting==2
SDEALLOCATE(FluxX)
SDEALLOCATE(FluxY)
SDEALLOCATE(FluxZ)
#endif
SDEALLOCATE(diffFluxX_L)
SDEALLOCATE(diffFluxX_R)
SDEALLOCATE(diffFluxY_L)
SDEALLOCATE(diffFluxY_R)
SDEALLOCATE(diffFluxZ_L)
SDEALLOCATE(diffFluxZ_R)
LiftingInitIsDone = .FALSE.
END SUBROUTINE FinalizeLifting

END MODULE MOD_Lifting
#endif /* PARABOLIC */
