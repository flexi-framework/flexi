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
!> \brief Contains the surface integral routine for the BR1 lifting operation.
!>
!> This module contains the subroutine that performs the surface integral operation for the BR1 lifting. The BR1 lifting is
!> impelemented in strong or weak form. The fluxes for the strong form are different for the slave and the master sides.
!==================================================================================================================================
MODULE MOD_Lifting_SurfInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE Lifting_SurfInt
  MODULE PROCEDURE Lifting_SurfInt
END INTERFACE


PUBLIC::Lifting_SurfInt
!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> \brief Surface integral in the BR1 scheme optimized for performance, for weak or strong formulation.
!>
!> Performs the surface integral for the BR1 routine. Uses the DoSurfInt routine from the DG operator to perform actual
!> integration.
!> If we use the strong formulation, the inner solution is substracted from the numerical flux. This is done in the FillFlux
!> routines for the master side. The flux on the master side is \f$ \frac{1}{2} (U^+ + U^-) - U^- = \frac{1}{2} (U^+ - U^-) \f$
!> since \f$ \frac{1}{2} (U^+ + U^-) \f$ is the numerical flux in the BR1 scheme. On the slave side, the flux becomes
!> \f$ \frac{1}{2} (U^+ + U^-) - U^+ = \frac{1}{2} (-U^+ + U^-) \f$ which is simply the master side flux multiplied by
!> \f$ -1 \f$. This means we don't have to flip the sign on the flux for the slave side in strong form as we normally do to
!> get the flux on the slave side.
!==================================================================================================================================
SUBROUTINE Lifting_SurfInt(Nloc,Flux,gradU,doMPISides,L_HatMinus,L_HatPlus,weak)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SurfintPrim,        ONLY: DoSurfIntPrim
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides,nElems
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_MINE
USE MOD_Mesh_Vars,          ONLY: S2V3,CS2V2
USE MOD_Mesh_Vars,          ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems_master,FV_Elems_slave
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                             !< polynomial degree
LOGICAL,INTENT(IN) :: doMPISides                                       !< = .TRUE. only MPISides_YOUR+MPIMortar are filled
                                                                       !< =.FALSE. BCSides+(Mortar-)InnerSides+MPISides_MINE
REAL,INTENT(IN)    :: Flux(1:PP_nVarPrim,0:Nloc,0:PP_NlocZ,nSides)         !< flux to be filled
REAL,INTENT(IN)    :: L_HatPlus(0:Nloc)                                !< lagrange polynomials at xi=+1 and pre-divided by
                                                                       !< integration weight
REAL,INTENT(IN)    :: L_HatMinus(0:Nloc)                               !< lagrange polynomials at xi=-1 and pre-divided by
                                                                       !< integration weight
REAL,INTENT(INOUT) :: gradU(PP_nVarPrim,0:Nloc,0:Nloc,0:PP_NlocZ,1:nElems) !< time derivative of solution
LOGICAL,INTENT(IN) :: weak                                             !< switch for weak or strong formulation
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,nbElemID,locSideID,nblocSideID,SideID,p,q,flip
INTEGER            :: firstSideID,lastSideID
REAL               :: FluxTmp(1:PP_nVarPrim,0:Nloc,0:PP_NlocZ)
!==================================================================================================================================
IF(doMPISides)THEN
  ! MPI YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  ! inner sides and MPI mine
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ElemID      = SideToElem(S2E_ELEM_ID,   SideID)
  nbElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)

  ! master sides
  IF(ElemID.GT.0)THEN
    IF (FV_Elems_master(SideID).EQ.0) THEN ! DG element
      locSideID   = SideToElem(S2E_LOC_SIDE_ID,SideID)
      ! orient flux to fit flip and locSide to element local system
      DO q=0,PP_NlocZ; DO p=0,Nloc
        FluxTmp(:,p,q)=Flux(:,CS2V2(1,p,q,locSideID),CS2V2(2,p,q,locSideID),SideID)
      END DO; END DO ! p,q
#if   (PP_NodeType==1)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus,   L_HatPlus,      locSideID,gradU(:,:,:,:,ElemID))
#elif (PP_NodeType==2)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),locSideID,gradU(:,:,:,:,ElemID))
#endif
    END IF
  END IF

  ! slave sides
  IF(nbElemID.GT.0)THEN
    IF (FV_Elems_slave(SideID).EQ.0) THEN ! DG element
      nblocSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
      flip        = SideToElem(S2E_FLIP,SideID)
      ! orient flux to fit flip and locSide to element local system
      IF(weak)THEN
        DO q=0,PP_NlocZ; DO p=0,Nloc
          FluxTmp(:,p,q)=-Flux(:,S2V3(1,p,q,flip,nblocSideID),S2V3(2,p,q,flip,nblocSideID),SideID)
        END DO; END DO ! p,q
      ELSE
        ! In strong form, don't flip the sign since the slave flux is the negative of the master flux
        DO q=0,PP_NlocZ; DO p=0,Nloc
          FluxTmp(:,p,q)= Flux(:,S2V3(1,p,q,flip,nblocSideID),S2V3(2,p,q,flip,nblocSideID),SideID)
        END DO; END DO ! p,q
      END IF
#if   (PP_NodeType==1)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus,   L_HatPlus,      nblocSideID,gradU(:,:,:,:,nbElemID))
#elif (PP_NodeType==2)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),nblocSideID,gradU(:,:,:,:,nbElemID))
#endif
    END IF
  END IF
END DO ! SideID=1,nSides
END SUBROUTINE Lifting_SurfInt

END MODULE MOD_Lifting_SurfInt
#endif /*PARABOLIC*/
