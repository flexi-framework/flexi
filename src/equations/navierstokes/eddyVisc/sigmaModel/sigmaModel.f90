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

!===================================================================================================================================
!> Subroutines necessary for calculating SigmaModel Eddy-Viscosity
!===================================================================================================================================
MODULE MOD_SigmaModel
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitSigmaModel
   MODULE PROCEDURE InitSigmaModel
END INTERFACE

INTERFACE FinalizeSigmaModel
   MODULE PROCEDURE FinalizeSigmaModel
END INTERFACE

PUBLIC::InitSigmaModel,SigmaModel,SigmaModel_surf,FinalizeSigmaModel
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get some parameters needed by SigmaModel modules and initialize SigmaModel
!===================================================================================================================================
SUBROUTINE InitSigmaModel()
! MODULES
USE MOD_Globals
USE MOD_PreProc                
USE MOD_EddyVisc_Vars              
USE MOD_ReadInTools       ,   ONLY:GETREAL,GETLOGICAL
USE MOD_Interpolation_Vars,   ONLY:InterpolationInitIsDone
USE MOD_Mesh_Vars         ,   ONLY:MeshInitIsDone,nElems
USE MOD_Interpolation_Vars,   ONLY:wGP
USE MOD_Mesh_Vars,            ONLY:sJ,nSides
USE MOD_Mesh_Vars,            ONLY:ElemToSide
USE MOD_Testcase_Vars,        ONLY:testcase
#if USE_MPI
USE MOD_MPI,                  ONLY:StartReceiveMPIData,FinishExchangeMPIData,StartSendMPIData
USE MOD_MPI_Vars,             ONLY:MPIRequest_DeltaS,nNbProcs
#endif /*USE_MPI*/ 
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iElem,j,k
REAL    :: CellVol
INTEGER :: iLocSide,SideID,FlipID
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.SigmaModelInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitSigmaModel not ready to be called or already called.")
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SigmaModel...'

! Read the variables used for LES model
! SigmaModel model
CS     = GETREAL('CS')
IF(testcase.EQ."channel") THEN
  ! Do Van Driest style damping or not
  VanDriest = GETLOGICAL('VanDriest','.FALSE.')
END IF

! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )

DO iElem=1,nElems                                        
  CellVol = 0.
  DO i=0,PP_N
    DO j=0,PP_N
      DO k=0,PP_N
        CellVol = CellVol +wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
      END DO
    END DO
  END DO
  DeltaS(iElem) = ( CellVol)**(1./3.)  / (REAL(PP_N)+1.)
  DO iLocSide=1,6
     SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
     FlipID=ElemToSide(E2S_FLIP,iLocSide,iElem) 
     IF(FlipID.EQ.0) THEN
       DeltaS_master(SideID)=DeltaS(iElem)
     ELSE
       DeltaS_slave(SideID)=DeltaS(iElem)
     END IF
  END DO
END DO
#if USE_MPI
! Send YOUR - receive MINE
CALL StartReceiveMPIData(DeltaS_slave, 1, 1,nSides,MPIRequest_DeltaS( :,SEND),SendID=1)
CALL StartSendMPIData(   DeltaS_slave, 1, 1,nSides,MPIRequest_DeltaS( :,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_DeltaS ) !Send MINE -receive YOUR
#endif /*USE_MPI*/

SigmaModelInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SigmaModel DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitSigmaModel

!===================================================================================================================================
!> Compute SigmaModel Eddy-Visosity at a given point in the volume
!===================================================================================================================================
SUBROUTINE SigmaModel(iElem,i,j,k,muSGS)
! MODULES
USE MOD_PreProc
USE MOD_EddyVisc_Vars,     ONLY:deltaS,CS,muSGSmax,SGS_Ind
USE MOD_Lifting_Vars,      ONLY:gradUx,gradUy,gradUz
USE MOD_DG_Vars,           ONLY:U
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem             !< index of current element
!> indices of the current volume point
INTEGER,INTENT(IN)                        :: i,j,k
!> gradients of the velocities w.r.t. all directions
REAL,INTENT(INOUT)                        :: muSGS             !< local SGS viscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYEV
! LOCAL VARIABLES
REAL               :: sigma1,sigma2,sigma3
REAL               :: pi
REAL               :: G_mat(3,3)
REAL               :: lambda(3)
REAL               :: work(9)  !lapack work array
INTEGER            :: info
REAL               :: d_model
!===================================================================================================================================
pi=acos(-1.)
G_Mat(1,1) = gradUx(2,i,j,k,iElem)*gradUx(2,i,j,k,iElem) + gradUx(3,i,j,k,iElem)*gradUx(3,i,j,k,iElem) &
           + gradUx(4,i,j,k,iElem)*gradUx(4,i,j,k,iElem)  
G_Mat(1,2) = gradUx(2,i,j,k,iElem)*gradUy(2,i,j,k,iElem) + gradUx(3,i,j,k,iElem)*gradUy(3,i,j,k,iElem) &
           + gradUx(4,i,j,k,iElem)*gradUy(4,i,j,k,iElem)  
G_Mat(1,3) = gradUx(2,i,j,k,iElem)*gradUz(2,i,j,k,iElem) + gradUx(3,i,j,k,iElem)*gradUz(3,i,j,k,iElem) &
           + gradUx(4,i,j,k,iElem)*gradUz(4,i,j,k,iElem)  
G_Mat(2,1) = G_Mat(1,2)  
G_Mat(2,2) = gradUy(2,i,j,k,iElem)*gradUy(2,i,j,k,iElem) + gradUy(3,i,j,k,iElem)*gradUy(3,i,j,k,iElem) &
           + gradUy(4,i,j,k,iElem)*gradUy(4,i,j,k,iElem)  
G_Mat(2,3) = gradUy(2,i,j,k,iElem)*gradUz(2,i,j,k,iElem) + gradUy(3,i,j,k,iElem)*gradUz(3,i,j,k,iElem) &
           + gradUy(4,i,j,k,iElem)*gradUz(4,i,j,k,iElem)  
G_Mat(3,1) = G_Mat(1,3)  
G_Mat(3,2) = G_Mat(2,3)  
G_Mat(3,3) = gradUz(2,i,j,k,iElem)*gradUz(2,i,j,k,iElem) + gradUz(3,i,j,k,iElem)*gradUz(3,i,j,k,iElem) &
           + gradUz(4,i,j,k,iElem)*gradUz(4,i,j,k,iElem)  
!LAPACK
CALL DSYEV('N','U',3,G_Mat,3,lambda,work,9,info)
IF(info .NE. 0) THEN
  WRITE(*,*)'Eigenvalue Computation failed 3D',info
  d_model=0.
ELSEIF(ANY(lambda.LE.0.))THEN
  sigma1 = SQRT(MAX(0.,lambda(3)))
  sigma2 = SQRT(MAX(0.,lambda(2)))
  sigma3 = SQRT(MAX(0.,lambda(1)))
  d_model = (sigma3*(sigma1-sigma2)*(sigma2-sigma3))/(sigma1**2)
ELSE 
  sigma1 = SQRT(lambda(3))
  sigma2 = SQRT(lambda(2))
  sigma3 = SQRT(lambda(1))
  d_model = (sigma3*(sigma1-sigma2)*(sigma2-sigma3))/(sigma1**2)
END IF
! SigmaModel model
muSGS= (CS*deltaS(iElem))**2. * d_model*U(1,i,j,k,iElem)
SGS_Ind(2,i,j,k,iElem) = muSGS
muSGSmax(iElem) = MAX(muSGS,muSGSmax(iElem))
END SUBROUTINE SigmaModel

!===================================================================================================================================
!> Compute SigmaModel Eddy-Visosity at a given point at the surface
!===================================================================================================================================
SUBROUTINE SigmaModel_surf(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,DeltaSS,SGS_Ind,muSGS,Face_xGP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!> gradients of the velocities w.r.t. all directions
REAL,INTENT(IN)                           :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32
REAL,INTENT(IN)                           :: rho               !< Density
REAL,INTENT(IN)                           :: DeltaSS           !< Filter width
REAL,INTENT(IN)                           :: SGS_Ind           !< Indicator for SGS model
REAL,INTENT(IN)                           :: Face_xGP          !< Coordinate for van-Driest damping
REAL,INTENT(OUT)                          :: muSGS             !< local SGS viscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
RETURN!use prolonged from volume
END SUBROUTINE SigmaModel_surf

!===============================================================================================================================
!> Deallocate arrays and finalize variables used by SigmaModel SGS model
!===============================================================================================================================
SUBROUTINE FinalizeSigmaModel()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!===============================================================================================================================
SigmaModelInitIsDone = .FALSE.
END SUBROUTINE FinalizeSigmaModel

END MODULE MOD_SigmaModel
