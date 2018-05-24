!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
!> Containes the routines that will initialize and finalize the walldistance routine as well as the main routine used in the wall
!> distance calculation.
!===================================================================================================================================
MODULE MOD_Walldistance
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitWalldistance
  MODULE PROCEDURE InitWalldistance
END INTERFACE

INTERFACE CalcWalldistance
  MODULE PROCEDURE CalcWalldistance
END INTERFACE

INTERFACE FinalizeWalldistance
  MODULE PROCEDURE FinalizeWalldistance
END INTERFACE

PUBLIC:: InitWalldistance,CalcWalldistance,FinalizeWalldistance

CONTAINS

!===================================================================================================================================
!> Read in user defined parameters and prepare data for walldistance.
!===================================================================================================================================
SUBROUTINE InitWalldistance()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PreProc
USE MOD_Globals
USE MOD_Walldistance_Vars
USE MOD_ReadInTools
USE MOD_StringTools,        ONLY: INTTOSTR
USE MOD_Mesh_Vars,          ONLY: nBCSides,Face_xGP,nElems
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisSurf
USE MOD_Interpolation,      ONLY: GetVandermonde,GetDerivativeMatrix
USE MOD_Interpolation_Vars, ONLY: NodeType,NodeTypeVisuInner
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER    :: iSide
!===================================================================================================================================
! Use NSuper+1 equidistant points for initial supersampling
NSuper = GETINT('NSuper',INTTOSTR(3*PP_N))
DebugVisu = GETLOGICAL('DebugVisu','.FALSE.')
IF (DebugVisu) NVisu  = GETINT('NVisu', INTTOSTR(2*PP_N))

! Vandermonde to interpolate the face coordinates to the supersampling points
ALLOCATE(Vdm_GaussN_EquiNSuper(0:NSuper,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NSuper,NodeTypeVisuInner,Vdm_GaussN_EquiNSuper)

! Get supersampling points
ALLOCATE(xSuper_Face(PP_dim,0:NSuper,0:ZDIM(NSuper),nBCSides))
DO iSide = 1,nBCSides
  CALL ChangeBasisSurf(PP_dim,PP_N,NSuper,Vdm_GaussN_EquiNSuper,Face_xGP(1:PP_dim,:,:,0,iSide),xSuper_Face(1:PP_dim,:,:,iSide))
  WRITE(UNIT_stdOut,'(A,F7.2,A35)',ADVANCE='NO') CHAR(13),REAL(iSide)/REAL(nBCSides)*100., '% of super sampling points created '
END DO
WRITE(UNIT_stdOut,'(A)') 'SUPER SAMPLING POINTS CREATED!'

! Array used to store the nearest supersampling point for each volume point
ALLOCATE(nearestFace(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))

! Output array of walldistance
ALLOCATE(distance(0:PP_N,0:PP_N,0:PP_NZ,nElems))
distance = 0.

! Derivative matrix needed to get gradient
ALLOCATE(D(0:PP_N,0:PP_N))
CALL GetDerivativeMatrix(PP_N,NodeType,D)
END SUBROUTINE InitWalldistance


!===================================================================================================================================
!> Main walldistance routine.
!===================================================================================================================================
SUBROUTINE CalcWalldistance()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PreProc
USE MOD_Globals
USE MOD_Walldistance_Vars
USE MOD_Interpolation_Vars
USE MOD_Mesh_Vars,          ONLY: Elem_xGP,nBCSides,nElems,nGlobalElems,offsetElem,Face_xGP
USE MOD_Mesh_Vars,          ONLY: BoundaryType,BC,MeshFile
USE MOD_HDF5_Output,        ONLY: WriteArray
USE MOD_IO_HDF5
USE MOD_VTK,                ONLY: WriteDataToVTK
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_Interpolation_Vars, ONLY: NodeType,NodeTypeVisu
USE MOD_Interpolation,      ONLY: GetVandermonde
USE MOD_Basis,              ONLY: LagrangeInterpolationPolys
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER    :: i,j,k,iElem,iSide,p,q,l
REAL       :: xi_i(1:PP_dim-1),xi_old(1:PP_dim-1),xVol(PP_dim),xCur(PP_dim)
REAL       :: dist,best
REAL       :: Jac(1:PP_dim,1:PP_dim-1)
REAL       :: dX(1:PP_dim,1:PP_dim-1,0:PP_N,0:PP_NZ)
CHARACTER(LEN=255) :: FileName,FileName_visu
CHARACTER(LEN=255),ALLOCATABLE:: StrVarNames_loc(:)
REAL,ALLOCATABLE,TARGET       :: Coords_NVisu(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET       :: distance_NVisu(:,:,:,:,:)
REAL,POINTER                  :: Coords_NVisu_p(:,:,:,:,:)
REAL,POINTER                  :: distance_NVisu_p(:,:,:,:,:)
REAL                          :: distanceTmp(1,0:PP_N,0:PP_N,0:PP_NZ)
REAL,ALLOCATABLE              :: Vdm_GaussN_NVisu(:,:)
REAL,EXTERNAL                 :: distanceFunc
INTEGER                       :: iter,maxIter=1000
REAL                          :: g(1:PP_dim-1)
REAL,PARAMETER                :: alpha = 0.1
REAL,PARAMETER                :: beta   = 0.1
REAL                          :: eta
REAL                          :: LagXi(0:PP_N),LagEta(0:PP_N)
REAL                          :: percentDone
!===================================================================================================================================
! First step: Coarse search using the supersampled points
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    best = HUGE(1.)
    xVol = Elem_xGP(1:PP_dim,i,j,k,iElem)
    DO iSide = 1, nBCSides
    IF ((BoundaryType(BC(iSide),BC_TYPE).NE.3).AND.(BoundaryType(BC(iSide),BC_TYPE).NE.4)) CYCLE
      DO q=0,ZDIM(NSuper); DO p=0,NSuper
        dist = NORM2(xVol-xSuper_Face(1:PP_dim,p,q,iSide))
        IF (dist.LT.best) THEN
          best = dist
          nearestFace(1,i,j,k,iElem) = p
          nearestFace(2,i,j,k,iElem) = q
          nearestFace(3,i,j,k,iElem) = iSide
        END IF
      END DO; END DO ! p,q=0,PP_N
    END DO ! iSide = 1, nBCSides
  END DO; END DO; END DO! i,j,k=0,PP_N
  WRITE(UNIT_stdOut,'(A,F7.2,A25)',ADVANCE='NO') CHAR(13),REAL(iElem)/REAL(nElems)*100., '% of coarse search done '
END DO ! iElem
WRITE(UNIT_stdOut,'(A)') 'COARSE SEARCH DONE!'

! Second step: fine search on the side that has been found using the coarse search, use a gradient descent with restrictions
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    ! Save the current volume coordinates
    xVol = Elem_xGP(1:PP_dim,i,j,k,iElem)
    ! Calculate the derivative of the face coordinates
    iSide = nearestFace(3,i,j,k,iElem)
    dX=0.
    DO q=0,PP_NZ; DO p=0,PP_N
      DO l=0,PP_N
        ! Matrix-vector multiplication
        dX(:,1,p,q)=dX(:,1,p,q) + D(p,l)*Face_xGP(1:PP_dim,l,q,0,iSide)
#if PP_dim == 3
        dX(:,2,p,q)=dX(:,2,p,q) + D(q,l)*Face_xGP(1:PP_dim,p,l,0,iSide)
#endif
      END DO
    END DO; END DO
    ! Start value from coarse search
    xi_i(1) = 1./REAL(NSuper+1)+2.*REAL(nearestFace(1,i,j,k,iElem))/REAL(NSuper+1) -1.
#if PP_dim == 3
    xi_i(2) = 1./REAL(NSuper+1)+2.*REAL(nearestFace(2,i,j,k,iElem))/REAL(NSuper+1) -1.
#endif
    iter = 0
    DO WHILE  (iter.LT.maxIter)
      ! Evaluate physical positions as well as the Jacobian for all the derivatives
      CALL LagrangeInterpolationPolys(xi_i(1),PP_N,xGP,wBary,LagXi)
#if PP_dim == 3
      CALL LagrangeInterpolationPolys(xi_i(2),PP_N,xGP,wBary,LagEta)
#else
      LagEta = 1.
#endif
      xCur = 0.
      Jac = 0.
      DO q=0,PP_NZ; DO p=0,PP_N
        xCur(:) = xCur(:) + LagXi(p)*LagEta(q)*Face_xGP(:,p,q,0,iSide)
        Jac(:,:) = Jac(:,:) + LagXi(p)*LagEta(q)*dX(:,:,p,q)
      END DO; END DO ! p,q=0,PP_N
      ! Compute the derivatives of the square of the distance function
      g(1) = 2.*(xCur(1)-xVol(1))*Jac(1,1) &
#if PP_dim == 3
           + 2.*(xCur(3)-xVol(3))*Jac(3,1) &
#endif
           + 2.*(xCur(2)-xVol(2))*Jac(2,1)
#if PP_dim == 3
      g(2) = 2.*(xCur(1)-xVol(1))*Jac(1,2) &
           + 2.*(xCur(2)-xVol(2))*Jac(2,2) &
           + 2.*(xCur(3)-xVol(3))*Jac(3,2)
#endif
      ! Step-size calculation by heuristic line search approach
      eta = 1.
      DO WHILE (distSquare(xi_i-eta*g,xVol,Face_xGP(1:PP_dim,:,:,0,iSide)).GT.&
               (distSquare(xi_i,      xVol,Face_xGP(1:PP_dim,:,:,0,iSide))-alpha*eta*(NORM2(g))**2))
        eta = eta*beta
      END DO
      ! Save old xi for abort criterion
      xi_old = xi_i
      ! Perform a gradient descent step
      xi_i = xi_i - eta*g
      ! Restrict the result to [-1,1] (projection onto allowed region)
      xi_i(1)=MIN(1.,MAX(-1.,xi_i(1)))
#if PP_dim == 3
      xi_i(2)=MIN(1.,MAX(-1.,xi_i(2)))
#endif
      iter = iter + 1
      ! Abort criterion
      IF (NORM2(xi_old-xi_i).LT.NORM2(xi_old)*1.E-08) EXIT
    END DO ! iter < maxIter
    ! Save result
    CALL LagrangeInterpolationPolys(xi_i(1),PP_N,xGP,wBary,LagXi)
#if PP_dim == 3
    CALL LagrangeInterpolationPolys(xi_i(2),PP_N,xGP,wBary,LagEta)
#else
    LagEta = 1.
#endif
    xCur = 0.
    DO q=0,PP_NZ; DO p=0,PP_N
      xCur(:) = xCur(:) + LagXi(p)*LagEta(q)*Face_xGP(:,p,q,0,iSide)
    END DO; END DO ! p,q=0,PP_N
#if PP_dim == 3
    distance(i,j,k,iElem) = SQRT((xVol(1)-xCur(1))**2+(xVol(2)-xCur(2))**2+(xVol(3)-xCur(3))**2)
#else
    distance(i,j,k,iElem) = SQRT((xVol(1)-xCur(1))**2+(xVol(2)-xCur(2))**2)
#endif
  END DO; END DO; END DO! i,j,k=0,PP_N
  WRITE(UNIT_stdOut,'(A,F7.2,A19)',ADVANCE='NO') CHAR(13),REAL(iElem)/REAL(nElems)*100., '% of elements done '
END DO ! iElem
WRITE(UNIT_stdOut,'(A)') 'FINE SEARCH DONE!'

! Third step: output of the result
FileName = MeshFile(1:INDEX(MeshFile,'_mesh.h5')-1)//'_walldistance.h5'
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
CALL WriteArray('walldistance',4,&
                (/PP_N+1,PP_N+1,PP_NZ+1,nGlobalElems/),&
                (/PP_N+1,PP_N+1,PP_NZ+1,nElems/),&
                (/0,     0,     0,      offsetElem/),.FALSE.,RealArray=distance)

IF (DebugVisu) THEN
  ! VTU visualization
  ALLOCATE(Vdm_GaussN_NVisu(0:NVisu,0:PP_N))
  CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVISU,Vdm_GaussN_NVisu)

  ALLOCATE(distance_NVisu(1,0:NVisu,0:NVisu,0:ZDIM(NVisu),1:(nElems)))
  distance_NVisu = 0.
  ALLOCATE(Coords_NVisu(1:3,0:NVisu,0:NVisu,0:ZDIM(NVisu),1:(nElems)))

  DO iElem=1,nElems
    ! Create coordinates of visualization points
    CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
    ! Interpolate distance onto visu grid
    distanceTmp(1,:,:,:) = distance(:,:,:,iElem)
    CALL ChangeBasisVolume(1,PP_N,NVisu,Vdm_GaussN_NVisu,distanceTmp,distance_NVisu(1:1,:,:,:,iElem))
  END DO !iElem

  ALLOCATE(StrVarNames_loc(1))
  StrVarNames_loc(1) = 'Walldistance'

  ! Visualize data
  FileName_visu = MeshFile(1:INDEX(MeshFile,'_mesh.h5')-1)//'_walldistance.vtu'
  Coords_NVisu_p => Coords_NVisu
  distance_NVisu_p => distance_NVisu
  CALL WriteDataToVTK(1,NVisu,nElems,StrVarNames_loc,Coords_NVisu_p,distance_NVisu_p,TRIM(FileName_visu),dim=PP_dim,DGFV=0)

  DEALLOCATE(distance_NVisu)
  DEALLOCATE(Coords_NVisu)
  DEALLOCATE(StrVarNames_loc)
  DEALLOCATE(Vdm_GaussN_NVisu)
END IF

END SUBROUTINE CalcWalldistance

!===================================================================================================================================
!> Function to return the square of the distance
!===================================================================================================================================
FUNCTION distSquare(xi_i,xVol,Face_xGP)
! MODULES                                                                                                                          !
USE MOD_PreProc
USE MOD_Globals
USE MOD_Interpolation_Vars
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: xi_i(  1:PP_dim-1)
REAL,INTENT(IN) :: xVol(1:PP_dim)
REAL,INTENT(IN) :: Face_xGP(1:PP_dim,0:PP_N,0:PP_NZ)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: distSquare
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: xCur(1:PP_dim)
INTEGER :: p,q
REAL    :: LagXi(0:PP_N),LagEta(0:PP_N)
!===================================================================================================================================
! Calculate the current position based on the input xi_i
CALL LagrangeInterpolationPolys(xi_i(1),PP_N,xGP,wBary,LagXi)
#if PP_dim == 3
CALL LagrangeInterpolationPolys(xi_i(2),PP_N,xGP,wBary,LagEta)
#else
LagEta = 1.
#endif
xCur = 0.
DO q=0,PP_NZ; DO p=0,PP_N
  xCur(:) = xCur(:) + LagXi(p)*LagEta(q)*Face_xGP(:,p,q)
END DO; END DO ! p,q=0,PP_N
#if PP_dim == 3
distSquare = (xVol(1)-xCur(1))**2+(xVol(2)-xCur(2))**2+(xVol(3)-xCur(3))**2
#else
distSquare = (xVol(1)-xCur(1))**2+(xVol(2)-xCur(2))**2
#endif

END FUNCTION distSquare

!===================================================================================================================================
!> Finalize walldistance variables
!===================================================================================================================================
SUBROUTINE FinalizeWalldistance()
! MODULES                                                                                                                          !
USE MOD_Walldistance_Vars
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(xSuper_Face)
SDEALLOCATE(Vdm_GaussN_EquiNSuper)
SDEALLOCATE(nearestFace)
SDEALLOCATE(distance)
SDEALLOCATE(D)

END SUBROUTINE Finalizewalldistance

END MODULE MOD_Walldistance
