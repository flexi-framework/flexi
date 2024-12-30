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
!> Module for the Finite Volume sub-cells shock capturing.
!>
!> DG elements, that are detected to contain a shock/high gradients/oscillations/..., can be switched to a Finite Volume scheme.
!> A DG element of polynomial degree N is subdivided into (N+1)^dim sub-cells (to each Gauss Point/DOF one FV sub-cell).
!> The FV sub-cells of such an element are updated using FV method with 2nd order TVD reconstruction (slope limiters).
!==================================================================================================================================
MODULE MOD_FV_Switching
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

#if FV_ENABLED == 1
PUBLIC::FV_Switch
PUBLIC::FV_ProlongFVElemsToFace
PUBLIC::FV_FillIni
PUBLIC::FV_Info
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Perform switching between DG element and FV sub-cells element (and vise versa) depending on the indicator value.
!> Optionally, the switching process can be done in the reference element to guarantee conservation on non-cartesian elements.
!==================================================================================================================================
SUBROUTINE FV_Switch(U,U2,U3,AllowToDG)
! MODULES
USE MOD_PreProc
USE MOD_Indicator_Vars  ,ONLY: IndValue
USE MOD_Indicator       ,ONLY: IndPersson
USE MOD_FV_Vars
USE MOD_Analyze
USE MOD_Mesh_Vars       ,ONLY: nElems, sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< state vector to be switched
REAL,INTENT(INOUT),OPTIONAL :: U2(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< optional additional state vector to be switched
REAL,INTENT(INOUT),OPTIONAL :: U3(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< optional additional state vector to be switched
LOGICAL,INTENT(IN)          :: AllowToDG                                  !< if .TRUE. FV element is allowed to switch to DG
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: U_DG(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL    :: ind
INTEGER :: iElem
!==================================================================================================================================
DO iElem=1,nElems
  IF (FV_Elems(iElem).EQ.0) THEN ! DG Element
    ! Switch DG to FV Element, if Indicator is higher then IndMin
    IF (IndValue(iElem).GT.FV_IndUpperThreshold) THEN
      ! switch Element to FV
      FV_Elems(iElem) = 1
      CALL FV_InterpolateDG2FV(U(:,:,:,:,iElem),sJ(:,:,:,iElem,0:FV_SIZE))
      IF (PRESENT(U2)) CALL FV_InterpolateDG2FV(U2(:,:,:,:,iElem),sJ(:,:,:,iElem,0:FV_SIZE))
      IF (PRESENT(U3)) CALL FV_InterpolateDG2FV(U3(:,:,:,:,iElem),sJ(:,:,:,iElem,0:FV_SIZE))
    END IF
  ELSE ! FV Element
    ! Switch FV to DG Element, if Indicator is lower then IndMax
    IF ((IndValue(iElem).LT.FV_IndLowerThreshold).AND.AllowToDG) THEN
      U_DG = U(:,:,:,:,iElem)
      CALL FV_InterpolateFV2DG(U_DG(:,:,:,:),sJ(:,:,:,iElem,0:FV_SIZE))
      IF (FV_toDG_indicator) THEN
        ind = IndPersson(U_DG(:,:,:,:))
        IF (ind.GT.FV_toDG_limit) CYCLE
      END IF
      ! switch Element to DG
      FV_Elems(iElem)  = 0
      U(:,:,:,:,iElem) = U_DG
      IF (PRESENT(U2)) CALL FV_InterpolateFV2DG(U2(:,:,:,:,iElem),sJ(:,:,:,iElem,0:FV_SIZE))
      IF (PRESENT(U3)) CALL FV_InterpolateFV2DG(U3(:,:,:,:,iElem),sJ(:,:,:,iElem,0:FV_SIZE))
    END IF
  END IF
END DO !iElem
! collect statistics
FV_Elems_counter  = FV_Elems_counter  + FV_Elems
FV_Switch_counter = FV_Switch_counter + 1
FV_Elems_Amount   = REAL(FV_Elems_Counter)/FV_Switch_counter

CALL FV_ProlongFVElemsToFace()
END SUBROUTINE FV_Switch


!==================================================================================================================================
!> Set FV_Elems_slave and FV_Elems_master information
!==================================================================================================================================
SUBROUTINE FV_ProlongFVElemsToFace()
! MODULES
USE MOD_FV_Vars         ,ONLY: FV_Elems,FV_Elems_master,FV_Elems_slave
USE MOD_Mesh_Vars       ,ONLY: SideToElem,nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iSide,ElemID,nbElemID
!==================================================================================================================================
! array not allocated in postiMode
IF (.NOT.ALLOCATED(SideToElem)) RETURN

! set information whether elements adjacent to a side are DG or FV elements
DO iSide = 1,nSides
  ElemID    = SideToElem(S2E_ELEM_ID   ,iSide)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID,iSide)
  !master sides
  IF(ElemID  .GT.0) FV_Elems_master(iSide) = FV_Elems(ElemID)
  !slave side (ElemID,locSide and flip =-1 if not existing)
  IF(nbElemID.GT.0) FV_Elems_slave( iSide) = FV_Elems(nbElemID)
END DO
END SUBROUTINE FV_ProlongFVElemsToFace


!==================================================================================================================================
!> Interpolate solution from DG representation to FV subcells.
!> Interpolation is done either conservatively in reference space or non-conservatively in phyiscal space.
!==================================================================================================================================
PPURE SUBROUTINE FV_InterpolateDG2FV(U_In,sJ_In)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars          ,ONLY: switchConservative,FV_Vdm
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisVolume
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT) ::  U_In(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)      !< state vector to be switched from DG to FV representation
REAL,INTENT(IN)    :: sJ_In(0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE)    !< inverse of Jacobian determinant at each Gauss point
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k
!==================================================================================================================================
IF (switchConservative) THEN
  ! Transform the DG solution into the reference element
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    U_In(:,i,j,k)=U_In(:,i,j,k)/sJ_In(i,j,k,0)
  END DO; END DO; END DO
  ! Perform interpolation from DG to FV
  CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,U_In)
  ! Transform back to physical space
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    U_In(:,i,j,k)=U_In(:,i,j,k)*sJ_In(i,j,k,1)
  END DO; END DO; END DO
ELSE
  CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,U_In)
END IF
END SUBROUTINE FV_InterpolateDG2FV


!==================================================================================================================================
!> Interpolate solution from FV subcell representation to DG.
!> Interpolation is done either conservatively in reference space or non-conservatively in phyiscal space.
!==================================================================================================================================
PPURE SUBROUTINE FV_InterpolateFV2DG(U_In,sJ_In)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars          ,ONLY: switchConservative,FV_sVdm
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisVolume
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT) ::  U_In(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)      !< state vector to be switched from FV to DG representation
REAL,INTENT(IN)    :: sJ_In(0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE)    !< inverse of Jacobian determinant at each Gauss point
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k
!==================================================================================================================================
IF (switchConservative) THEN
  ! Transform the FV solution into the reference element
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    U_In(:,i,j,k)=U_In(:,i,j,k)/sJ_In(i,j,k,1)
  END DO; END DO; END DO
  ! Perform interpolation from FV to DG
  CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_sVdm,U_In)
  ! Transform back to physical space
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    U_In(:,i,j,k)=U_In(:,i,j,k)*sJ_In(i,j,k,0)
  END DO; END DO; END DO
ELSE
  CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_sVdm,U_In)
END IF
END SUBROUTINE FV_InterpolateFV2DG


!==================================================================================================================================
!> Print information on the amount of FV subcells
!==================================================================================================================================
SUBROUTINE FV_Info(iter)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars    ,ONLY: nGlobalElems
USE MOD_Analyze_Vars ,ONLY: totalFV_nElems
USE MOD_FV_Vars      ,ONLY: FV_Elems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN) :: iter !< number of iterations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF (iter.EQ.1_8) totalFV_nElems = SUM(FV_Elems) ! counter for output of FV amount during analyze
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,totalFV_nElems,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  ! totalFV_nElems is counted in PrintStatusLine
ELSE
  CALL MPI_REDUCE(totalFV_nElems,0           ,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif /*USE_MPI*/
SWRITE(UNIT_stdOut,'(A,F8.3,A)')' FV amount %: ', REAL(totalFV_nElems) / REAL(nGlobalElems) / iter*100
totalFV_nElems = 0
END SUBROUTINE FV_Info


!==================================================================================================================================
!> Initialize all FV elements and overwrite data of DG FillIni. Each subcell is supersampled with PP_N points in each space
!> dimension and the mean value is taken as value for this subcell.
!==================================================================================================================================
SUBROUTINE FV_FillIni()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis             ,ONLY: InitializeVandermonde
USE MOD_ChangeBasis       ,ONLY: ChangeBasis2D_XYZ, ChangeBasis3D_XYZ
USE MOD_ChangeBasisByDim  ,ONLY: ChangeBasisVolume
USE MOD_DG_Vars           ,ONLY: U
USE MOD_Equation_Vars     ,ONLY: IniExactFunc
USE MOD_Exactfunc         ,ONLY: ExactFunc
USE MOD_FV_Vars           ,ONLY: FV_Elems,FV_Vdm,FV_CellType,FV_IniSharp,FV_IniSupersample
USE MOD_FV_Basis          ,ONLY: FV_Build_X_w_BdryX
USE MOD_Interpolation     ,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars,ONLY: NodeType,NodeTypeVISUInner
USE MOD_Mesh_Vars         ,ONLY: nElems
USE MOD_Mesh_Vars         ,ONLY: Elem_xGP
USE MOD_ReadInTools       ,ONLY: GETLOGICAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: i,iElem, j,k,ii,jj,kk,iVar
REAL                   :: FV_w(0:PP_N),FV_BdryX(0:PP_N+1)
REAL,DIMENSION(0:PP_N) :: FV_X,xGP,wGP,wBary,SubxGP
REAL                   :: VDM(0:PP_N,0:PP_N,0:PP_N)
REAL,ALLOCATABLE       :: xx(:,:,:,:)
REAL                   :: tmp(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                   :: Elem_xFV(1:3,0:PP_N,0:PP_N,0:PP_NZ)
!===================================================================================================================================
! initial call of indicator
FV_Elems = 0
! Switch DG elements to FV if necessary (converts initial DG solution to FV solution)
CALL FV_Switch(U,AllowToDG=.FALSE.)

IF (.NOT.FV_IniSharp) THEN
  ! Super sample initial solution of all FV elements. Necessary if already initial DG solution contains oscillations, which
  ! may lead to non valid solutions inside a sub-cell.!
  !!! THIS IS EXPENSIVE !!!
  IF (FV_IniSupersample) THEN

    CALL GetNodesAndWeights(PP_N,NodeType,xGP,wGP,wBary)
    CALL FV_Build_X_w_BdryX(PP_N,FV_X,FV_w,FV_BdryX,FV_CellType)
    DO i=0,PP_N
      ! compute equidistant supersampling points inside FV sub-cell
      DO j=0,PP_N
        SubxGP(j) = FV_BdryX(i) + (j+0.5)/(PP_N+1)*FV_w(i)
      END DO
      ! build Vandermonde for mapping the whole interval [-1,1] to the i-th FV subcell
      CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,SubxGP,VDM(:,:,i))
    END DO


    ALLOCATE(xx(1:3,0:PP_N,0:PP_N,0:PP_NZ)) ! coordinates supersampled to FV subcell
    DO iElem=1,nElems
      IF (FV_Elems(iElem).EQ.0) CYCLE ! DG element
      DO k=0,PP_NZ
        DO j=0,PP_N
          DO i=0,PP_N
            ! supersample coordinates to i,j,k-th subcells
#if PP_dim == 3
            CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm(:,:,i),Vdm(:,:,j),Vdm(:,:,k),Elem_xGP(1:3,:,:,:,iElem),xx)
#else
            CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm(:,:,i),Vdm(:,:,j),Elem_xGP(1:3,:,:,0,iElem),xx(:,:,:,0))
#endif
            ! evaluate ExactFunc for all supersampled points of subcell (i,j,k)
            DO kk=0,PP_NZ; DO jj=0,PP_N; DO ii=0,PP_N
              CALL ExactFunc(IniExactFunc,0.,xx(1:3,ii,jj,kk),tmp(:,ii,jj,kk))
            END DO; END DO; END DO
            ! mean value
            DO iVar=1,PP_nVar
              U(iVar,i,j,k,iElem) = SUM(tmp(iVar,:,:,:)) / ((PP_N+1)**2*(PP_NZ+1))
            END DO
          END DO ! i
        END DO ! j
      END DO !k
    END DO ! iElem=1,nElems
    DEALLOCATE(xx)
  END IF
ELSE
  ! maintain a sharp interface in the FV region
  DO iElem=1,nElems
    IF (FV_Elems(iElem).EQ.0) CYCLE ! DG element
    ! get coordinates of the FV elements
    CALL ChangeBasisVolume(3,PP_N,PP_N,FV_Vdm,Elem_xGP(1:3,:,:,:,iElem),Elem_xFV(1:3,:,:,:))
    DO k=0,PP_NZ
      DO j=0,PP_N
        DO i=0,PP_N
          CALL ExactFunc(IniExactFunc,0.,Elem_xFV(:,i,j,k),U(1:PP_nVar,i,j,k,iElem))
        END DO ! i
      END DO ! j
    END DO !k
  END DO ! iElem=1,nElems
END IF

END SUBROUTINE FV_FillIni

#endif /*FV_ENABLED == 1*/
END MODULE MOD_FV_Switching
