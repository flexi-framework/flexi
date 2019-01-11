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
!> Module to handle overintegration/de-aliasing  operations
!==================================================================================================================================
MODULE MOD_Overintegration
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: OVERINTEGRATIONTYPE_NONE       = 0
INTEGER,PARAMETER :: OVERINTEGRATIONTYPE_CUTOFF     = 1
INTEGER,PARAMETER :: OVERINTEGRATIONTYPE_CONSCUTOFF = 2

!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitOverintegration
  MODULE PROCEDURE InitOverintegration
END INTERFACE

INTERFACE Overintegration
  MODULE PROCEDURE Overintegration
END INTERFACE

INTERFACE FinalizeOverintegration
  MODULE PROCEDURE FinalizeOverintegration
END INTERFACE

PUBLIC :: InitOverintegration
PUBLIC :: Overintegration
PUBLIC :: FinalizeOverintegration
PUBLIC :: DefineParametersOverintegration
!==================================================================================================================================



CONTAINS

!==================================================================================================================================
!> Define parameters of overintegration module
!==================================================================================================================================
SUBROUTINE DefineParametersOverintegration()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Overintegration")
CALL prms%CreateIntFromStringOption('OverintegrationType', "Type of overintegration. None, CutOff, ConsCutOff","none")
CALL addStrListEntry('OverintegrationType','none',      OVERINTEGRATIONTYPE_NONE)
CALL addStrListEntry('OverintegrationType','cutoff',    OVERINTEGRATIONTYPE_CUTOFF)
CALL addStrListEntry('OverintegrationType','conscutoff',OVERINTEGRATIONTYPE_CONSCUTOFF)
CALL prms%CreateIntOption('NUnder',             "Polynomial degree to which solution is filtered (OverintegrationType == 1 or 2")
END SUBROUTINE DefineParametersOverintegration



!==================================================================================================================================
!> Initialize all necessary information to perform overintegration
!==================================================================================================================================
SUBROUTINE InitOverintegration()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Overintegration_Vars
USE MOD_Interpolation       ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars  ,ONLY:InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType
USE MOD_ChangeBasisByDim    ,ONLY:ChangeBasisVolume
USE MOD_ReadInTools         ,ONLY:GETINT,GETINTFROMSTR
USE MOD_Mesh_Vars           ,ONLY:DetJac_Ref,NGeoRef,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iDeg,iElem,i,j,k                             !< Loop counters
REAL,ALLOCATABLE    :: DetJac_NUnder(:,:,:,:)                       !< Determinant of the Jacobian on Nunder,
                                                                    !< size [1,0..Nunder,0..Nunder,0..Nunder]
REAL,ALLOCATABLE    :: Vdm_NGeoRef_NUnder(:,:)
!==================================================================================================================================
! Check if the necessary prerequisites are met
IF(OverintegrationInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
  CALL CollectiveStop(__STAMP__,&
    'InitOverintegration not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OVERINTEGRATION...'

!Set default values
NUnder=PP_N

OverintegrationType = GETINTFROMSTR('OverintegrationType')
SELECT CASE(OverintegrationType)
CASE (OVERINTEGRATIONTYPE_NONE) ! no overintegration, collocation DGSEM
  OverintegrationType = 0
CASE (OVERINTEGRATIONTYPE_CUTOFF) ! modal cut-off filter on Ju_t
  ! Prepare filter matrix for modal filtering: Identity matrix up to Nunder, then zero diagonal entries
  ALLOCATE(OverintegrationMat(0:PP_N,0:PP_N))
  OverintegrationMat = 0.
  NUnder = GETINT('NUnder')
  DO iDeg=0,NUnder
    OverintegrationMat(iDeg,iDeg) = 1.
  END DO
  ! Assemble filter matrix in nodal space: transform nodal solution to modal (Legendre) polynomials, apply filter matrix and
  ! transfer back to nodal (Lagrange) basis. Combine these operations into single ready-to-use matrix
  OverintegrationMat=MATMUL(MATMUL(Vdm_Leg,OverintegrationMat),sVdm_Leg)
  SWRITE(UNIT_stdOut,'(A)') ' Method of overintegration: cut-off filter'

CASE (OVERINTEGRATIONTYPE_CONSCUTOFF) ! conservative modal cut-off filter: Here, Ju_t is projection filtered from N to Nunder,
                        ! then sJ is applied on Nunder u_t is then interpolated back onto N grid: ensures that the modal content
                        ! of u_t between NUnder and N is zero. Here, the sJ on Nunder is prepared for later use from projection
                        ! of J from N to Under and inverting.
  NUnder = GETINT('NUnder')
  IF(NUnder.LT.PP_N)THEN
    !global
    ALLOCATE(Vdm_N_NUnder(0:NUnder,0:PP_N),Vdm_NUnder_N(0:PP_N,0:NUnder))
    ALLOCATE(sJNUnder(0:NUnder,0:NUnder,0:ZDIM(NUnder),nElems))
    CALL GetVandermonde(PP_N,NodeType,NUnder,NodeType,Vdm_N_NUnder,Vdm_NUnder_N,modal=.TRUE.)
    !local
    ALLOCATE(Vdm_NGeoRef_NUnder(0:NUnder,0:NGeoRef),DetJac_NUnder(1,0:NUnder,0:NUnder,0:ZDIM(NUnder)))
    CALL GetVandermonde(NGeoRef,NodeType,NUnder,NodeType,Vdm_NGeoRef_NUnder,modal=.TRUE.)

    DO iElem=1,nElems
      CALL ChangeBasisVolume(1,NGeoRef,NUnder,Vdm_NGeoRef_NUnder,DetJac_Ref(:,:,:,:,iElem),DetJac_NUnder)
      DO k=0,ZDIM(NUnder); DO j=0,NUnder; DO i=0,NUnder
        sJNUnder(i,j,k,iElem)=1./DetJac_NUnder(1,i,j,k)
      END DO; END DO; END DO !i,j,k=0,PP_N
    END DO
    DEALLOCATE(Vdm_NGeoRef_NUnder,DetJac_NUnder)
  ELSE
    SWRITE(UNIT_stdOut,'(A)') ' WARNING: Overintegration is disabled for NUnder >= N !!!'
    OverintegrationType = OVERINTEGRATIONTYPE_NONE
  END IF
  SWRITE(UNIT_stdOut,'(A)') ' Method of overintegration: cut-off filter (conservative)'

CASE DEFAULT
  CALL Abort(__STAMP__, &
      "Unknown OverintegrationType!")
END SELECT

OverintegrationInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OVERINTEGRATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOverintegration



!==================================================================================================================================
!> \brief Performs the overintegration according to the selected type.
!> The Overintegration routine will call the specfic functions needed to perform the selected overintegration type.
!> For the cut off version (option 1): Call the filter routine to apply the modal cut off filter.
!> For the conservative cut off version (option 2): Call the special conservative filter routine.
!==================================================================================================================================
SUBROUTINE Overintegration(U_in)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars            ,ONLY: nElems
USE MOD_Overintegration_Vars ,ONLY: OverintegrationType,OverintegrationMat
USE MOD_Filter               ,ONLY: Filter_Pointer
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< Time derivative vector to be filtered
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SELECT CASE (OverintegrationType)
CASE (OVERINTEGRATIONTYPE_CUTOFF)
  CALL Filter_Pointer(U_in, OverintegrationMat)
CASE (OVERINTEGRATIONTYPE_CONSCUTOFF)
  CALL FilterConservative(U_in)
CASE DEFAULT
  CALL Abort(__STAMP__, &
    "OverintegrationType unknown or not allowed!", IntInfo=OverintegrationType)
END SELECT
END SUBROUTINE Overintegration

!TODO Implement 2d

!==================================================================================================================================
!> Modal cutoff filter conserving both JU and U
!> project JU_t down from degree N to NUnder, then divide by Jacobian built on NUnder, interpolate resulting U up to N again
!> input : JU_t on degree N, output: filtered Ut on degree N
!==================================================================================================================================
SUBROUTINE FilterConservative(U_in)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Overintegration_Vars,  ONLY: Vdm_N_NUnder,Vdm_NUnder_N,NUnder,sJNUnder
USE MOD_ChangeBasisByDim,      ONLY: ChangeBasisVolume
USE MOD_Vector,                ONLY: VNullify
USE MOD_Mesh_Vars,             ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars,               ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)    !< Time derivative / JU_t to be filtered
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if PP_dim==2
INTEGER             :: i,j,iElem
#else
INTEGER             :: i,j,k,iElem
REAL                :: X3D_Buf1(PP_nVar,0:NUnder,0:PP_N,0:PP_N)     ! intermediate results from 1D interpolations
REAL                :: X3D_Buf2(PP_nVar,0:NUnder,0:NUnder,0:PP_N)   !
REAL                :: X3D_Buf3(PP_nVar,0:PP_N,0:NUnder,0:NUnder)   !
REAL                :: X3D_Buf4(PP_nVar,0:PP_N,0:PP_N,0:NUnder)     !
INTEGER             :: iU,jU,kU,nDOF_N,nDOF_NUnder
#endif
REAL                :: U_loc(   PP_nVar,0:NUnder,0:NUnder,0:NUnder) ! U_t / JU_t on NUnder
!==================================================================================================================================
#if PP_dim==2
! The following 5 lines are original lines of code for the conservative filtering.
! The below code does the same, but optimized for performance. For 2D, we use the non-optimized code.
! TODO: Use optimized code for 2D!
! === BEGIN ORIGINAL CODE ===
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) CYCLE ! Do only, when DG element
#endif
  CALL ChangeBasisVolume(PP_nVar,PP_N,NUnder,Vdm_N_NUnder,U_in(:,:,:,:,iElem),U_loc)
  DO j=0,NUnder; DO i=0,NUnder
    U_loc(:,i,j,0)=U_loc(:,i,j,0)*sJNUnder(i,j,0,iElem)
  END DO; END DO
  CALL ChangeBasisVolume(PP_nVar,NUnder,PP_N,Vdm_NUnder_N,U_loc,U_in(:,:,:,:,iElem))
END DO
! === END ORIGINAL CODE ===
#else
nDOF_N     =PP_nVar*(PP_N+1  )**3
nDOF_NUnder=PP_nVar*(NUnder+1)**3
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) CYCLE ! Do only, when DG element
#endif
  ! First transform JU_N to JU_NUnder
  ! first direction i
  DO k=0,PP_N; DO j=0,PP_N
    DO iU=0,NUnder
      X3D_Buf1(:,iU,j,k)=Vdm_N_NUnder(iU,0)*U_In(:,0,j,k,iElem)
    END DO ! iU
    DO i=1,PP_N
      DO iU=0,NUnder
        X3D_Buf1(:,iU,j,k)=X3D_Buf1(:,iU,j,k)+Vdm_N_NUnder(iU,i)*U_In(:,i,j,k,iElem)
      END DO ! iU
    END DO ! i
  END DO; END DO ! k,j

  CALL VNullify(nDOF_NUnder,U_loc)
  CALL VNullify(nDOF_N,U_in(:,:,:,:,iElem))

  ! second direction j
  DO k=0,PP_N
    DO jU=0,NUnder; DO iU=0,NUnder
      X3D_Buf2(:,iU,jU,k)=Vdm_N_NUnder(jU,0)*X3D_Buf1(:,iU,0,k)
    END DO; END DO ! iU, jU
    DO j=1,PP_N
      DO jU=0,NUnder; DO iU=0,NUnder
        X3D_Buf2(:,iU,jU,k)=X3D_Buf2(:,iU,jU,k)+Vdm_N_NUnder(jU,j)*X3D_Buf1(:,iU,j,k)
      END DO; END DO ! iU, jU
    END DO ! j
  END DO ! k
  ! last direction k
  DO k=0,PP_N
    DO kU=0,NUnder; DO jU=0,NUnder; DO iU=0,NUnder
      U_loc(:,iU,jU,kU)=U_loc(:,iU,jU,kU)+Vdm_N_NUnder(kU,k)*X3D_Buf2(:,iU,jU,k)
    END DO; END DO; END DO ! iU, jU, kU
  END DO ! k


  ! Apply Jacobian (JU_NUnder -> U_NUnder)
  DO k=0,NUnder; DO j=0,NUnder; DO i=0,NUnder
    U_loc(:,i,j,k)=U_loc(:,i,j,k)*sjNUnder(i,j,k,iElem)
  END DO; END DO; END DO


  ! Now transform U_NUnder back to U_N
  ! First direction iU
  DO kU=0,NUnder; DO jU=0,NUnder
    DO i=0,PP_N
      X3D_Buf3(:,i,jU,kU)=Vdm_NUnder_N(i,0)*U_loc(:,0,jU,kU)
    END DO
    DO iU=1,NUnder
      DO i=0,PP_N
        X3D_Buf3(:,i,jU,kU)=X3D_Buf3(:,i,jU,kU)+Vdm_NUnder_N(i,iU)*U_loc(:,iU,jU,kU)
      END DO
    END DO
  END DO; END DO ! jU, jU
  ! second direction jU
  DO kU=0,NUnder
    DO j=0,PP_N; DO i=0,PP_N
      X3D_Buf4(:,i,j,kU)=Vdm_NUnder_N(j,0)*X3D_Buf3(:,i,0,kU)
    END DO; END DO ! i,j
    DO jU=1,NUnder
      DO j=0,PP_N; DO i=0,PP_N
        X3D_Buf4(:,i,j,kU)=X3D_Buf4(:,i,j,kU)+Vdm_NUnder_N(j,jU)*X3D_Buf3(:,i,jU,kU)
      END DO; END DO ! i,j
    END DO ! jU
  END DO ! kU
  ! last direction kU
  DO kU=0,NUnder
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      U_in(:,i,j,k,iElem)=U_in(:,i,j,k,iElem)+Vdm_NUnder_N(k,kU)*X3D_Buf4(:,i,j,kU)
    END DO; END DO; END DO ! i,j,k
  END DO ! kU

END DO
#endif

END SUBROUTINE FilterConservative



!==================================================================================================================================
!> Deallocate Overintegration arrays
!==================================================================================================================================
SUBROUTINE FinalizeOverintegration()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Overintegration_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(OverintegrationMat)
SDEALLOCATE(sJNUnder)
SDEALLOCATE(Vdm_NUnder_N)
SDEALLOCATE(Vdm_N_NUnder)
OverintegrationInitIsDone = .FALSE.
END SUBROUTINE FinalizeOverintegration



END MODULE MOD_Overintegration
