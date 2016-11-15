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


!==================================================================================================================================
!> Changes a 2D or 3D tensor product polynomial with Lagrange Basis of degree NIn to
!> 2D or 3D tensor product polynomial of a Lagrange Basis NOut, using two
!> arbitrary point distributions xi_In(0:NIn) and xi_Out(0:NOut) and a series of 1D operations 
!> \f[ \tilde{u}_{:,j} = \mathcal{V}_{1D,(Nout+1)x(Nin+1)}^{-1} \hat{u}_{:,j} \f]
!> \f[ \hat{p}_{i,:} = \mathcal{V}_{1D,(Nout+1)x(Nin+1)}^{-1} \tilde{u}_{i,:} \f]
!==================================================================================================================================
MODULE MOD_ChangeBasis
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ChangeBasis3D
  MODULE PROCEDURE ChangeBasis3D_Single
  MODULE PROCEDURE ChangeBasis3D_Mult
END INTERFACE

INTERFACE ChangeBasis3D_XYZ
  MODULE PROCEDURE ChangeBasis3D_XYZ
END INTERFACE

INTERFACE ChangeBasis2D
  MODULE PROCEDURE ChangeBasis2D
END INTERFACE

INTERFACE ChangeBasis2D_selective
  MODULE PROCEDURE ChangeBasis2D_selective
  MODULE PROCEDURE ChangeBasis2D_selective_overwrite
END INTERFACE

INTERFACE ChangeBasis1D 
  MODULE PROCEDURE ChangeBasis1D
END INTERFACE

PUBLIC :: ChangeBasis3D
PUBLIC :: ChangeBasis3D_XYZ
PUBLIC :: ChangeBasis2D
PUBLIC :: ChangeBasis2D_selective
PUBLIC :: ChangeBasis1D
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Interpolate a 3D tensor product Lagrange polynomial defined by (NIn+1) 1D Lagrange basis functions of order (Nin) and node 
!> positions xi_In(0:Nin) to another 3D tensor product Lagrange basis defined by (NOut+1) 1D interpolation points on the node
!> positions xi_out(0:NOut).
!> xi is defined in the 1D referent element \f$ \xi \in [-1,1] \f$.
!> _Mult means that the routine is not restricted to one element, but the data fields of several elements (nElems) is
!> processed.
!==================================================================================================================================
SUBROUTINE ChangeBasis3D_Mult(nVar,nElems,NIn,NOut,Vdm,UIn,UOut,addToOutput)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: nVar                                    !< Number of variables
INTEGER,INTENT(IN)  :: nElems                                  !< Number of elements
INTEGER,INTENT(IN)  :: NIn                                     !< Input polynomial degree, no. of points = NIn+1
INTEGER,INTENT(IN)  :: NOut                                    !< Output polynomial degree, no. of points = NOut+1
                                   
REAL,INTENT(IN)     :: UIn(nVar,0:NIn,0:NIn,0:NIn,nElems)      !< Input field, dimensions must match nVar,NIn and nElems
REAL,INTENT(IN)     :: Vdm(0:NOut,0:NIn)                       !< 1D Vandermonde In -> Out
REAL,INTENT(INOUT)  :: UOut(nVar,0:NOut,0:NOut,0:NOut,nElems)  !< Output field
LOGICAL,INTENT(IN)  :: addToOutput                             !< TRUE: add the result to 'in' state of Uout, FALSE: overwrite Uout
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iI,jI,kI,iO,jO,kO,iElem,a,b,nVar2
REAL,ALLOCATABLE    :: UBuf1(:,:,:,:),UBuf2(:,:,:,:)           ! first and second intermediate results from 1D interpolations
!==================================================================================================================================
nVar2=nVar*nElems
IF(nVar2.GT.2*nVar)THEN
  ALLOCATE(UBuf2(nVar2,0:NIn,0:NIn,0:NIn))
  ALLOCATE(UBuf1(nVar2,0:NOut,0:NIn,0:NIn))

  ! pack solution
  DO iElem=1,nElems
    a=nVar*(iElem-1)+1
    b=nVar*iElem
    DO kI=0,NIn; DO jI=0,NIn; DO iI=0,NIn
      Ubuf2(a:b,iI,jI,kI)=UIn(:,iI,jI,kI,iElem)
    END DO; END DO; END DO
  END DO

  ! first direction iI
  DO kI=0,NIn; DO jI=0,NIn
    DO iO=0,NOut
      UBuf1(:,iO,jI,kI)=Vdm(iO,0)*Ubuf2(:,0,jI,kI)
    END DO
    DO iI=1,NIn
      DO iO=0,NOut
        UBuf1(:,iO,jI,kI)=UBuf1(:,iO,jI,kI)+Vdm(iO,iI)*Ubuf2(:,iI,jI,kI)
      END DO
    END DO
  END DO; END DO

  DEALLOCATE(Ubuf2)
  ALLOCATE(UBuf2(nVar2,0:NOut,0:NOut,0:NIn))

  ! second direction jI
  DO kI=0,NIn
    DO jO=0,NOut; DO iO=0,NOut
      UBuf2(:,iO,jO,kI)=Vdm(jO,0)*UBuf1(:,iO,0,kI)
    END DO; END DO
    DO jI=1,NIn
      DO jO=0,NOut; DO iO=0,NOut
        UBuf2(:,iO,jO,kI)=UBuf2(:,iO,jO,kI)+Vdm(jO,jI)*UBuf1(:,iO,jI,kI)
      END DO; END DO
    END DO
  END DO

  DEALLOCATE(Ubuf1)
  ALLOCATE(UBuf1(nVar2,0:NOut,0:NOut,0:NOut))

  ! last direction kI
  DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
    Ubuf1(:,iO,jO,kO)=Vdm(kO,0)*UBuf2(:,iO,jO,0)
  END DO; END DO; END DO
  DO kI=1,NIn
    DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
      Ubuf1(:,iO,jO,kO)=Ubuf1(:,iO,jO,kO)+Vdm(kO,kI)*UBuf2(:,iO,jO,kI)
    END DO; END DO; END DO
  END DO

  ! unpack solution
  IF(addToOutput)THEN
    DO iElem=1,nElems
      a=nVar*(iElem-1)+1
      b=nVar*iElem
      DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,kO,iElem)=UOut(:,iO,jO,kO,iElem)+Ubuf1(a:b,iO,jO,kO)
      END DO; END DO; END DO
    END DO
  ELSE
    DO iElem=1,nElems
      a=nVar*(iElem-1)+1
      b=nVar*iElem
      DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,kO,iElem)=Ubuf1(a:b,iO,jO,kO)
      END DO; END DO; END DO
    END DO
  END IF
  DEALLOCATE(UBuf1,Ubuf2)

ELSE

  ALLOCATE(UBuf1(nVar,0:NOut,0:NIn,0:NIn))
  ALLOCATE(UBuf2(nVar,0:NOut,0:NOut,0:NIn))
  DO iElem=1,nElems
    ! first direction iI
    DO kI=0,NIn; DO jI=0,NIn
      DO iO=0,NOut
        UBuf1(:,iO,jI,kI)=Vdm(iO,0)*UIn(:,0,jI,kI,iElem)
      END DO
      DO iI=1,NIn
        DO iO=0,NOut
          UBuf1(:,iO,jI,kI)=UBuf1(:,iO,jI,kI)+Vdm(iO,iI)*UIn(:,iI,jI,kI,iElem)
        END DO
      END DO
    END DO; END DO
    ! second direction jI
    DO kI=0,NIn
      DO jO=0,NOut; DO iO=0,NOut
        UBuf2(:,iO,jO,kI)=Vdm(jO,0)*UBuf1(:,iO,0,kI)
      END DO; END DO
      DO jI=1,NIn
        DO jO=0,NOut; DO iO=0,NOut
          UBuf2(:,iO,jO,kI)=UBuf2(:,iO,jO,kI)+Vdm(jO,jI)*UBuf1(:,iO,jI,kI)
        END DO; END DO
      END DO
    END DO
    ! last direction kI
    IF(addToOutput)THEN
      DO kI=0,NIn
        DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
          UOut(:,iO,jO,kO,iElem)=UOut(:,iO,jO,kO,iElem)+Vdm(kO,kI)*UBuf2(:,iO,jO,kI)
        END DO; END DO; END DO
      END DO
    ELSE
      DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,kO,iElem)=Vdm(kO,0)*UBuf2(:,iO,jO,0)
      END DO; END DO; END DO
      DO kI=1,NIn
        DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
          UOut(:,iO,jO,kO,iElem)=UOut(:,iO,jO,kO,iElem)+Vdm(kO,kI)*UBuf2(:,iO,jO,kI)
        END DO; END DO; END DO
      END DO
    END IF
  END DO
  DEALLOCATE(UBuf1,Ubuf2)

END IF
END SUBROUTINE ChangeBasis3D_Mult



!==================================================================================================================================
!> Interpolate a 3D tensor product Lagrange polynomial defined by (NIn+1) 1D Lagrange basis functions of order (Nin) and node 
!> positions xi_In(0:Nin) to another 3D tensor product Lagrange basis defined by (NOut+1) 1D interpolation points on the node
!> positions xi_out(0:NOut).
!> xi is defined in the 1D referent element \f$ \xi \in [-1,1] \f$.
!>  _Single is only suitable for one tensor product element
!==================================================================================================================================
SUBROUTINE ChangeBasis3D_Single(Dim1,NIn,NOut,Vdm,X3D_In,X3D_Out)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1                                    !< Number of variables
INTEGER,INTENT(IN)  :: NIn                                     !< Input polynomial degree, no. of points = NIn+1
INTEGER,INTENT(IN)  :: NOut                                    !< Output polynomial degree, no. of points = NOut+1
REAL,INTENT(IN)     :: X3D_In(1:Dim1,0:NIn,0:NIn,0:NIn)        !< Input field, dimensions must match Dim1,NIn
REAL,INTENT(OUT)    :: X3D_Out(1:Dim1,0:NOut,0:NOut,0:NOut)    !< Output field, dimensions must match Dim1,NOut
REAL,INTENT(IN)     :: Vdm(0:NOut,0:NIn)                       !< 1D Vandermonde In -> Out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNIn,jNIn,kNIn,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(1:Dim1,0:NOut,0:NIn,0:NIn)     ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:Dim1,0:NOut,0:NOut,0:NIn)    ! second intermediate results from 1D interpolations
!==================================================================================================================================
X3D_buf1=0.
! first direction iNIn
DO kNIn=0,NIn
  DO jNIn=0,NIn
    DO iNIn=0,NIn
      DO iN_Out=0,NOut
        X3D_Buf1(:,iN_Out,jNIn,kNIn)=X3D_Buf1(:,iN_Out,jNIn,kNIn)+Vdm(iN_Out,iNIn)*X3D_In(:,iNIn,jNIn,kNIn)
      END DO
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jNIn
DO kNIn=0,NIn
  DO jNIn=0,NIn
    DO jN_Out=0,NOut
      DO iN_Out=0,NOut
        X3D_Buf2(:,iN_Out,jN_Out,kNIn)=X3D_Buf2(:,iN_Out,jN_Out,kNIn)+Vdm(jN_Out,jNIn)*X3D_Buf1(:,iN_Out,jNIn,kNIn)
      END DO
    END DO
  END DO
END DO
X3D_Out=0.
! last direction kNIn
DO kNIn=0,NIn
  DO kN_Out=0,NOut
    DO jN_Out=0,NOut
      DO iN_Out=0,NOut
        X3D_Out(:,iN_Out,jN_Out,kN_Out)=X3D_Out(:,iN_Out,jN_Out,kN_Out)+Vdm(kN_Out,kNIn)*X3D_Buf2(:,iN_Out,jN_Out,kNIn)
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis3D_Single

!==================================================================================================================================
!> Interpolate a 3D tensor product Lagrange polynomial defined by (NIn+1) 1D Lagrange basis functions of order (Nin) and node 
!> positions xi_In(0:Nin) to another 3D tensor product Lagrange basis defined by (NOut+1) 1D interpolation points on the node
!> positions xi_out(0:NOut) using DIFFERENT 1D Vdm matrices in the xi,eta and zeta directions.
!> xi is defined in the 1D referent element \f$ \xi \in [-1,1] \f$.
!==================================================================================================================================
SUBROUTINE ChangeBasis3D_XYZ(Dim1,NIn,NOut,Vdm_xi,Vdm_eta,Vdm_zeta,X3D_In,X3D_Out)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1                                    !< Number of variables
INTEGER,INTENT(IN)  :: NIn                                     !< Input polynomial degree, no. of points = NIn+1
INTEGER,INTENT(IN)  :: NOut                                    !< Output polynomial degree, no. of points = NOut+1
REAL,INTENT(IN)     :: X3D_In(1:Dim1,0:NIn,0:NIn,0:NIn)        !< Input field, dimensions must match Dim1,NIn
REAL,INTENT(OUT)    :: X3D_Out(1:Dim1,0:NOut,0:NOut,0:NOut)    !< Output field, dimensions must match Dim1,NOut
REAL,INTENT(IN)     :: Vdm_xi(0:NOut,0:NIn)                    !< 1D Vandermonde In -> Out xi direction
REAL,INTENT(IN)     :: Vdm_eta(0:NOut,0:NIn)                   !< 1D Vandermonde In -> Out eta direction
REAL,INTENT(IN)     :: Vdm_zeta(0:NOut,0:NIn)                  !< 1D Vandermonde In -> Out zeta direction

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNIn,jNIn,kNIn,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(1:Dim1,0:NOut,0:NIn,0:NIn)     ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:Dim1,0:NOut,0:NOut,0:NIn)    ! second intermediate results from 1D interpolations
!==================================================================================================================================
X3D_buf1=0.
! first direction iNIn
DO kNIn=0,NIn
  DO jNIn=0,NIn
    DO iNIn=0,NIn
      DO iN_Out=0,NOut
        X3D_Buf1(:,iN_Out,jNIn,kNIn)=X3D_Buf1(:,iN_Out,jNIn,kNIn)+Vdm_xi(iN_Out,iNIn)*X3D_In(:,iNIn,jNIn,kNIn)
      END DO
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jNIn
DO kNIn=0,NIn
  DO jNIn=0,NIn
    DO jN_Out=0,NOut
      DO iN_Out=0,NOut
        X3D_Buf2(:,iN_Out,jN_Out,kNIn)=X3D_Buf2(:,iN_Out,jN_Out,kNIn)+Vdm_eta(jN_Out,jNIn)*X3D_Buf1(:,iN_Out,jNIn,kNIn)
      END DO
    END DO
  END DO
END DO
X3D_Out=0.
! last direction kNIn
DO kNIn=0,NIn
  DO kN_Out=0,NOut
    DO jN_Out=0,NOut
      DO iN_Out=0,NOut
        X3D_Out(:,iN_Out,jN_Out,kN_Out)=X3D_Out(:,iN_Out,jN_Out,kN_Out)+Vdm_zeta(kN_Out,kNIn)*X3D_Buf2(:,iN_Out,jN_Out,kNIn)
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis3D_XYZ


!==================================================================================================================================
!> Interpolate a 2D tensor product Lagrange polynomial defined by (NIn+1) 1D Lagrange basis functions of order (Nin) and node 
!> positions xi_In(0:Nin) to another 2D tensor product Lagrange basis defined by (NOut+1) 1D interpolation points on the node
!> positions xi_out(0:NOut).
!> xi is defined in the 1D referent element \f$ \xi \in [-1,1] \f$.
!>  _Single is only suitable for one tensor product element
!==================================================================================================================================
SUBROUTINE ChangeBasis2D(Dim1,NIn,NOut,Vdm,X2D_In,X2D_Out)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1                                    !< Number of variables
INTEGER,INTENT(IN)  :: NIn                                     !< Input polynomial degree, no. of points = NIn+1
INTEGER,INTENT(IN)  :: NOut                                    !< Output polynomial degree, no. of points = NOut+1
REAL,INTENT(IN)     :: X2D_In(1:Dim1,0:NIn,0:NIn)              !< Input field, dimensions must match Dim1,NIn
REAL,INTENT(OUT)    :: X2D_Out(1:Dim1,0:NOut,0:NOut)           !< Output field, dimensions must match Dim1,NOut
REAL,INTENT(IN)     :: Vdm(0:NOut,0:NIn)                       !< 1D Vandermonde In -> Out

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNIn,jNIn,iN_Out,jN_Out
REAL                :: X2D_Buf1(1:Dim1,0:NOut,0:NIn)           ! first intermediate results from 1D interpolations
!==================================================================================================================================
X2D_buf1=0.
! first direction iNIn
DO jNIn=0,NIn
  DO iNIn=0,NIn
    DO iN_Out=0,NOut
      X2D_Buf1(:,iN_Out,jNIn)=X2D_Buf1(:,iN_Out,jNIn)+Vdm(iN_Out,iNIn)*X2D_In(:,iNIn,jNIn)
    END DO
  END DO
END DO
X2D_Out=0.
! second direction jNIn
DO jNIn=0,NIn
  DO jN_Out=0,NOut
    DO iN_Out=0,NOut
      X2D_Out(:,iN_Out,jN_Out)=X2D_Out(:,iN_Out,jN_Out)+Vdm(jN_Out,jNIn)*X2D_Buf1(:,iN_Out,jNIn)
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis2D


! TODO: documentation, variable bezeichnung im gleichen Stil wie oben
SUBROUTINE ChangeBasis2D_selective(Dim1,N_In,N_Out,SideID_1,SideID_N,firstSideID,lastSideID,Vdm,X2D_In,X2D_Out,mask,mask_ref,&
        addToOutput)
!==================================================================================================================================
!==================================================================================================================================
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: Dim1,N_In,N_Out,firstSideID,lastSideID,SideID_1,SideID_N
REAL,INTENT(IN)              :: X2D_In(1:Dim1,0:N_In,0:N_In,SideID_1:SideID_N)
REAL,INTENT(IN)              :: Vdm(0:N_Out,0:N_In)
REAL,INTENT(INOUT)           :: X2D_Out(1:Dim1,0:N_Out,0:N_Out,SideID_1:SideID_N)
INTEGER,INTENT(IN),OPTIONAL  :: mask(SideID_1:SideID_N)
INTEGER,INTENT(IN),OPTIONAL  :: mask_ref
LOGICAL,INTENT(IN),OPTIONAL  :: addToOutput
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iN_In,jN_In,iN_Out,jN_Out,SideID,SideID2,iVec
REAL                :: X2D_Buf1(Dim1*PP_VEC,0:N_Out,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: VECIN (Dim1*PP_VEC,0:N_In,0:N_In)
REAL                :: VECOUT(Dim1*PP_VEC,0:N_Out,0:N_Out)
INTEGER             :: mask_loc(firstSideID:lastSideID), mask_ref_loc
LOGICAL             :: addToOutput_loc
!==================================================================================================================================
IF (PRESENT(mask)) THEN
  mask_loc = mask(firstSideID:lastSideID)
  mask_ref_loc = mask_ref
ELSE
  mask_loc = 0
  mask_ref_loc = 0
END IF
addToOutput_loc=.FALSE.
IF (PRESENT(addToOutput)) THEN
  addToOutput_loc = addToOutput
END IF
SideID = firstSideID
SideID2 = firstSideID
DO WHILE (SideID.LE.lastSideID)
  ! pack solution into vector VECIN
  iVec = 1
  DO WHILE (iVec < Dim1*PP_VEC)
    IF (mask_loc(SideID).EQ.mask_ref_loc) THEN
      VECIN(iVec:iVec+Dim1-1,:,:) = X2D_In(:,:,:,SideID)
      IF (addToOutput_loc) THEN
        ! initialize VECOUT with X2D_Out
        VECOUT(iVec:iVec+Dim1-1,:,:) = X2D_Out(:,:,:,SideID)
      END IF
      iVec = iVec + Dim1
    END IF
    SideID = SideID + 1
    IF (SideID.GT.lastSideID) EXIT
  END DO
  ! nullify VECOUT if not initialized with data from X2D_Out
  IF (.NOT.addToOutput_loc) THEN
    VECOUT=0.
  END IF

  X2D_buf1=0.
  ! first direction iN_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      DO iN_Out=0,N_Out
        X2D_Buf1(1:iVec-1,iN_Out,jN_In)=X2D_Buf1(1:iVec-1,iN_Out,jN_In)+Vdm(iN_Out,iN_In)*VECIN(1:iVec-1,iN_In,jN_In)
      END DO
    END DO
  END DO
  ! second direction jN_In
  DO jN_In=0,N_In
    DO jN_Out=0,N_Out
      DO iN_Out=0,N_Out
        VECOUT(1:iVec-1,iN_Out,jN_Out)=VECOUT(1:iVec-1,iN_Out,jN_Out)+Vdm(jN_Out,jN_In)*X2D_Buf1(1:iVec-1,iN_Out,jN_In)
      END DO
    END DO
  END DO

  ! unpack solution from vector VECOUT
  iVec = 1
  DO WHILE (iVec < Dim1*PP_VEC)
    IF (mask_loc(SideID2).EQ.mask_ref_loc) THEN
      X2D_Out(:,:,:,SideID2) = VECOUT(iVec:iVec+Dim1-1,:,:)
      iVec = iVec + Dim1
    END IF
    SideID2 = SideID2 + 1
    IF (SideID2.GT.lastSideID) EXIT
  END DO
END DO  
END SUBROUTINE ChangeBasis2D_selective

! TODO: documentation, variable bezeichnung im gleichen Stil wie oben
SUBROUTINE ChangeBasis2D_selective_overwrite(Dim1,N_In,SideID_1,SideID_N,firstSideID,lastSideID,Vdm,X2D_InOut,mask,mask_ref)
!==================================================================================================================================
!==================================================================================================================================
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: Dim1,N_In,firstSideID,lastSideID,SideID_1,SideID_N
REAL,INTENT(INOUT)           :: X2D_InOut(1:Dim1,0:N_In,0:N_In,SideID_1:SideID_N)
REAL,INTENT(IN)              :: Vdm(0:N_In,0:N_In)
INTEGER,INTENT(IN),OPTIONAL  :: mask(SideID_1:SideID_N)
INTEGER,INTENT(IN),OPTIONAL  :: mask_ref
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iN_In,jN_In,iN_Out,jN_Out,SideID,SideID2,iVec
REAL                :: X2D_Buf1(Dim1*PP_VEC,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: VEC     (Dim1*PP_VEC,0:N_In,0:N_In)
INTEGER             :: mask_loc(firstSideID:lastSideID), mask_ref_loc
!==================================================================================================================================
IF (PRESENT(mask)) THEN
  mask_loc = mask(firstSideID:lastSideID)
  mask_ref_loc = mask_ref
ELSE
  mask_loc = 0
  mask_ref_loc = 0
END IF    
SideID = firstSideID
SideID2 = firstSideID
DO WHILE (SideID < lastSideID)
  ! pack solution into vector VECIN
  iVec = 1
  DO WHILE (iVec < Dim1*PP_VEC)
    IF (mask_loc(SideID).EQ.mask_ref_loc) THEN
      VEC(iVec:iVec+Dim1-1,:,:) = X2D_InOut(:,:,:,SideID)
      iVec = iVec + Dim1
    END IF
    IF (SideID.EQ.lastSideID) EXIT
    SideID = SideID + 1
  END DO

  X2D_buf1=0.
  ! first direction iN_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      DO iN_Out=0,N_In
        X2D_Buf1(1:iVec-1,iN_Out,jN_In)=X2D_Buf1(1:iVec-1,iN_Out,jN_In)+Vdm(iN_Out,iN_In)*VEC(1:iVec-1,iN_In,jN_In)
      END DO
    END DO
  END DO
  VEC=0.
  ! second direction jN_In
  DO jN_In=0,N_In
    DO jN_Out=0,N_In
      DO iN_Out=0,N_In
        VEC(1:iVec-1,iN_Out,jN_Out)=VEC(1:iVec-1,iN_Out,jN_Out)+Vdm(jN_Out,jN_In)*X2D_Buf1(1:iVec-1,iN_Out,jN_In)
      END DO
    END DO
  END DO

  ! unpack solution from vector VEC
  iVec = 1
  DO WHILE (iVec < Dim1*PP_VEC)
    IF (mask_loc(SideID2).EQ.mask_ref_loc) THEN
      X2D_InOut(:,:,:,SideID2) = VEC(iVec:iVec+Dim1-1,:,:)
      iVec = iVec + Dim1
    END IF
    IF (SideID2.EQ.lastSideID) EXIT
    SideID2 = SideID2 + 1
  END DO
END DO  
END SUBROUTINE ChangeBasis2D_selective_overwrite

! TODO: documentation, variable bezeichnung im gleichen Stil wie oben
SUBROUTINE ChangeBasis1D(Dim1,N_In,N_Out,Vdm,X1D_In,X1D_Out)
!==================================================================================================================================
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Dim1,N_In,N_Out
REAL,INTENT(IN)    :: X1D_In(1:Dim1,0:N_In)
REAL,INTENT(IN)    :: Vdm(0:N_Out,0:N_In)
REAL,INTENT(OUT)   :: X1D_Out(1:Dim1,0:N_Out)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: iN_In,iN_Out
!==================================================================================================================================
X1D_Out=0.
! first direction iN_In
DO iN_In=0,N_In
  DO iN_Out=0,N_Out
    X1D_Out(:,iN_Out)=X1D_Out(:,iN_Out)+Vdm(iN_Out,iN_In)*X1D_In(:,iN_In)
  END DO
END DO
END SUBROUTINE ChangeBasis1D

END MODULE MOD_ChangeBasis
