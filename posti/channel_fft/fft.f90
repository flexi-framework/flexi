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

MODULE MOD_FFT
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitFFT
  MODULE PROCEDURE InitFFT
END INTERFACE

INTERFACE FinalizeFFT
  MODULE PROCEDURE FinalizeFFT
END INTERFACE

INTERFACE PerformFFT
  MODULE PROCEDURE PerformFFT
END INTERFACE

INTERFACE FFTOutput
  MODULE PROCEDURE FFTOutput
END INTERFACE

INTERFACE PrimStateAtFFTCoords
  MODULE PROCEDURE PrimStateAtFFTCoords
END INTERFACE

PUBLIC:: InitFFT,FinalizeFFT,PerformFFT,FFTOutput,PrimStateAtFFTCoords

CONTAINS

!===================================================================================================================================
!> Initialize FFT. Read in user-defined parameters, allocate arrays for solution and auxilliary FFT variables.
!> Prepare the Vandermonde to interpolate from the state to the FFT grid. Prepare coordinates of FFT grind.
!===================================================================================================================================
SUBROUTINE InitFFT()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE FFTW3
USE MOD_FFT_Vars
USE MOD_Commandline_Arguments
USE MOD_Interpolation,           ONLY: GetVandermonde
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_ReadInTools,             ONLY: GETREAL,GETINT
USE MOD_Mesh_Vars,               ONLY: nElems_IJK,nElems
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Interpolation_Vars ,     ONLY: NodeType,NodeTypeVISUInner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample,iVar,iRP,iStart,iEnd
!===================================================================================================================================
! Read in user-defined parameters
NCalc  = GETINT('NCalc')         ! Polynomial degree to perfrom DFFT on
Re_tau = GETREAL('Re_tau')       ! Reynolds number bases on friction velocity and channel half height

! Allocate solution array in DG vars, used to store state
ALLOCATE(U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))

! Allocate auxilliary arrays for FFT
N_FFT(1)=(NCalc+1)*nElems_IJK(1)
N_FFT(2)=(NCalc+1)*nElems_IJK(2)
N_FFT(3)=(NCalc+1)*nElems_IJK(3)
nSamples = (NCalc+1)*nElems_IJK(1)
nSamples_spec  = INT((nSamples)/2)+1
ALLOCATE(in(nSamples))
ALLOCATE(out(nSamples))
ALLOCATE(Ex_uu(N_FFT(2),nSamples_spec))
ALLOCATE(Ex_vv(N_FFT(2),nSamples_spec))
ALLOCATE(Ex_ww(N_FFT(2),nSamples_spec))
ALLOCATE(Ex_pp(N_FFT(2),nSamples_spec))
ALLOCATE(Ez_uu(N_FFT(2),nSamples_spec))
ALLOCATE(Ez_vv(N_FFT(2),nSamples_spec))
ALLOCATE(Ez_ww(N_FFT(2),nSamples_spec))
ALLOCATE(Ez_pp(N_FFT(2),nSamples_spec))
Ex_uu=0.;Ex_vv=0.;Ex_ww=0.;Ex_pp=0.
Ez_uu=0.;Ez_vv=0.;Ez_ww=0.;Ez_pp=0.
ALLOCATE(U_FFT(1:5,1:N_FFT(1),1:N_FFT(2),1:N_FFT(3)))
ALLOCATE(MS_PSD(N_FFT(2),4))
ALLOCATE(MS_t(N_FFT(2)))
ALLOCATE(M_t(N_FFT(2),4))
MS_t=0.
MS_PSD=0.
M_t=0.

! Get vandermonde from computation N on Gauss or Gauss-Lobatto points to equidistand FFT points
ALLOCATE(VdmGaussEqui(0:NCalc,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NCalc,NodeTypeVISUInner,VdmGaussEqui)

CALL DFFTW_PLAN_DFT_1D(plan,nSamples,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
CALL FFTCoords()

END SUBROUTINE InitFFT

!===================================================================================================================================
!> Interpolate the coordinates of the mesh to the FFT grid.
!===================================================================================================================================
SUBROUTINE FFTCoords
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FFT_Vars
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_Mesh_Vars             ,ONLY: Elem_IJK,nElems_IJK,Elem_xGP,nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: ii,jj,kk,i,j,k,iElem
INTEGER                             :: iGlob,jGlob,kGlob
REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: X_aux
!===================================================================================================================================
ALLOCATE(X_aux(1:3,0:NCalc,0:NCalc,0:NCalc))
ALLOCATE(X_FFT(1:3,1:N_FFT(1),1:N_FFT(2),1:N_FFT(3)))
X_FFT=0.
DO iElem=1,nGlobalElems
  ii=Elem_IJK(1,iElem)
  jj=Elem_IJK(2,iElem)
  kk=Elem_IJK(3,iElem)
  CALL ChangeBasis3D(3,PP_N,NCalc,VdmGaussEqui,Elem_xGP(:,:,:,:,iElem),X_aux)
  DO k=0,NCalc
    kGlob=(kk-1)*(NCalc+1)+k+1
    DO j=0,NCalc
      jGlob=(jj-1)*(NCalc+1)+j+1
      DO i=0,NCalc
        iGlob=(ii-1)*(NCalc+1)+i+1
        X_FFT(:,iGlob,jGlob,kGlob)=X_aux(:,i,j,k)
      END DO
    END DO
  END DO
END DO

DEALLOCATE(X_aux)
END SUBROUTINE FFTCoords

!===================================================================================================================================
!> Interpolate the state to the FFT grid and convert solution to primite variables.
!===================================================================================================================================
SUBROUTINE PrimStateAtFFTCoords
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FFT_Vars
USE MOD_ChangeBasis                 ,ONLY: ChangeBasis3D
USE MOD_Mesh_Vars                   ,ONLY: Elem_IJK,nElems_IJK,Elem_xGP,nGlobalElems
USE MOD_DG_Vars                     ,ONLY: U
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: ii,jj,kk,i,j,k,iElem
INTEGER                             :: iGlob,jGlob,kGlob
REAL                                :: kappaM1
REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: U_aux
!===================================================================================================================================
ALLOCATE(U_aux(1:5,0:NCalc,0:NCalc,0:NCalc))

KappaM1=0.4
U_FFT=0.
DO iElem=1,nGlobalElems
  ii=Elem_IJK(1,iElem)
  jj=Elem_IJK(2,iElem)
  kk=Elem_IJK(3,iElem)
  CALL ChangeBasis3D(PP_nVar,PP_N,NCalc,VdmGaussEqui,U(:,:,:,:,iElem),U_aux)
  DO k=0,NCalc
    kGlob=(kk-1)*(NCalc+1)+k+1
    DO j=0,NCalc
      jGlob=(jj-1)*(NCalc+1)+j+1
      DO i=0,NCalc
        iGlob=(ii-1)*(NCalc+1)+i+1
        U_FFT(1,iGlob,jGlob,kGlob)=U_aux(1,i,j,k)
        U_FFT(2:4,iGlob,jGlob,kGlob)=U_aux(2:4,i,j,k)/U_aux(1,i,j,k)
        U_FFT(5,iGlob,jGlob,kGlob)=KappaM1*(U_aux(5,i,j,k)-U_aux(1,i,j,k)*SUM(U_FFT(2:4,iGlob,jGlob,kGlob)**2)*0.5)
      END DO
    END DO
  END DO
END DO

DEALLOCATE(U_Aux)
END SUBROUTINE PrimStateAtFFTCoords

!===================================================================================================================================
!> Calls the actual FFT routines, called for each statefile.
!===================================================================================================================================
SUBROUTINE PerformFFT()
! MODULES                                                                                                                          !
USE FFTW3
USE MOD_FFT_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: j,k
!===================================================================================================================================
DO j=1,N_FFT(2)
  DO k=1,N_FFT(3)
    CALL FFT(U_FFT(2,1:nSamples,j,k),Ex_uu(j,1:nSamples_spec))
    CALL FFT(U_FFT(3,1:nSamples,j,k),Ex_vv(j,1:nSamples_spec))
    CALL FFT(U_FFT(4,1:nSamples,j,k),Ex_ww(j,1:nSamples_spec))
    CALL FFT(U_FFT(5,1:nSamples,j,k),Ex_pp(j,1:nSamples_spec))
  END DO
END DO
DO j=1,N_FFT(2)
  DO k=1,N_FFT(3)
    CALL FFT(U_FFT(2,k,j,1:nSamples),Ez_uu(j,1:nSamples_spec))
    CALL FFT(U_FFT(3,k,j,1:nSamples),Ez_vv(j,1:nSamples_spec))
    CALL FFT(U_FFT(4,k,j,1:nSamples),Ez_ww(j,1:nSamples_spec))
    CALL FFT(U_FFT(5,k,j,1:nSamples),Ez_pp(j,1:nSamples_spec))
  END DO
END DO
DO j=1,N_FFT(2)
  M_t(j,1)=M_t(j,1)+SUM(U_FFT(2,1:N_FFT(1),j,1:N_FFT(3)))
  M_t(j,2)=M_t(j,2)+SUM(U_FFT(3,1:N_FFT(1),j,1:N_FFT(3)))
  M_t(j,3)=M_t(j,3)+SUM(U_FFT(4,1:N_FFT(1),j,1:N_FFT(3)))
  M_t(j,4)=M_t(j,4)+SUM(U_FFT(5,1:N_FFT(1),j,1:N_FFT(3)))
  ! uv
  MS_t(j)=MS_t(j)+SUM(U_FFT(2,1:N_FFT(1),j,1:N_FFT(3))*U_FFT(3,1:N_FFT(1),j,1:N_FFT(3)))
END DO

END SUBROUTINE PerformFFT

!===================================================================================================================================
!> Low-level wrapper routine to DFFTW call.
!===================================================================================================================================
SUBROUTINE FFT(U_in,U_hat)
! MODULES
USE MOD_Globals
USE FFTW3
USE MOD_FFT_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                  ::U_in(1:N_FFT(1))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)               ::U_hat(nSamples_spec)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
in=U_in
CALL DFFTW_EXECUTE_DFT(plan, in, out)
out(2:nSamples_spec-1)=2*(1./REAL(nSamples)*ABS(out(2:nSamples_spec-1)))**2
out(1)=ABS(out(1)/REAL(nSamples))**2 !mean value 
IF(MOD(nSamples,2).EQ.0.)THEN
  out(nSamples_spec)=ABS(out(nSamples_spec)/REAL(nSamples))**2
!  WRITE(*,*)"Even number of samples"
ELSE
!  WRITE(*,*)"Odd number of samples"
  out(nSamples_spec)=2.*(ABS(out(nSamples_spec))/REAL(nSamples))**2
END IF
U_hat(2:nSamples_spec)=U_hat(2:nSamples_spec)+out(2:nSamples_spec)
!sum of mean square into first index
U_hat(1)=U_hat(1)+SUM(out(1:nSamples_spec))

END SUBROUTINE FFT

!===================================================================================================================================
!> Average the results and write the output files
!===================================================================================================================================
SUBROUTINE FFTOutput()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_FFT_Vars
USE MOD_Commandline_Arguments
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: j,k,Fileunit_EK
CHARACTER(LEN=255) :: FileName_EK
LOGICAL            :: connected
!===================================================================================================================================
!Average over the lower and upper channel halfes
DO j=N_FFT(2)/2+1,N_FFT(2)
  Ex_uu(j,:)=(Ex_uu(N_FFT(2)-j+1,:)+Ex_uu(j,:))
  Ex_vv(j,:)=(Ex_vv(N_FFT(2)-j+1,:)+Ex_vv(j,:))
  Ex_ww(j,:)=(Ex_ww(N_FFT(2)-j+1,:)+Ex_ww(j,:))
  Ex_pp(j,:)=(Ex_pp(N_FFT(2)-j+1,:)+Ex_pp(j,:))
  Ez_uu(j,:)=(Ez_uu(N_FFT(2)-j+1,:)+Ez_uu(j,:))
  Ez_vv(j,:)=(Ez_vv(N_FFT(2)-j+1,:)+Ez_vv(j,:))
  Ez_ww(j,:)=(Ez_ww(N_FFT(2)-j+1,:)+Ez_ww(j,:))
  Ez_pp(j,:)=(Ez_pp(N_FFT(2)-j+1,:)+Ez_pp(j,:))
  MS_t(j)=(MS_t(N_FFT(2)-j+1)-MS_t(j))
  M_t(j,:)=(M_t(N_FFT(2)-j+1,:)+M_t(j,:))
END DO
M_t=M_t/((nArgs-1)*N_FFT(1)*N_FFT(3)*2)
MS_PSD(:,1)=Ex_uu(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,1))**2
MS_PSD(:,2)=Ex_vv(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,2))**2
MS_PSD(:,3)=Ex_ww(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,3))**2
MS_PSD(:,4)=Ex_pp(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,4))**2
MS_t(:)=  MS_t(:)  /((nArgs-1)*N_FFT(1)*N_FFT(3)*2)-(M_t(:,1)*M_t(:,2))
!TWO Times meansquare=amplituden quadrate=>Mittelung entlang z und files
Ex_uu(:,2:nSamples_Spec)=Ex_uu(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)
Ex_vv(:,2:nSamples_Spec)=Ex_vv(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)
Ex_ww(:,2:nSamples_Spec)=Ex_ww(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)
Ex_pp(:,2:nSamples_Spec)=Ex_pp(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)
Ez_uu(:,2:nSamples_Spec)=Ez_uu(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)
Ez_vv(:,2:nSamples_Spec)=Ez_vv(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)
Ez_ww(:,2:nSamples_Spec)=Ez_ww(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)
Ez_pp(:,2:nSamples_Spec)=Ez_pp(:,2:nSamples_Spec)/(N_FFT(3)*REAL(nArgs-1)*2)

!Write Spectra
FileUnit_EK=155
INQUIRE(UNIT=FileUnit_EK, OPENED=connected)
IF (Connected) THEN
  DO
    FileUnit_EK=Fileunit_EK+1
    INQUIRE(UNIT=FileUnit_EK, OPENED=connected)
    IF (.NOT. connected) EXIT
  END DO
END IF
!SWRITE(UNIT_stdOut,*)'------------------------------'
!SWRITE(UNIT_stdOut,*)' Max. relative error in RMS: ', SQRT(maxdev),RMS_PSD,RMS_t
!SWRITE(UNIT_stdOut,*)' MEANS: ', E_uu(N(2),1),M_t**2
FileName_EK=TIMESTAMP(TRIM(ProjectName)//'MS',time)
FileName_EK=TRIM(Filename_EK)//'.dat'
OPEN(FileUnit_Ek,FILE=Filename_EK,STATUS="REPLACE")
WRITE(FileUnit_EK,*)'TITLE     = "MeanSquares "'
WRITE(FileUnit_EK,'(a)')'VARIABLES ="yPlus"'
WRITE(FileUnit_EK,'(a)')'"uu"'
WRITE(FileUnit_EK,'(a)')'"vv"'
WRITE(FileUnit_EK,'(a)')'"ww"'
WRITE(FileUnit_EK,'(a)')'"uv"'
WRITE(FileUnit_EK,'(a)')'"u_mean"'
WRITE(FileUnit_EK,'(a)')'ZONE T="ZONE1"'
WRITE(FileUnit_EK,'(a)')' STRANDID=0, SOLUTIONTIME=0'
WRITE(FileUnit_EK,*)' I=',N_FFT(2)/2,', J=1, K=1, ZONETYPE=Ordered'
WRITE(FileUnit_EK,'(a)')' DATAPACKING=POINT'
WRITE(FileUnit_EK,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
DO j=N_FFT(2)/2+1,N_FFT(2)
  WRITE(FileUnit_EK,'(6(E20.12,X))')(1-(X_FFT(2,1,j,1)))*Re_tau,MS_PSD(j,1),MS_PSD(j,2),MS_PSD(j,3),MS_t(j),m_t(j,1)
END DO
CLOSE(FILEUnit_EK)
!-------------------------------------------------
FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_EnergySpectra_x',time)
FileName_EK=TRIM(Filename_EK)//'.dat'
OPEN(FileUnit_Ek,FILE=Filename_EK,STATUS="REPLACE")
WRITE(FileUnit_EK,*)'TITLE     = "Energy Spectra_x "'
WRITE(FileUnit_EK,'(a)')'VARIABLES ="k"'
WRITE(FileUnit_EK,'(a)')'"E_uu_x"'
WRITE(FileUnit_EK,'(a)')'"E_vv_x"'
WRITE(FileUnit_EK,'(a)')'"E_ww_x"'
WRITE(FileUnit_EK,'(a)')'"E_pp_x"'
DO j=N_FFT(2)/2+1,N_FFT(2)
  WRITE(FileUnit_EK,'(A14,I4,A1)')'ZONE T="E_xx,yPlus=',INT((1-ABS(X_FFT(2,1,j,1)))*Re_tau),'"'
  WRITE(FileUnit_EK,'(a)')' STRANDID=0, SOLUTIONTIME=0'
  WRITE(FileUnit_EK,*)' I=',nSamples_spec,', J=1, K=1, ZONETYPE=Ordered'
  WRITE(FileUnit_EK,'(a)')' DATAPACKING=POINT'
  WRITE(FileUnit_EK,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
  DO k=1,nSamples_spec
    WRITE(FileUnit_EK,'((I5.5,X),4(E20.12,X))')k-1,Ex_uu(j,k),Ex_vv(j,k),Ex_ww(j,k),Ex_pp(j,k)
  END DO
END DO
CLOSE(FILEUnit_EK)
!-------------------------------------------------
FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_EnergySpectra_z',time)
FileName_EK=TRIM(Filename_EK)//'.dat'
OPEN(FileUnit_Ek,FILE=Filename_EK,STATUS="REPLACE")
WRITE(FileUnit_EK,*)'TITLE     = "Energy Spectra_z "'
WRITE(FileUnit_EK,'(a)')'VARIABLES ="k"'
WRITE(FileUnit_EK,'(a)')'"E_uu_z"'
WRITE(FileUnit_EK,'(a)')'"E_vv_z"'
WRITE(FileUnit_EK,'(a)')'"E_ww_z"'
WRITE(FileUnit_EK,'(a)')'"E_pp_z"'
DO j=N_FFT(2)/2+1,N_FFT(2)
  WRITE(FileUnit_EK,'(A14,I4,A1)')'ZONE T="E_zz,yPlus=',INT((1-ABS(X_FFT(2,1,j,1)))*Re_tau),'"'
  WRITE(FileUnit_EK,'(a)')' STRANDID=0, SOLUTIONTIME=0'
  WRITE(FileUnit_EK,*)' I=',nSamples_spec,', J=1, K=1, ZONETYPE=Ordered'
  WRITE(FileUnit_EK,'(a)')' DATAPACKING=POINT'
  WRITE(FileUnit_EK,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
  DO k=1,nSamples_spec
    WRITE(FileUnit_EK,'((I5.5,X),4(E20.12,X))')2*(k-1),Ez_uu(j,k),Ez_vv(j,k),Ez_ww(j,k),Ez_pp(j,k)
  END DO
END DO
CLOSE(FILEUnit_EK)

END SUBROUTINE FFTOutput


!===================================================================================================================================
!> Finalize FFT
!===================================================================================================================================
SUBROUTINE FinalizeFFT()
! MODULES                                                                                                                          !
USE MOD_FFT_Vars
USE MOD_DG_Vars,     ONLY: U
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(U)
SDEALLOCATE(U_FFT)
SDEALLOCATE(X_FFT)
SDEALLOCATE(Ex_uu)
SDEALLOCATE(Ex_vv)
SDEALLOCATE(Ex_pp)
SDEALLOCATE(Ez_uu)
SDEALLOCATE(Ez_vv)
SDEALLOCATE(Ez_pp)
SDEALLOCATE(in)
SDEALLOCATE(out)
SDEALLOCATE(VdmGaussEqui)
SDEALLOCATE(MS_t)
SDEALLOCATE(MS_PSD)
SDEALLOCATE(M_t)
END SUBROUTINE FinalizeFFT

END MODULE MOD_FFT
