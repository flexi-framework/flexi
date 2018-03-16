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

INTERFACE Rogallo
  MODULE PROCEDURE Rogallo
END INTERFACE

INTERFACE FinalizeFFT
  MODULE PROCEDURE FinalizeFFT
END INTERFACE

PUBLIC:: InitFFT,Rogallo,FinalizeFFT

CONTAINS

!===================================================================================================================================
!> Initialize FFT. Define necessary parameters, allocate arrays for solution and auxiliary FFT variables.
!===================================================================================================================================
SUBROUTINE InitFFT()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Init_Hit_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: i,j,k
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT BASIS...'

! No Spectral code without pi and 2pi
Pi     = ACOS(-1.)
TwoPi  = 2.*Pi
II     = cmplx(0.,1.0)
kmax   = NINT(sqrt(REAL(3*N_FFT**2)))+1
scalefactor = 1/N_FFT**3

! ALLOCATE the solution array Uloc
ALLOCATE(Uloc(1:PP_nVar,1:N_FFT,1:N_FFT,1:N_FFT))
! allocate an array for the physical position of the nodes per proc
ALLOCATE(Localxyz(1:3,1:N_FFT,1:N_FFT,1:N_FFT))
Nc=FLOOR(REAL(N_FFT)/2.)
Endw(1)=Nc+1
Endw(2:3)=N_FFT
ALLOCATE(Localk(1:4,1:Endw(1),1:Endw(2),1:Endw(3)))


! fill the wave space 
DO i=1,Endw(1)
  DO j=1,Endw(2)
    DO k=1,Endw(3)
      ! due to symmetry, only nx/2+1 kx are stored, ranging from 0..1..2... to Nx/2 (Nx/2+1 waves) e.g. Nx=7: k1=:0,1,2,3; Nx=4:0,1,2 
      localk(1,i,j,k)=i-1
      ! ny waves are stored, order: 0,1,2,..,Nyquist,..,..,-2,-1, e.g Ny=4: 0,1,2,-1, Ny=7: 0,1,2,3,-3,-2,-1
      IF (j.LE.Nc+1) THEN
        localk(2,i,j,k)=j-1
      ELSE 
        localk(2,i,j,k)=-(N_FFT+1-j)
      END IF
      IF (k.LE.Nc+1) THEN
        localk(3,i,j,k)=k-1
      ELSE 
        localk(3,i,j,k)=-(N_FFT+1-k)
      END IF
      localk(4,i,j,k)=NINT(sqrt(REAL(localk(1,i,j,k)*localk(1,i,j,k)&
                                    +localk(2,i,j,k)*localk(2,i,j,k)+localk(3,i,j,k)*localk(3,i,j,k))))
    END DO 
  END DO 
END DO

DO i=1,N_FFT
  DO j=1,N_FFT
    DO k=1,N_FFT
      localxyz(1,i,j,k)=REAL((i-1))/REAL(N_FFT)
      localxyz(2,i,j,k)=REAL((j-1))/REAL(N_FFT)
      localxyz(3,i,j,k)=REAL((k-1))/REAL(N_FFT)
    END DO
  END DO
END DO


! Prepare local Array for result of Forward FFT (so for wavespace results) 
ALLOCATE(U_FFT(1:PP_nVar,1:Endw(1),1:Endw(2),1:Endw(3)))
ALLOCATE(F_vv(1:3,1:3,1:endw(1),1:endw(2),1:endw(3)))
ALLOCATE(fhat(1:3,1:endw(1),1:endw(2),1:endw(3)))
ALLOCATE(phat(1:endw(1),1:endw(2),1:endw(3)))
GlobalMeshOffset(1:3)=0

SWRITE(UNIT_stdOut,'(A)')' INIT BASIS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitFFT

!===================================================================================================================================
! Computes the initial energy spectrum and creates a velocity field with random phase matching the energy spectrum 
!===================================================================================================================================
SUBROUTINE Rogallo
! MODULES
USE MOD_Globals
USE MOD_Init_Hit_Vars,      ONLY: endw,localk,U_FFT,TwoPi,II,Uloc,N_FFT,kmax,Nc,InitSpec,plan
USE FFTW3
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k
REAL    :: Theta1,Theta2,Phi,kw,kw2,kk(1:3),r,Nenner,a1,a2,a3,a4,a0,E11,specscale,k0,kp,u0
REAL    :: a,e,kn,ke,nu,l
COMPLEX :: alpha,beta
REAL,ALLOCATABLE:: E_P(:)
INTEGER,ALLOCATABLE:: kperShell(:)
integer, allocatable :: seed(:)
integer :: nn,seed1
!===================================================================================================================================
ALLOCATE(kperShell(0:kmax+1))
kperShell=0
DO i=-N_FFT,N_FFT; DO j=-N_FFT,N_FFT; DO k=-N_FFT,N_FFT   
  r=SQRT(REAL(i**2+j**2+k**2))
  kperShell(NINT(r))=kperShell(NINT(r))+1
END DO; END DO; END DO

ALLOCATE(E_P(0:kmax+3))
E_P=0.
SELECT CASE(InitSpec)
 CASE(1)
   a0=4.7935398
   a1=-1.3284141
   a2=-0.2146974
   a3=-0.0314604
   a4=-0.0169870
   DO i=1,kmax+3
     r=REAL(i)
     ! adding a scaling of 0.01...
     E11=0.001*EXP(a0+a1*LOG(r)+a2*(LOG(r))**2+a3*(LOG(r))**3+a4*LOG(r)**4)
     E_P(i)=E11*(0.5*(a1+2*a2*LOG(r)+3*a3*LOG(r)**2+4*a4*LOG(r)**3)**2 + a2 - a1 +(3*a3-2*a2)*LOG(r) +(6*a4-3*a3)*(LOG(r))**2&
    -4*a4*(LOG(r))**4)
   END DO
 CASE(2) ! blaisdell
   specscale=0.01; k0=6
   DO i=1,kmax+3
     E_p(i)= specscale*(i**4*EXP(-2.*(i/k0)**2))
   END DO
 CASE(3) !Chasnov
   write(*,*)"Chasnov"
   a0=5.319230405352436e-01 !for s=4 according to Batchelor-Proudman flow
   kp=4  !to be chosen, scaling Re
   u0=5. !scaling the total Energy
   a1=a0*u0**2.*kp**(-5.)  !Batchelor-Proudman again
   print*, 'specscale',a1
   DO i=1,kmax+3
     E_p(i)= a1*i**4.*EXP(-2.*(i/kp)**2.)
   END DO
 CASE(4) !inf inertial range
   write(*,*)'Infinite inertial range spectrum k^(-5/3)'
   DO i=1,kmax+3
     E_p(i)= i**(-5/3.) 
  END DO
 CASE(5) ! karman-pao
   a = 1.453 !scaling const Bailly 99
   u0= 0.3 !rms of u
   ke= 2.!related to peak of E wavenumber, w ~~ sqrt(12/5) ke
   nu= 5e-4
   L=  0.746834/ke
   e=  u0**3/L
   kn= e**0.25*nu**(-0.75) !kolmogorov wavenumber e^1/4 nu^3/4 : zB nu=5e-4; e~~u0/L L~~0.746834/ke ke=2. u0=1. e=
   DO i=1,kmax+3
     E_p(i) = a * u0**2/ke *(i/ke)**4/(1+(i/ke)**2)**(17/6.)*EXP(-2*(i/kn)**2.)
   END DO
END SELECT


call random_seed(size = nn)
allocate(seed(nn))
call system_clock(count=seed1)
seed = seed1
call random_seed(put=seed)
call random_seed(get=seed)
write (*, *)'seed',seed

DO i=1,Endw(1); DO j=1,Endw(2); DO k=1,Endw(3)
  IF ((((localk(4,i,j,k)+1).LE.Kmax+3).AND.(localk(4,i,j,k).GT.0)) .AND.&
              ((localk(1,i,j,k).GT.0 .OR. localk(2,i,j,k).GT.0))) THEN
    CALL RANDOM_NUMBER(Theta1)
    CALL RANDOM_NUMBER(Theta2)
    CALL RANDOM_NUMBER(Phi)
    Theta1=Theta1*TwoPi
    Theta2=Theta2*TwoPi
    Phi=Phi*TwoPi
    kk(1:3)=REAL(localk(1:3,i,j,k))
    kw2=kk(1)**2+kk(2)**2+kk(3)**2
    kw=sqrt(kw2)
    nenner=(Kpershell(NINT(kw))-0)/2.!2*TwoPi*kw2
    Alpha =SQRT(E_p(NINT(kw))/(Nenner))*EXP(II*Theta1)*COS(PHI)
    Beta  =SQRT(E_p(NINT(kw))/(Nenner))*EXP(II*Theta2)*SIN(PHI)

    U_FFT(2,i,j,k)=(alpha*kw*kk(2)+beta*kk(1)*kk(3))/(kw*SQRT(kk(1)**2+kk(2)**2))
    U_FFT(3,i,j,k)=(beta*kk(2)*kk(3)-alpha*kw*kk(1))/(kw*SQRT(kk(1)**2+kk(2)**2))
    U_FFT(4,i,j,k)=-1./kw*(beta*SQRT(kk(1)**2+kk(2)**2))
  ELSE
    U_FFT(2:4,i,j,k)=0.0
  END IF
END DO; END DO; END DO

DEALLOCATE(E_P,kperShell,seed)
!2/3 Filter for clean incompressible data
DO i=1,Endw(1); DO j=1,Endw(2); DO k=1,Endw(3) 
  IF(localk(4,i,j,k).GE.INT(2/3.*Nc)) U_FFT(:,i,j,k) = 0.
END DO; END DO; END DO

Uloc=1.
DO i=2,4
  CALL DFFTW_PLAN_DFT_C2R_3D(plan,N_FFT,N_FFT,N_FFT,U_FFT(i,:,:,:),Uloc(i,:,:,:),FFTW_ESTIMATE)
  CALL DFFTW_Execute(plan,U_FFT(i,:,:,:),Uloc(i,:,:,:))
END DO
! set constant density
Uloc(1,:,:,:) = 1.0
! compute rho*v
Uloc(2:4,:,:,:)= 1.0*Uloc(2:4,:,:,:)
CALL Compute_incompressible_P()
END SUBROUTINE Rogallo


SUBROUTINE Compute_incompressible_P
!===================================================================================================================================
! Transformation from primitive to conservative variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Init_Hit_Vars,      ONLY: endw,localk,TwoPi,II,Uloc,N_FFT,plan,scalefactor,F_vv,fhat,phat
USE FFTW3
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                 :: sKappaM1, Kappa, Mach,ksquared,vmax,p0
REAL                 :: vv(1:3,1:3,1:N_FFT,1:N_FFT,1:N_FFT)
REAL                 :: v(1:3,1:N_FFT,1:N_FFT,1:N_FFT)
REAL                 :: p(1:N_FFT,1:N_FFT,1:N_FFT)
INTEGER              :: i,j,k
!===================================================================================================================================
sKappaM1=1/0.4 
Kappa = 1.4 
Mach = 0.1

phat=0.;vmax=0.
DO i=1,N_FFT; DO j=1,N_FFT; DO k=1,N_FFT
  V(1:3,i,j,k) = Uloc(2:4,i,j,k)/Uloc(1,i,j,k)
  vmax=MAX(vmax,SQRT(V(1,i,j,k)**2+V(2,i,j,k)**2+V(3,i,j,k)**2))
END DO; END DO; END DO! i,j,k=1,N_FFT

DO j=1,3; DO k=1,3
    vv(j,k,:,:,:)=v(j,:,:,:)*v(k,:,:,:)
END DO; END DO

DO j=1,3; DO k=1,3
  CALL DFFTW_PLAN_DFT_R2C_3D(plan,N_FFT,N_FFT,N_FFT,vv(j,k,:,:,:),F_vv(j,k,:,:,:),FFTW_ESTIMATE)
  CALL DFFTW_Execute(plan,vv(j,k,:,:,:),F_vv(j,k,:,:,:))
END DO; END DO

DO j=1,3; DO k=1,3
  F_vv(j,k,:,:,:) = F_vv(j,k,:,:,:) *Scalefactor*TwoPi*II*localk(k,:,:,:)
END DO; END DO

fhat(1,:,:,:)= F_vv(1,1,:,:,:)+F_vv(1,2,:,:,:)+F_vv(1,3,:,:,:)
fhat(2,:,:,:)= F_vv(2,1,:,:,:)+F_vv(2,2,:,:,:)+F_vv(2,3,:,:,:)
fhat(3,:,:,:)= F_vv(3,1,:,:,:)+F_vv(3,2,:,:,:)+F_vv(3,3,:,:,:)

DO i=1,Endw(1); DO j=1,Endw(2); DO k=1,Endw(3)
  ksquared=localk(1,i,j,k)**2+localk(2,i,j,k)**2+localk(3,i,j,k)**2
  IF (ksquared.eq.0) THEN
     Phat(i,j,k)=0
  ELSE
     phat(i,j,k)=1./ksquared*(localk(1,i,j,k)*fhat(1,i,j,k)+&
                          localk(2,i,j,k)*fhat(2,i,j,k)+&
                          localk(3,i,j,k)*fhat(3,i,j,k))
  END IF
END DO; END DO; END DO! i,j,k=1,End(1,2,3)

CALL DFFTW_PLAN_DFT_C2R_3D(plan,N_FFT,N_FFT,N_FFT,phat(:,:,:),p(:,:,:),FFTW_ESTIMATE)
CALL DFFTW_Execute(plan,phat(:,:,:),p(:,:,:))

SWRITE(*,*) "Vmax in field is",vmax
p0=1.*vmax**2/(Kappa*Mach**2)
SWRITE(*,*) "For a selected Machnumber of ",Mach,", the mean pressure is",p0
! add mean pressure to fluctuations
p=p+p0
  
DO i=1,N_FFT; DO j=1,N_FFT; DO k=1,N_FFT
  Uloc(5,i,j,k)=sKappaM1*p(i,j,k)+0.5*SUM(Uloc(2:4,i,j,k)*v(1:3,i,j,k))
END DO; END DO; END DO! i,j,k=1,N_FFT
    
END SUBROUTINE Compute_incompressible_P

!===================================================================================================================================
!> Finalize FFT
!===================================================================================================================================
SUBROUTINE FinalizeFFT()
! MODULES                                                                                                                          !
USE MOD_Init_Hit_Vars
USE MOD_DG_Vars,       ONLY: U
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(U)
SDEALLOCATE(Uloc)
SDEALLOCATE(Uloc_c)
SDEALLOCATE(U_j)
SDEALLOCATE(U_k)
SDEALLOCATE(U_FFT)
SDEALLOCATE(LocalXYZ)
SDEALLOCATE(LocalK)
SDEALLOCATE(phat)
SDEALLOCATE(fhat)
SDEALLOCATE(F_vv)
END SUBROUTINE FinalizeFFT

END MODULE MOD_FFT
