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

MODULE MOD_INIT_HIT
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE Init_InitHIT
  MODULE PROCEDURE Init_InitHIT
END INTERFACE

INTERFACE Rogallo
  MODULE PROCEDURE Rogallo
END INTERFACE

INTERFACE Finalize_InitHIT
  MODULE PROCEDURE Finalize_InitHIT
END INTERFACE

PUBLIC:: Init_InitHIT,Rogallo,Finalize_InitHIT

CONTAINS

!===================================================================================================================================
!> Initialize FFT. Define necessary parameters, allocate arrays for solution and auxiliary FFT variables.
!===================================================================================================================================
SUBROUTINE Init_InitHit()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Init_Hit_Vars
USE MOD_FFT_Vars,      ONLY: N_FFT,Endw
USE MOD_DG_Vars,       ONLY: U
USE MOD_Mesh_Vars,     ONLY: nElems
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT InitHit...'

! Allocate the solution array Uloc
ALLOCATE(U_FFT(   1:PP_nVar,1:Endw(1),1:Endw(2),1:Endw(3)))
ALLOCATE(U_Global(1:PP_nVar,1:N_FFT,1:N_FFT,1:N_FFT))
ALLOCATE(U(       1:PP_nVar,0:N,0:N,0:N,nElems))
U_FFT    = 0.
U_Global = 0.
U        = 0.

SWRITE(UNIT_stdOut,'(A)')' INIT BASIS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE Init_InitHit

!===================================================================================================================================
! Computes the initial energy spectrum and creates a velocity field with random phase matching the energy spectrum
!===================================================================================================================================
SUBROUTINE Rogallo()
! MODULES
USE MOD_Globals
USE MOD_PreProc,            ONLY: PP_PI
USE MOD_Init_Hit_Vars,      ONLY: U_FFT,U_Global,InitSpec,seed
USE MOD_FFT,                ONLY: ComputeFFT_C2R
USE MOD_FFT_Vars,           ONLY: Endw,localk,N_FFT,Nc,II,kmax
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k
INTEGER :: nn
REAL    :: Theta1,Theta2,Phi,kw,kw2,kk(1:3),r,Nenner,a1,a2,a3,a4,a0,E11,specscale,k0,kp,u0
REAL    :: a,e,kn,ke,nu,l
COMPLEX :: alpha,beta
REAL    :: E_P(0:kmax+3)
INTEGER :: kperShell(0:kmax+1)
INTEGER,ALLOCATABLE :: seed_loc(:)
!===================================================================================================================================
! Choose energy spectrum
SELECT CASE(InitSpec)
 CASE(1) ! Rogallo
   SWRITE(*,*) "SPECTRUM: Rogallo"
   a0 =  4.7935398
   a1 = -1.3284141
   a2 = -0.2146974
   a3 = -0.0314604
   a4 = -0.0169870
   DO i=1,kmax+3
     r=REAL(i)
     ! adding a scaling of 0.01...
     E11=0.001*EXP(a0+a1*LOG(r)+a2*(LOG(r))**2+a3*(LOG(r))**3+a4*LOG(r)**4)
     E_P(i)=E11*(0.5*(a1+2*a2*LOG(r)+3*a3*LOG(r)**2+4*a4*LOG(r)**3)**2 + a2 - a1 +(3*a3-2*a2)*LOG(r) +(6*a4-3*a3)*(LOG(r))**2&
    -4*a4*(LOG(r))**4)
   END DO

 CASE(2) ! Blaisdell
   SWRITE(*,*) "SPECTRUM: Blaisdell"
   specscale=0.01
   k0=6
   DO i=1,kmax+3
     E_p(i)= specscale*(i**4*EXP(-2.*(i/k0)**2))
   END DO

 CASE(3) ! Chasnov
   SWRITE(*,*) "SPECTRUM: Chasnov"
   a0=5.319230405352436e-01 ! for s=4 according to Batchelor-Proudman flow
   kp=4  ! to be chosen, scaling Re
   u0=5. ! scaling the total Energy
   a1=a0*u0**2.*kp**(-5.)  ! Batchelor-Proudman again
   DO i=1,kmax+3
     E_p(i)= a1*i**4.*EXP(-2.*(i/kp)**2.)
   END DO

 CASE(4) ! inf inertial range
   SWRITE(*,*) "SPECTRUM: Infinite inertial range spectrum k^(-5/3)"
   DO i=1,kmax+3
     E_p(i)= i**(-5/3.)
  END DO

 CASE(5) ! karman-pao
   SWRITE(*,*) "SPECTRUM: Karman-Pao"
   a  = 1.453 ! scaling const Bailly 99
   u0 = 0.3   ! rms of u
   ke = 2.    ! related to peak of E wavenumber, w ~~ sqrt(12/5) ke
   nu = 5e-4
   L  = 0.746834/ke
   e  = u0**3/L
   kn = e**0.25*nu**(-0.75) ! kolmogorov wavenumber e^1/4 nu^3/4 : zB nu=5e-4; e~~u0/L L~~0.746834/ke ke=2. u0=1. e=
   DO i=1,kmax+3
     E_p(i) = a * u0**2/ke *(i/ke)**4/(1+(i/ke)**2)**(17/6.)*EXP(-2*(i/kn)**2.)
   END DO

END SELECT

! Get random seed from system clock or from parameter file if specified
CALL RANDOM_SEED(SIZE = nn)
ALLOCATE(seed_loc(nn))
IF (seed.EQ.0) CALL SYSTEM_CLOCK(COUNT=seed)
seed_loc = seed
CALL RANDOM_SEED(put=seed_loc)
CALL RANDOM_SEED(get=seed_loc)
SWRITE (*,*) 'Seed:',seed
DEALLOCATE(seed_loc)

kperShell=0
DO i=-N_FFT,N_FFT; DO j=-N_FFT,N_FFT; DO k=-N_FFT,N_FFT
  r=SQRT(REAL(i**2+j**2+k**2))
  kperShell(NINT(r))=kperShell(NINT(r))+1
END DO; END DO; END DO

DO i=1,Endw(1); DO j=1,Endw(2); DO k=1,Endw(3)
  IF ((((localk(4,i,j,k)+1).LE.Kmax+3).AND.(localk(4,i,j,k).GT.0)) .AND.&
              ((localk(1,i,j,k).GT.0 .OR. localk(2,i,j,k).GT.0))) THEN
    CALL RANDOM_NUMBER(Theta1)
    CALL RANDOM_NUMBER(Theta2)
    CALL RANDOM_NUMBER(Phi)
    Theta1=Theta1*2.*PP_PI
    Theta2=Theta2*2.*PP_PI
    Phi=Phi*2.*PP_PI
    kk(1:3)=REAL(localk(1:3,i,j,k))
    kw2=kk(1)**2+kk(2)**2+kk(3)**2
    kw=SQRT(kw2)
    Nenner = (Kpershell(NINT(kw))-0)/2.!2*TwoPi*kw2
    Alpha  = SQRT(E_p(NINT(kw))/(Nenner))*EXP(II*Theta1)*COS(PHI)
    Beta   = SQRT(E_p(NINT(kw))/(Nenner))*EXP(II*Theta2)*SIN(PHI)

    U_FFT(2,i,j,k)=(alpha*kw*kk(2)+beta*kk(1)*kk(3))/(kw*SQRT(kk(1)**2+kk(2)**2))
    U_FFT(3,i,j,k)=(beta*kk(2)*kk(3)-alpha*kw*kk(1))/(kw*SQRT(kk(1)**2+kk(2)**2))
    U_FFT(4,i,j,k)=-1./kw*(beta*SQRT(kk(1)**2+kk(2)**2))
  ELSE
    U_FFT(2:4,i,j,k)=0.0
  END IF
END DO; END DO; END DO

! 2/3 Filter for clean incompressible data
DO k=1,Endw(3); DO j=1,Endw(2); DO i=1,Endw(1)
  IF(localk(4,i,j,k).GE.INT(2/3.*Nc)) U_FFT(:,i,j,k) = 0.
END DO; END DO; END DO

! Pull back to phyiscal space and normalize
CALL ComputeFFT_C2R(3,U_FFT(2:4,:,:,:),U_Global(2:4,:,:,:))

! Compute compressible state from incompressible Rogallo data
CALL Compute_Compressible_State(U_Global)

END SUBROUTINE Rogallo


!===================================================================================================================================
! Compute thermodymamically consistent compressible state for incompressible input state
!===================================================================================================================================
SUBROUTINE Compute_Compressible_State(U_In)
! MODULES
USE MOD_Globals
USE MOD_PreProc,            ONLY: PP_PI
USE MOD_FFT_Vars,           ONLY: Endw,localk,II,N_FFT
USE MOD_FFT,                ONLY: ComputeFFT_R2C,ComputeFFT_C2R
#if USE_OPENMP
USE OMP_Lib
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: U_In(1:5,1:N_FFT,1:N_FFT,1:N_FFT)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER       :: Ma0   = 0.1
REAL,PARAMETER       :: rho0  = 1.0
REAL,PARAMETER       :: Kappa = 1.4
REAL                 :: ksquared,vmax,p0
REAL                 :: p(       1,1:N_FFT,1:N_FFT,1:N_FFT)
REAL                 :: v(     1:3,1:N_FFT,1:N_FFT,1:N_FFT)
REAL                 :: vv(1:3,1:3,1:N_FFT,1:N_FFT,1:N_FFT)
COMPLEX              :: p_c(       1,1:Endw(1),1:Endw(2),1:Endw(3))
COMPLEX              :: f_c(     1:3,1:Endw(1),1:Endw(2),1:Endw(3))
COMPLEX              :: vv_c(1:3,1:3,1:Endw(1),1:Endw(2),1:Endw(3))
INTEGER              :: i,j,k
!===================================================================================================================================
! Save primitive velocity and compute maximum velocity in field
vmax=0.
DO i=1,N_FFT; DO j=1,N_FFT; DO k=1,N_FFT
  v(1:3,i,j,k) = U_In(2:4,i,j,k)
  vmax=MAX(vmax,SQRT(v(1,i,j,k)**2+v(2,i,j,k)**2+v(3,i,j,k)**2))
END DO; END DO; END DO! i,j,k=1,N_FFT

! Compute tensor product of velocities, perform FFT and normalize
DO j=1,3
  DO k=1,3
    vv(j,k,:,:,:)=v(j,:,:,:)*v(k,:,:,:)
  END DO
  CALL ComputeFFT_R2C(3,vv(j,:,:,:,:),vv_c(j,:,:,:,:))
  DO k=1,3
    vv_c(j,k,:,:,:) = vv_c(j,k,:,:,:)*2.*PP_PI*II*localk(k,:,:,:)
  END DO
END DO

! Sum up velocity correlations in all three directions
f_c(1,:,:,:) = vv_c(1,1,:,:,:)+vv_c(1,2,:,:,:)+vv_c(1,3,:,:,:)
f_c(2,:,:,:) = vv_c(2,1,:,:,:)+vv_c(2,2,:,:,:)+vv_c(2,3,:,:,:)
f_c(3,:,:,:) = vv_c(3,1,:,:,:)+vv_c(3,2,:,:,:)+vv_c(3,3,:,:,:)

! Compute pressure fluctuations from velocity fluctuations in Fourier space
DO i=1,Endw(1); DO j=1,Endw(2); DO k=1,Endw(3)
  ksquared=localk(1,i,j,k)**2+localk(2,i,j,k)**2+localk(3,i,j,k)**2
  IF (ksquared.EQ.0) THEN
    p_c(1,i,j,k)=0.
  ELSE
    p_c(1,i,j,k)=1./ksquared*(  localk(1,i,j,k)*f_c(1,i,j,k) &
                              + localk(2,i,j,k)*f_c(2,i,j,k) &
                              + localk(3,i,j,k)*f_c(3,i,j,k) )
  END IF
END DO; END DO; END DO! i,j,k=1,Endw(1,2,3)

! Transform pressure fluctuations from Fourier to physical space
CALL ComputeFFT_C2R(1,p_c(:,:,:,:),p(:,:,:,:))
p=p/REAL(N_FFT**3)**3   ! Correct normalization necessary

! Compute mean pressure from Mach number and add to fluctuations
p0 = vmax**2/(Kappa*Ma0**2)
p=p+p0
SWRITE(Unit_StdOut,'(A,F4.2,A,F7.1,A,F4.2)') 'For the selected Mach number ',Ma0, &
                                             ', the mean pressure is ',p0,', with Vmax in field ',vmax

! Compute compressible state with ideal gas EOS
DO i=1,N_FFT; DO j=1,N_FFT; DO k=1,N_FFT
  U_In(  1,i,j,k) = rho0
  U_In(2:4,i,j,k) = U_In(1,i,j,k)*U_In(2:4,i,j,k)
  U_In(  5,i,j,k) = p(1,i,j,k)/(Kappa-1.)+0.5*SUM(U_In(2:4,i,j,k)*v(1:3,i,j,k))
END DO; END DO; END DO

END SUBROUTINE Compute_Compressible_State

!===================================================================================================================================
!> Finalize FFT
!===================================================================================================================================
SUBROUTINE Finalize_InitHIT
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
SDEALLOCATE(U_FFT)
SDEALLOCATE(U_Global)
END SUBROUTINE Finalize_InitHIT

END MODULE MOD_INIT_HIT
