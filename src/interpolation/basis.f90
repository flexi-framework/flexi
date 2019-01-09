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
!> Routines to provide and evaluate basis function coefficients, or provide fast interpolation coefficients
!==================================================================================================================================
MODULE MOD_Basis
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE BuildLegendreVdm
   MODULE PROCEDURE BuildLegendreVdm
END INTERFACE

INTERFACE InitializeVandermonde
   MODULE PROCEDURE InitializeVandermonde
END INTERFACE

INTERFACE ChebyshevGaussNodesAndWeights
   MODULE PROCEDURE ChebyshevGaussNodesAndWeights
END INTERFACE

INTERFACE ChebyGaussLobNodesAndWeights
   MODULE PROCEDURE ChebyGaussLobNodesAndWeights
END INTERFACE

INTERFACE ClenshawCurtisNodesAndWeights
   MODULE PROCEDURE ClenshawCurtisNodesAndWeights
END INTERFACE

INTERFACE LegendreGaussNodesAndWeights
   MODULE PROCEDURE LegendreGaussNodesAndWeights
END INTERFACE

INTERFACE LegGaussLobNodesAndWeights
   MODULE PROCEDURE LegGaussLobNodesAndWeights
END INTERFACE

INTERFACE LegendrePolynomialAndDerivative
   MODULE PROCEDURE LegendrePolynomialAndDerivative
END INTERFACE

INTERFACE PolynomialDerivativeMatrix
   MODULE PROCEDURE PolynomialDerivativeMatrix
END INTERFACE

INTERFACE BarycentricWeights
   MODULE PROCEDURE BarycentricWeights
END INTERFACE

INTERFACE LagrangeInterpolationPolys
   MODULE PROCEDURE LagrangeInterpolationPolys
END INTERFACE

PUBLIC::BuildLegendreVdm
PUBLIC::InitializeVandermonde
PUBLIC::LegGaussLobNodesAndWeights
PUBLIC::LegendreGaussNodesAndWeights
PUBLIC::ChebyshevGaussNodesAndWeights
PUBLIC::ChebyGaussLobNodesAndWeights
PUBLIC::ClenshawCurtisNodesAndWeights
PUBLIC::LegendrePolynomialAndDerivative
PUBLIC::PolynomialDerivativeMatrix
PUBLIC::BarycentricWeights
PUBLIC::LagrangeInterpolationPolys
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Build a 1D Vandermonde matrix from an orthonormal Legendre basis to a nodal basis and reverse
!==================================================================================================================================
SUBROUTINE buildLegendreVdm(N_In,xi_In,Vdm_Leg,sVdm_Leg)
! MODULES
USE MOD_Globals,ONLY:abort
USE MOD_PreProc,ONLY:PP_RealTolerance
USE MOD_Mathtools,ONLY:INVERSE
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_In                    !< input polynomial degree
REAL,INTENT(IN)    :: xi_In(0:N_In)           !< nodal positions [-1,1]
REAL,INTENT(OUT)   ::  Vdm_Leg(0:N_In,0:N_In) !< Vandermonde from Legendre to nodal basis
REAL,INTENT(OUT)   :: sVdm_Leg(0:N_In,0:N_In) !< Vandermonde from nodal basis to Legendre
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j
REAL               :: dummy
!REAL               :: wBary_Loc(0:N_In)
!REAL               :: xGauss(0:N_In),wGauss(0:N_In)
!==================================================================================================================================
! Alternative to matrix inversion: Compute inverse Vandermonde directly
! Direct inversion seems to be more accurate

!CALL BarycentricWeights(N_In,xi_in,wBary_loc)
!! Compute first the inverse (by projection)
!CALL LegendreGaussNodesAndWeights(N_In,xGauss,wGauss)
!!Vandermonde on xGauss
!DO i=0,N_In
!  DO j=0,N_In
!    CALL LegendrePolynomialAndDerivative(j,xGauss(i),Vdm_Leg(i,j),dummy)
!  END DO !i
!END DO !j
!Vdm_Leg=TRANSPOSE(Vdm_Leg)
!DO j=0,N_In
!  Vdm_Leg(:,j)=Vdm_Leg(:,j)*wGauss(j)
!END DO
!!evaluate nodal basis (depends on NodeType, for Gauss: unity matrix)
!CALL InitializeVandermonde(N_In,N_In,wBary_Loc,xi_In,xGauss,sVdm_Leg)
!sVdm_Leg=MATMUL(Vdm_Leg,sVdm_Leg)

!compute the Vandermonde on xGP (Depends on NodeType)
DO i=0,N_In; DO j=0,N_In
  CALL LegendrePolynomialAndDerivative(j,xi_In(i),Vdm_Leg(i,j),dummy)
END DO; END DO !j
sVdm_Leg=INVERSE(Vdm_Leg)
!check (Vdm_Leg)^(-1)*Vdm_Leg := I
dummy=ABS(SUM(ABS(MATMUL(sVdm_Leg,Vdm_Leg)))/(N_In+1.)-1.)
IF(dummy.GT.10.*PP_RealTolerance) CALL abort(__STAMP__,&
                                         'problems in MODAL<->NODAL Vandermonde ',999,dummy)
END SUBROUTINE buildLegendreVdm



!===================================================================================================================================
!> Build a 1D Vandermonde matrix using the Lagrange basis functions of degree
!> N_In, evaluated at the interpolation points xi_Out
!===================================================================================================================================
PPURE SUBROUTINE InitializeVandermonde(N_In,N_Out,wBary_In,xi_In,xi_Out,Vdm)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_In                       !< (IN)  input polynomial degree
INTEGER,INTENT(IN) :: N_Out                      !< (IN)  output polynomial degree
REAL,INTENT(IN)    :: xi_In(0:N_In)              !< (IN)  input nodal positions [-1,1]
REAL,INTENT(IN)    :: xi_Out(0:N_Out)            !< (IN)  outout nodal positions [-1,1]
REAL,INTENT(IN)    :: wBary_In(0:N_In)           !< (IN)  input interpolation weights
REAL,INTENT(OUT)   :: Vdm(0:N_Out,0:N_In)        !< (OUT) nodal Vandermonde from N_In to N_out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iXi
!==================================================================================================================================
DO iXi=0,N_Out
  CALL LagrangeInterpolationPolys(xi_Out(iXi),N_In,xi_In,wBary_In,Vdm(iXi,:)) !l(0:N_In)
END DO
END SUBROUTINE InitializeVandermonde



!===================================================================================================================================
!> Evaluate the Legendre polynomial L_N and its derivative at position x[-1,1]
!> recursive algorithm using the N_in-1 N_in-2 Legendre polynomials
!> algorithm 22, Kopriva book
!===================================================================================================================================
ELEMENTAL SUBROUTINE LegendrePolynomialAndDerivative(N_in,x,L,Lder)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in   !< (IN)  polynomial degree, (N+1) CLpoints
REAL,INTENT(IN)    :: x      !< (IN)  coordinate value in the interval [-1,1]
REAL,INTENT(OUT)   :: L      !< (OUT) Legedre polynomial evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
REAL,INTENT(OUT)   :: Lder   !< (OUT) Legedre polynomial deriv. evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iLegendre
REAL    :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
REAL    :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
IF(N_in .EQ. 0)THEN
  L=1.
  Lder=0.
ELSEIF(N_in .EQ. 1) THEN
  L=x
  Lder=1.
ELSE ! N_in > 1
  L_Nm2=1.
  L_Nm1=x
  Lder_Nm2=0.
  Lder_Nm1=1.
  DO iLegendre=2,N_in
    L=(REAL(2*iLegendre-1)*x*L_Nm1 - REAL(iLegendre-1)*L_Nm2)/REAL(iLegendre)
    Lder=Lder_Nm2 + REAL(2*iLegendre-1)*L_Nm1
    L_Nm2=L_Nm1
    L_Nm1=L
    Lder_Nm2=Lder_Nm1
    Lder_Nm1=Lder
  END DO !iLegendre=2,N_in
END IF ! N_in
!normalize
L=L*SQRT(REAL(N_in)+0.5)
Lder=Lder*SQRT(REAL(N_in)+0.5)
END SUBROUTINE LegendrePolynomialAndDerivative




!==================================================================================================================================
!> Compute Chebychev-Gauss nodes and integration weights (algorithm 27, Kopriva book)
!==================================================================================================================================
PPURE SUBROUTINE ChebyshevGaussNodesAndWeights(N_in,xGP,wGP)
USE MOD_Preproc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in         !< polynomial degree, (N_in+1) CLpoints
REAL,INTENT(OUT)          :: xGP(0:N_in)  !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(OUT),OPTIONAL :: wGP(0:N_in)  !< Gauss point integration weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iGP
!==================================================================================================================================
DO iGP=0,N_in
  xGP(iGP)=-cos((2*iGP+1)/(2*REAL(N_in)+2)*PP_PI)
END DO
IF(PRESENT(wGP))THEN
  DO iGP=0,N_in
    wGP(iGP)=PP_PI/REAL(N_in+1)
  END DO
END IF
END SUBROUTINE ChebyshevGaussNodesAndWeights



!==================================================================================================================================
!> Compute Chebychev-Gauss-Lobatto nodes and integration weights (algorithm 27, Kopriva book)
!==================================================================================================================================
PPURE SUBROUTINE ChebyGaussLobNodesAndWeights(N_in,xGP,wGP)
USE MOD_Preproc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in         !< polynomial degree, (N_in+1) CLpoints
REAL,INTENT(OUT)          :: xGP(0:N_in)  !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(OUT),OPTIONAL :: wGP(0:N_in)  !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP
!==================================================================================================================================
DO iGP=0,N_in
  xGP(iGP)=-COS(iGP/REAL(N_in)*PP_PI)
END DO
IF(PRESENT(wGP))THEN
  DO iGP=0,N_in
    wGP(iGP)=PP_PI/REAL(N_in)
  END DO
  wGP(0)=wGP(0)*0.5
  wGP(N_in)=wGP(N_in)*0.5
END IF
END SUBROUTINE ChebyGaussLobNodesAndWeights


!==================================================================================================================================
!> Compute Clenshaw-Curtis nodes and integration weights
!==================================================================================================================================
PPURE SUBROUTINE ClenshawCurtisNodesAndWeights(N_in,xGP,wGP)
USE MOD_Preproc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in         !< polynomial degree, (N_in+1) CLpoints
REAL,INTENT(OUT)          :: xGP(0:N_in)  !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(OUT),OPTIONAL :: wGP(0:N_in)  !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iGP,j
REAL                      :: b,theta
!==================================================================================================================================
IF(N_in.EQ.0)THEN
  xGP(0) = 0.
  wGP(0) = 2.
ELSE
  DO iGP=0,N_in
    xGP(iGP) = -COS(iGP*PP_PI/REAL(N_in))
  END DO
  xGP(0)   =-1.
  IF(MOD(N_in+1,2).EQ.1)THEN
    xGP(N_in/2)=0.
  END IF
  xGP(N_in)= 1.
  IF(PRESENT(wGP))THEN
    wGP=1.
    DO iGP=0,N_in
      theta = REAL(iGP)*PP_PI/REAL(N_in)
      DO j=1,N_in/2
        b=MERGE(1.,2.,2*j.EQ.N_in)
        wGP(iGP) = wGP(iGP) - b * COS(2.*REAL(j)*theta) / REAL(4*j*j-1)
      END DO
    END DO
    wGP(1:N_in-1)=2.*wGP(1:N_in-1)
    wGP=wGP/REAL(N_in)
  END IF
END IF
END SUBROUTINE ClenshawCurtisNodesAndWeights

!==================================================================================================================================
!> @brief Compute Legendre-Gauss nodes and integration weights (algorithm 23, Kopriva book)
!>
!> Starting with Chebychev point positions, a Newton method is used to find the roots
!> of the Legendre Polynomial L_(N_in+1), which are the positions of Gausspoints
!> uses LegendrePolynomialAndDerivative subroutine
!==================================================================================================================================
SUBROUTINE LegendreGaussNodesAndWeights(N_in,xGP,wGP)
!MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in              !< polynomial degree, (N_in+1) Gausspoints
REAL,INTENT(OUT)          :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(OUT),OPTIONAL :: wGP(0:N_in)       !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: nIter = 10        ! max. number of newton iterations
REAL                      :: Tol   = 1.E-15    ! tolerance for Newton iteration: TODO: use variable tolerance here!
INTEGER                   :: iGP,iter
REAL                      :: L_Np1,Lder_Np1    ! L_{N_in+1},Lder_{N_in+1}
REAL                      :: dx                ! Newton step
REAL                      :: cheb_tmp          ! temporary variable for evaluation of chebychev node positions
!==================================================================================================================================
IF(N_in .EQ. 0) THEN
  xGP=0.
  IF(PRESENT(wGP))wGP=2.
  RETURN
ELSEIF(N_in.EQ.1)THEN
  xGP(0)=-sqrt(1./3.)
  xGP(N_in)=-xGP(0)
  IF(PRESENT(wGP))wGP=1.
  RETURN
ELSE ! N_in>1
  cheb_tmp=2.*atan(1.)/REAL(N_in+1) ! pi/(2N+2)
  DO iGP=0,(N_in+1)/2-1 !since points are symmetric, only left side is computed
    xGP(iGP)=-cos(cheb_tmp*REAL(2*iGP+1)) !initial guess
    ! Newton iteration
    DO iter=0,nIter
      CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
      dx=-L_Np1/Lder_Np1
      xGP(iGP)=xGP(iGP)+dx
      IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
    END DO ! iter
    IF(iter.GT.nIter) THEN
      SWRITE(*,*) 'maximum iteration steps >10 in Newton iteration for Legendre Gausspoint'
      xGP(iGP)=-cos(cheb_tmp*REAL(2*iGP+1)) !initial guess
      ! Newton iteration
      DO iter=0,nIter
        SWRITE(*,*)iter,xGP(iGP)    !DEBUG
        CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
        dx=-L_Np1/Lder_Np1
        xGP(iGP)=xGP(iGP)+dx
        IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
      END DO !iter
      CALL abort(__STAMP__,&
                 'ERROR: Legendre Gauss nodes could not be computed up to desired precision. Code stopped!')
    END IF ! (iter.GT.nIter)
    CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
    xGP(N_in-iGP)=-xGP(iGP)
    IF(PRESENT(wGP))THEN
      !wGP(iGP)=2./((1.-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1) !if Legendre not normalized
      wGP(iGP)=(2.*N_in+3)/((1.-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1)
      wGP(N_in-iGP)=wGP(iGP)
    END IF
  END DO !iGP
END IF ! N_in
IF(mod(N_in,2) .EQ. 0) THEN
  xGP(N_in/2)=0.
  CALL LegendrePolynomialAndDerivative(N_in+1,xGP(N_in/2),L_Np1,Lder_Np1)
  !IF(PRESENT(wGP))wGP(N_in/2)=2./(Lder_Np1*Lder_Np1) !if Legendre not normalized
  IF(PRESENT(wGP))wGP(N_in/2)=(2.*N_in+3)/(Lder_Np1*Lder_Np1)
END IF ! (mod(N_in,2) .EQ. 0)
END SUBROUTINE LegendreGaussNodesAndWeights



!==================================================================================================================================
!> Evaluate the polynomial q=L_{N_in+1}-L_{N_in-1} and its derivative at position x in [-1,1]
!> Recursive algorithm using the N_in-1 N_in-2 Legendre polynomials. (Algorithm 24, Kopriva book)
!==================================================================================================================================
ELEMENTAL SUBROUTINE qAndLEvaluation(N_in,x,q,qder,L)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in                            !< polynomial degree
REAL,INTENT(IN)    :: x                               !< coordinate value in the interval [-1,1]
REAL,INTENT(OUT)   :: L                               !< \f$ L_N(\xi) \f$
REAL,INTENT(OUT)   :: q                               !< \f$ q_N(\xi) \f$
REAL,INTENT(OUT)   :: qder                            !< \f$ \partial/\partial\xi \; L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iLegendre
REAL               :: L_Nm1,L_Nm2                     ! L_{N_in-2},L_{N_in-1}
REAL               :: Lder,Lder_Nm1,Lder_Nm2          ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
L_Nm2=1.
L_Nm1=x
Lder_Nm2=0.
Lder_Nm1=1.
DO iLegendre=2,N_in
  L=(REAL(2*iLegendre-1)*x*L_Nm1 - REAL(iLegendre-1)*L_Nm2)/REAL(iLegendre)
  Lder=Lder_Nm2 + REAL(2*iLegendre-1)*L_Nm1
  L_Nm2=L_Nm1
  L_Nm1=L
  Lder_Nm2=Lder_Nm1
  Lder_Nm1=Lder
END DO ! iLegendre
q=REAL(2*N_in+1)/REAL(N_in+1)*(x*L -L_Nm2) !L_{N_in+1}-L_{N_in-1} !L_Nm2 is L_Nm1, L_Nm1 was overwritten!
qder= REAL(2*N_in+1)*L             !Lder_{N_in+1}-Lder_{N_in-1}
END SUBROUTINE qAndLEvaluation



!==================================================================================================================================
!> Starting with initial guess by Parter Relation, a Newton method is used to find the roots
!> of the Legendre Polynomial Lder_(N_in), which are the positions of Gauss-Lobatto points.
!> Uses qAndLEvaluation subroutine.
!> algorithm 25, Kopriva
!==================================================================================================================================
SUBROUTINE LegGaussLobNodesAndWeights(N_in,xGP,wGP)
! MODULES
USE MOD_Preproc
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in             !< polynomial degree (N_in+1) Gausspoints
REAL,INTENT(OUT)          :: xGP(0:N_in)      !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(OUT),OPTIONAL :: wGP(0:N_in)      !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: nIter = 10       ! max. number of newton iterations
REAL                      :: Tol   = 1.E-15   ! tolerance for Newton iteration : TODO: use variable tolerance here!
INTEGER                   :: iGP,iter
REAL                      :: q,qder,L         ! \f$ q=L_{N_in+1}-L_{N_in-1} \f$ ,qder is derivative, \f$ L=L_{N_in} \f$
REAL                      :: dx               ! Newton step
REAL                      :: cont1,cont2      !temporary variable for evaluation of parter nodes positions
!==================================================================================================================================
xGP(0)=-1.
xGP(N_in)= 1.
IF(PRESENT(wGP))THEN
  wGP(0)= 2./REAL(N_in*(N_in+1))
  wGP(N_in)=wGP(0)
END IF
IF(N_in.GT.1)THEN
  cont1=PP_Pi/REAL(N_in) ! pi/N_in
  cont2=3./(REAL(8*N_in)*PP_Pi) ! 3/(8*N_in*pi)
  DO iGP=1,(N_in+1)/2-1 !since points are symmetric, only left side is computed
    xGP(iGP)=-cos(cont1*(REAL(iGP)+0.25)-cont2/(REAL(iGP)+0.25)) !initial guess
    ! Newton iteration
    DO iter=0,nIter
      CALL qAndLEvaluation(N_in,xGP(iGP),q,qder,L)
      dx=-q/qder
      xGP(iGP)=xGP(iGP)+dx
      IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
    END DO ! iter
    IF(iter.GT.nIter) THEN
      SWRITE(*,*) 'maximum iteration steps >10 in Newton iteration for LGL point:'
      xGP(iGP)=-cos(cont1*(REAL(iGP)+0.25)-cont2/(REAL(iGP)+0.25)) !initial guess
      ! Newton iteration
      DO iter=0,nIter
        SWRITE(*,*)'iter,x^i',iter,xGP(iGP)     !DEBUG
        CALL qAndLEvaluation(N_in,xGP(iGP),q,qder,L)
        dx=-q/qder
        xGP(iGP)=xGP(iGP)+dx
        IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
      END DO ! iter
      CALL abort(__STAMP__,&
                 'ERROR: Legendre Gauss Lobatto nodes could not be computed up to desired precision. Code stopped!')
    END IF ! (iter.GT.nIter)
    CALL qAndLEvaluation(N_in,xGP(iGP),q,qder,L)
    xGP(N_in-iGP)=-xGP(iGP)
    IF(PRESENT(wGP))THEN
      wGP(iGP)=wGP(0)/(L*L)
      wGP(N_in-iGP)=wGP(iGP)
    END IF
  END DO ! iGP
END IF !(N_in.GT.1)
IF(mod(N_in,2) .EQ. 0) THEN
  xGP(N_in/2)=0.
  CALL qAndLEvaluation(N_in,xGP(N_in/2),q,qder,L)
  IF(PRESENT(wGP))wGP(N_in/2)=wGP(0)/(L*L)
END IF ! (mod(N_in,2) .EQ. 0)
END SUBROUTINE LegGaussLobNodesAndWeights



!==================================================================================================================================
!> Computes barycentric (interpolation) weights for interpolation polynomial given by set of nodes. (Algorithm 30, Kopriva book)
!==================================================================================================================================
PPURE SUBROUTINE BarycentricWeights(N_in,xGP,wBary)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in               !< polynomial degree
REAL,INTENT(IN)    :: xGP(0:N_in)        !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(OUT)   :: wBary(0:N_in)      !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,jGP
!==================================================================================================================================
wBary(:)=1.
DO iGP=1,N_in
  DO jGP=0,iGP-1
    wBary(jGP)=wBary(jGP)*(xGP(jGP)-xGP(iGP))
    wBary(iGP)=wBary(iGP)*(xGP(iGP)-xGP(jGP))
  END DO ! jGP
END DO ! iGP
wBary(:)=1./wBary(:)
END SUBROUTINE BarycentricWeights



!==================================================================================================================================
!> Computes polynomial differentiation matrix for interpolation polynomial given by set of nodes. (Algorithm 37, Kopriva book)
!==================================================================================================================================
PPURE SUBROUTINE PolynomialDerivativeMatrix(N_in,xGP,D)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in              !< polynomial degree
REAL,INTENT(IN)    :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(OUT)   :: D(0:N_in,0:N_in)  !< differentiation Matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,iLagrange
REAL               :: wBary(0:N_in)
!==================================================================================================================================
CALL BarycentricWeights(N_in,xGP,wBary)
D(:,:)=0.
DO iLagrange=0,N_in
  DO iGP=0,N_in
    IF(iLagrange.NE.iGP)THEN
      D(iGP,iLagrange)=wBary(iLagrange)/(wBary(iGP)*(xGP(iGP)-xGP(iLagrange)))
      D(iGP,iGP)=D(iGP,iGP)-D(iGP,iLagrange)
    END IF ! (iLagrange.NE.iGP)
  END DO ! iGP
END DO ! iLagrange
END SUBROUTINE PolynomialDerivativeMatrix


!==================================================================================================================================
!> Determines if two real numbers are equal up to a specified tolerance (=PP_RealTolerance, normaly set to machine precision)
!> Takes into account that x,y are located in-between [-1;1] for additional accuracy
!> Based on Algorithm 139, Kopriva
!==================================================================================================================================
ELEMENTAL FUNCTION ALMOSTEQUAL(x,y)
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: x                !< (IN)  first scalar to be compared
REAL,INTENT(IN) :: y                !< (IN)  second scalar to be compared
LOGICAL         :: AlmostEqual      !< (OUT) TRUE if |x-y| < 2*PP_RealTolerance
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
AlmostEqual=.FALSE.
IF((x.EQ.0.).OR.(y.EQ.0.)) THEN
  IF(ABS(x-y).LE.2.*PP_RealTolerance) AlmostEqual=.TRUE.
ELSE ! x, y not zero
  IF((ABS(x-y).LE.PP_RealTolerance*ABS(x)).AND.((ABS(x-y).LE.PP_RealTolerance*ABS(y)))) AlmostEqual=.TRUE.
END IF ! x,y zero
END FUNCTION ALMOSTEQUAL


!============================================================================================================================
!> Computes all Lagrange functions evaluated at position x in [-1,1]
!> For details see paper Barycentric Lagrange Interpolation by Berrut and Trefethen (SIAM 2004)
!> Uses function ALMOSTEQUAL
!> Algorithm 34, Kopriva book
!============================================================================================================================
PPURE SUBROUTINE LagrangeInterpolationPolys(x,N_in,xGP,wBary,L)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)   :: x                !< Coordinate
INTEGER,INTENT(IN) :: N_in             !< polynomial degree
REAL,INTENT(IN)    :: xGP(0:N_in)      !< Gauss point positions for the reference interval [-1,1]
REAL,INTENT(IN)    :: wBary(0:N_in)    !< Barycentric weights
REAL,INTENT(OUT)   :: L(0:N_in)        !< Lagrange basis functions evaluated at x
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP
LOGICAL            :: xEqualGP         ! is x equal to a Gauss Point
REAL               :: DummySum
!============================================================================================================================
xEqualGP=.FALSE.
DO iGP=0,N_in
  L(iGP)=0.
  IF(ALMOSTEQUAL(x,xGP(iGP))) THEN
    L(iGP)=1.
    xEqualGP=.TRUE.
  END IF
END DO

! if x is equal to a Gauss point, L=(0,....,1,....0)
IF(xEqualGP) RETURN
DummySum=0.
DO iGP=0, N_in
  L(iGP)=wBary(iGP)/(x-xGP(iGP))
  DummySum=DummySum+L(iGP)
END DO

DO iGP=0,N_in
  L(iGP)=L(iGP)/DummySum
END DO
END SUBROUTINE LagrangeInterpolationPolys


END MODULE MOD_Basis
