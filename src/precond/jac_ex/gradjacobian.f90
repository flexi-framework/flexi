#include "flexi.h"
#if EQNSYSNR==2
#include "eos.h"
#endif

MODULE MOD_GradJacobian
#if PARABOLIC   
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!INTERFACE EvalFluxGradJacobian
!  MODULE PROCEDURE EvalFluxGradJacobian
!END INTERFACE


PUBLIC::EvalFluxGradJacobian
!===================================================================================================================================

CONTAINS


#if EQNSYSNR==1
SUBROUTINE EvalFluxGradJacobian(Nloc,U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz)
!===================================================================================================================================
! Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                          :: Nloc
REAL,DIMENSION(PP_nVar,Nloc),INTENT(IN)     :: U
REAL,DIMENSION(PP_nVarPrim,Nloc),INTENT(IN) :: UPrim
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
                                                        ! Gradient of the diffusive Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,PP_nVar,Nloc),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
fJacQx = -DiffC
fJacQy = 0.
fJacQz = 0.

gJacQx = 0.
gJacQy = -DiffC
gJacQz = 0.

#if PP_dim==3
hJacQx = 0.
hJacQy = 0.
hJacQz = -DiffC
#endif
END SUBROUTINE EvalFluxGradJacobian
#endif /*linearscalaradvection*/

#if EQNSYSNR==2
SUBROUTINE EvalFluxGradJacobian(Nloc,U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz)
!===================================================================================================================================
! Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:s43,s23
USE MOD_Viscosity
USE MOD_EOS_Vars,     ONLY:cp,Pr,KappasPr
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                   :: Nloc
REAL,DIMENSION(PP_nVar,Nloc),INTENT(IN)              :: U
REAL,DIMENSION(PP_nVarPrim,Nloc),INTENT(IN)          :: UPrim
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
                                                        ! Gradient of the diffusive Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,PP_nVar,Nloc),INTENT(OUT)     :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: srho,e                                  ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2                                   ! auxiliary variables
REAL                :: v1v1,v2v2,v1v2
REAL                :: KappasPrEmVabs
REAL                :: s43mKappasPr                            ! (4/3-kappa/Pr)
REAL                :: mKappasPr                               ! (1-kappa/Pr)
REAL                :: s13
REAL                :: muS                                     ! viscosity
REAL                :: lambda                                  
INTEGER             :: i
#if PP_dim==3
REAL                :: v3                              ! auxiliary variables
REAL                :: v3v3,v2v3,v1v3
#endif
!===================================================================================================================================
s13          = REAL(1./3.)
s43mKappasPr = s43-KappasPr
mKappasPr    = 1.-KappasPr

DO i=1,Nloc
  srho = 1. / U(1,i) ! 1/rho
  v1   = UPrim(2,i)
  v2   = UPrim(3,i)
  v1v1 = v1*v1
  v2v2 = v2*v2
  v1v2 = v1*v2

  e=U(5,i)*srho                              ! e...specific energy
  ! Viscous part
  ! ideal gas law
  muS=VISCOSITY_PRIM(UPrim)
  lambda=THERMAL_CONDUCTIVITY_H(muS)
  srho=srho*muS !add viscosity  ! = mu / rho
#if PP_dim==3
  v3   = UPrim(4,i)
  v3v3 = v3*v3
  v1v3 = v1*v3
  v2v3 = v2*v3

  KappasPrEmVabs = KappasPr*(e-v1v1-v2v2-v3v3)
  ! compute the derivative of F with respect to Qx
  fJacQx(1,1:5,i) = 0.
  fJacQx(2,1:5,i) = (/ srho*s43*v1,-srho*s43,      0.,      0.,0. /) 
  fJacQx(3,1:5,i) = (/     srho*v2,       0.,   -srho,      0.,0. /) 
  fJacQx(4,1:5,i) = (/     srho*v3,       0.,      0.,   -srho,0. /) 
  ! for each element
  fJacQx(5,1  ,i) =  srho*(s43*v1v1+v2v2+v3v3 + KappasPrEmVabs )
  fJacQx(5,2  ,i) = -srho*s43mKappasPr*v1
  fJacQx(5,3  ,i) = -srho*mKappasPr*v2
  fJacQx(5,4  ,i) = -srho*mKappasPr*v3
  fJacQx(5,5  ,i) = -srho*KappasPr
  ! compute the derivative of G with respect to Qx
  gJacQx(1,1:5,i) = 0.
  gJacQx(2,1:5,i) = (/     srho*v2,       0.,   -srho,      0.,0. /) 
  gJacQx(3,1:5,i) = (/-srho*s23*v1, srho*s23,      0.,      0.,0. /) 
  gJacQx(4,1:5,i) = 0.
  ! for each element
  gJacQx(5,1  ,i) =  srho*s13*v1v2
  gJacQx(5,2  ,i) =  srho*s23*v2
  gJacQx(5,3  ,i) = -srho*v1
  gJacQx(5,4  ,i) =     0.
  gJacQx(5,5  ,i) =  0.
  ! compute the derivative of H with respect to Qx
  hJacQx(1,1:5,i) = 0.
  hJacQx(2,1:5,i) = (/     srho*v3,       0.,      0.,   -srho,0. /) 
  hJacQx(3,1:5,i) = 0.
  hJacQx(4,1:5,i) = (/-srho*s23*v1, srho*s23,      0.,      0.,0. /)
  ! for each element
  hJacQx(5,1  ,i) =  srho*s13*v1v3
  hJacQx(5,2  ,i) =  srho*s23*v3
  hJacQx(5,3  ,i) =  0.
  hJacQx(5,4  ,i) =  -srho*v1
  hJacQx(5,5  ,i) =  0.
  ! compute the derivative of F with respect to Qy
  fJacQy(1,1:5,i) = 0.
  fJacQy(2,1:5,i) = (/-srho*s23*v2,       0.,srho*s23,      0.,0. /) 
  fJacQy(3,1:5,i) = (/     srho*v1,    -srho,      0.,      0.,0. /) 
  fJacQy(4,1:5,i) = 0.
  ! for each element
  fJacQy(5,1  ,i) =  srho*s13*v1v2
  fJacQy(5,2  ,i) = -srho*v2
  fJacQy(5,3  ,i) =  srho*s23*v1
  fJacQy(5,4  ,i) =  0.
  fJacQy(5,5  ,i) =  0.
  ! compute the derivative of G with respect to Qy
  gJacQy(1,1:5,i) = 0.
  gJacQy(2,1:5,i) = (/     srho*v1,    -srho,        0.,      0.,0. /) 
  gJacQy(3,1:5,i) = (/ srho*s43*v2,        0.,-srho*s43,      0.,0. /) 
  gJacQy(4,1:5,i) = (/     srho*v3,        0.,       0.,   -srho,0. /)
  ! for each element
  gJacQy(5,1  ,i) =  srho*(KappasPrEmVabs+v1v1+s43*v2v2+v3v3)
  gJacQy(5,2  ,i) = -srho*v1*mKappasPr
  gJacQy(5,3  ,i) = -srho*v2*s43mKappasPr
  gJacQy(5,4  ,i) = -srho*mKappasPr*v3
  gJacQy(5,5  ,i) = -srho*KappasPr
  ! compute the derivative of H with respect to Qy
  hJacQy(1,1:5,i) = 0.
  hJacQy(2,1:5,i) = 0.
  hJacQy(3,1:5,i) = (/     srho*v3,       0.,      0.,   -srho,0. /)
  hJacQy(4,1:5,i) = (/-srho*s23*v2,       0.,srho*s23,      0.,0. /)
  ! for each element
  hJacQy(5,1  ,i) =  srho*s13*v2v3
  hJacQy(5,2  ,i) =  0.
  hJacQy(5,3  ,i) =  srho*s23*v3
  hJacQy(5,4  ,i) = -srho*v2
  hJacQy(5,5  ,i) =  0.
  ! compute the derivative of F with respect to Qz
  fJacQz(1,1:5,i) = 0.
  fJacQz(2,1:5,i) = (/-srho*s23*v3,       0.,       0.,srho*s23,0. /) 
  fJacQz(3,1:5,i) = 0.
  fJacQz(4,1:5,i) = (/     srho*v1,    -srho,       0.,      0.,0.  /)
  ! for each element
  fJacQz(5,1  ,i) =  srho*s13*v1v3
  fJacQz(5,2  ,i) = -srho*v3
  fJacQz(5,3  ,i) =  0.
  fJacQz(5,4  ,i) =  srho*s23*v1
  fJacQz(5,5  ,i) =  0.
  ! compute the derivative of G with respect to Qz
  gJacQz(1,1:5,i) = 0.
  gJacQz(2,1:5,i) = 0.
  gJacQz(3,1:5,i) = (/-srho*s23*v3,        0.,       0.,srho*s23,0. /) 
  gJacQz(4,1:5,i) = (/     srho*v2,        0.,    -srho,      0.,0. /)
  ! for each element
  gJacQz(5,1  ,i) =  srho*s13*v2v3
  gJacQz(5,2  ,i) =  0.
  gJacQz(5,3  ,i) = -srho*v3
  gJacQz(5,4  ,i) =  srho*s23*v2
  gJacQz(5,5  ,i) =  0.
  ! compute the derivative of H with respect to Qz
  hJacQz(1,1:5,i) = 0.
  hJacQz(2,1:5,i) = (/     srho*v1,    -srho,       0.,        0.,0. /)
  hJacQz(3,1:5,i) = (/     srho*v2,       0.,    -srho,        0.,0. /)
  hJacQz(4,1:5,i) = (/ srho*s43*v3,       0.,       0., -srho*s43,0. /)
  ! for each element
  hJacQz(5,1  ,i) =  srho*(v1v1+v2v2+s43*v3v3+KappasPrEmVabs)
  hJacQz(5,2  ,i) = -srho*mKappasPr*v1 
  hJacQz(5,3  ,i) = -srho*mKappasPr*v2 
  hJacQz(5,4  ,i) = -srho*s43mKappasPr*v3
  hJacQz(5,5  ,i) = -srho*KappasPr
#else

  KappasPrEmVabs = KappasPr*(e-v1v1-v2v2)
  ! compute the derivative of F with respect to Qx
  fJacQx(1,1:5,i) = 0.
  fJacQx(2,1:5,i) = (/ srho*s43*v1,-srho*s43,      0.,      0.,0. /) 
  fJacQx(3,1:5,i) = (/     srho*v2,       0.,   -srho,      0.,0. /) 
  fJacQx(4,1:5,i) = 0.
  !fJacQx(4,1:5,i) = (/          0.,       0.,      0.,   -srho,0. /) 
  ! for each element
  fJacQx(5,1  ,i) =  srho*(s43*v1v1+v2v2 + KappasPrEmVabs )
  fJacQx(5,2  ,i) = -srho*s43mKappasPr*v1
  fJacQx(5,3  ,i) = -srho*mKappasPr*v2
  fJacQx(5,4  ,i) = 0. 
  fJacQx(5,5  ,i) = -srho*KappasPr
  ! compute the derivative of G with respect to Qx
  gJacQx(1,1:5,i) = 0.
  gJacQx(2,1:5,i) = (/     srho*v2,       0.,   -srho,      0.,0. /) 
  gJacQx(3,1:5,i) = (/-srho*s23*v1, srho*s23,      0.,      0.,0. /) 
  gJacQx(4,1:5,i) = 0.
  ! for each element
  gJacQx(5,1  ,i) =  srho*s13*v1v2
  gJacQx(5,2  ,i) =  srho*s23*v2
  gJacQx(5,3  ,i) = -srho*v1
  gJacQx(5,4  ,i) =  0.
  gJacQx(5,5  ,i) =  0.
  ! compute the derivative of H with respect to Qx
  hJacQx(:,:,i) = 0.
  ! compute the derivative of F with respect to Qy
  fJacQy(1,1:5,i) = 0.
  fJacQy(2,1:5,i) = (/-srho*s23*v2,       0.,srho*s23,      0.,0. /) 
  fJacQy(3,1:5,i) = (/     srho*v1,    -srho,      0.,      0.,0. /) 
  fJacQy(4,1:5,i) = 0.
  ! for each element
  fJacQy(5,1  ,i) =  srho*s13*v1v2
  fJacQy(5,2  ,i) = -srho*v2
  fJacQy(5,3  ,i) =  srho*s23*v1
  fJacQy(5,4  ,i) =  0.
  fJacQy(5,5  ,i) =  0.
  ! compute the derivative of G with respect to Qy
  gJacQy(1,1:5,i) = 0.
  gJacQy(2,1:5,i) = (/     srho*v1,    -srho,        0.,      0.,0. /) 
  gJacQy(3,1:5,i) = (/ srho*s43*v2,        0.,-srho*s43,      0.,0. /) 
  gJacQy(4,1:5,i) = 0.
  !gJacQy(4,1:5,i) = (/          0.,        0.,       0.,   -srho,0. /)
  ! for each element
  gJacQy(5,1  ,i) =  srho*(KappasPrEmVabs+v1v1+s43*v2v2)
  gJacQy(5,2  ,i) = -srho*v1*mKappasPr
  gJacQy(5,3  ,i) = -srho*v2*s43mKappasPr
  gJacQy(5,4  ,i) = 0. 
  gJacQy(5,5  ,i) = -srho*KappasPr
  ! compute the derivative of H with respect to Qy
  hJacQy(:,:,i) = 0.
  ! compute the derivative of F with respect to Qz
  fJacQz(:,:,i) = 0.
  ! compute the derivative of G with respect to Qz
  gJacQz(:,:,i) = 0.
  ! compute the derivative of H with respect to Qz
  hJacQz(:,:,i) = 0.
#endif
END DO ! i

END SUBROUTINE EvalFluxGradJacobian

#endif /*navierstokes*/

#endif /*PARABOLIC*/
END MODULE MOD_GradJacobian
