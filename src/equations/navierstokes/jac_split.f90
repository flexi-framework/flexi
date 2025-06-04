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

!===================================================================================================================================
!> This module contains the routine used to calculate the point-wise flux jacobian for the split flux formulations. They are either
!> implemented analytically, or a finite difference approach is used to get an approximate jacobian.
!===================================================================================================================================
MODULE MOD_Jac_Split
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
#ifdef SPLIT_DG
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: Jac_Split
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Calculate jacobian of split flux
!===================================================================================================================================
PPURE SUBROUTINE Jac_Split(U,UPrim,URef,UPrimRef,Metric,MetricRef,dfdu)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_DG_Vars, ONLY: SplitDG
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: Metric        !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MetricRef     !< metric terms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT) :: dfdu      !< dof local jacobian of volume split flux
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(SplitDG)
CASE(0)
  ! SD flux, retains normal DG formulation for debugging
  CALL Jac_Split_SD(U,UPrim,URef,UPrimRef,Metric,MetricRef,dfdu)
CASE(4)
  ! KEP flux after Pirozzoli
  CALL Jac_Split_PI(U,UPrim,URef,UPrimRef,Metric,MetricRef,dfdu)
CASE DEFAULT
  ! If no analytical derivative is implemented, use a FD approach
  CALL Jac_Split_FD(U,UPrim,URef,UPrimRef,Metric,MetricRef,dfdu)
END SELECT

END SUBROUTINE Jac_Split


!===================================================================================================================================
!> Calculate analytical jacobian of split flux (SD flux), which corresponds to the normal flux Jacobi
!===================================================================================================================================
PPURE SUBROUTINE Jac_Split_SD(U,UPrim,URef,UPrimRef,Metric,MetricRef,dfdu)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Jacobian,   ONLY: EvalAdvFluxJacobianPoint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: Metric        !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MetricRef     !< metric terms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT) :: dfdu      !< dof local jacobian of volume split flux
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar)         :: fJac,gJac,hJac
!===================================================================================================================================
CALL EvalAdvFluxJacobianPoint(U,UPrim,fJac,gJac,hJac)

dfdu = 0.5*(MetricRef(1)+Metric(1))*fJac(:,:) + &
#if PP_dim == 3
       0.5*(MetricRef(3)+Metric(3))*hJac(:,:) + &
#endif
       0.5*(MetricRef(2)+Metric(2))*gJac(:,:)

END SUBROUTINE Jac_Split_SD


!===================================================================================================================================
!> Calculate analytical jacobian of split flux (PI flux)
!===================================================================================================================================
PPURE SUBROUTINE Jac_Split_PI(U,UPrim,URef,UPrimRef,Metric,MetricRef,dfdu)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_EOS_Vars  ,ONLY: KappaM1
USE MOD_SplitFlux ,ONLY: SplitDGVolume_pointer
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: Metric        !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MetricRef     !< metric terms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT) :: dfdu      !< dof local jacobian of volume split flux
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: dpdU(5),df1dU(5),dg1dU(5),absu2,f1,g1,Hmean
#if PP_dim == 3
REAL                                    :: dh1dU(5),h1
REAL,DIMENSION(PP_nVar,PP_nVar)         :: hJac
#endif
REAL,DIMENSION(PP_nVar,PP_nVar)         :: fJac,gJac
!===================================================================================================================================
#if PP_dim==3
absu2=UPrim(2)**2+UPrim(3)**2+UPrim(4)**2
#else
absu2=UPrim(2)**2+UPrim(3)**2
#endif
Hmean     = 0.5*((URef(5)+UPrimRef(5))/URef(1)+(U(5)+UPrim(5))/U(1))
dpdU(1)   =  KappaM1*0.5*absu2
dpdU(2:4) = -KappaM1*UPrim(2:4)
dpdU(5)   =  KappaM1

f1 = 0.5 *(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))    ! {rho}*{u}

df1dU(1)   = 0.5*(UPrimRef(2)-URef(1)*UPrim(2)/U(1))
df1dU(2)   = 0.5*(1.+URef(1)/U(1))
df1dU(3:5) = 0.

fJac(1,1:5)= (/ df1dU(1), df1dU(2), df1dU(3),  df1dU(4), df1dU(5)/)

fJac(2,1)= df1dU(1)*0.5*(UPrimRef(2)+UPrim(2)) - 0.5*f1*UPrim(2)/U(1)+dpdU(1)
fJac(2,2)= df1dU(2)*0.5*(UPrimRef(2)+UPrim(2)) + 0.5*f1/U(1)         +dpdU(2)
fJac(2,3)= df1dU(3)*0.5*(UPrimRef(2)+UPrim(2))                       +dpdU(3)
#if PP_dim==3
fJac(2,4)= df1dU(4)*0.5*(UPrimRef(2)+UPrim(2))                       +dpdU(4)
#endif
fJac(2,5)= df1dU(5)*0.5*(UPrimRef(2)+UPrim(2))                       +dpdU(5)

fJac(3,1)= df1dU(1)*0.5*(UPrimRef(3)+UPrim(3)) - 0.5*f1*UPrim(3)/U(1)
fJac(3,2)= df1dU(2)*0.5*(UPrimRef(3)+UPrim(3))
fJac(3,3)= df1dU(3)*0.5*(UPrimRef(3)+UPrim(3)) + 0.5*f1/U(1)
#if PP_dim==3
fJac(3,4)= df1dU(4)*0.5*(UPrimRef(3)+UPrim(3))
#endif
fJac(3,5)= df1dU(5)*0.5*(UPrimRef(3)+UPrim(3))

#if PP_dim==3
fJac(4,1)= df1dU(1)*0.5*(UPrimRef(4)+UPrim(4)) - 0.5*f1*UPrim(4)/U(1)
fJac(4,2)= df1dU(2)*0.5*(UPrimRef(4)+UPrim(4))
fJac(4,3)= df1dU(3)*0.5*(UPrimRef(4)+UPrim(4))
fJac(4,4)= df1dU(4)*0.5*(UPrimRef(4)+UPrim(4)) + 0.5*f1/U(1)
fJac(4,5)= df1dU(5)*0.5*(UPrimRef(4)+UPrim(4))
#endif

fJac(5,1)= df1dU(1)*Hmean + 0.5*f1*((-U(5)-UPrim(5))/(U(1)**2)+dpdU(1)/U(1))
fJac(5,2)= df1dU(2)*Hmean + 0.5*f1*                            dpdU(2)/U(1)
fJac(5,3)= df1dU(3)*Hmean + 0.5*f1*                            dpdU(3)/U(1)
#if PP_dim==3
fJac(5,4)= df1dU(4)*Hmean + 0.5*f1*                            dpdU(4)/U(1)
#endif
fJac(5,5)= df1dU(5)*Hmean + 0.5*f1*(                   1./U(1)+dpdU(5)/U(1))

#if PP_dim==2
fJac(4,:)= 0.
fJac(:,4)= 0.
#endif

g1 = 0.5 *(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))    ! {rho}*{v}

dg1dU(1)   = 0.5*(UPrimRef(3)-URef(1)*UPrim(3)/U(1))
dg1dU(2)   = 0.
dg1dU(3)   = 0.5*(1.+URef(1)/U(1))
dg1dU(4:5) = 0.

gJac(1,1:5)= (/ dg1dU(1), dg1dU(2), dg1dU(3),  dg1dU(4), dg1dU(5)/)

gJac(2,1)= dg1dU(1)*0.5*(UPrimRef(2)+UPrim(2)) - 0.5*g1*UPrim(2)/U(1)
gJac(2,2)= dg1dU(2)*0.5*(UPrimRef(2)+UPrim(2)) + 0.5*g1/U(1)
gJac(2,3)= dg1dU(3)*0.5*(UPrimRef(2)+UPrim(2))
#if PP_dim==3
gJac(2,4)= dg1dU(4)*0.5*(UPrimRef(2)+UPrim(2))
#endif
gJac(2,5)= dg1dU(5)*0.5*(UPrimRef(2)+UPrim(2))

gJac(3,1)= dg1dU(1)*0.5*(UPrimRef(3)+UPrim(3)) - 0.5*g1*UPrim(3)/U(1)+dpdU(1)
gJac(3,2)= dg1dU(2)*0.5*(UPrimRef(3)+UPrim(3))                       +dpdU(2)
gJac(3,3)= dg1dU(3)*0.5*(UPrimRef(3)+UPrim(3)) + 0.5*g1/U(1)         +dpdU(3)
#if PP_dim==3
gJac(3,4)= dg1dU(4)*0.5*(UPrimRef(3)+UPrim(3))                       +dpdU(4)
#endif
gJac(3,5)= dg1dU(5)*0.5*(UPrimRef(3)+UPrim(3))                       +dpdU(5)

#if PP_dim==3
gJac(4,1)= dg1dU(1)*0.5*(UPrimRef(4)+UPrim(4)) - 0.5*g1*UPrim(4)/U(1)
gJac(4,2)= dg1dU(2)*0.5*(UPrimRef(4)+UPrim(4))
gJac(4,3)= dg1dU(3)*0.5*(UPrimRef(4)+UPrim(4))
gJac(4,4)= dg1dU(4)*0.5*(UPrimRef(4)+UPrim(4)) + 0.5*g1/U(1)
gJac(4,5)= dg1dU(5)*0.5*(UPrimRef(4)+UPrim(4))
#endif

gJac(5,1)= dg1dU(1)*Hmean + 0.5*g1*((-U(5)-UPrim(5))/(U(1)**2)+dpdU(1)/U(1))
gJac(5,2)= dg1dU(2)*Hmean + 0.5*g1*                            dpdU(2)/U(1)
gJac(5,3)= dg1dU(3)*Hmean + 0.5*g1*                            dpdU(3)/U(1)
#if PP_dim==3
gJac(5,4)= dg1dU(4)*Hmean + 0.5*g1*                            dpdU(4)/U(1)
#endif
gJac(5,5)= dg1dU(5)*Hmean + 0.5*g1*(                   1./U(1)+dpdU(5)/U(1))

#if PP_dim==2
gJac(4,:)= 0.
gJac(:,4)= 0.
#endif

#if PP_dim==3
h1 = 0.5 *(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))    ! {rho}*{w}

dh1dU(1)   = 0.5*(UPrimRef(4)-URef(1)*UPrim(4)/U(1))
dh1dU(2:3) = 0.
dh1dU(4)   = 0.5*(1.+URef(1)/U(1))
dh1dU(5)   = 0.

hJac(1,1:5)= (/ dh1dU(1), dh1dU(2), dh1dU(3),  dh1dU(4), dh1dU(5)/)

hJac(2,1)= dh1dU(1)*0.5*(UPrimRef(2)+UPrim(2)) - 0.5*h1*UPrim(2)/U(1)
hJac(2,2)= dh1dU(2)*0.5*(UPrimRef(2)+UPrim(2)) + 0.5*h1/U(1)
hJac(2,3)= dh1dU(3)*0.5*(UPrimRef(2)+UPrim(2))
hJac(2,4)= dh1dU(4)*0.5*(UPrimRef(2)+UPrim(2))
hJac(2,5)= dh1dU(5)*0.5*(UPrimRef(2)+UPrim(2))

hJac(3,1)= dh1dU(1)*0.5*(UPrimRef(3)+UPrim(3)) - 0.5*h1*UPrim(3)/U(1)
hJac(3,2)= dh1dU(2)*0.5*(UPrimRef(3)+UPrim(3))
hJac(3,3)= dh1dU(3)*0.5*(UPrimRef(3)+UPrim(3)) + 0.5*h1/U(1)
hJac(3,4)= dh1dU(4)*0.5*(UPrimRef(3)+UPrim(3))
hJac(3,5)= dh1dU(5)*0.5*(UPrimRef(3)+UPrim(3))

hJac(4,1)= dh1dU(1)*0.5*(UPrimRef(4)+UPrim(4)) - 0.5*h1*UPrim(4)/U(1)+dpdU(1)
hJac(4,2)= dh1dU(2)*0.5*(UPrimRef(4)+UPrim(4))                       +dpdU(2)
hJac(4,3)= dh1dU(3)*0.5*(UPrimRef(4)+UPrim(4))                       +dpdU(3)
hJac(4,4)= dh1dU(4)*0.5*(UPrimRef(4)+UPrim(4)) + 0.5*h1/U(1)         +dpdU(4)
hJac(4,5)= dh1dU(5)*0.5*(UPrimRef(4)+UPrim(4))                       +dpdU(5)

hJac(5,1)= dh1dU(1)*Hmean + 0.5*h1*((-U(5)-UPrim(5))/(U(1)**2)+dpdU(1)/U(1))
hJac(5,2)= dh1dU(2)*Hmean + 0.5*h1*                            dpdU(2)/U(1)
hJac(5,3)= dh1dU(3)*Hmean + 0.5*h1*                            dpdU(3)/U(1)
hJac(5,4)= dh1dU(4)*Hmean + 0.5*h1*                            dpdU(4)/U(1)
hJac(5,5)= dh1dU(5)*Hmean + 0.5*h1*(                   1./U(1)+dpdU(5)/U(1))
#endif

dfdu = 0.5*(MetricRef(1)+Metric(1))*fJac(:,:) + &
#if PP_dim == 3
       0.5*(MetricRef(3)+Metric(3))*hJac(:,:) + &
#endif
       0.5*(MetricRef(2)+Metric(2))*gJac(:,:)

END SUBROUTINE Jac_Split_PI


!===================================================================================================================================
!> Calculate jacobian of an arbitrary split flux using the finite difference approach (EXPENSIVE BUT FOR ALL FLUXES)
!===================================================================================================================================
PPURE SUBROUTINE Jac_Split_FD(U,UPrim,URef,UPrimRef,Metric,MetricRef,dfdu)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_SplitFlux     ,ONLY: SplitDGVolume_pointer
USE MOD_Implicit_Vars ,ONLY: reps0_O1,sreps0_O1
USE MOD_EOS           ,ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: Metric        !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MetricRef     !< metric terms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT) :: dfdu      !< dof local jacobian of volume split flux
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: FluxTilde(PP_nVar),UTilde(PP_nVar),UPrimTilde(PP_nVarPrim),Flux(PP_nVar)
INTEGER                                 :: iVar, jVar
!===================================================================================================================================

CALL SplitDGVolume_pointer(Uref,UPrimRef,U,UPrim,MetricRef,Metric,Flux)
UTilde = U
! Pertubate each variable and use FD (first order) to approximate the derivative
DO jVar=1,PP_nVar
  UTilde(jVar) = UTilde(jVar) + reps0_O1
  CALL ConsToPrim(UPrimTilde,UTilde)
  CALL SplitDGVolume_pointer(Uref,UPrimRef,UTilde,UPrimTilde,MetricRef,Metric,FluxTilde)
  DO iVar=1,PP_nVar
    dfdu(iVar,jVar) = (FluxTilde(iVar)-Flux(iVar))*sreps0_O1
  END DO ! iVar
  UTilde = U
END DO
END SUBROUTINE Jac_Split_FD

#endif /*SPLIT_DG*/
END MODULE MOD_Jac_Split
