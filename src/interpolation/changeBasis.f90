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
INTERFACE ChangeBasis3D_XYZ
  MODULE PROCEDURE ChangeBasis3D_XYZ
END INTERFACE

INTERFACE ChangeBasis3D
  MODULE PROCEDURE ChangeBasis3D
  MODULE PROCEDURE ChangeBasis3D_singleVar
END INTERFACE

INTERFACE ChangeBasis2D_XYZ
  MODULE PROCEDURE ChangeBasis2D_XYZ
END INTERFACE

INTERFACE ChangeBasis2D
  MODULE PROCEDURE ChangeBasis2D
  MODULE PROCEDURE ChangeBasis2D_singleVar
END INTERFACE

INTERFACE ChangeBasis1D
  MODULE PROCEDURE ChangeBasis1D
  MODULE PROCEDURE ChangeBasis1D_singleVar
END INTERFACE

PUBLIC :: ChangeBasis3D_XYZ
PUBLIC :: ChangeBasis3D
PUBLIC :: ChangeBasis2D_XYZ
PUBLIC :: ChangeBasis2D
PUBLIC :: ChangeBasis1D
!==================================================================================================================================
CONTAINS

#define _ADD_DIM
#include "changeBasis.t90"
#undef _ADD_DIM
END MODULE MOD_ChangeBasis

!==================================================================================================================================
!> Changes surface or volume data present as polynomial with Lagrange Basis of degree NIn to
!> a representation of polynomial of a Lagrange Basis NOut, using two
!> arbitrary point distributions xi_In(0:NIn) and xi_Out(0:NOut) and a series of 1D operations.
!==================================================================================================================================
MODULE MOD_ChangeBasisByDim
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
#if PP_dim == 3
INTERFACE ChangeBasisVolume
  MODULE PROCEDURE ChangeBasis3D
  MODULE PROCEDURE ChangeBasis3D_singleVar
END INTERFACE
#endif

#if PP_dim == 3
INTERFACE ChangeBasisSurf
#else
INTERFACE ChangeBasisVolume
#endif
  MODULE PROCEDURE ChangeBasis2D
  MODULE PROCEDURE ChangeBasis2D_singleVar
END INTERFACE

#if PP_dim == 2
INTERFACE ChangeBasisSurf
  MODULE PROCEDURE ChangeBasis1D
  MODULE PROCEDURE ChangeBasis1D_singleVar
END INTERFACE
#endif

PUBLIC :: ChangeBasisVolume
PUBLIC :: ChangeBasisSurf
!==================================================================================================================================
CONTAINS

#if PP_dim == 2
#  define _ADD_DIM ,1
#else
#  define _ADD_DIM
#endif
#include "changeBasis.t90"
#undef _ADD_DIM
END MODULE MOD_ChangeBasisByDim
