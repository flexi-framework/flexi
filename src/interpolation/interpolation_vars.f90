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
!> Contains global variables used by the DG modules.
!==================================================================================================================================
MODULE MOD_Interpolation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! reserved for Gauss Points with polynomial degree N, all allocated (0:N)
REAL,ALLOCATABLE   :: L_Plus(:)                   !< L for boundary flux computation at plus side  (1)
REAL,ALLOCATABLE   :: L_Minus(:)                  !< L for boundary flux computation at minus side (-1)
REAL,ALLOCATABLE   :: xGP(:)                      !< Gauss point coordinates
REAL,ALLOCATABLE   :: wGP(:)                      !< GP integration weights
REAL,ALLOCATABLE   :: wBary(:)                    !< barycentric weights
REAL,ALLOCATABLE   :: Vdm_Leg(:,:)                !< Legendre Vandermonde matrix
REAL,ALLOCATABLE   :: sVdm_Leg(:,:)               !< inverse of Legendre Vandermonde matrix

!==================================================================================================================================
!@{ Named nodetype parameters
!==================================================================================================================================
CHARACTER(LEN=255),PARAMETER :: NodeTypeG    = 'GAUSS'                    !< Gauss nodes (-1,1)
CHARACTER(LEN=255),PARAMETER :: NodeTypeGL   = 'GAUSS-LOBATTO'            !< Gauss-Lobatto nodes [-1,1]
CHARACTER(LEN=255),PARAMETER :: NodeTypeGR   = 'GAUSS-RADAU'              !< Gauss-Radau nodes [-1,1)
CHARACTER(LEN=255),PARAMETER :: NodeTypeCL   = 'CHEBYSHEV-GAUSS-LOBATTO'  !< Chabyshev-Gauss-Lobatto nodes [-1,1]
CHARACTER(LEN=255),PARAMETER :: NodeTypeVISU = 'VISU'                     !< equidistant nodes [-1,1]
CHARACTER(LEN=255),PARAMETER :: NodeTypeVISUInner  = 'VISU_INNER'         !< inner visualization nodes
CHARACTER(LEN=255),PARAMETER :: NodeTypeVISUFVEqui = 'VISU_FVEQUI'        !< equidistant visualization nodes for FV
!@}
#if (PP_NodeType==1)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS'
#elif (PP_NodeType==2)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS-LOBATTO'
#elif (PP_NodeType==3)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'CHEBYSHEV-GAUSS-LOBATTO'
#endif

LOGICAL           :: InterpolationInitIsDone = .FALSE.
END MODULE MOD_Interpolation_Vars
