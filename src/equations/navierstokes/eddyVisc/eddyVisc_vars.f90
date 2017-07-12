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

!===================================================================================================================================
!> Contains the parameters needed for eddy viscosity models
!===================================================================================================================================
MODULE MOD_EddyVisc_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

ABSTRACT INTERFACE
  SUBROUTINE EddyViscInt(iElem,i,j,k,muSGS)
  INTEGER,INTENT(IN)  :: iElem  !< index of current element
  !> indices of the c
  INTEGER,INTENT(IN)  :: i,j,k
  !> gradients of the directions
  REAL,INTENT(INOUT)  :: muSGS  !< local SGS viscosity
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE EddyVisc_surfInt(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,DeltaSS,SGS_Ind,muSGS,Face_xGP)
  !> gradients of the velocities w.r.t. all directions
  REAL,INTENT(IN)   :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32
  REAL,INTENT(IN)   :: rho      !< Density
  REAL,INTENT(IN)   :: DeltaSS  !< Filter width
  REAL,INTENT(IN)   :: SGS_Ind  !< Indicator for SGS model
  REAL,INTENT(IN)   :: Face_xGP !< Coordinate for van-Driest damping
  REAL,INTENT(OUT)  :: muSGS    !< local SGS viscosity
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE testfilterInt(U_in)
    USE MOD_PreProc
    USE MOD_Mesh_Vars,              ONLY:nElems
    REAL,INTENT(IN)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE FinalizeEddyViscosityInt()
  END SUBROUTINE
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                             :: eddyViscType           !< type of eddy viscosity
CHARACTER(LEN=255)                  :: WallDistFile
PROCEDURE(EddyViscInt)     ,POINTER :: eddyViscosity          !< pointer to routine for computing volume eddy viscosity
PROCEDURE(EddyVisc_surfInt),POINTER :: eddyViscosity_surf     !< pointer to routine for computing surface eddy viscosity
PROCEDURE(testfilterInt),POINTER       :: testfilter     !< pointer to routine for computing test filter
PROCEDURE(FinalizeEddyViscosityInt),POINTER       :: FinalizeEddyViscosity     !< pointer tofinalize routine 

!Smagosinsky Standard
REAL,ALLOCATABLE  :: DeltaS(:)         !< filter width, used by Smagorinsky modell
REAL,ALLOCATABLE  :: DeltaS_m(:,:)     !< filter width per direction, used by Vreman modell, apr. by cartesian cells
REAL,ALLOCATABLE  :: DeltaS_master(:)
REAL,ALLOCATABLE  :: DeltaS_slave(:)
REAL,ALLOCATABLE  :: muSGS(:,:,:,:,:)  !< Viscosity for the sub-grid
REAL,ALLOCATABLE  :: muSGSmax(:)       !< Viscosity for the sub-grid
REAL              :: CS                !< Smagorinsky constant, LES
REAL              :: PrSGS             !< Prandtl number for the sub-grid scales
! Dynamic Smagorinsky
REAL,ALLOCATABLE     :: FilterMat_Testfilter(:,:) 
REAL,ALLOCATABLE     :: SGS_Ind(:,:,:,:,:) 
REAL,ALLOCATABLE     :: SGS_Ind_master(:,:,:,:) 
REAL,ALLOCATABLE     :: SGS_Ind_slave(:,:,:,:) 
REAL,ALLOCATABLE     :: IntElem(:,:,:,:) 
LOGICAL,ALLOCATABLE  :: filter_ind(:,:) !< Do filter along i,j,k index?
LOGICAL,ALLOCATABLE  :: average_ind(:,:) !< Do average along i,j,k index?
INTEGER,ALLOCATABLE  :: average_type(:)  !< Type of average_ind for select case
!MATTEO:debug output
!REAL,ALLOCATABLE  :: S_en_out(:,:,:,:,:)  !< Debug output of |S|
!REAL,ALLOCATABLE  :: filtdir_out(:)  !< Debug output of filtering directions
!REAL,ALLOCATABLE  :: walldist_out(:)  !< Debug output of wall distance
!REAL,ALLOCATABLE  :: walldist_x(:)  !< Debug output of wall distance
!REAL,ALLOCATABLE  :: walldist_y(:)  !< Debug output of wall distance
!REAL,ALLOCATABLE  :: walldist_z(:)  !< Debug output of wall distance


LOGICAL           :: VanDriest=.FALSE.
LOGICAL           :: SmagorinskyInitIsDone=.FALSE.
LOGICAL           :: DynSmagInitIsDone=.FALSE.
LOGICAL           :: VremanInitIsDone=.FALSE.
LOGICAL           :: SigmaModelInitIsDone=.FALSE.

END MODULE MOD_EddyVisc_Vars
