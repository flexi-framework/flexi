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
  SUBROUTINE EddyViscInt()
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
!===================================================================================================================================
!> Finalizes the parameters needed for the eddy viscosity models
!===================================================================================================================================
  SUBROUTINE FinalizeEddyViscosityInt()
  END SUBROUTINE
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                      :: eddyViscType          !< type of eddy viscosity
PROCEDURE(EddyViscInt),POINTER               :: ComputeEddyViscosity  !< pointer to routine for computing volume eddy viscosity
PROCEDURE(FinalizeEddyViscosityInt),POINTER  :: FinalizeEddyViscosity !< pointer tofinalize routine 

!Smagosinsky Standard
REAL,ALLOCATABLE  :: Damp(:,:,:,:,:)       !< damping factor
REAL,ALLOCATABLE  :: DeltaS(:)             !< filter width, used by Smagorinsky modell
REAL,ALLOCATABLE  :: muSGS(:,:,:,:,:)      !< Viscosity for the sub-grid
REAL,ALLOCATABLE  :: muSGS_master(:,:,:,:) !< Viscosity for the sub-grid on master sides
REAL,ALLOCATABLE  :: muSGS_slave (:,:,:,:) !< Viscosity for the sub-grid on slave sides
REAL,ALLOCATABLE  :: muSGSmax(:)           !< Viscosity for the sub-grid
REAL              :: CS                    !< Smagorinsky constant, LES
REAL              :: PrSGS                 !< Prandtl number for the sub-grid scales

LOGICAL           :: VanDriest=.FALSE.     !< Logical indicating if Van Driest damping is activated (only use for channel flow)
LOGICAL           :: SmagorinskyInitIsDone=.FALSE. !< Logical indicating if smagorinsky model has been initialized
LOGICAL           :: SigmaModelInitIsDone=.FALSE.  !< Logical indicating if sigma model has been initialized

END MODULE MOD_EddyVisc_Vars
