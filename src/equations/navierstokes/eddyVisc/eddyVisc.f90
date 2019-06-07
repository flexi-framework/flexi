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

!===================================================================================================================================
!> In this module the parameters for all the eddy viscosity models will be set, so we only need to call this single routine
!> from flexi.f90.
!===================================================================================================================================
MODULE MOD_EddyVisc
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE DefineParametersEddyVisc
  MODULE PROCEDURE DefineParametersEddyVisc
END INTERFACE

PUBLIC:: DefineParametersEddyVisc,InitEddyVisc,FinalizeEddyVisc
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Define the parameters for all the eddy viscosity models.
!===================================================================================================================================
SUBROUTINE DefineParametersEddyVisc()
! MODULES
USE MOD_ReadInTools,        ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("EddyViscParameters")
! Smagorinksy model parameters
CALL prms%CreateIntFromStringOption('eddyViscType', "(0) none: No eddy viscosity, (1) Smagorinsky",'none')
CALL addStrListEntry(               'eddyViscType','none',0)
CALL addStrListEntry(               'eddyViscType','smagorinsky',1)
CALL prms%CreateRealOption(         'CS',          "EddyViscParameters constant")
CALL prms%CreateRealOption(         'PrSGS',       "Turbulent Prandtl number",'0.7')
CALL prms%CreateLogicalOption(      'VanDriest',   "Van Driest damping, only for channel flow!", '.FALSE.')
END SUBROUTINE DefineParametersEddyVisc

!===================================================================================================================================
!> Initialize eddy viscosity routines
!===================================================================================================================================
SUBROUTINE InitEddyVisc()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_EddyVisc_Vars
USE MOD_Smagorinsky
USE MOD_DefaultEddyVisc
USE MOD_SigmaModel
USE MOD_Mesh_Vars  ,ONLY: nElems,nSides
USE MOD_ReadInTools,ONLY: GETINTFROMSTR, GETREAL
USE MOD_IO_HDF5    ,ONLY: AddToFieldData,FieldOut
USE MOD_EOS_Vars   ,ONLY: mu0
!===================================================================================================================================
eddyViscType = GETINTFROMSTR('eddyViscType')

! Allocate arrays needed by all SGS models
ALLOCATE(DeltaS(nElems))
DeltaS=0.
ALLOCATE(muSGS(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
muSGS = 0.
ALLOCATE(muSGS_master(1,0:PP_N,0:PP_NZ,nSides))
ALLOCATE(muSGS_slave (1,0:PP_N,0:PP_NZ,nSides))
muSGS_master=0.
muSGS_slave =0.
! Turbulent Prandtl number
PrSGS  = GETREAL('PrSGS','0.7')

SELECT CASE(eddyViscType)
  CASE(0) ! No eddy viscosity model, set function pointers to dummy subroutines which do nothing
    ! Nothing to init
    ComputeEddyViscosity  => DefaultEddyVisc
    FinalizeEddyViscosity => FinalizeDefaultEddyViscosity
  CASE(1) !Smagorinsky with optional Van Driest damping for channel flow
    CALL InitSmagorinsky()
    ComputeEddyViscosity  => Smagorinsky_Volume
    FinalizeEddyViscosity => Finalizesmagorinsky
  CASE(4) !sigma Model 2015
    CALL InitSigmaModel()
    ComputeEddyViscosity  => SigmaModel_Volume
    FinalizeEddyViscosity => FinalizeSigmaModel
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,&
      'Eddy Viscosity Type not specified!')
END SELECT
CALL AddToFieldData(FieldOut,(/1,PP_N+1,PP_N+1,PP_NZ+1/),'VMSData',(/'muSGS'/),RealArray=muSGS)
END SUBROUTINE

!===================================================================================================================================
!> Finalize eddy viscosity routines
!===================================================================================================================================
SUBROUTINE FinalizeEddyVisc()
! MODULES
USE MOD_EddyVisc_Vars
!===================================================================================================================================
SDEALLOCATE(DeltaS)
SDEALLOCATE(muSGS)
SDEALLOCATE(muSGS_master)
SDEALLOCATE(muSGS_slave)
IF (ASSOCIATED(FinalizeEddyViscosity)) CALL FinalizeEddyViscosity()
END SUBROUTINE

END MODULE MOD_EddyVisc
