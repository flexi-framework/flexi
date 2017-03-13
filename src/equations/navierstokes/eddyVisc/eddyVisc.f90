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
USE MOD_ReadInTools,        ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("EddyViscParameters")
! Smagorinksy model parameters
CALL prms%CreateRealOption(   'CS',       "EddyViscParameters constant")
CALL prms%CreateRealOption(   'PrSGS',    "Turbulent Prandtl number",'0.7')
CALL prms%CreateLogicalOption('VanDriest',"Van Driest damping, only for channel flow!", '.FALSE.')
CALL prms%CreateIntOption('N_Testfilter',"Number of basis for the test filter")
CALL prms%CreateStringOption(   'WallDistFile',    "File containing the wall distance",'')
CALL prms%CreateRealOption(   'WallDistanceTreshold',    "Treshold for zonal model filtering",'0.0')
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
USE MOD_DynSmag
USE MOD_DefaultEddyVisc
USE MOD_Mesh_Vars  ,ONLY:nElems,nSides
USE MOD_ReadInTools,ONLY: GETINTFROMSTR
USE MOD_IO_HDF5    ,ONLY: AddToFieldData, AddToElemData
USE MOD_EOS_Vars, ONLY: mu0
!===================================================================================================================================
eddyViscType = GETINTFROMSTR('eddyViscType')

! Allocate arrays needed by all SGS models
ALLOCATE(DeltaS(nElems))
ALLOCATE(DeltaS_master(1:nSides))
ALLOCATE(DeltaS_slave (1:nSides))
DeltaS_master=0.
DeltaS_slave=0.
DeltaS=0.
ALLOCATE(SGS_Ind(2,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(SGS_Ind_master(2,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(SGS_Ind_slave (2,0:PP_N,0:PP_N,1:nSides))
SGS_Ind=0.
SGS_Ind_master=0.
SGS_Ind_slave=0.
ALLOCATE(muSGS(1,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(muSGSmax(nElems))
muSGS = 0.
muSGSmax=8.*mu0

IF(eddyViscType.EQ.2) THEN !dynamic Smagorinsky
  !MATTEO: debug output
  ALLOCATE(S_en_out(1,0:PP_N,0:PP_N,0:PP_N,nElems))
  S_en_out = 0.
  ALLOCATE(filtdir_out(nElems))
  ALLOCATE(walldist_out(nElems))
  ALLOCATE(walldist_x(nElems))
  ALLOCATE(walldist_y(nElems))
  ALLOCATE(walldist_z(nElems))
  !!!!!!!!!!
  ALLOCATE(FilterMat_Testfilter(0:PP_N,0:PP_N))
  FilterMat_Testfilter = 0.
  ALLOCATE(filter_ind(3,nElems))
  ALLOCATE(average_ind(3,nElems))
  ALLOCATE(average_type(nElems))
  ALLOCATE(IntELem(0:PP_N,0:PP_N,0:PP_N,nElems))
  PrSGS = 0.7
END IF

SELECT CASE(eddyViscType)
  CASE(0) ! No eddy viscosity model, set function pointers to dummy subroutines which do nothing
    ! Nothing to init
    eddyViscosity          => DefaultEddyVisc
    eddyViscosity_surf     => DefaultEddyVisc_surf
    FinalizeEddyViscosity => FinalizeDefaultEddyViscosity
  CASE(1) !Smagorinsky with optional Van Driest damping for channel flow
    CALL InitSmagorinsky()
    eddyViscosity          => Smagorinsky
    eddyViscosity_surf     => Smagorinsky_surf
    FinalizeEddyViscosity => Finalizesmagorinsky
  CASE(2) !Smagorinsky*DynSmag indicator with optional Van Driest damping for channel flow
    CALL InitDynSmag
    eddyViscosity      => DynSmag 
    eddyViscosity_surf => DefaultEddyVisc_surf !NOT used, surface muSGS from volume prolongated 
    FinalizeEddyViscosity => Finalizedynsmag
    testfilter         => Compute_cd
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,&
      'Eddy Viscosity Type not specified!')
END SELECT
CALL AddToFieldData((/1,PP_N+1,PP_N+1,PP_N+1/),'VMSData',(/'muSGS'/),RealArray=muSGS)
!MATTEO: debug output
CALL AddToFieldData((/2,PP_N+1,PP_N+1,PP_N+1/),'VMSData',(/'Csmag   ','muSgsInd'/),RealArray=SGS_Ind)
CALL AddToFieldData((/1,PP_N+1,PP_N+1,PP_N+1/),'VMSData',(/'S_norm'/),RealArray=S_en_out)
CALL AddToElemData('FilterInd',RealArray=filtdir_out(:))
CALL AddToElemData('WallDist',RealArray=walldist_out(:))
CALL AddToElemData('WallDist_x',RealArray=walldist_x(:))
CALL AddToElemData('WallDist_y',RealArray=walldist_y(:))
CALL AddToElemData('WallDist_z',RealArray=walldist_z(:))

END SUBROUTINE

!===================================================================================================================================
!> Finalize eddy viscosity routines
!===================================================================================================================================
SUBROUTINE FinalizeEddyVisc()
! MODULES
USE MOD_EddyVisc_Vars
USE MOD_Smagorinsky   ,ONLY: FinalizeSmagorinsky
!===================================================================================================================================
SDEALLOCATE(DeltaS)
SDEALLOCATE(DeltaS_master)
SDEALLOCATE(DeltaS_slave)
SDEALLOCATE(muSGS)
SDEALLOCATE(muSGSmax)
SDEALLOCATE(SGS_Ind)
SDEALLOCATE(SGS_Ind_master)
SDEALLOCATE(SGS_Ind_slave)
SELECT CASE(eddyViscType)
  CASE(1)
    CALL FinalizeSmagorinsky()
END SELECT
END SUBROUTINE

END MODULE MOD_EddyVisc
