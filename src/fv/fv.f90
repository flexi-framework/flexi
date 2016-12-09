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
#if FV_ENABLED
#include "flexi.h"

!==================================================================================================================================
!> Module for the Finite Volume sub-cells shock capturing.
!>
!> DG elements, that are detected to contain a shock/high gradients/oscillations/..., can be switched to a Finite Volume scheme.
!> A DG element of polynomial degree N is subdivided into (N+1)^3 sub-cells (to each Gauss Point/DOF one FV sub-cell).
!> The FV sub-cells of such an element are updated using FV method with 2nd order TVD reconstruction (slope limiters).
!==================================================================================================================================
MODULE MOD_FV
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DefineParametersFV
  MODULE PROCEDURE DefineParametersFV
END INTERFACE

INTERFACE InitFV
  MODULE PROCEDURE InitFV
END INTERFACE

INTERFACE FV_Switch
  MODULE PROCEDURE FV_Switch
END INTERFACE

INTERFACE FV_Info
  MODULE PROCEDURE FV_Info
END INTERFACE

INTERFACE FV_FillIni
  MODULE PROCEDURE FV_FillIni
END INTERFACE

INTERFACE FV_DGtoFV
  MODULE PROCEDURE FV_DGtoFV
END INTERFACE

INTERFACE FinalizeFV
  MODULE PROCEDURE FinalizeFV
END INTERFACE

PUBLIC::DefineParametersFV
PUBLIC::InitFV
PUBLIC::FV_Switch
PUBLIC::FV_Info
PUBLIC::FV_FillIni
PUBLIC::FV_DGtoFV
PUBLIC::FinalizeFV
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FV
!==================================================================================================================================
SUBROUTINE DefineParametersFV()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
#if FV_RECONSTRUCT
USE MOD_FV_Limiter  ,ONLY: DefineParametersFV_Limiter
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('FV')
CALL prms%CreateRealOption('FV_IndUpperThreshold',"Upper threshold: Element is switched from DG to FV if indicator \n"//&
                                                        "rises above this value" )
CALL prms%CreateRealOption('FV_IndLowerThreshold',"Lower threshold: Element is switched from FV to DG if indicator \n"//& 
                                                        "falls below this value")
CALL prms%CreateLogicalOption('FV_toDG_indicator',"Apply additional Persson indicator to check if DG solution after switch \n"//&
                                                  "from FV to DG is valid.", '.FALSE.')
CALL prms%CreateRealOption   ('FV_toDG_limit',"Threshold for FV_toDG_indicator")
CALL prms%CreateLogicalOption('FV_toDGinRK',  "Allow switching of FV elements to DG during Runge Kutta stages. \n"//&
                                              "This may violated the DG timestep restriction of the element.", '.FALSE.')
#if FV_RECONSTRUCT
CALL DefineParametersFV_Limiter()
#endif
END SUBROUTINE DefineParametersFV

!==================================================================================================================================
!> Read in parameters needed for FV sub-cells (indicator min/max and type of limiter) and allocate several arrays.
!> Build metrics for FV sub-cells and performe initial switch from DG to FV sub-cells for all troubled cells.
!==================================================================================================================================
SUBROUTINE InitFV()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_FV_Basis
USE MOD_Basis        ,ONLY: InitializeVandermonde
USE MOD_Indicator    ,ONLY: doCalcIndicator
USE MOD_Mesh_Vars    ,ONLY: nElems,nSides
USE MOD_FV_Limiter
USE MOD_ReadInTools
USE MOD_IO_HDF5      ,ONLY: AddToElemData
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
IF(.NOT.FVInitBasisIsDone)THEN
   CALL CollectiveStop(__STAMP__,&
     'InitFV not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FV...'

doCalcIndicator=.TRUE.
! read minimal and maximal threshold for the indicator
FV_IndLowerThreshold = GETREAL('FV_IndLowerThreshold','-99.')
FV_IndUpperThreshold = GETREAL('FV_IndUpperThreshold', '99.')


#if FV_RECONSTRUCT
CALL InitFV_Limiter()
#endif

! read flag indicating, if an additional Persson indicator should check if a FV sub-cells element really contains no oscillations
! anymore. 
FV_toDG_indicator = GETLOGICAL('FV_toDG_indicator')
IF (FV_toDG_indicator) FV_toDG_limit = GETREAL('FV_toDG_limit')
FV_toDGinRK = GETLOGICAL("FV_toDGinRK")


! allocate array for indicators
ALLOCATE(FV_Elems(nElems))
! all cells are initially DG cells
FV_Elems = 0
CALL AddToElemData('FV_Elems',IntArray=FV_Elems)

ALLOCATE(FV_Elems_counter(nElems))
ALLOCATE(FV_Elems_Amount(nElems))
FV_Elems_counter  = 0
FV_Switch_counter = 0
FV_Elems_Amount = 0
CALL AddToElemData('FV_Elems_Amount',RealArray=FV_Elems_Amount)

! allocate arrays for indicators at faces
!ALLOCATE(FV_Elems_master(1:nSides)) ! moved to InitFV_Metrics, since needed there for U_Mortar routine
ALLOCATE(FV_Elems_slave(1:nSides))
ALLOCATE(FV_Elems_Sum(1:nSides))
FV_Elems_master = 0
FV_Elems_slave = 0
FV_Elems_Sum = 0

#if FV_RECONSTRUCT
! allocate array for multipurpose 
ALLOCATE(FV_multi_master(1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FV_multi_slave (1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides))
FV_multi_slave = 0.0
FV_multi_master = 0.0

! allocate array for FD-gradient over faces
ALLOCATE(FV_surf_gradU_master(1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)) 
ALLOCATE(FV_surf_gradU_slave (1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)) 
FV_surf_gradU_master = 0.0
FV_surf_gradU_slave = 0.0

! The gradients of the conservative variables are stored at each volume integration point
ALLOCATE(gradUxi(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUeta(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUzeta(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems))
gradUxi=0.
gradUeta=0.
gradUzeta=0.
ALLOCATE(gradUxi_central  (PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUeta_central (PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(gradUzeta_central(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems))
gradUxi_central  =0.
gradUeta_central =0.
gradUzeta_central=0.
#endif

FVInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT FV DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitFV

!==================================================================================================================================
!> Performe switching between DG element and FV sub-cells element (and vise versa) depending on the indicator value
!==================================================================================================================================
SUBROUTINE FV_Switch(AllowToDG)
! MODULES
USE MOD_PreProc
USE MOD_ChangeBasis    ,ONLY: ChangeBasis3D
USE MOD_DG_Vars        ,ONLY: U
USE MOD_Indicator_Vars ,ONLY: IndValue
USE MOD_Indicator      ,ONLY: IndPersson
USE MOD_FV_Vars
USE MOD_Analyze
USE MOD_Mesh_Vars      ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: AllowToDG
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: U_DG(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL    :: U_FV(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL    :: ind
INTEGER :: iElem
!==================================================================================================================================
DO iElem=1,nElems
  IF (FV_Elems(iElem).EQ.0) THEN ! DG Element
    ! Switch DG to FV Element, if Indicator is higher then IndMin
    IF (IndValue(iElem).GT.FV_IndUpperThreshold) THEN
      ! switch Element to FV
      FV_Elems(iElem) = 1
      CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_Vdm,U(:,:,:,:,iElem),U_FV)
      U(:,:,:,:,iElem) = U_FV
    END IF
  ELSE ! FV Element
    ! Switch FV to DG Element, if Indicator is lower then IndMax
    IF ((IndValue(iElem).LT.FV_IndLowerThreshold).AND.AllowToDG) THEN
      CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_sVdm,U(:,:,:,:,iElem),U_DG)
      IF (FV_toDG_indicator) THEN
        ind = IndPersson(U_DG(1,:,:,:))
        IF (ind.GT.FV_toDG_limit) CYCLE
      END IF
      U(:,:,:,:,iElem) = U_DG
      FV_Elems(iElem) = 0  ! switch Elemnent to DG
    END IF
  END IF
END DO !iElem
FV_Elems_counter  = FV_Elems_counter  + FV_Elems
FV_Switch_counter = FV_Switch_counter + 1
FV_Elems_Amount   = REAL(FV_Elems_Counter)/FV_Switch_counter
END SUBROUTINE FV_Switch

!==================================================================================================================================
!> Print information on the amount of FV subcells
!==================================================================================================================================
SUBROUTINE FV_Info(iter)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nGlobalElems
USE MOD_Analyze_Vars,ONLY: totalFV_nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN) :: iter !< number of iterations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,totalFV_nElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
  ! totalFV_nElems is counted in PrintStatusLine
ELSE
  CALL MPI_REDUCE(totalFV_nElems,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(A,F8.3,A)')' FV amount %: ', totalFV_nElems / REAL(nGlobalElems) / iter*100 
totalFV_nElems = 0
END SUBROUTINE FV_Info


!==================================================================================================================================
!> Initialize all FV elements and overwrite data of DG FillIni. Each subcell is supersampled with PP_N points in each space 
!> dimension and the mean value is taken as value for this subcell.
!==================================================================================================================================
SUBROUTINE FV_FillIni()
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars           ,ONLY: U
USE MOD_Mesh_Vars         ,ONLY: nElems
USE MOD_FV_Vars           ,ONLY: FV_Elems
USE MOD_Indicator         ,ONLY: CalcIndicator
USE MOD_Interpolation     ,ONLY: GetVandermonde 
USE MOD_ChangeBasis       ,ONLY: ChangeBasis3D
USE MOD_Mesh_Vars         ,ONLY: Elem_xGP
USE MOD_Equation_Vars     ,ONLY: IniExactFunc
USE MOD_Exactfunc         ,ONLY: ExactFunc
USE MOD_Interpolation_Vars,ONLY: NodeTypeG, NodeTypeVISUInner
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,iElem, j,k,ii,jj,kk,iVar
REAL              :: Vdm(0:(PP_N+1)**2-1,0:PP_N), xx(1:3,0:(PP_N+1)**2-1,0:(PP_N+1)**2-1,0:(PP_N+1)**2-1)
REAL              :: tmp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
!===================================================================================================================================
! initial call of indicator
CALL CalcIndicator(U,0.)
FV_Elems = 0
CALL FV_Switch(AllowToDG=.FALSE.)

! build vandermonde to supersample each subcell with PP_N points per direction
CALL GetVandermonde(PP_N,NodetypeG,(PP_N+1)**2-1,NodeTypeVISUInner,Vdm)
DO iElem=1,nElems
  IF (FV_Elems(iElem).EQ.0) CYCLE ! DG element
  ! supersample all subcells
  CALL ChangeBasis3D(3,PP_N,(PP_N+1)**2-1,Vdm,Elem_xGP(1:3,:,:,:,iElem),xx)
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        ! evaluate ExactFunc for all supersampled points of subcell (i,j,k)
        DO kk=0,PP_N; DO jj=0,PP_N; DO ii=0,PP_N
          CALL ExactFunc(IniExactFunc,0.,xx(1:3,i*(PP_N+1)+ii,j*(PP_N+1)+jj,k*(PP_N+1)+kk),tmp(:,ii,jj,kk))
        END DO; END DO; END DO
        ! mean value 
        DO iVar=1,PP_nVar
          U(iVar,i,j,k,iElem) = SUM(tmp(iVar,:,:,:)) / (PP_N+1)**3
        END DO
      END DO ! i
    END DO ! j
  END DO !k
END DO ! iElem=1,nElems
END SUBROUTINE FV_FillIni

!==================================================================================================================================
!> Switch DG solution at faces between a DG element and a FV sub-cells element to Finite Volume.
!==================================================================================================================================
SUBROUTINE FV_DGtoFV(nVar,U_master,U_slave)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ChangeBasis ,ONLY: ChangeBasis2D_selective
USE MOD_FV_Vars
USE MOD_Mesh_Vars   ,ONLY: firstInnerSide,lastMPISide_MINE,nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVar
REAL,INTENT(INOUT) :: U_master(nVar,0:PP_N,0:PP_N,1:nSides) !< Solution on master side
REAL,INTENT(INOUT) :: U_slave (nVar,0:PP_N,0:PP_N,1:nSides) !< Solution on slave side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: firstSideID,lastSideID
!==================================================================================================================================
firstSideID = firstInnerSide
lastSideID  = lastMPISide_MINE

CALL ChangeBasis2D_selective(nVar,PP_N,1,nSides,firstSideID,lastSideID,FV_Vdm,U_master,FV_Elems_Sum,2)
CALL ChangeBasis2D_selective(nVar,PP_N,1,nSides,firstSideID,lastSideID,FV_Vdm,U_slave ,FV_Elems_Sum,1)
END SUBROUTINE FV_DGtoFV


!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeFV()
! MODULES
USE MOD_FV_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(FV_Elems)
!SDEALLOCATE(FV_Elems_master) ! moved to mesh.f90
SDEALLOCATE(FV_Elems_slave)
SDEALLOCATE(FV_Elems_Counter)
SDEALLOCATE(FV_Elems_Amount)
SDEALLOCATE(FV_Elems_Sum)
#if FV_RECONSTRUCT
SDEALLOCATE(FV_surf_gradU_master)
SDEALLOCATE(FV_surf_gradU_slave)
SDEALLOCATE(FV_multi_master)
SDEALLOCATE(FV_multi_slave)
SDEALLOCATE(gradUxi)
SDEALLOCATE(gradUeta)
SDEALLOCATE(gradUzeta)
SDEALLOCATE(gradUxi_central)
SDEALLOCATE(gradUeta_central)
SDEALLOCATE(gradUzeta_central)
#endif

FVInitIsDone=.FALSE.
END SUBROUTINE FinalizeFV


END MODULE MOD_FV
#endif /* FV_ENABLED */
