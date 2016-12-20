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
!> \brief Routines that handle restart capabilities.
!> 
!> With this feature a simulation can be resumed from a state file that has been created during a previous
!> simulation (restart file). The restart file is passed to FLEXI as a second command line argument.
!> The restart can also be performed from a file with a different polynomial degree or node type than the current simulation.
!==================================================================================================================================
MODULE MOD_Restart
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitRestart
  MODULE PROCEDURE InitRestart
END INTERFACE

INTERFACE Restart
  MODULE PROCEDURE Restart
END INTERFACE

INTERFACE FinalizeRestart
  MODULE PROCEDURE FinalizeRestart
END INTERFACE

PUBLIC :: InitRestart,FinalizeRestart
PUBLIC :: Restart
!==================================================================================================================================

PUBLIC::DefineParametersRestart
CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE DefineParametersRestart()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Restart")
CALL prms%CreateLogicalOption('ResetTime', "Override solution time to t=0 on restart.", '.FALSE.')
END SUBROUTINE DefineParametersRestart

!==================================================================================================================================
!> \brief Initialize all necessary information to perform the restart.
!>
!> The routine checks if two arguments have been passed to FLEXI on the command line. If so, the second one is supposed
!> to be the restart state. If only one argument has been passed, no restart will be performed.
!> - In the restart case, it is checked if the restart file exists at all. If so, the properties of the restart file will be
!>   read (polynomial degree and node type are needed) and stored to use later.The flag DoRestart indicating the restart is set
!>   to be used by other routines. Also the simulation time of the restart is read.
!>   A optional parameter ResetTime can be used to set the restart time to 0.
!> - If no restart is performed, the RestartTime is set to 0.
!>
!> The routine also checks if the node type and polynomial degree of the restart file is the same than in the current simulation.
!> If not, a flag InterpolateSolution is set. This will be used by the actual Restart routine.
!==================================================================================================================================
SUBROUTINE InitRestart()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Restart_Vars
USE MOD_HDF5_Input,         ONLY: ISVALIDHDF5FILE
USE MOD_Interpolation_Vars, ONLY: InterpolationInitIsDone,NodeType
USE MOD_HDF5_Input,         ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_ReadInTools,        ONLY: GETLOGICAL,GETREALARRAY
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: ResetTime,validHDF5
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.RestartInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,'InitRestart not ready to be called or already called.')
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'

! Check if we want to perform a restart
IF (LEN_TRIM(RestartFile).GT.0) THEN
  SWRITE(UNIT_StdOut,'(A,A,A)')' | Restarting from file "',TRIM(RestartFile),'":'
  ! Check if restart file is a valid state
  validHDF5 = ISVALIDHDF5FILE(RestartFile)
  IF(.NOT.validHDF5) &
      CALL CollectiveStop(__STAMP__,'ERROR - Restart file not a valid state file.')
  ! Set flag indicating a restart to other routines
  DoRestart = .TRUE.
  ! Read in parameters of restart solution
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart)
  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  ! Option to set the calculation time to 0 even tho performing a restart
  ResetTime=GETLOGICAL('ResetTime','.FALSE.')
  IF(ResetTime) RestartTime=0.
  CALL CloseDataFile()
ELSE
  ! No restart
  RestartTime = 0.
  SWRITE(UNIT_StdOut,'(A)')' | No restart wanted, doing a fresh computation!'
END IF

! Check if we need to interpolate the restart file to our current polynomial degree and node type
IF(DoRestart .AND. ((N_Restart.NE.PP_N) .OR. (TRIM(NodeType_Restart).NE.TRIM(NodeType))))THEN
  InterpolateSolution=.TRUE.
#if FV_ENABLED
  CALL CollectiveStop(__STAMP__,'ERROR: The restart to a different polynomial degree is not available for FV.') 
#endif
ELSE
  InterpolateSolution=.FALSE.
END IF

RestartInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RESTART DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRestart


!==================================================================================================================================
!> \brief This routine performs the actual restart. It is called from the FLEXI main program just before the TimeDisc routine.
!>
!> For the restart there are 2 cases, depending on the flag InterpolateSolution set by InitRestart:
!> - No interpolation is necesary. Then simply read in the DG_Solution array from the restart file and store it in the solution
!>   array U.
!> - We need to interpolate the restart solution. If the polynomial degree of our computation is lower than in the restart file,
!>   a simple change basis is used to get the current solution U. If the polynomial degree is higher than in the restart file,
!>   special care is taken to ensure a conservative projection of the restart solution. To do this, the restart solution is
!>   transformed to reference space using the Jacobian build on a polynomial degree of 3*NGeo (so it can be represented exactly),
!>   then the change basis is applied. The resulting solution U is then transformed back to physical space.
!>
!> All state files that would be re-written by the simulation (with the same project name and a time stamp after the restart time
!> or all files with the same project name if no restart is performed) are deleted at the end of the routine.
!==================================================================================================================================
SUBROUTINE Restart()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Restart_Vars
#if FV_ENABLED
USE MOD_HDF5_input         ,ONLY: ReadAttribute
USE MOD_FV_Vars            ,ONLY: FV_Elems
USE MOD_IO_HDF5            ,ONLY: File_ID,nDims,HSize
USE MOD_HDF5_Input         ,ONLY: GetDataSize,DatasetExists
USE MOD_Indicator_Vars     ,ONLY: IndValue
#endif
USE MOD_DG_Vars,            ONLY: U
USE MOD_Mesh_Vars,          ONLY: offsetElem,detJac_Ref,Ngeo
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_HDF5_input,         ONLY: OpenDataFile,CloseDataFile,ReadArray
USE MOD_HDF5_Output,        ONLY: FlushFiles
USE MOD_Interpolation,      ONLY: GetVandermonde
USE MOD_ApplyJacobianCons,  ONLY: ApplyJacobianCons
USE MOD_Interpolation_Vars, ONLY: NodeType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE   :: U_local(:,:,:,:,:)
INTEGER            :: iElem,i,j,k
REAL               :: JNR(1,0:N_Restart,0:N_Restart,0:N_Restart)
REAL               :: Vdm_NRestart_N(0:PP_N,0:N_Restart)
REAL               :: Vdm_3Ngeo_NRestart(0:N_Restart,0:3*NGeo)
#if FV_ENABLED
LOGICAL                        :: ElemDataFound,  FV_ElemsFound
INTEGER                        :: nVarAdd_HDF5,iVarAdd
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAdd_HDF5(:)
REAL,ALLOCATABLE               :: ElemData(:,:)
#endif

!==================================================================================================================================
IF(DoRestart)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#if FV_ENABLED    
  FV_ElemsFound = .FALSE.
  CALL DatasetExists(File_ID, 'ElemData', ElemDataFound)
  IF (ElemDataFound) THEN
    ! get size of ElemData array
    CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
    nVarAdd_HDF5=INT(HSize(1),4)
    ! read ElemData
    ALLOCATE(ElemData(nVarAdd_HDF5,nElems))
    CALL ReadArray('ElemData',2,(/nVarAdd_HDF5,nElems/),OffsetElem,2,RealArray=ElemData)
    ! read variable names of additional data in ElemData
    ALLOCATE(VarNamesAdd_HDF5(nVarAdd_HDF5))
    CALL ReadAttribute(File_ID,'VarNamesAdd',nVarAdd_HDF5,StrArray=VarNamesAdd_HDF5)
    ! search for FV_Elems and IndValue
    DO iVarAdd=1,nVarAdd_HDF5
      IF (TRIM(VarNamesAdd_HDF5(iVarAdd)).EQ."FV_Elems") THEN
        FV_Elems = INT(ElemData(iVarAdd,:))
        FV_ElemsFound = .TRUE.
      END IF
      IF (TRIM(VarNamesAdd_HDF5(iVarAdd)).EQ."IndValue") THEN
        IndValue = ElemData(iVarAdd,:)
      END IF
    END DO
  END IF
#endif
  ! Read in state
  IF(.NOT. InterpolateSolution)THEN
    ! No interpolation needed, read solution directly from file
    CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),OffsetElem,5,RealArray=U)
  ELSE
    ! We need to interpolate the solution to the new computational grid
    SWRITE(UNIT_stdOut,*)'Interpolating solution from restart grid with N=',N_restart,' to computational grid with N=',PP_N

    CALL GetVandermonde(N_Restart, NodeType_Restart,PP_N,      NodeType,         &
                        Vdm_NRestart_N,     modal=.TRUE.)
    CALL GetVandermonde(3*Ngeo,    NodeType,        N_Restart, NodeType_Restart, &
                        Vdm_3Ngeo_NRestart, modal=.TRUE.)

    ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,nElems))
    CALL ReadArray('DG_Solution',5,&
                   (/PP_nVar,N_Restart+1,N_Restart+1,N_Restart+1,nElems/),&
                   OffsetElem,5,RealArray=U_local)

    ! Transform solution to refspace and project solution to N
    ! For conservativity deg of detJac should be identical to EFFECTIVE polynomial deg of solution
    ! (e.g. beware when filtering the jacobian )
    IF(N_Restart.GT.PP_N)THEN
      DO iElem=1,nElems
        IF (FV_Elems(iElem).EQ.0) THEN ! DG element
          CALL ChangeBasisVolume(1,3*Ngeo,N_Restart,Vdm_3Ngeo_NRestart,detJac_Ref(:,:,:,:,iElem),JNR)
          DO k=0,N_Restart; DO j=0,N_Restart; DO i=0,N_Restart
            U_local(:,i,j,k,iElem)=U_local(:,i,j,k,iElem)*JNR(1,i,j,k)
          END DO; END DO; END DO
          CALL ChangeBasisVolume(PP_nVar,N_Restart,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
#if FV_ENABLED          
        ELSE ! FV element
          STOP 'Not implemented yet'
#endif
        END IF
      END DO
      ! Transform back
      CALL ApplyJacobianCons(U,toPhysical=.TRUE.)
    ELSE
      DO iElem=1,nElems
        IF (FV_Elems(iElem).EQ.0) THEN ! DG element
          CALL ChangeBasisVolume(PP_nVar,N_Restart,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
#if FV_ENABLED          
        ELSE ! FV element
          STOP 'Not implemented yet'
#endif
        END IF
      END DO
    END IF

    DEALLOCATE(U_local)
    SWRITE(UNIT_stdOut,*)'DONE!'
  END IF
  CALL CloseDataFile()
  ! Delete all files that will be rewritten
  CALL FlushFiles(RestartTime)
ELSE
  ! Delete all files since we are doing a fresh start
  CALL FlushFiles()
END IF
END SUBROUTINE Restart



!==================================================================================================================================
!> Finalizes variables necessary for restart subroutines
!==================================================================================================================================
SUBROUTINE FinalizeRestart()
! MODULES
USE MOD_Restart_Vars
IMPLICIT NONE
!==================================================================================================================================
RestartInitIsDone = .FALSE.
END SUBROUTINE FinalizeRestart

END MODULE MOD_Restart
