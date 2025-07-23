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

PUBLIC:: DefineParametersRestart
PUBLIC:: InitRestartFile
PUBLIC:: InitRestart
PUBLIC:: Restart
#if FV_ENABLED
PUBLIC:: SupersampleFVCell
#endif
PUBLIC:: FinalizeRestart
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE DefineParametersRestart()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Restart")
CALL prms%CreateLogicalOption('ResetTime', "Override solution time to t=0 on restart.", '.FALSE.')
#if FV_ENABLED
CALL prms%CreateIntOption(    'NFVRestartSuper', "Polynomial degree for equidistant supersampling of FV subcells when restarting&
                                                  &on a different polynomial degree. Default 2*MAX(N,NRestart).")
#endif
CALL prms%CreateLogicalOption('FlushInitialState',"Check whether (during restart) the statefile from which the restart is performed&
                                                  &should be deleted.", '.FALSE.')
END SUBROUTINE DefineParametersRestart


!==================================================================================================================================
!> \brief Check the presence and file type of the restart file (state, timeAvg)
!>
!> The routine checks if two arguments have been passed to FLEXI on the command line. If so, the second one is supposed
!> to be the restart state. If only one argument has been passed, no restart will be performed.
!> - In the restart case, it is checked if the restart file exists at all. If so, the type of the restart file will be assessed
!>   (state or timeAvg file). If the latter one is detected, the presence of all conservative or primitive variables is checked
!>   and PrimToCons performed if necessary. The type of restart file is indicated with the variable RestartMode.
!==================================================================================================================================
SUBROUTINE InitRestartFile(RestartFile_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars      ,ONLY: StrVarNames,StrVarNamesPrim
USE MOD_HDF5_Input         ,ONLY: ISVALIDHDF5FILE,GetVarNames,DatasetExists
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_ReadInTools        ,ONLY: GETLOGICAL,GETREAL
USE MOD_Restart_Vars
#if FV_ENABLED
USE MOD_StringTools        ,ONLY: INTTOSTR
USE MOD_ReadInTools        ,ONLY: GETINT
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: RestartFile_in !< state file to restart from
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: validHDF5
LOGICAl            :: RestartMean,VarNamesExist
INTEGER            :: iVar
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames_tmp(:)
REAL                            :: StartT,EndT
!==================================================================================================================================
RestartFile = RestartFile_in

! Check if we want to perform a restart
IF (LEN_TRIM(RestartFile).LE.0) RETURN

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' CHECK RESTART FILE...'
SWRITE(UNIT_stdOut,'(A,A,A)')' | Checking restart from file "',TRIM(RestartFile),'":'
GETTIME(StartT)

! Check if restart file is a valid state. This routine requires the file to be closed.
validHDF5 = ISVALIDHDF5FILE(RestartFile)
IF(.NOT.validHDF5) &
    CALL CollectiveStop(__STAMP__,'ERROR - Restart file not a valid state file.')

! Read in parameters of restart solution
CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#if EQNSYSNR != 1 && EQNSYSNR != 4
! Check if the file is a time-averaged file
CALL DatasetExists(File_ID,'Mean',RestartMean)
! Read in attributes
IF (.NOT.RestartMean) THEN
#endif /* EQNSYSNR != 1 && EQNSYSNR != 4 */
  CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart)
  RestartMode = 1
  SWRITE(UNIT_StdOut,'(A)') ' | Restarting from state file ...'
#if EQNSYSNR != 1 && EQNSYSNR != 4
ELSE
  CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart,'Mean')
  ! Get the VarNames to compare later
  VarNamesExist = .FALSE.
  SDEALLOCATE(VarNames_tmp)
  CALL GetVarNames('VarNames_Mean',VarNames_tmp,VarNamesExist)

  ! Check if variables match
  RestartCons = -1
  RestartPrim = -1
  DO iVar = 1,nVar_Restart
    SELECT CASE(TRIM(VarNames_tmp(iVar)))
      CASE(TRIM(StrVarNames(1)))     ! Density
        RestartCons(1) = iVar
        RestartPrim(1) = iVar
      CASE(TRIM(StrVarNames(2)))     ! MomentumX
        RestartCons(2) = iVar
      CASE(TRIM(StrVarNames(3)))     ! MomentumY
        RestartCons(3) = iVar
      CASE(TRIM(StrVarNames(4)))     ! MomentumZ
        RestartCons(4) = iVar
      CASE(TRIM(StrVarNames(5)))     ! EnergyStagnationDensity
        RestartCons(5) = iVar
      CASE(TRIM(StrVarNamesPrim(2))) ! VelocityX
        RestartPrim(2) = iVar
      CASE(TRIM(StrVarNamesPrim(3))) ! VelocityY
        RestartPrim(3) = iVar
      CASE(TRIM(StrVarNamesPrim(4))) ! VelocityZ
        RestartPrim(4) = iVar
      CASE(TRIM(StrVarNamesPrim(5))) ! Pressure
        RestartPrim(5) = iVar
      CASE(TRIM(StrVarNamesPrim(6))) ! Temperature
        RestartPrim(6) = iVar
#if EQNSYSNR == 3
      CASE(TRIM(StrVarNames(6)))     ! muTilde
        RestartCons(6) = iVar
      CASE(TRIM(StrVarNamesPrim(7))) ! nuTilde
        RestartPrim(7) = iVar
#endif /* EQNSYSNR == 3 */
    END SELECT
  END DO
  ! Use conservative variables available
  IF (ALL(RestartCons.NE.-1)) THEN
    RestartMode = 2
    SWRITE(UNIT_stdOut,'(A)') ' | Restarting from time-averaged file using conservative variables...'
  ELSE IF (ALL(RestartPrim.NE.-1)) THEN
    RestartMode = 3
    SWRITE(UNIT_stdOut,'(A)') ' | Restarting from time-averaged file using primitive variables...'
  ELSE
    RestartMode = 0
  END IF
END IF
#endif /* EQNSYSNR != 1 && EQNSYSNR != 4 */

CALL CloseDataFile()

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'CHECK RESTART FILE DONE!', DisplayLine=.TRUE.)

END SUBROUTINE InitRestartFile


!==================================================================================================================================
!> \brief Initialize all necessary information to perform the restart.
!>
!> The routine checks if two arguments have been passed to FLEXI on the command line. If so, the second one is supposed
!> to be the restart state. If only one argument has been passed, no restart will be performed.
!> - In the restart case, it is checked if the restart file exists at all. If so, the properties of the restart file will be
!>   read (polynomial degree and node type are needed) and stored to use later.The flag DoRestart indicating the restart is set
!>   to be used by other routines. Also the simulation time of the restart is read.
!>   An optional parameter ResetTime can be used to set the restart time to 0.
!> - If no restart is performed, the RestartTime is set to 0.
!>
!> The routine also checks if the node type and polynomial degree of the restart file is the same than in the current simulation.
!> If not, a flag InterpolateSolution is set. This will be used by the actual Restart routine.
!==================================================================================================================================
SUBROUTINE InitRestart(RestartFile_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input,         ONLY: ISVALIDHDF5FILE,GetVarNames,DatasetExists
USE MOD_HDF5_Input,         ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_Interpolation_Vars, ONLY: InterpolationInitIsDone,NodeType
USE MOD_Mesh_Vars,          ONLY: nGlobalElems,NGeo
USE MOD_ReadInTools,        ONLY: GETLOGICAL,GETREAL!,ExtractParameterFile,CompareParameterFile
USE MOD_Restart_Vars
#if FV_ENABLED
USE MOD_ReadInTools,        ONLY: GETINT
USE MOD_StringTools,        ONLY: INTTOSTR
#endif /*FV_ENABLED*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: RestartFile_in !< state file to restart from
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: ResetTime,validHDF5,WriteSuccessful!,prmChanged,userblockFound
REAL               :: StartT,EndT
! CHARACTER(LEN=255) :: ParameterFileOld
!==================================================================================================================================
IF(.NOT.InterpolationInitIsDone .OR. RestartInitIsDone) &
  CALL CollectiveStop(__STAMP__,'InitRestart not ready to be called or already called.')

! If not done previously, check the restart file
IF (RestartMode.EQ.-1) CALL InitRestartFile(RestartFile_in)

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'
GETTIME(StartT)

! Check if we want to perform a restart
IF (LEN_TRIM(RestartFile).GT.0) THEN
  ! Restart not possible, some variables might be missing
  IF (RestartMode.LT. 1 .AND. .NOT. postiMode) CALL CollectiveStop(__STAMP__, &
    'Provided file for restart has not all conservative/primitive variables available!')

  SWRITE(UNIT_stdOut,'(A,A,A)')' | Restarting from file "',TRIM(RestartFile),'":'
  ! Check if restart file is a valid state
  validHDF5 = ISVALIDHDF5FILE(RestartFile)
  IF(.NOT.validHDF5) CALL CollectiveStop(__STAMP__,'ERROR - Restart file not a valid state file.')
  ! Set flag indicating a restart to other routines
  DoRestart = .TRUE.
  ! Read in parameters of restart solution
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  ! Check if restart file was written successfully
  CALL DatasetExists(File_ID,'TIME',WriteSuccessful,attrib=.TRUE.)
  IF (.NOT.WriteSuccessful) &
    CALL Abort(__STAMP__,'Restart file missing WriteSuccessful marker. Aborting...')

  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  ! Option to set the calculation time to 0 even tho performing a restart
  ResetTime = GETLOGICAL('ResetTime')
  IF (postiMode) ResetTime = .FALSE.
  IF (ResetTime) RestartTime = 0.
  CALL CloseDataFile()

  ! ! Ensure this is not the same run starting over with ResetTime=T
  ! IF (ResetTime) THEN
  !   IF (.NOT.GETLOGICAL('ResetTimeOverride')) THEN
  !     ! Extract the old parameter file
  !     IF (MPIRoot) THEN
  !       ParameterFileOld = ".flexi.old.ini"
  !       CALL ExtractParameterFile(RestartFile,ParameterFileOld,userblockFound)
  !
  !       ! Compare it against the current file
  !       IF (userblockFound) THEN
  !         CALL CompareParameterFile(ParameterFile,ParameterFileOld,prmChanged)
  !         IF (.NOT.prmChanged) &
  !           CALL Abort(__STAMP__,'Running simulation with ResetTime=T, same parameter file and ResetTimeOverride=F!')
  !       END IF ! userblockFound
  !     END IF ! MPIRoot
  !   END IF ! .NOT.ResetTimeOverride
  ! END IF ! ResetTime

  ! Check if number of elements match
  IF (nElems_Restart.NE.nGlobalElems) THEN
    CALL CollectiveStop(__STAMP__, "Restart File has different number of elements!")
  END IF
! No restart
ELSE
  RestartTime = 0.
  SWRITE(UNIT_stdOut,'(A)')' | No restart wanted, doing a fresh computation!'
END IF

! Check if we need to interpolate the restart file to our current polynomial degree and node type
IF(DoRestart .AND. ((N_Restart.NE.PP_N) .OR. (TRIM(NodeType_Restart).NE.TRIM(NodeType))))THEN
  InterpolateSolution = .TRUE.
  IF(MIN(N_Restart,PP_N).LT.NGeo) &
    CALL PrintWarning('The geometry is or was underresolved and will potentially change on restart!')
#if FV_ENABLED
  NFVRestartSuper = GETINT('NFVRestartSuper',INTTOSTR(2*MAX(PP_N,N_Restart)))
#endif
ELSE
  InterpolateSolution=.FALSE.
END IF

! Check whether (during restart) the statefile from which the restart is performed should be deleted
FlushInitialState = GETLOGICAL('FlushInitialState')

RestartWallTime   = FLEXITIME()
RestartInitIsDone = .TRUE.
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'INIT RESTART DONE!', DisplayLine=.TRUE.)

END SUBROUTINE InitRestart


!==================================================================================================================================
!> \brief This routine performs the actual restart. It is called from the FLEXI main program just before the TimeDisc routine.
!>
!> For the restart there are 2 cases, depending on the flag InterpolateSolution set by InitRestart:
!> - No interpolation is necessary. Then simply read in the DG_Solution array from the restart file and store it in the solution
!>   array U.
!> - We need to interpolate the restart solution. If the polynomial degree of our computation is lower than in the restart file,
!>   a simple change basis is used to get the current solution U. If the polynomial degree is higher than in the restart file,
!>   special care is taken to ensure a conservative projection of the restart solution. To do this, the restart solution is
!>   transformed to reference space using the Jacobian built on a polynomial degree of 3*NGeo (so it can be represented exactly),
!>   then the change basis is applied. The resulting solution U is then transformed back to physical space.
!>
!> All state files that would be re-written by the simulation (with the same project name and a time stamp after the restart time
!> or all files with the same project name if no restart is performed) are deleted at the end of the routine.
!==================================================================================================================================
SUBROUTINE Restart(doFlushFiles)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ApplyJacobianCons,  ONLY: ApplyJacobianCons
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_DG_Vars,            ONLY: U
USE MOD_EOS,                ONLY: PrimToCons
USE MOD_HDF5_Input,         ONLY: GetDataSize
USE MOD_HDF5_Input,         ONLY: OpenDataFile,CloseDataFile,ReadArray,GetArrayAndName
USE MOD_HDF5_Output,        ONLY: FlushFiles
USE MOD_Interpolation,      ONLY: GetVandermonde
USE MOD_Interpolation_Vars, ONLY: NodeType
USE MOD_IO_HDF5
USE MOD_Mesh_Vars,          ONLY: offsetElem,detJac_Ref,Ngeo
USE MOD_Mesh_Vars,          ONLY: nElems,nGlobalElems
USE MOD_Restart_Vars
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems
#endif /*FV_ENABLED*/
#if FV_ENABLED == 1
USE MOD_FV_Switching,       ONLY: FV_ProlongFVElemsToFace
USE MOD_Indicator_Vars,     ONLY: IndValue
USE MOD_StringTools,        ONLY: STRICMP
#endif /*FV_ENABLED==1*/
#if PP_dim == 3
USE MOD_2D,                 ONLY: ExpandArrayTo3D
#else
USE MOD_2D,                 ONLY: to2D_rank5
#endif /*PP_dim*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: doFlushFiles !< flag to delete old state files
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE   :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE   :: U_localNVar(:,:,:,:,:)
REAL,ALLOCATABLE   :: U_local2(:,:,:,:,:)
INTEGER            :: iElem,i,j,k
INTEGER            :: iVar
INTEGER            :: HSize_proc(5)
REAL,ALLOCATABLE   :: JNR(:,:,:,:)
REAL               :: Vdm_NRestart_N(0:PP_N,0:N_Restart)
REAL               :: Vdm_3Ngeo_NRestart(0:N_Restart,0:3*NGeo)
LOGICAL            :: doFlushFiles_loc
#if FV_ENABLED == 1
INTEGER             :: nVal(15)
REAL,ALLOCATABLE    :: ElemData(:,:),tmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesElemData(:)
#endif /*FV_ENABLED==1*/
! Timers
REAL                            :: StartT,EndT
!==================================================================================================================================

IF (PRESENT(doFlushFiles)) THEN; doFlushFiles_loc = doFlushFiles
ELSE                           ; doFlushFiles_loc = .TRUE.
END IF

IF (DoRestart) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' PERFORMING RESTART...'
  GETTIME(StartT)

  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#if FV_ENABLED == 1
  ! Read FV element distribution and indicator values from elem data array if possible
  CALL GetArrayAndName('ElemData','VarNamesAdd',nVal,tmp,VarNamesElemData)
  IF (ALLOCATED(VarNamesElemData)) THEN
    ALLOCATE(ElemData(nVal(1),nVal(2)))
    ElemData = RESHAPE(tmp,(/nVal(1),nVal(2)/))
    ! search for FV_Elems and IndValue
    FV_Elems=0
    IndValue=0.
    DO iVar=1,nVal(1)
      IF (STRICMP(VarNamesElemData(iVar),"FV_Elems")) THEN
        FV_Elems = INT(ElemData(iVar,:))
      END IF
      IF (STRICMP(VarNamesElemData(iVar),"IndValue")) THEN
        IndValue = ElemData(iVar,:)
      END IF
    END DO
  ELSE
    IndValue=0.
    FV_Elems=0.
  END IF
  SDEALLOCATE(ElemData)
  SDEALLOCATE(VarNamesElemData)
  SDEALLOCATE(tmp)
  CALL FV_ProlongFVElemsToFace()
#endif /*FV_ENABLED==1*/

  ! Mean files only have a dummy DG_Solution, we have to pick the "Mean" array in this case
  IF (RestartMode.GT.1) THEN
    CALL GetDataSize(File_ID,'Mean'       ,nDims,HSize)
  ELSE
    CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
  END IF

  ! Sanity check, number of elements
  IF ((HSize(2).NE.N_Restart+1).OR.(HSize(3).NE.N_Restart+1).OR.(HSize(5).NE.nGlobalElems)) &
    CALL Abort(__STAMP__,"Dimensions of restart file do not match!")

  HSize_proc    = INT(HSize)
  HSize_proc(5) = nElems
  ! Allocate array to hold the restart data
  ALLOCATE(U_local(nVar_Restart,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
  DEALLOCATE(HSize)
  ! Mean files only have a dummy DG_Solution, we have to pick the "Mean" array in this case
  IF (RestartMode.GT.1) THEN
    CALL ReadArray('Mean'       ,5,HSize_proc,OffsetElem,5,RealArray=U_local)
  ELSE
    CALL ReadArray('DG_Solution',5,HSize_proc,OffsetElem,5,RealArray=U_local)
  END IF

  ! Truncate the solution if we read a restart file from a different equation system or from a time-averaged file
  SELECT CASE(RestartMode)
    ! Conservative Variables, truncated
    CASE(1)
      IF (PP_nVar.LT.nVar_Restart) THEN
        ALLOCATE(U_localNVar(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
        ! Pass to truncated variables
        U_localNVar(1:PP_nVar,:,:,:,:) = U_local(1:PP_nVar,:,:,:,:)
      END IF
    ! Conservative Variables, time-averaged
    ! Pass to corresponding variables from the time-averaged file
    CASE(2)
      ! Variables might not be continuous
      ALLOCATE(U_localNVar(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
      DO iVar = 1,PP_nVar
        U_localNVar(iVar,:,:,:,:) = U_local(RestartCons(iVar),:,:,:,:)
      END DO
    ! Primitive Variables, time-averaged
    ! Pass to corresponding variables from the time-averaged file
    CASE(3)
      ALLOCATE(U_localNVar(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
      ! Variables might not be continuous
      ALLOCATE(U_local2(1:PP_nVarPrim,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
      DO iVar = 1,PP_nVarPrim
        U_local2(iVar,:,:,:,:) = U_local(RestartPrim(iVar),:,:,:,:)
      END DO
      CALL PrimToCons(HSize_proc(2)-1,U_local2(:,:,:,:,:),U_localNVar(RestartCons(1):RestartCons(1)+PP_nVar-1,:,:,:,:))
      DEALLOCATE(U_local2)
  END SELECT

  ! Write truncated array back to U_local
  IF (ALLOCATED(U_localNVar)) THEN
    CALL MOVE_ALLOC(U_localNVar,U_local)
  END IF

  ! Read in state
  IF(.NOT. InterpolateSolution)THEN
    ! No interpolation needed, read solution directly from file
#if PP_dim == 3
    IF (HSize_proc(4).EQ.1) THEN
      ! FLEXI compiled 3D, but data is 2D => expand third space dimension
      CALL ExpandArrayTo3D(5,(/PP_nVar,PP_N+1,PP_N+1,1,nElems/),4,PP_N+1,U_local,U)
    ELSE
      ! FLEXI compiled 3D + data 3D
      U = U_local
    END IF
#else
    IF (HSize_proc(4).EQ.1) THEN
      ! FLEXI compiled 2D + data 2D
      U = U_local
    ELSE
      ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
      CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVar,PP_N,PP_N,PP_N,nElems/),4,U_local)
      U = U_local
    END IF
#endif
  ELSE ! InterpolateSolution
    ! We need to interpolate the solution to the new computational grid
    SWRITE(UNIT_stdOut,'(A,I0,3A,I0,3A)') ' | Interpolating solution from restart grid with N=',N_restart,' (',TRIM(NodeType_Restart), &
                                          ') to computational grid with N='                    ,PP_N     ,' (',TRIM(NodeType),')'

    CALL GetVandermonde(N_Restart, NodeType_Restart,PP_N,      NodeType,         &
                        Vdm_NRestart_N,     modal=.TRUE.)
    CALL GetVandermonde(3*Ngeo,    NodeType,        N_Restart, NodeType_Restart, &
                        Vdm_3Ngeo_NRestart, modal=.TRUE.)

#if PP_dim == 3
    IF (HSize_proc(4).EQ.1) THEN
      ! FLEXI compiled 3D, but data is 2D => expand third space dimension
      ! use temporary array 'U_local2' to store 3D data
      ALLOCATE(U_local2(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,nElems))
      CALL ExpandArrayTo3D(5,HSize_proc,4,N_Restart,U_local,U_local2)
      ! Reallocate 'U_local' to 3D and mv data from U_local2 to U_local
      CALL MOVE_ALLOC(U_local2,U_local)
    END IF
#else
    IF (HSize_proc(4).NE.1) THEN
      ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
      CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVar,N_Restart,N_Restart,N_Restart,nElems/),4,U_local)
    END IF
#endif /*PP_dim == 3*/

    ! Transform solution to refspace and project solution to N
    ! For conservativity deg of detJac should be identical to EFFECTIVE polynomial deg of solution
    ! (e.g. beware when filtering the jacobian )
    IF(N_Restart.GT.PP_N)THEN
      ALLOCATE(JNR(1,0:N_Restart,0:N_Restart,0:N_Restart*(PP_dim-2)))
      DO iElem=1,nElems
        IF (FV_Elems(iElem).EQ.0) THEN ! DG element
          CALL ChangeBasisVolume(1,3*Ngeo,N_Restart,Vdm_3Ngeo_NRestart,detJac_Ref(:,:,:,:,iElem),JNR)
          DO k=0,N_Restart*(PP_dim-2); DO j=0,N_Restart; DO i=0,N_Restart
            U_local(:,i,j,k,iElem)=U_local(:,i,j,k,iElem)*JNR(1,i,j,k)
          END DO; END DO; END DO
          CALL ChangeBasisVolume(PP_nVar,N_Restart,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
#if FV_ENABLED
        ELSE ! FV element
          CALL SupersampleFVCell(U_local(:,:,:,:,iElem),U(:,:,:,:,iElem),N_Restart,PP_N,NFVRestartSuper)
#endif
        END IF
      END DO
      DEALLOCATE(JNR)
      ! Transform back
      CALL ApplyJacobianCons(U,toPhysical=.TRUE.,FVE=0)
    ELSE
      DO iElem=1,nElems
        IF (FV_Elems(iElem).EQ.0) THEN ! DG element
          CALL ChangeBasisVolume(PP_nVar,N_Restart,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
#if FV_ENABLED
        ELSE ! FV element
          CALL SupersampleFVCell(U_local(:,:,:,:,iElem),U(:,:,:,:,iElem),N_Restart,PP_N,NFVRestartSuper)
#endif
        END IF
      END DO
    END IF ! N_Restart.GT.PP_N

    DEALLOCATE(U_local)
  END IF ! InterpolateSolution
  CALL CloseDataFile()

  IF (RestartMode.GT.1) THEN
    SWRITE(UNIT_stdOut,'(A,ES14.7)') ' | Restart from time-averaged file successful t = ',RestartTime
  ELSE
    SWRITE(UNIT_stdOut,'(A,ES14.7)') ' | Restart from state file successful at t = '     ,RestartTime
  END IF

  ! Delete all files that will be rewritten
  IF (doFlushFiles_loc) CALL FlushFiles(RestartTime)
  GETTIME(EndT)
  CALL DisplayMessageAndTime(EndT-StartT, 'PERFORMING RESTART DONE!', DisplayLine=.TRUE.)
ELSE
  ! Delete all files since we are doing a fresh start
  IF (doFlushFiles_loc) CALL FlushFiles()
END IF

END SUBROUTINE Restart


#if FV_ENABLED
!==================================================================================================================================
!> This routine will take a FV element with a certain number of subcells and convert it to a different number of subcells.
!> Is used during the restart from one polynomial degree to another.
!> The procedure is as follows: The old solution will be supersampled with an adjustable number of equidistant points.
!> The new solution is then simply calculated by taking the mean value of the supersampling points inside of the new subcell.
!> Attention: This procedure is not conservative! We do not take the Jacobian into account.
!==================================================================================================================================
SUBROUTINE SupersampleFVCell(UOld,UNew,NOld,NNew,NSuper)
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: UOld(1:PP_nVar,0:NOld,0:NOld,0:ZDIM(NOld)) !< One FV element on NOld
REAL,INTENT(OUT)   :: UNew(1:PP_nVar,0:NNew,0:NNew,0:ZDIM(NNew)) !< FV Element on NNew
INTEGER,INTENT(IN) :: NOld                                       !< Old polynomial degree
INTEGER,INTENT(IN) :: NNew                                       !< New polynomial degree
INTEGER,INTENT(IN) :: NSuper                                     !< Polynomial degree for supersampling
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: deltaXiOld,deltaXi,deltaXiSuper,xiSuper
REAL                :: U_mean(PP_nVar)
INTEGER             :: iOld,jOld,kOld,iSuper,jSuper,kSuper
INTEGER             :: i,j,k
!==================================================================================================================================
! Supersample the old FV solution (with (NSuper+1)**dim superampling points in each new
! sub cell), then take the mean value
deltaXiOld = 2.0/(REAL(NOld)+1.)       ! Length (in reference space) of a FV element in the old element
deltaXi          = 2.0/(REAL(NNew)+1.)       ! Length (in reference space) of a FV element in the new element
deltaXiSuper     = deltaXi/(REAL(NSuper)+1.) ! Spacing (in reference space) between supersampling points
DO k=0,ZDIM(NNew); DO j=0,NNew; DO i=0,NNew
  U_mean = 0.
  DO kSuper=0,ZDIM(NSuper); DO jSuper=0,NSuper; DO iSuper=0,NSuper
    ! Calculate the index that the current super sampling point has in the old solution
    xiSuper = (REAL(i)*deltaXi) + (REAL(iSuper)+0.5)*deltaXiSuper
    iOld = INT(xiSuper/deltaXiOld)
    xiSuper = (REAL(j)*deltaXi) + (REAL(jSuper)+0.5)*deltaXiSuper
    jOld = INT(xiSuper/deltaXiOld)
#if PP_dim == 3
    xiSuper = (REAL(k)*deltaXi) + (REAL(kSuper)+0.5)*deltaXiSuper
    kOld = INT(xiSuper/deltaXiOld)
#else
    kOld = 0
#endif
    ! Calculate sum for mean value
    U_mean(:) = U_mean(:)+UOld(:,iOld,jOld,kOld)
  END DO; END DO; END DO! iSuper,jSuper,kSuper=0,NOld
  UNew(:,i,j,k) = U_mean(:)/REAL((NSuper+1)**PP_dim)
END DO; END DO; END DO! i,j,k=0,NNew
END SUBROUTINE SupersampleFVCell
#endif

!==================================================================================================================================
!> Finalizes variables necessary for restart subroutines
!==================================================================================================================================
SUBROUTINE FinalizeRestart()
! MODULES
USE MOD_Restart_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
RestartMode       = -1
RestartInitIsDone = .FALSE.

END SUBROUTINE FinalizeRestart

END MODULE MOD_Restart
