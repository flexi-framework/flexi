!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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

!===================================================================================================================================
!> Contains routines to read in the state (with or without gradients) - if we want the gradients, DGTimeDerivative_weakForm is
!> called once.
!===================================================================================================================================
MODULE MOD_Posti_ReadState
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: ReadState
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> This routine will read in the current state from the statefile. Will call one of two routines: ReadStateWithoutGradients if no
!> gradients have to be visualized or calculated and so no DG operator call is necessary, or ReadStateAndGradients if gradients
!> are needed and we need to calculate the DG operator once.
!> If the DG operator has to be called, we need some parameters. Either a seperate parameter file is passed,  then this one will
!> be used, or we try to extract the parameter file from the userblock.
!> If both fails and we need to compute the DG operator, the program will abort.
!> If the DG operator should not be called and no parameter file (seperate or from userblock) can be found,
!> we specify PP_N from the state file as our polynomial degree (later needed by InitInterpolation).
!===================================================================================================================================
SUBROUTINE ReadState(prmfile,statefile)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars   ,ONLY:withDGOperator
USE MOD_ReadInTools ,ONLY:ExtractParameterFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile      !< FLEXI parameter file, used if DG operator is called
CHARACTER(LEN=255),INTENT(IN)    :: statefile    !< HDF5 state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: userblockFound
!===================================================================================================================================
userblockFound = .TRUE. ! Set to true to later test for existing parameters either form userblock or from seperate file
IF (LEN_TRIM(prmfile).EQ.0) THEN ! No separate parameter file has been given
  ! Try to extract parameter file
  prmfile = ".flexi.ini"
  CALL ExtractParameterFile(statefile,prmfile,userblockFound)
  ! Only abort if we need some parameters to call the DG operator
  IF (.NOT.userblockFound.AND.withDGOperator) THEN
    CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(statefile)//"' and no parameter file specified.")
  END IF
END IF

SWRITE(UNIT_stdOut,'(A,L1)') " [ALL] get solution withDGOperator = ", withDGOperator
IF (withDGOperator) THEN
  CALL ReadStateAndGradients(prmfile,statefile)
ELSE
  ! If no parameters have been specified, use PP_N to later initialize the interpolation routines.
  IF (.NOT.userblockFound) THEN
    CALL ReadStateWithoutGradients(prmfile,statefile,NIN=PP_N)
  ELSE
    CALL ReadStateWithoutGradients(prmfile,statefile)
  END IF
END IF

END SUBROUTINE ReadState


!===================================================================================================================================
!> Read a state file via Restart routine and preforms one DGTimeDerivative_weakForm.
!> This fill at least 'U', 'UPrim', and the gradients 'gradUx/y/z'.
!===================================================================================================================================
SUBROUTINE ReadStateAndGradients(prmfile,statefile)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars          ,ONLY: doAvg2D
USE MOD_Baseflow            ,ONLY: DefineParametersBaseflow,InitBaseflow,FinalizeBaseflow
USE MOD_DG                  ,ONLY: InitDG,DGTimeDerivative_weakForm,FinalizeDG
USE MOD_EOS                 ,ONLY: DefineParametersEos
USE MOD_Equation            ,ONLY: DefineParametersEquation,InitEquation,FinalizeEquation
USE MOD_TestCase            ,ONLY: DefineParametersTestcase
USE MOD_Exactfunc           ,ONLY: DefineParametersExactFunc
USE MOD_Filter              ,ONLY: DefineParametersFilter,InitFilter,FinalizeFilter
USE MOD_Interpolation       ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5             ,ONLY: DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Mesh                ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mortar              ,ONLY: InitMortar,FinalizeMortar
USE MOD_MPI                 ,ONLY: DefineParametersMPI
USE MOD_Overintegration     ,ONLY: DefineParametersOverintegration,InitOverintegration,FinalizeOverintegration
USE MOD_ReadInTools         ,ONLY: prms
USE MOD_ReadInTools         ,ONLY: FinalizeParameters
USE MOD_Restart             ,ONLY: DefineParametersRestart,InitRestart,Restart,FinalizeRestart
USE MOD_Restart_Vars        ,ONLY: RestartTime
USE MOD_Visu_Vars           ,ONLY: changedMeshFile,changedWithDGOperator
USE MOD_Visu_Vars           ,ONLY: MeshFile
#if USE_MPI
USE MOD_MPI                 ,ONLY: InitMPIvars,FinalizeMPI
#endif /*USE_MPI*/
#if FV_ENABLED
USE MOD_FV                  ,ONLY: DefineParametersFV,InitFV,FinalizeFV
USE MOD_FV_Basis            ,ONLY: InitFV_Basis,FinalizeFV_Basis
USE MOD_Indicator           ,ONLY: DefineParametersIndicator,InitIndicator,FinalizeIndicator
#endif /*FV_ENABLED*/
#if PARABOLIC
USE MOD_Lifting             ,ONLY: DefineParametersLifting,InitLifting,FinalizeLifting
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN):: prmfile       !< FLEXI parameter file, used if DG operator is called
CHARACTER(LEN=255),INTENT(IN):: statefile     !< HDF5 state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL FinalizeInterpolation()
#if FV_ENABLED
CALL FinalizeIndicator()
CALL FinalizeFV_Basis()
CALL FinalizeFV()
#endif
CALL FinalizeMortar()
CALL FinalizeRestart()
#if USE_MPI
IF (changedMeshFile.OR.changedWithDGOperator) THEN
  CALL FinalizeMPI()
END IF
#endif
CALL FinalizeEquation()
CALL FinalizeDG()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
#if PARABOLIC
CALL FinalizeLifting()
#endif
CALL FinalizeBaseflow()

! read options from parameter file
CALL FinalizeParameters()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
CALL DefineParametersMesh()
CALL DefineParametersFilter()
CALL DefineParametersOverintegration()
#if FV_ENABLED
CALL DefineParametersIndicator()
CALL DefineParametersFV()
#endif
CALL DefineParametersEos()
CALL DefineParametersEquation()
CALL DefineParametersExactFunc()
CALL DefineParametersTestcase()
#if PARABOLIC
CALL DefineParametersLifting()
#endif
CALL DefineParametersBaseflow()
CALL prms%read_options(prmfile)

! Initialization of I/O routines
CALL InitIOHDF5()

CALL InitInterpolation()
#if FV_ENABLED
CALL InitFV_Basis()
#endif
CALL InitMortar()
CALL InitRestart(statefile)

IF (changedMeshFile.OR.changedWithDGOperator) THEN
  CALL FinalizeMesh()
  CALL InitMesh(meshMode=2,MeshFile_IN=MeshFile)
END IF

CALL InitFilter()
CALL InitOverintegration()
#if USE_MPI
IF (changedMeshFile.OR.changedWithDGOperator) THEN
  CALL InitMPIvars()
END IF
#endif
CALL InitEquation()
CALL InitBaseflow()
doAvg2D = .FALSE.

CALL InitDG()
#if FV_ENABLED
CALL InitIndicator()
CALL InitFV()
#endif
#if PARABOLIC
CALL InitLifting()
#endif /*PARABOLIC*/
CALL Restart(doFlushFiles=.FALSE.)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Call DGTimeDerivative_weakForm..."
CALL DGTimeDerivative_weakForm(RestartTime)
SWRITE(UNIT_stdOut,'(A)')             "DONE"

CALL FinalizeParameters()

END SUBROUTINE ReadStateAndGradients


!===================================================================================================================================
!> Read 'U' directly from a state file.
!===================================================================================================================================
SUBROUTINE ReadStateWithoutGradients(prmfile,statefile,Nin)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ApplyJacobian       ,ONLY: ApplyJacobian
USE MOD_ChangeBasisByDim    ,ONLY: ChangeBasisVolume
USE MOD_DG                  ,ONLY: FinalizeDG
USE MOD_DG_Vars             ,ONLY: U
USE MOD_EOS                 ,ONLY: DefineParametersEos,InitEOS,PrimToCons
USE MOD_HDF5_Input          ,ONLY: OpenDataFile,ReadAttribute,ReadArray,CloseDataFile
USE MOD_HDF5_Input          ,ONLY: GetVarnames,GetDataSize,DatasetExists
USE MOD_Interpolation       ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Interpolation       ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars  ,ONLY: NodeType
USE MOD_IO_HDF5             ,ONLY: DefineParametersIO_HDF5,InitIOHDF5
USE MOD_IO_HDF5             ,ONLY: File_ID,nDims,HSize
USE MOD_Mesh                ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars           ,ONLY: nElems,nGlobalElems,offsetElem
USE MOD_Mesh_Vars           ,ONLY: detJac_Ref,Ngeo
USE MOD_MPI                 ,ONLY: DefineParametersMPI
USE MOD_ReadInTools         ,ONLY: prms
USE MOD_ReadInTools         ,ONLY: FinalizeParameters
USE MOD_Restart             ,ONLY: DefineParametersRestart,InitRestart,Restart,FinalizeRestart
USE MOD_Restart_Vars        ,ONLY: N_Restart,RestartTime
USE MOD_Restart_Vars        ,ONLY: NodeType_Restart,InterpolateSolution
USE MOD_Visu_Vars           ,ONLY: changedMeshFile,doSurfVisu,hasFV_Elems
USE MOD_Visu_Vars           ,ONLY: MeshFile,meshMode_old,nVar_State
#if USE_MPI
USE MOD_MPI                 ,ONLY: FinalizeMPI
#endif /*USE_MPI*/
#if FV_ENABLED
USE MOD_FV_Basis            ,ONLY: InitFV_Basis,FinalizeFV_Basis,DefineParametersFV_Basis
USE MOD_FV_Vars             ,ONLY: FV_Elems
USE MOD_Indicator_Vars      ,ONLY: IndValue
USE MOD_Mortar              ,ONLY: InitMortar,FinalizeMortar
USE MOD_Restart             ,ONLY: SupersampleFVCell
USE MOD_Restart_Vars        ,ONLY: NFVRestartSuper
#endif /*FV_ENABLED*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN):: prmfile       !< FLEXI parameter file, used if DG operator is called
CHARACTER(LEN=255),INTENT(IN):: statefile     !< HDF5 state file
INTEGER,INTENT(IN),OPTIONAL  :: Nin           !< Polynomial degree used in InitInterpolation (OPTIONAL)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: meshMode_loc
LOGICAL                        :: changedMeshMode
! Interpolation from different equation systems
INTEGER                        :: iElem,i,j,k
REAL,ALLOCATABLE               :: JNR(    :,:,:,:)
REAL,ALLOCATABLE               :: U_local(:,:,:,:,:)
INTEGER                        :: HSize_proc(5)
REAL                           :: Vdm_NRestart_N(    0:PP_N     ,0:N_Restart)
REAL                           :: Vdm_3Ngeo_NRestart(0:N_Restart,0:3*NGeo)
! Var names
CHARACTER(LEN=255),ALLOCATABLE :: VarNames_tmp(:)
LOGICAl                        :: VarNamesExist
! Timers
REAL                           :: StartT,EndT
!===================================================================================================================================
CALL FinalizeInterpolation()

! Some features require normal vectors or metrics
meshMode_loc = 0 ! Minimal mesh init
IF (doSurfVisu)  meshMode_loc = MAX(meshMode_loc,2)
IF (hasFV_Elems) meshMode_loc = MAX(meshMode_loc,2)

#if FV_ENABLED
! For FV and higher mesh modes the FV basis is needed
IF (meshMode_loc.EQ.2)THEN
  CALL FinalizeFV_Basis()
  CALL FinalizeMortar()
END IF
#endif
CALL FinalizeRestart()
CALL FinalizeDG()

! check if the mesh mode has changed from the last time
changedMeshMode = (meshMode_loc.NE.meshMode_old)

#if USE_MPI
IF ((changedMeshFile).OR.(changedMeshMode)) THEN
  CALL FinalizeMPI()
END IF
#endif

! read options from parameter file
CALL FinalizeParameters()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
CALL DefineParametersMesh()
#if FV_ENABLED
CALL DefineParametersFV_Basis()
#endif
CALL DefineParametersEOS()
CALL prms%read_options(prmfile)

! Initialization of I/O routines
CALL InitIOHDF5()

IF (PRESENT(NIn)) THEN
  CALL InitInterpolation(Nin)
ELSE
  CALL InitInterpolation()
END IF

#if FV_ENABLED
! We need to call the FV basis init to allocate some arrays needed in mesh init
IF (meshMode_loc.EQ.2)THEN
  CALL InitFV_Basis()
  CALL InitMortar()
END IF

IF (.NOT.ALLOCATED(FV_Elems)) THEN
  ALLOCATE(FV_Elems(nElems)) ! holds information if element is DG (0) or FV (1)
  FV_Elems = 0.
END IF
IF (.NOT.ALLOCATED(IndValue)) THEN
  ALLOCATE(IndValue(nElems)) ! holds information if element is DG (0) or FV (1)
  IndValue = 0.
END IF
#endif

CALL InitRestart(statefile)

! Call mesh init if the mesh file changed or we need a different mesh mode
IF ((changedMeshFile).OR.(changedMeshMode)) THEN
  CALL FinalizeMesh()
  CALL InitMesh(meshMode=meshMode_loc,MeshFile_IN=MeshFile)
END IF

! Initialize EOS since some quantities need gas properties like R and kappa
CALL InitEOS()

! save old mesh mode for future comparison
meshMode_old = meshMode_loc

! Allocate the local DG solution: element-based
SDEALLOCATE(U)
ALLOCATE(U(1:nVar_State,0:PP_N,0:PP_N,0:PP_NZ,nElems))

! Restart subroutine assumes PP_nVar. Only call it if the number of variables in the file to be visualized matches
IF (nVar_State.EQ.PP_nVar) THEN
  CALL Restart(doFlushFiles=.FALSE.)
! Otherwise, call our own subroutine
ELSE
  IF(MPIRoot)THEN
    WRITE(UNIT_stdOut,'(132("-"))')
    WRITE(UNIT_stdOut,'(A)') ' PERFORMING READIN...'
    WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='YES')' | Reading field from data file "',TRIM(statefile),'"...'
    GETTIME(StartT)
  END IF

  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  ! Attempt to set the restart time
  CALL DatasetExists(File_ID,'Time',VarNamesExist,attrib=.TRUE.)
  IF (VarNamesExist) CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)

  ! Check if the file is a State or a TimeAvg file
  CALL GetVarNames('VarNames',VarNames_tmp,VarNamesExist)
  ! > State file with at least SOME variables
  IF (VarNamesExist) THEN
    ! Sanity check
    CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
    IF ((HSize(2).NE.PP_N+1).OR.(HSize(3).NE.PP_N+1).OR.(HSize(4).NE.PP_NZ+1).OR.(HSize(5).NE.nGlobalElems)) &
      CALL Abort(__STAMP__,'Dimensions of restart file do not match! Check mesh file and 2D/3D mode!')

    CALL ReadArray('DG_Solution',5,(/nVar_State,PP_N+1,PP_N+1,PP_NZ+1,nElems/),offsetElem,5,RealArray=U)
  ! > Hopefully a timeAvg file with at least SOME variable available
  ELSE
    ! Sanity check
    CALL GetDataSize(File_ID,'Mean',nDims,HSize)
    IF ((HSize(2).NE.PP_N+1).OR.(HSize(3).NE.PP_N+1).OR.(HSize(4).NE.PP_NZ+1).OR.(HSize(5).NE.nGlobalElems)) &
      CALL Abort(__STAMP__,'Dimensions of restart file do not match! Check mesh file and 2D/3D mode!')

    HSize_proc    = INT(HSize)
    HSize_proc(5) = nElems
    ALLOCATE(U_local(nVar_State,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))

    CALL ReadArray('Mean',5,HSize_proc,OffsetElem,5,RealArray=U_local)

    ! Copied from restart.f90 but adapted to conform to any number of variables
    IF (InterpolateSolution) THEN
      ! We need to interpolate the solution to the new computational grid
      SWRITE(UNIT_stdOut,'(A,I0,3A,I0,3A)') ' | Interpolating solution from restart grid with node type "',TRIM(NodeType_Restart), &
                                            '" to computational grid with node type"'                     ,TRIM(NodeType),'"'

      CALL GetVandermonde(N_Restart, NodeType_Restart,PP_N,      NodeType,         &
                          Vdm_NRestart_N,     modal=.TRUE.)
      CALL GetVandermonde(3*Ngeo,    NodeType,        N_Restart, NodeType_Restart, &
                          Vdm_3Ngeo_NRestart, modal=.TRUE.)

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
        CALL ApplyJacobian(nVar_State,U,toPhysical=.TRUE.,FVE=0)
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
    ELSE
      CALL MOVE_ALLOC(U_local,U)
    END IF

    SDEALLOCATE(U_local)
  END IF
  CALL CloseDataFile()

  IF (MPIRoot) THEN
    GETTIME(EndT)
    CALL DisplayMessageAndTime(EndT-StartT, 'PERFORMING READIN DONE!', DisplayLine=.TRUE.)
  END IF
END IF

CALL FinalizeParameters()

END SUBROUTINE ReadStateWithoutGradients

END MODULE MOD_Posti_ReadState
