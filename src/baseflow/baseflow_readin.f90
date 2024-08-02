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
#include "eos.h"

!==================================================================================================================================
!> Subroutines needed for the general base flow based on a moving time average of the instationary flow field, also known as Pruett
!> damping. See "The temporally filtered Navier–Stokes equations: Properties of the residual stress" for details.
!==================================================================================================================================
MODULE MOD_Baseflow_Readin
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ReadBaseFlow
  MODULE PROCEDURE ReadBaseFlow
END INTERFACE

#if EQNSYSNR == 2 /* NAVIER-STOKES */
INTERFACE ReadBaseFlowRMS
  MODULE PROCEDURE ReadBaseFlowRMS
END INTERFACE
#endif /* NAVIER-STOKES */

PUBLIC :: ReadBaseFlow
#if EQNSYSNR == 2 /* NAVIER-STOKES */
PUBLIC :: ReadBaseflowRMS
#endif /* NAVIER-STOKES */
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> \brief Read baseflow from HDF5 file
!> This routine reads the base flow from a .h5 file. This is used for the pruett damping or for the base
!> flow type three for sponge or for the BCs 204,205.
!> Since the base flow is stored on the whole domain, there are no problems if we e.g. change the shape of the sponge region.
!==================================================================================================================================
SUBROUTINE ReadBaseFlow(FileName)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Baseflow_Vars
USE MOD_ApplyJacobianCons  ,ONLY: ApplyJacobianCons
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_Equation_Vars      ,ONLY: StrVarNames
USE MOD_HDF5_Input         ,ONLY: CloseDataFile,ReadAttribute,ReadArray,DatasetExists,GetDataSize,GetDataProps
USE MOD_HDF5_Output        ,ONLY: WriteBaseFlow
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_IO_HDF5            ,ONLY: OpenDataFile
USE MOD_IO_HDF5            ,ONLY: File_ID,nDims,HSize
USE MOD_Mesh_Vars          ,ONLY: detJac_Ref,Ngeo,MeshFile
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nElems,nGlobalElems
USE MOD_Output_Vars        ,ONLY: ProjectName
USE MOD_StringTools        ,ONLY: STRICMP
#if PP_dim == 3
USE MOD_2D                 ,ONLY: ExpandArrayTo3D
#else
USE MOD_2D                 ,ONLY: to2D_rank5
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName                 !< HDF5 filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nVar_BaseFlow
INTEGER                        :: N_BaseFlow
INTEGER                        :: nElems_BaseFlow
CHARACTER(LEN=255)             :: NodeType_BaseFlow
LOGICAL                        :: InterpolateSolution
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesTmp(:)
LOGICAL                        :: foundField,foundVar(PP_nVar)
REAL,ALLOCATABLE               :: BaseFlow_local(:,:,:,:,:)
REAL,ALLOCATABLE               :: BaseFlow_localNVar(:,:,:,:,:)
REAL,ALLOCATABLE               :: BaseFlow_local2(:,:,:,:,:)
INTEGER                        :: iElem,i,j,k
INTEGER                        :: HSize_proc(5)
REAL,ALLOCATABLE               :: JNR(:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_NBaseFlow_N(:,:)
REAL,ALLOCATABLE               :: Vdm_3Ngeo_NBaseFlow(:,:)
REAL                           :: StartT,EndT
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' READIN MEAN BASEFLOW ...'
#if USE_MPI
StartT=MPI_WTIME()
#else
CALL CPU_TIME(StartT)
#endif

CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

CALL DatasetExists(File_ID, 'DG_Solution', foundField)
IF (.NOT.foundField) CALL Abort(__STAMP__,'There is no mean data in baseflow file')

CALL GetDataProps(nVar_BaseFlow,N_BaseFlow,nElems_BaseFlow,NodeType_BaseFlow)
CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)

! Sanity check, number of elements
IF ((HSize(2).NE.N_BaseFlow+1).OR.(HSize(3).NE.N_BaseFlow+1).OR.(HSize(5).NE.nGlobalElems)) &
  CALL Abort(__STAMP__,"Dimensions of baseflow file do not match!")

HSize_proc    = INT(HSize)
HSize_proc(5) = nElems
! Allocate array to hold the BaseFlow data
ALLOCATE(BaseFlow_local(nVar_BaseFlow,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))

! Read array
CALL ReadArray('DG_Solution',5,HSize_proc,OffsetElem,5,RealArray=BaseFlow_local)

! Check if all arrays are available
ALLOCATE(BaseFlow_localNVar(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
ALLOCATE(VarNamesTmp(nVar_BaseFlow))
! Read variable names
CALL ReadAttribute(File_ID,'VarNames',nVar_BaseFlow,StrArray=VarNamesTmp)
! Read in time from restart file
CALL ReadAttribute(File_ID,'Time',1,RealScalar=BaseFlowTime)

foundVar = .FALSE.
DO i = 1,nVar_BaseFlow
  DO j = 1,PP_nVar
    IF(STRICMP(TRIM(VarNamesTmp(i)),StrVarNames(j))) THEN
      BaseFlow_localNVar(j,:,:,:,:)=BaseFlow_local(i,:,:,:,:)
      foundVar(j)=.TRUE.
    END IF
  END DO ! j = 1,PP_nVar
END DO !  i = 1,nVar_BaseFlow

IF(.NOT. ANY(foundVar)) CALL Abort(__STAMP__,&
  'Not a valid baseflow file, missing conservative mean values of DG solution!')

! Write truncated array back to BaseFlow_local
DEALLOCATE(BaseFlow_local)
ALLOCATE(BaseFlow_local(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
BaseFlow_local = BaseFlow_localNVar
DEALLOCATE(BaseFlow_localNVar)

! Check if we need to interpolate the BaseFlow file to our current polynomial degree and node type
IF((N_BaseFlow.NE.PP_N) .OR. (TRIM(NodeType_BaseFlow).NE.TRIM(NodeType)))THEN
  InterpolateSolution = .TRUE.
  IF(MIN(N_BaseFlow,PP_N).LT.NGeo) &
    CALL PrintWarning('The geometry is or was underresolved and will potentially change on baseflow!')
ELSE
  InterpolateSolution=.FALSE.
END IF

! Read in state
IF(.NOT. InterpolateSolution)THEN
  ! No interpolation needed, read solution directly from file
#if PP_dim == 3
  IF (HSize_proc(4).EQ.1) THEN
    ! FLEXI compiled 3D, but data is 2D => expand third space dimension
    CALL ExpandArrayTo3D(5,(/PP_nVar,PP_N+1,PP_N+1,1,nElems/),4,PP_N+1,BaseFlow_local,BaseFlow)
  ELSE
    ! FLEXI compiled 3D + data 3D
    BaseFlow = BaseFlow_local
  END IF
#else
  IF (HSize_proc(4).EQ.1) THEN
    ! FLEXI compiled 2D + data 2D
    BaseFlow = BaseFlow_local
  ELSE
    ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
    CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVar,PP_N,PP_N,PP_N,nElems/),4,BaseFlow_local)
    BaseFlow = BaseFlow_local
  END IF
#endif
ELSE ! InterpolateSolution
  ! We need to interpolate the solution to the new computational grid
  SWRITE(UNIT_stdOut,'(A,I0,3A,I0,3A)') ' | Interpolating solution from baseflow grid with N=',N_BaseFlow,' (',TRIM(NodeType_BaseFlow), &
                                        ') to computational grid with N='                   ,PP_N     ,' (',TRIM(NodeType),')'

  ALLOCATE(Vdm_NBaseFlow_N(0:PP_N,0:N_BaseFlow))
  ALLOCATE(Vdm_3Ngeo_NBaseFlow(0:N_BaseFlow,0:3*NGeo))
  CALL GetVandermonde(N_BaseFlow, NodeType_BaseFlow,PP_N,      NodeType,         &
                      Vdm_NBaseFlow_N,     modal=.TRUE.)
  CALL GetVandermonde(3*Ngeo,    NodeType,        N_BaseFlow, NodeType_BaseFlow, &
                      Vdm_3Ngeo_NBaseFlow, modal=.TRUE.)

#if PP_dim == 3
  IF (HSize_proc(4).EQ.1) THEN
    ! FLEXI compiled 3D, but data is 2D => expand third space dimension
    ! use temporary array 'BaseFlow_local2' to store 3D data
    ALLOCATE(BaseFlow_local2(PP_nVar,0:N_BaseFlow,0:N_BaseFlow,0:N_BaseFlow,nElems))
    CALL ExpandArrayTo3D(5,HSize_proc,4,N_BaseFlow,BaseFlow_local,BaseFlow_local2)
    ! Reallocate 'BaseFlow_local' to 3D and mv data from BaseFlow_local2 to BaseFlow_local
    DEALLOCATE(BaseFlow_local)
    ALLOCATE(BaseFlow_local(PP_nVar,0:N_BaseFlow,0:N_BaseFlow,0:N_BaseFlow,nElems))
    BaseFlow_local = BaseFlow_local2
    DEALLOCATE(BaseFlow_local2)
  END IF
#else
  IF (HSize_proc(4).NE.1) THEN
    ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
    CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVar,N_BaseFlow,N_BaseFlow,N_BaseFlow,nElems/),4,BaseFlow_local)
  END IF
#endif
  ! Transform solution to refspace and project solution to N
  ! For conservativity deg of detJac should be identical to EFFECTIVE polynomial deg of solution
  ! (e.g. beware when filtering the jacobian )
  IF(N_BaseFlow.GT.PP_N)THEN
    ALLOCATE(JNR(1,0:N_BaseFlow,0:N_BaseFlow,0:N_BaseFlow*(PP_dim-2)))
    DO iElem=1,nElems
      CALL ChangeBasisVolume(1,3*Ngeo,N_BaseFlow,Vdm_3Ngeo_NBaseFlow,detJac_Ref(:,:,:,:,iElem),JNR)
      DO k=0,N_BaseFlow*(PP_dim-2); DO j=0,N_BaseFlow; DO i=0,N_BaseFlow
        BaseFlow_local(:,i,j,k,iElem)=BaseFlow_local(:,i,j,k,iElem)*JNR(1,i,j,k)
      END DO; END DO; END DO
      CALL ChangeBasisVolume(PP_nVar,N_BaseFlow,PP_N,Vdm_NBaseFlow_N,BaseFlow_local(:,:,:,:,iElem),BaseFlow(:,:,:,:,iElem))
    END DO
    DEALLOCATE(JNR)
    ! Transform back
    CALL ApplyJacobianCons(BaseFlow,toPhysical=.TRUE.,FVE=0)
  ELSE
    DO iElem=1,nElems
      CALL ChangeBasisVolume(PP_nVar,N_BaseFlow,PP_N,Vdm_NBaseFlow_N,BaseFlow_local(:,:,:,:,iElem),BaseFlow(:,:,:,:,iElem))
    END DO
  END IF

  DEALLOCATE(BaseFlow_local)
  DEALLOCATE(Vdm_NBaseFlow_N)
  DEALLOCATE(Vdm_3Ngeo_NBaseFlow)
END IF

DEALLOCATE(HSize)
DEALLOCATE(VarNamesTmp)
CALL CloseDataFile()

IF(.NOT.doBaseFlowRMS)THEN
  IF(InterpolateSolution)THEN
    SWRITE(UNIT_stdOut,'(A,A,A)') ' | Writing interpolated baseflow to file "',TRIM(TIMESTAMP(TRIM(ProjectName)//'_interpolated_BaseFlow',BaseFlowTime))//'.h5','"'
    CALL WriteBaseFlow(ProjectName=TRIM(ProjectName)//'_interpolated',MeshFileName=TRIM(MeshFile),OutputTime=BaseFlowTime,FutureTime=BaseFlowTime)
  END IF
  SWRITE(UNIT_stdOut,'(A,ES13.7)')' | Readin of baseflow file successful at t = ',BaseFlowTime
END IF

EndT = FLEXITIME()
SWRITE(UNIT_stdOut,'(A,F0.3,A)')' READIN MEAN BASEFLOW DONE  [',EndT-StartT,'s]'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE ReadBaseFlow


#if EQNSYSNR == 2 /* NAVIER-STOKES */
!==================================================================================================================================
!> Read in the current RMS values from a previous calculation to BaseFlow from
!==================================================================================================================================
SUBROUTINE ReadBaseFlowRMS(FileName)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Baseflow_Vars
USE MOD_ApplyJacobianCons  ,ONLY: ApplyJacobianCons
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_Equation_Vars      ,ONLY: StrVarNamesFluc
USE MOD_HDF5_Input         ,ONLY: CloseDataFile,ReadAttribute,ReadArray,DatasetExists,GetDataSize,GetDataProps
USE MOD_HDF5_Output        ,ONLY: WriteBaseFlow
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_IO_HDF5            ,ONLY: OpenDataFile
USE MOD_IO_HDF5            ,ONLY: File_ID,nDims,HSize
USE MOD_Mesh_Vars          ,ONLY: detJac_Ref,Ngeo,MeshFile
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nElems,nGlobalElems
USE MOD_Output_Vars        ,ONLY: ProjectName
USE MOD_StringTools        ,ONLY: STRICMP
#if PP_dim == 3
USE MOD_2D                 ,ONLY: ExpandArrayTo3D
#else
USE MOD_2D                 ,ONLY: to2D_rank5
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName                 !< HDF5 filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nVar_BaseFlowRMS
INTEGER                        :: N_BaseFlowRMS
INTEGER                        :: nElems_BaseFlowRMS
REAL                           :: BaseFlowRMSTime
CHARACTER(LEN=255)             :: NodeType_BaseFlowRMS
LOGICAL                        :: InterpolateSolution
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesTmp(:)
LOGICAL                        :: foundField
LOGICAL,ALLOCATABLE            :: foundVar(:)
REAL,ALLOCATABLE               :: BaseFlowRMS_local(:,:,:,:,:)
REAL,ALLOCATABLE               :: BaseFlowRMS_localNVar(:,:,:,:,:)
REAL,ALLOCATABLE               :: BaseFlowRMS_local2(:,:,:,:,:)
INTEGER                        :: iElem,i,j,k
INTEGER                        :: HSize_proc(5)
REAL,ALLOCATABLE               :: JNR(:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_NBaseFlowRMS_N(:,:)
REAL,ALLOCATABLE               :: Vdm_3Ngeo_NBaseFlowRMS(:,:)
REAL                           :: StartT,EndT
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' READIN RMS BASEFLOW ...'
#if USE_MPI
StartT=MPI_WTIME()
#else
CALL CPU_TIME(StartT)
#endif

CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

CALL DatasetExists(File_ID, 'Fluc', foundField)
IF (.NOT.foundField) CALL Abort(__STAMP__,'There is no RMS data in baseflow file')

CALL GetDataProps(nVar_BaseFlowRMS,N_BaseFlowRMS,nElems_BaseFlowRMS,NodeType_BaseFlowRMS,'Fluc')
CALL GetDataSize(File_ID,'Fluc',nDims,HSize)

! Sanity check, number of elements
IF ((HSize(2).NE.N_BaseFlowRMS+1).OR.(HSize(3).NE.N_BaseFlowRMS+1).OR.(HSize(5).NE.nGlobalElems)) &
  CALL Abort(__STAMP__,"Dimensions of baseflow file do not match!")

HSize_proc    = INT(HSize)
HSize_proc(5) = nElems
! Allocate array to hold the BaseFlow data
ALLOCATE(BaseFlowRMS_local(nVar_BaseFlowRMS,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))

! Read array
CALL ReadArray('Fluc',5,HSize_proc,OffsetElem,5,RealArray=BaseFlowRMS_local)

! Check if all arrays are available
ALLOCATE(BaseFlowRMS_localNVar(PP_nVarRMS,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
ALLOCATE(VarNamesTmp(nVar_BaseFlowRMS))
! Read variable names
CALL ReadAttribute(File_ID,'VarNames_Fluc',nVar_BaseFlowRMS,StrArray=VarNamesTmp)
! Read in time from restart file
CALL ReadAttribute(File_ID,'Time',1,RealScalar=BaseFlowRMSTime)

ALLOCATE(foundVar(nVar_BaseFlowRMS))
foundVar = .FALSE.
DO i = 1,nVar_BaseFlowRMS
  DO j = 1,PP_nVarRMS
    IF(STRICMP(TRIM(VarNamesTmp(i)),StrVarNamesFluc(j))) THEN
      BaseFlowRMS_localNVar(j,:,:,:,:)=BaseFlowRMS_local(i,:,:,:,:)
      foundVar(j)=.TRUE.
    END IF
  END DO ! j = 1,PP_nVarRMS
END DO !  i = 1,nVal(1)

IF(.NOT. ANY(foundVar)) CALL Abort(__STAMP__,&
  'Not a valid baseflow file, missing Reynolds stress component!')
DEALLOCATE(foundVar)

! Write truncated array back to BaseFlowRMS_local
DEALLOCATE(BaseFlowRMS_local)
ALLOCATE(BaseFlowRMS_local(PP_nVarRMS,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
BaseFlowRMS_local = BaseFlowRMS_localNVar
DEALLOCATE(BaseFlowRMS_localNVar)

! Check if we need to interpolate the BaseFlowRMS file to our current polynomial degree and node type
IF((N_BaseFlowRMS.NE.PP_N) .OR. (TRIM(NodeType_BaseFlowRMS).NE.TRIM(NodeType)))THEN
  InterpolateSolution = .TRUE.
  IF(MIN(N_BaseFlowRMS,PP_N).LT.NGeo) &
    CALL PrintWarning('The geometry is or was underresolved and will potentially change on baseflow!')
ELSE
  InterpolateSolution=.FALSE.
END IF

! Read in state
IF(.NOT. InterpolateSolution)THEN
  ! No interpolation needed, read solution directly from file
#if PP_dim == 3
  IF (HSize_proc(4).EQ.1) THEN
    ! FLEXI compiled 3D, but data is 2D => expand third space dimension
    CALL ExpandArrayTo3D(5,(/PP_nVarRMS,PP_N+1,PP_N+1,1,nElems/),4,PP_N+1,BaseFlowRMS_local,BaseFlowRMS)
  ELSE
    ! FLEXI compiled 3D + data 3D
    BaseFlowRMS = BaseFlowRMS_local
  END IF
#else
  IF (HSize_proc(4).EQ.1) THEN
    ! FLEXI compiled 2D + data 2D
    BaseFlowRMS = BaseFlowRMS_local
  ELSE
    ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
    CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVarRMS,PP_N,PP_N,PP_N,nElems/),4,BaseFlowRMS_local)
    BaseFlowRMS = BaseFlowRMS_local
  END IF
#endif
ELSE ! InterpolateSolution
  ! We need to interpolate the solution to the new computational grid
  SWRITE(UNIT_stdOut,'(A,I0,3A,I0,3A)') ' | Interpolating solution from baseflow grid with N=',N_BaseFlowRMS,' (',TRIM(NodeType_BaseFlowRMS), &
                                        ') to computational grid with N='                   ,PP_N     ,' (',TRIM(NodeType),')'

  ALLOCATE(Vdm_NBaseFlowRMS_N(0:PP_N,0:N_BaseFlowRMS))
  ALLOCATE(Vdm_3Ngeo_NBaseFlowRMS(0:N_BaseFlowRMS,0:3*NGeo))
  CALL GetVandermonde(N_BaseFlowRMS, NodeType_BaseFlowRMS,PP_N,      NodeType,         &
                      Vdm_NBaseFlowRMS_N,     modal=.TRUE.)
  CALL GetVandermonde(3*Ngeo,    NodeType,        N_BaseFlowRMS, NodeType_BaseFlowRMS, &
                      Vdm_3Ngeo_NBaseFlowRMS, modal=.TRUE.)

#if PP_dim == 3
  IF (HSize_proc(4).EQ.1) THEN
    ! FLEXI compiled 3D, but data is 2D => expand third space dimension
    ! use temporary array 'BaseFlowRMS_local2' to store 3D data
    ALLOCATE(BaseFlowRMS_local2(PP_nVarRMS,0:N_BaseFlowRMS,0:N_BaseFlowRMS,0:N_BaseFlowRMS,nElems))
    CALL ExpandArrayTo3D(5,HSize_proc,4,N_BaseFlowRMS,BaseFlowRMS_local,BaseFlowRMS_local2)
    ! Reallocate 'BaseFlowRMS_local' to 3D and mv data from BaseFlowRMS_local2 to BaseFlowRMS_local
    DEALLOCATE(BaseFlowRMS_local)
    ALLOCATE(BaseFlowRMS_local(PP_nVarRMS,0:N_BaseFlowRMS,0:N_BaseFlowRMS,0:N_BaseFlowRMS,nElems))
    BaseFlowRMS_local = BaseFlowRMS_local2
    DEALLOCATE(BaseFlowRMS_local2)
  END IF
#else
  IF (HSize_proc(4).NE.1) THEN
    ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
    CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVarRMS,N_BaseFlowRMS,N_BaseFlowRMS,N_BaseFlowRMS,nElems/),4,BaseFlowRMS_local)
  END IF
#endif
  ! Transform solution to refspace and project solution to N
  ! For conservativity deg of detJac should be identical to EFFECTIVE polynomial deg of solution
  ! (e.g. beware when filtering the jacobian )
  IF(N_BaseFlowRMS.GT.PP_N)THEN
    ALLOCATE(JNR(1,0:N_BaseFlowRMS,0:N_BaseFlowRMS,0:N_BaseFlowRMS*(PP_dim-2)))
    DO iElem=1,nElems
      CALL ChangeBasisVolume(1,3*Ngeo,N_BaseFlowRMS,Vdm_3Ngeo_NBaseFlowRMS,detJac_Ref(:,:,:,:,iElem),JNR)
      DO k=0,N_BaseFlowRMS*(PP_dim-2); DO j=0,N_BaseFlowRMS; DO i=0,N_BaseFlowRMS
        BaseFlowRMS_local(:,i,j,k,iElem)=BaseFlowRMS_local(:,i,j,k,iElem)*JNR(1,i,j,k)
      END DO; END DO; END DO
      CALL ChangeBasisVolume(PP_nVarRMS,N_BaseFlowRMS,PP_N,Vdm_NBaseFlowRMS_N,BaseFlowRMS_local(:,:,:,:,iElem),BaseFlowRMS(:,:,:,:,iElem))
    END DO
    DEALLOCATE(JNR)
    ! Transform back
    CALL ApplyJacobianCons(BaseFlowRMS,toPhysical=.TRUE.,FVE=0)
  ELSE
    DO iElem=1,nElems
      CALL ChangeBasisVolume(PP_nVarRMS,N_BaseFlowRMS,PP_N,Vdm_NBaseFlowRMS_N,BaseFlowRMS_local(:,:,:,:,iElem),BaseFlowRMS(:,:,:,:,iElem))
    END DO
  END IF

  DEALLOCATE(BaseFlowRMS_local)
  DEALLOCATE(Vdm_NBaseFlowRMS_N)
  DEALLOCATE(Vdm_3Ngeo_NBaseFlowRMS)
END IF

DEALLOCATE(HSize)
DEALLOCATE(VarNamesTmp)
CALL CloseDataFile()

EndT = FLEXITIME()
SWRITE(UNIT_stdOut,'(A,ES13.7)')' | Readin of baseflow file successful at t = ',BaseFlowTime
SWRITE(UNIT_stdOut,'(A,F0.3,A)')' READIN RMS BASEFLOW DONE  [',EndT-StartT,'s]'
SWRITE(UNIT_stdOut,'(132("-"))')


IF(InterpolateSolution)THEN
  SWRITE(UNIT_stdOut,'(A,A,A)') ' | Writing interpolated baseflow to file "',TRIM(TIMESTAMP(TRIM(ProjectName)//'_interpolated_BaseFlow',BaseFlowTime))//'.h5','"'
  CALL WriteBaseFlow(ProjectName=TRIM(ProjectName)//'_interpolated',MeshFileName=TRIM(MeshFile),OutputTime=BaseFlowTime,FutureTime=BaseFlowTime)
  SWRITE(UNIT_stdOut,'(132("-"))')
END IF

END SUBROUTINE ReadBaseFlowRMS
#endif /* NAVIER-STOKES */

END MODULE MOD_Baseflow_Readin
