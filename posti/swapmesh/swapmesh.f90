!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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

MODULE MOD_SwapMesh
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitSwapmesh
  MODULE PROCEDURE InitSwapmesh
END INTERFACE

INTERFACE ReadOldStateFile
  MODULE PROCEDURE ReadOldStateFile
END INTERFACE

INTERFACE WriteNewStateFile
  MODULE PROCEDURE WriteNewStateFile
END INTERFACE

INTERFACE FinalizeSwapmesh
  MODULE PROCEDURE FinalizeSwapmesh
END INTERFACE

PUBLIC:: InitSwapmesh,ReadOldStateFile,WriteNewStateFile,FinalizeSwapmesh

CONTAINS

!===================================================================================================================================
!> Read in user defined parameters and prepare data for swapmesh.
!> The old and new mesh will be read and stored, the necessary Vandermonde matrizes are build and the parametric coordinates
!> of the new gauss points in the old mesh are found.
!===================================================================================================================================
SUBROUTINE InitSwapmesh()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_SwapMesh_Vars
USE MOD_ReadInTools
USE MOD_Commandline_Arguments
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension,INTTOSTR
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Interpolation,           ONLY: GetVandermonde
USE MOD_Interpolation_Vars,      ONLY: NodeTypeVISU,NodeTypeCL
USE MOD_ChangeBasis,             ONLY: ChangeBasis3D
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute
USE MOD_IO_HDF5,                 ONLY: File_ID
USE MOD_SMParametricCoordinates, ONLY: GetParametricCoordinates
USE MOD_Interpolation,           ONLY: InitInterpolation
USE MOD_Output_Vars,             ONLY: NOut,ProjectName
USE MOD_Mesh_Vars,               ONLY: nElems,OffsetElem,nGlobalElems
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
CHARACTER(LEN=255)  :: MeshFile_state
INTEGER             :: nElems_State
REAL                :: Time
CHARACTER(LEN=255)  :: tmp
!===================================================================================================================================
! Open the first statefile to read necessary attributes
CALL OpenDataFile(Args(2),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=MeshFile_state)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL GetDataProps(nVar_State,NState,nElems_State,NodeTypeState)
CALL CloseDataFile()

! change projectname to not re-write the same state
ProjectName = TRIM(ProjectName)//'_newMesh'

! Read user-defined parameters
MeshFileOld = GETSTR('MeshFileOld',TRIM(MeshFile_state))

! New mesh file, the state will be interpolated to this one
MeshFileNew = GETSTR('MeshFileNew')

! Curved meshes or not
useCurvedsOld = GETLOGICAL("useCurvedsOld")
useCurvedsNew = GETLOGICAL("useCurvedsNew")

! Read in polynomial degrees used for interpolation, supersampling and for the new state
! If they have not been specified, use the polynomial degree of the old state
NInter      = GETINT('NInter',INTTOSTR(NState))
NSuper      = GETINT('NSuper',INTTOSTR(NState))
NNew        = GETINT('NNew',  INTTOSTR(NState))
! Set NOut for output state routine
NOut        = NNew

printTroublemakers = GETLOGICAL('printTroublemakers','.TRUE.')

! Tolerance used in coarse element search - standard is 5% tolerance
maxTol = 1. + GETREAL('maxTolerance','5.e-2')

IF (CountOption('RefState').GE.1) THEN
  ! If a reference state is given in the parameter file, it will be used whenever a gauss point can not be found.
  ! This means the program will never abort, so the abort tolerance will be set to a huge number!
  RefState = GETREALARRAY('RefState',nVar_State)
  abortTol = HUGE(1.)
ELSE
  ! Set the abort tolerance, standard is the same as coarse search tolerance
  WRITE(tmp,*) maxTol
  abortTol     =GETREAL('abortTolerance',TRIM(tmp))
END IF

IF (CountOption('displacement').GE.1) THEN
  ! Possibility to set a displacement between the old and new mesh by defining a vector
  displacement = GETREALARRAY('displacement',3)
END IF

! Init interpolation on new polynomial degree
CALL InitInterpolation(NNew)

! Initialize the old mesh, store the mesh coordinates (transformed to CL points) and the number of elements as well as the old NGeo
Time=FLEXITIME()
SWRITE(UNIT_stdOut,'(A)') ' INIT OLD MESH ...'
CALL ReadMeshCoords(MeshFileOld,useCurvedsOld,NGeoOld,nElemsOld,xCLOld)
SWRITE(UNIT_stdOut,*)'done in ',FLEXITIME()-Time

! Translate the old mesh along the displacement vector if needed
IF (CountOption('displacement').GE.1) THEN
  DO i=1,3
    xCLOld(i,:,:,:,:) = xCLOld(i,:,:,:,:)+displacement(i)
  END DO
END IF

! Initialize new mesh
Time=FLEXITIME()
SWRITE(UNIT_stdOut,'(A)') ' INIT NEW MESH ...'
CALL ReadMeshCoords(MeshFileNew,useCurvedsNew,NGeoNew,nElemsNew,xCLNew)
SWRITE(UNIT_stdOut,*)'done in ',FLEXITIME()-Time

! Set offset elem and local and global number of elements in mesh vars (later needed for output routine)
nGlobalElems = nElemsNew
nElems       = nElemsNew
OffsetElem   = 0

! Allocate array for new mesh on interpolation points
ALLOCATE(xCLInter(3,0:NInter,0:NInter,0:NInter,nElemsNew))

! Prepare the necessary Vandermonde matrizes
CALL prepareVandermonde()

! Evaluate parametric coordinates
SWRITE(UNIT_stdOut,'(A)') ' EVALUATING PARAMETRIC COORDINATES ...'
ALLOCATE(xiInter(3,  0:NInter,0:NInter,0:NInter,nElemsNew))
ALLOCATE(InterToElem(0:NInter,0:NInter,0:NInter,nElemsNew))
ALLOCATE(IPDone(     0:NInter,0:NInter,0:NInter,nElemsNew))
ALLOCATE(equalElem(nElemsNew))

! Find the parametric coordinates of the new gauss points in the old mesh
CALL GetParametricCoordinates()

! Allocate aray for new and the old solution - for the new solution, we use the U array from DG vars
! to be able to later use the WriteState routine from FLEXI
ALLOCATE(UOld(nVar_State,0:NState,0:NState,0:NState,nElemsOld))
ALLOCATE(U   (nVar_State,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
END SUBROUTINE InitSwapmesh

!===================================================================================================================================
!> This routine will read in the specified mesh file, convert the equidistant mesh coordinates to CL points and return them.
!> Additionally the number of elements in the mesh as well as NGeo will be returned.
!> The user can specify if curved meshes should be used or not.
!===================================================================================================================================
SUBROUTINE ReadMeshCoords(MeshFile,useCurveds,NGeo,nElems,XCL)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_HDF5_Input
USE MOD_Interpolation,         ONLY: GetVandermonde
USE MOD_Interpolation_Vars,    ONLY: NodeTypeVISU,NodeTypeCL
USE MOD_IO_HDF5,               ONLY: File_ID
USE MOD_ChangeBasis,           ONLY: ChangeBasis3D
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),INTENT(IN)  :: MeshFile
LOGICAL,INTENT(IN)             :: useCurveds
REAL,ALLOCATABLE,INTENT(OUT)   :: XCL(:,:,:,:,:)
INTEGER,INTENT(OUT)            :: NGeo
INTEGER,INTENT(OUT)            :: nElems
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: NodeCoords(:,:,:,:,:)
REAL,ALLOCATABLE               :: NodeCoordsTmp(:,:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_EQNgeo_CLNgeo(:,:)
INTEGER                        :: iElem
!===================================================================================================================================
! Open the mesh file
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
! Read NGeo from mesh file
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo)

! Get the number of elements in the mesh file by reading the size of the ElemInfo array
CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
IF(HSize(1).NE.6) THEN
  CALL abort(__STAMP__,&
    'ERROR: Wrong size of ElemInfo, should be 6')
END IF
CHECKSAFEINT(HSize(2),4)
nElems=INT(HSize(2),4)
DEALLOCATE(HSize)

! Now read in the equidistant coordinates of the mesh nodes
IF(useCurveds)THEN
  ALLOCATE(NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),0,2,RealArray=NodeCoords)
ELSE
  ! If the mesh should be interpreted as not curved, store only the corner nodes
  ALLOCATE(NodeCoords(   3,0:1,   0:1,   0:1,   nElems))
  ALLOCATE(NodeCoordsTmp(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),0,2,RealArray=NodeCoordsTmp)
  NodeCoords(:,0,0,0,:)=NodeCoordsTmp(:,0,   0,   0,   :)
  NodeCoords(:,1,0,0,:)=NodeCoordsTmp(:,NGeo,0,   0,   :)
  NodeCoords(:,0,1,0,:)=NodeCoordsTmp(:,0,   NGeo,0,   :)
  NodeCoords(:,1,1,0,:)=NodeCoordsTmp(:,NGeo,NGeo,0,   :)
  NodeCoords(:,0,0,1,:)=NodeCoordsTmp(:,0,   0,   NGeo,:)
  NodeCoords(:,1,0,1,:)=NodeCoordsTmp(:,NGeo,0,   NGeo,:)
  NodeCoords(:,0,1,1,:)=NodeCoordsTmp(:,0,   NGeo,NGeo,:)
  NodeCoords(:,1,1,1,:)=NodeCoordsTmp(:,NGeo,NGeo,NGeo,:)
  DEALLOCATE(NodeCoordsTmp)
  NGeo=1
END IF

! Convert the equidistant mesh nodes to CL points
ALLOCATE(Vdm_EQNgeo_CLNgeo(0:Ngeo,0:Ngeo))
CALL GetVandermonde(Ngeo,NodeTypeVISU,Ngeo,NodeTypeCL,Vdm_EQNgeo_CLNgeo,modal=.FALSE.)
ALLOCATE(xCL(3,0:NGeo,0:NGeo,0:NGeo,nElems))
DO iElem=1,nElems
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem),xCL(:,:,:,:,iElem))
END DO ! iElem
DEALLOCATE(Vdm_EQNgeo_CLNgeo,NodeCoords)

CALL CloseDataFile()
END SUBROUTINE ReadMeshCoords

!=================================================================================================================================
!> Prepare Vandemonde matrizes used in swapmesh. This includes:
!>   * Vandemonde CL    Ngeo   --> Equidistant NSuper
!>   * Vandemonde CL    NInter --> Gauss       NNew
!>   * Vandemonde Gauss NState --> Gauss       NNew
!>   * Vandemonde CL    Ngeo   --> CL          NInter
!=================================================================================================================================
SUBROUTINE prepareVandermonde()
! MODULES
USE MOD_Swapmesh_Vars,       ONLY: NInter,NNew,NState,NState,NGeoOld,NGeoNew,NSuper
USE MOD_Swapmesh_Vars,       ONLY: NodeTypeState,nElemsNew,xCLInter
USE MOD_Swapmesh_Vars,       ONLY: Vdm_CLNGeo_EquiNSuper,Vdm_CLNInter_GPNNew,Vdm_GPNState_GPNNew
USE MOD_Swapmesh_Vars,       ONLY: xCLNew
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_Interpolation_Vars,  ONLY: NodeType,NodeTypeCL
USE MOD_ChangeBasis,         ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElemNew
REAL               :: Vdm_CLNGeo_CLNInter(0:NInter,0:NGeoNew)
!=================================================================================================================================
! Build visu Vandermonde
ALLOCATE(Vdm_CLNGeo_EquiNSuper(0:NSuper,0:NGeoOld))
CALL GetVandermonde(NGeoOld,NodeTypeState,NSuper,NodeType,Vdm_CLNGeo_EquiNSuper)

! Vandermonde from interpolation CL to new solution G/GL
ALLOCATE(Vdm_CLNInter_GPNNew(0:NNew,0:NInter))
CALL GetVandermonde(NInter,NodeTypeCL,NNew,NodeTypeState,Vdm_CLNInter_GPNNew)

! Vandermonde for direct interpolation in equal elements
IF(NNew.NE.NState)THEN
  ALLOCATE(Vdm_GPNState_GPNNew(0:NNew,0:NState))
  CALL GetVandermonde(NState,NodeTypeState,NNew,NodeTypeState,Vdm_GPNState_GPNNew)
END IF

IF(NGeoNew.NE.NInter)THEN
  CALL GetVandermonde(NGeoNew,NodeTypeCL,NInter,NodeTypeCL,Vdm_CLNGeo_CLNInter)
  DO iElemNew=1,nElemsNew
    CALL ChangeBasis3D(3,NGeoNew,NInter,Vdm_CLNGeo_CLNInter,xCLNew(:,:,:,:,iElemNew),xCLInter(:,:,:,:,iElemNew))
  END DO
ELSE
  xCLInter = xCLNew
END IF

END SUBROUTINE prepareVandermonde

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadOldStateFile(StateFile)
! MODULES                                                                                                                          !
USE MOD_HDF5_Input,    ONLY: OpenDataFile,CloseDataFile,ReadArray,ReadAttribute
USE MOD_IO_HDF5,       ONLY: File_ID
USE MOD_Swapmesh_Vars, ONLY: nVar_State,NState,nElemsOld,Time_State,UOld
USE MOD_ReadInTools,   ONLY: ExtractParameterFile
USE MOD_Output_Vars,   ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Output,        ONLY: insert_userblock
USE ISO_C_BINDING,     ONLY: C_NULL_CHAR
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),INTENT(IN)      :: StateFile !< State file to be read
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: userblockFound
CHARACTER(LEN=255)               :: prmfile=".parameter.ini"
!===================================================================================================================================
! Open the data file
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Read the DG solution and store in UNew
CALL ReadArray('DG_Solution',5,&
               (/nVar_State,NState+1,NState+1,NState+1,nElemsOld/),0,5,RealArray=UOld)

! Read the current time 
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_State)

! Extract parameter file from userblock (if found)
CALL ExtractParameterFile(StateFile,TRIM(prmfile),userblockFound)
! prepare userblock file
CALL insert_userblock(TRIM(UserBlockTmpFile)//C_NULL_CHAR,TRIM(prmfile)//C_NULL_CHAR)
INQUIRE(FILE=TRIM(UserBlockTmpFile),SIZE=userblock_total_len)

! Close the data file
CALL CloseDataFile()
END SUBROUTINE ReadOldStateFile

!===================================================================================================================================
!> Write the new state file by calling the WriteState routine from FLEXI. All necessary variables must have been set correctly!
!===================================================================================================================================
SUBROUTINE WriteNewStateFile()
! MODULES                                                                                                                          !
USE MOD_HDF5_Output,        ONLY: WriteState
USE MOD_Swapmesh_Vars,      ONLY: Time_State,MeshFileNew
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL WriteState(TRIM(MeshFileNew),Time_State,Time_State,isErrorFile=.FALSE.)
END SUBROUTINE WriteNewStateFile

!===================================================================================================================================
!> Finalize swapmesh variables
!===================================================================================================================================
SUBROUTINE FinalizeSwapmesh()
! MODULES                                                                                                                          !
USE MOD_Swapmesh_Vars
USE MOD_DG_Vars,         ONLY: U
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Vdm_CLNGeo_EquiNSuper)
SDEALLOCATE(Vdm_CLNInter_GPNNew)
SDEALLOCATE(Vdm_GPNState_GPNNew)
SDEALLOCATE(xCLInter)
SDEALLOCATE(xCLOld)
SDEALLOCATE(xCLNew)
SDEALLOCATE(xiInter)
SDEALLOCATE(InterToElem)
SDEALLOCATE(equalElem)
SDEALLOCATE(IPDone)
SDEALLOCATE(UOld)
SDEALLOCATE(U)

END SUBROUTINE FinalizeSwapmesh

END MODULE MOD_SwapMesh
