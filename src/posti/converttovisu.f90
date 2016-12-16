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

!===================================================================================================================================
!> Contains routines that convert the calculated FV or DG quantities to the visualization grid. There are separate routines
!> to convert the ElemData and FieldData to the visualization grid.
!===================================================================================================================================
MODULE MOD_Posti_ConvertToVisu
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ConvertToVisu_DG
  MODULE PROCEDURE ConvertToVisu_DG
END INTERFACE
PUBLIC:: ConvertToVisu_DG

INTERFACE ConvertToVisu_GenericData
  MODULE PROCEDURE ConvertToVisu_GenericData
END INTERFACE
PUBLIC:: ConvertToVisu_GenericData

#if FV_ENABLED
INTERFACE ConvertToVisu_FV
  MODULE PROCEDURE ConvertToVisu_FV
END INTERFACE
PUBLIC:: ConvertToVisu_FV

#if FV_RECONSTRUCT
INTERFACE ConvertToVisu_FV_Reconstruct
  MODULE PROCEDURE ConvertToVisu_FV_Reconstruct
END INTERFACE
PUBLIC:: ConvertToVisu_FV_Reconstruct
#endif /* FV_RECONSTRUCT */
#endif /* FV_ENABLED */

CONTAINS

!===================================================================================================================================
!> Perform a ChangeBasis of the calculated DG quantities to the visualization grid.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_DG() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,iVar,iVarVisu,iVarCalc
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
!===================================================================================================================================

! compute UVisu_DG 
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVisuPosti,Vdm_N_NVisu,modal=.FALSE.)
! convert DG solution to UVisu_DG
SDEALLOCATE(UVisu_DG)
ALLOCATE(UVisu_DG(0:NVisu,0:NVisu,0:NVisu,nElems_DG,nVarVisuTotal))
DO iVar=1,nVarDep
  IF (mapVisu(iVar).GT.0) THEN
    iVarCalc = mapCalc(iVar) 
    iVarVisu = mapVisu(iVar) 
    DO iElem = 1,nElems_DG
      CALL ChangeBasis3D(PP_N,NVisu,Vdm_N_NVisu,UCalc_DG(:,:,:,iElem,iVarCalc),UVisu_DG(:,:,:,iElem,iVarVisu))
    END DO
  END IF
END DO 
SDEALLOCATE(Vdm_N_NVisu)
END SUBROUTINE ConvertToVisu_DG

#if FV_ENABLED        
!===================================================================================================================================
!> Convert the calculated FV quantities to the visualization grid.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_FV(mapCalc,maskVisu)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars         ,ONLY: nVarDep,VarNamesTotal
USE MOD_Posti_Vars         ,ONLY: mapVisu,UVisu_FV,nElems_FV,UCalc_FV
USE MOD_ReadInTools        ,ONLY: GETINT
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)          :: mapCalc(nVarDep)
INTEGER,INTENT(IN),OPTIONAL :: maskVisu(nVarDep)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVar,i,j,k,iElem
INTEGER            :: iVarVisu,iVarCalc
!===================================================================================================================================
! compute UVisu_FV
DO iVar=1,nVarDep
  iVarVisu = mapVisu(iVar) 
  IF (PRESENT(maskVisu)) iVarVisu = maskVisu(iVar)*iVarVisu
  IF (iVarVisu.GT.0) THEN
    SWRITE(*,*) "    ", TRIM(VarNamesTotal(iVar))
    iVarCalc = mapCalc(iVar) 
    DO iElem = 1,nElems_FV
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        UVisu_FV(i*2:i*2+1, j*2:j*2+1, k*2:k*2+1,iElem,iVarVisu) = UCalc_FV(i,j,k,iElem,iVarCalc)
      END DO; END DO; END DO
    END DO
  END IF
END DO 

END SUBROUTINE ConvertToVisu_FV


#if FV_RECONSTRUCT
!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE ConvertToVisu_FV_Reconstruct()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_ReadInTools        ,ONLY: GETINT
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
USE MOD_FV_Vars            ,ONLY: gradUxi,gradUeta,gradUzeta
USE MOD_FV_Vars            ,ONLY: FV_dx_XI_L,FV_dx_ETA_L,FV_dx_ZETA_L
USE MOD_FV_Vars            ,ONLY: FV_dx_XI_R,FV_dx_ETA_R,FV_dx_ZETA_R
USE MOD_EOS                ,ONLY: PrimToCons
USE MOD_DG_Vars            ,ONLY: UPrim
USE MOD_EOS_Posti          ,ONLY: GetMaskPrim
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar,i,j,k,iElem,iElem_FV
INTEGER             :: iVarCalc
INTEGER             :: nVarPrim,iVarPrim
INTEGER             :: mapUPrim(PP_nVarPrim)
INTEGER             :: mapUCalc(PP_nVarPrim)
INTEGER             :: maskPrim(nVarDep)
!===================================================================================================================================
! Build local maps of maximal size PP_nVarPrim:
! - mapUCalc(1:nVarPrim) = indices of the nVarPrim primitive quantities that should be visualized in the UCalc_FV array
! - mapUPrim(1:nVarPrim) = indices of the nVarPrim primitive quantities in the UPrim array
! Example: 
!   If only velocityX and pressure should be visualized then:
!     nVarPrim = 2 
!     mapUPrim(1) = 2     mapUCalc(1) = index of velocityX in UCalc_FV 
!     mapUPrim(2) = 5     mapUCalc(2) = index of pressure  in UCalc_FV 
nVarPrim = 0
iVarPrim = 0
maskPrim = GetMaskPrim()
DO iVar=1,nVarDep
  IF (maskPrim(iVar).GT.0) THEN
    iVarPrim = iVarPrim + 1
    IF (mapCalc_FV(iVar).GT.0) THEN
      nVarPrim = nVarPrim + 1
      mapUPrim(nVarPrim) = iVarPrim
      mapUCalc(nVarPrim) = mapCalc_FV(iVar)
    END IF
  END IF
END DO
SWRITE(*,*) "  nVarPrim", nVarPrim
SWRITE(*,*) "  mapUPrim", mapUPrim(1:nVarPrim)
SWRITE(*,*) "  mapUCalc", mapUCalc(1:nVarPrim)


DO iElem_FV=1,nElems_FV
  iElem = mapElems_FV(iElem_FV)  
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO iVar=1,nVarPrim
      iVarPrim = mapUPrim(iVar)
      iVarCalc = mapUCalc(iVar)
      UCalc_FV(i*2  ,j*2  ,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2  ,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2+1,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem)  &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2  ,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem)  &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2+1,k*2  ,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2  ,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2+1,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem) &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2+1,k*2+1,iElem_FV,iVarCalc) = UPrim(iVarPrim,i,j,k,iElem)  &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(i,j,k,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,j,k,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
    END DO
  END DO; END DO; END DO
END DO ! iElem_FV
END SUBROUTINE ConvertToVisu_FV_Reconstruct

#endif /* FV_RECONSTRUCT */

#endif /* FV_ENABLED */

!===================================================================================================================================
!> This routine will read all variables that are not conservative or derived quantities and convert the ones that should be 
!> visualized to the visu grid.
!> These variables include the additional data from the ElemData and FieldData datasetes as well as other datasets that are
!> present in the HDF5 file. The variables will be named DATASETNAME:VARIABLENAME if a attribute VarNames_DATASETNAME exist
!> where we can read the variable names from. If this  attribute does not exist, the name will be a generic DATASETNAME:1,2... .
!> For each dataset a new Vandermonde matrix is build to convert from the specific polynomial degree to the visu grid,
!> so the datasets are not limited to one polynomial degree. Either elementwise (2 dimensions) or pointwise (5 dimensions) datasets
!> are allowed.
!> The addtional variables will always be sorted AFTER the conservative or derived quantities.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_GenericData(statefile) 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_IO_HDF5            ,ONLY: HSize
USE MOD_HDF5_Input         ,ONLY: File_ID
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,ReadArray,CloseDataFile,DatasetExists,ReadAttribute,GetDataSize
USE MOD_Mesh_Vars          ,ONLY: nElems,offsetElem
USE MOD_StringTools        ,ONLY: STRICMP,split_string
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),INTENT(IN)  :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVarVisu,iElem_DG,iElem_FV,iElem,iVarDataset,iVar,iVar2
INTEGER                        :: substring_count,nDims,nVal,nSize,stat
CHARACTER(LEN=255)             :: substrings(2),DatasetName,VariableName,DataSetOld
LOGICAL                        :: datasetFound,varnamesExist,varFound,datasetChanged
REAL,ALLOCATABLE               :: ElemData(:,:),FieldData(:,:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_DG_Visu(:,:),Vdm_FV_Visu(:,:)
CHARACTER(LEN=255),ALLOCATABLE :: DatasetVarNames(:)
!===================================================================================================================================
SWRITE(*,*) "Convert generic datasets to Visu grid"
! Open HDF5 file
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

DataSetOld = ''  ! Used to decide if arrays and Vandermonde matrix should be re-allocated

! Loop over all generic variables that should be visualized - sorted after the dependant variables
DO iVar=nVarDep+1,nVarTotal
  ! Check if this variable should be visualized
  IF ((mapVisu(iVar)).GT.0) THEN
    ! The format of the generic data varnames is DATASETNAME:VARIABLENAME - split into DATASETNAME and VARIABLENAME
    CALL split_string(TRIM(VarNamesTotal(iVar)),':',substrings,substring_count)
    ! If we find more than one substring, the variable is additional data
    IF (substring_count.GT.1) THEN
      ! Store dataset and variable name
      DatasetName =  TRIM(substrings(1))
      VariableName = TRIM(substrings(2))

      ! Check if we have a new dataset
      datasetChanged = .NOT.STRICMP(TRIM(DataSetName),TRIM(DataSetOld))

      SWRITE(*,*) "Convert variable ",TRIM(VariableName)," from dataset ", TRIM(DatasetName)

      ! Get metadata if dataset changed
      IF (datasetChanged) THEN
        ! Try to open the dataset
        CALL DatasetExists(File_ID,TRIM(DatasetName),datasetFound)
        ! Abort if the dataset was not found
        IF (.NOT.datasetFound)  THEN
          CALL CloseDataFile()
          CALL Abort(__STAMP__,'Dataset '//TRIM(DatasetName)//' does not exist.')
        END IF
        ! Get dimensions of the dataset and store number of variables as well as size of array
        CALL GetDataSize(File_ID,TRIM(DatasetName),nDims,HSize)
        nVal   = INT(HSize(1))
        nSize  = INT(HSize(2))
        ! Check if an attribute with the variable names exists
        CALL DatasetExists(File_ID,"VarNames_"//TRIM(DatasetName),varnamesExist,attrib=.TRUE.)
        ! If so, read the varnames
        IF (varnamesExist) THEN
          SDEALLOCATE(DataSetVarNames)
          ALLOCATE(DataSetVarNames(nVal))
          CALL ReadAttribute(File_ID,"VarNames_"//TRIM(DatasetName),nVal,StrArray=DataSetVarNames)
        END IF
      END IF

      ! Compute the index of the variable in the dataset
      IF (varnamesExist) THEN
        ! To check which of these variables should be visualized, compare against the varnames from the dataset
        varFound = .FALSE.
        DO iVar2=1,nVal
          IF (STRICMP(TRIM(VariableName),DataSetVarNames(iVar2))) THEN
            ! This variable should be visualized
            iVarDataset = iVar2
            varFound=.TRUE.
            RETURN
          END IF
        END DO !iVar2=1,nVal
        ! Abort if the variable name has not been found in the variable names from the dataset
        IF (.NOT.varFound) THEN
          CALL CloseDataFile()
          CALL Abort(__STAMP__,'Variable '//TRIM(VariableName)//' not found.')
        END IF
      ELSE ! varnamesExist
        ! If no variable name attribute is present in the HDF5 file, the variable is named in the format 1,2,3,...
        READ(VariableName,'(I2)',IOSTAT=stat) iVarDataset
        ! Abort if read has failed
        IF (stat.NE.0) THEN
          CALL CloseDataFile()
          CALL Abort(__STAMP__,'Variable '//TRIM(VariableName)//' could not be converted to integer value.')
        END IF
      END IF ! varnamesExist

      ! Read in the data if we have a new dataset. Also allocate Vandermonde matrix used in conversion to visu grid.
      IF (datasetChanged) THEN
        SELECT CASE(nDims)
        CASE(2) ! Elementwise data
          ! Allocate array and read dataset
          SDEALLOCATE(ElemData)
          ALLOCATE(ElemData(nVal,nElems))
          CALL ReadArray(TRIM(DatasetName),2,(/nVal,nElems/),offsetElem,2,RealArray=ElemData)  
        CASE(5) ! Pointwise data
          ! Allocate array and read dataset
          SDEALLOCATE(FieldData)
          ALLOCATE(FieldData(nVal,nSize,nSize,nSize,nElems))
          CALL ReadArray(TRIM(DatasetName),5,(/nVal,nSize,nSize,nSize,nElems/),offsetElem,5,RealArray=FieldData)  
          ! Get Vandermonde matrix used to convert to the visu grid
          SDEALLOCATE(Vdm_DG_Visu)
          ALLOCATE(Vdm_DG_Visu(0:NVisu,0:nSize-1))
          CALL GetVandermonde(nSize-1,NodeType,NVisu,NodeTypeVisuPosti,Vdm_DG_Visu,modal=.FALSE.)
          SDEALLOCATE(Vdm_FV_Visu)
          ALLOCATE(Vdm_FV_Visu(0:NVisu_FV,0:nSize-1))
          CALL GetVandermonde(nSize-1,NodeType,NVisu_FV,NodeTypeVisuPosti,Vdm_FV_Visu,modal=.FALSE.)
        CASE DEFAULT
          CALL Abort(__STAMP__,'Dataset '//TRIM(DatasetName)//' does not have 2 or 5 dimensions.')
        END SELECT

        ! Store current name of dataset
        DataSetOld = TRIM(DatasetName)
      END IF ! New dataset

      ! Get index of visu array that we should write to
      iVarVisu= mapVisu(iVar)
      ! Convert the generic data to visu grid
      SELECT CASE(nDims)
      CASE(2) ! Elementwise data
        ! Simply write the elementwise data to all visu points
        DO iElem_DG=1,nElems_DG
          iElem = mapElems_DG(iElem_DG)
          UVisu_DG(:,:,:,iElem_DG,iVarVisu) = ElemData(iVarDataset,iElem)
        END DO
        DO iElem_FV=1,nElems_FV
          iElem = mapElems_FV(iElem_FV)
          UVisu_FV(:,:,:,iElem_FV,iVarVisu) = ElemData(iVarDataset,iElem)
        END DO
      CASE(5) ! Pointwise data
        ! Perform changebasis to visu grid
        DO iElem_DG=1,nElems_DG
          iElem = mapElems_DG(iElem_DG)
          CALL ChangeBasis3D(nSize-1,NVisu,Vdm_DG_Visu,FieldData(iVarDataset,:,:,:,iElem),&
                                                    UVisu_DG(:,:,:,iElem_DG,iVarVisu))
        END DO
        DO iElem_FV=1,nElems_FV
          iElem = mapElems_FV(iElem_FV)
          CALL ChangeBasis3D(nSize-1,NVisu,Vdm_DG_Visu,FieldData(iVarDataset,:,:,:,iElem),&
                                                       UVisu_FV(:,:,:,iElem_DG,iVarVisu))
        END DO
      END SELECT

    END IF ! substring_count.GT.1
  END IF ! mapVisu(iVar).GT.0

END DO !iVar=1,

! Close HDF5 file
CALL CloseDataFile()

! Cleanup of allocatable arrays
SDEALLOCATE(ElemData)
SDEALLOCATE(FieldData)
SDEALLOCATE(Vdm_DG_Visu)
SDEALLOCATE(Vdm_FV_Visu)
SDEALLOCATE(DataSetVarNames)

END SUBROUTINE ConvertToVisu_GenericData

END MODULE MOD_Posti_ConvertToVisu
