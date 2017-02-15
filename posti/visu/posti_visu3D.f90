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
!> Standalone version of the Visu3D tool. Read in parameter file, loop over all given State files and call the visu3D routine for
!> all of them.
!>
!> Usage: posti parameter_posti.ini [parameter_flexi.ini] State1.h5 State2.h5 ...
!> The optional parameter_flexi.ini is used for FLEXI parameters instead of the ones that are found in the userblock of the 
!> State file.
!===================================================================================================================================
PROGRAM Posti_Visu3D
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_Commandline_Arguments
USE MOD_Visu3D
USE MOD_ISO_VARYING_STRING
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_VTK                   ,ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet,CARRAY
USE MOD_Output_Vars           ,ONLY: ProjectName
USE MOD_StringTools           ,ONLY: STRICMP,GetFileExtension
impliCIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iArg
CHARACTER(LEN=255),TARGET      :: prmfile
CHARACTER(LEN=255),TARGET      :: postifile
CHARACTER(LEN=255),TARGET      :: statefile
TYPE (CARRAY)                  :: coordsDG_out
TYPE (CARRAY)                  :: valuesDG_out
TYPE (CARRAY)                  :: nodeidsDG_out
TYPE (CARRAY)                  :: coordsFV_out
TYPE (CARRAY)                  :: valuesFV_out
TYPE (CARRAY)                  :: nodeidsFV_out
TYPE (CARRAY)                  :: varnames_out
TYPE (CARRAY)                  :: components_out

INTEGER                        :: skipArgs
INTEGER                        :: nVal
CHARACTER(LEN=255),POINTER     :: VarNames_p(:)
REAL,POINTER                   :: Coords_p(:,:,:,:,:)
REAL,POINTER                   :: Values_p(:,:,:,:,:)
CHARACTER(LEN=255)             :: FileString_DG
#if FV_ENABLED                            
CHARACTER(LEN=255)             :: FileString_FV
CHARACTER(LEN=255)             :: FileString_multiblock
INTEGER                        :: NVisu_k_FV
#endif
INTEGER                        :: NVisu_k
#if !USE_MPI
INTEGER                        :: MPI_COMM_WORLD = 0
#endif
!==================================================================================================================================
CALL InitMPI()
CALL ParseCommandlineArguments()
IF (nArgs.LT.1) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti [posti-prm-file [flexi-prm-file]] statefile [statefiles]')
END IF

prmfile = ""
! check if parameter file is given
IF(STRICMP(GetFileExtension(Args(1)),'ini')) THEN
  skipArgs = 1 ! first argument is the parameter file
  postifile = Args(1)
  ! check if a second parameter file is given (this is used instead of the parameter file stored in the userblock of a state file)
  IF (nArgs.GT.2) THEN
    IF (STRICMP(GetFileExtension(Args(2)),'ini')) THEN
      prmfile = Args(2)
      skipArgs = 2
    END IF
  END IF
ELSE IF(STRICMP(GetFileExtension(Args(1)),'h5')) THEN
  skipArgs = 0 ! do not skip a argument. first argument is a h5 file
  postifile = ""
ELSE
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti [posti-prm-file [flexi-prm-file]] statefile [statefiles]')
END IF

DO iArg=1+skipArgs,nArgs
  statefile = TRIM(Args(iArg))
  SWRITE(*,*) "Processing state-file: ",TRIM(statefile)
  
  CALL visu3D(MPI_COMM_WORLD, prmfile, postifile, statefile, &
      coordsDG_out,valuesDG_out,nodeidsDG_out, &
      coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)

#if FV_ENABLED                            
  FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DG',OutputTime))//'.vtu'
#else
  FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
#endif
  nVal = varnames_out%len/255
  CALL C_F_POINTER(varnames_out%data, VarNames_p, [nVal])
  NVisu_k = NVisu * (VisuDimension-2)
  CALL C_F_POINTER(coordsDG_out%data, Coords_p, [3,NVisu+1,NVisu+1,NVisu_k+1,nElems_DG])
  CALL C_F_POINTER(valuesDG_out%data, Values_p, [  NVisu+1,NVisu+1,NVisu_k+1,nElems_DG,nVal])
  CALL WriteDataToVTK(nVal,NVisu   ,nElems_DG,VarNames_p,Coords_p,Values_p,FileString_DG,&
      dim=VisuDimension,DGFV=0,nValAtLastDimension=.TRUE.)
#if FV_ENABLED                            
  FileString_FV=TRIM(TIMESTAMP(TRIM(ProjectName)//'_FV',OutputTime))//'.vtu'
  NVisu_k_FV = NVisu_FV * (VisuDimension-2)
  CALL C_F_POINTER(coordsFV_out%data, Coords_p, [3,NVisu_FV+1,NVisu_FV+1,NVisu_k_FV+1,nElems_FV])
  CALL C_F_POINTER(valuesFV_out%data, Values_p, [  NVisu_FV+1,NVisu_FV+1,NVisu_k_FV+1,nElems_FV,nVal])
  CALL WriteDataToVTK(nVal,NVisu_FV,nElems_FV,VarNames_p,Coords_p,Values_p,FileString_FV,&
      dim=VisuDimension,DGFV=1,nValAtLastDimension=.TRUE.)

  IF (MPIRoot) THEN                   
    ! write multiblock file
    FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtm'
    CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_DG,FileString_FV)
  ENDIF
#endif
END DO

END PROGRAM 

