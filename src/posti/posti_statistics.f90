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
!> Module containing the main procedures for the POSTI tool: visu3d_requestInformation is called by ParaView to create a
!> list of available variables and visu3D is the main routine of POSTI.
!===================================================================================================================================
MODULE MOD_Statistics
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE postistat
  MODULE PROCEDURE postistat
END INTERFACE

PUBLIC:: postistat

CONTAINS

!===================================================================================================================================
!> Routine using the visu3D data to perform a statistic evaluation 
!===================================================================================================================================
SUBROUTINE postistat(doINIT,Values_p,Values_out,NVisu,Nvisu_k,nElems,nVal)
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL,INTENT(IN)          :: doINIT                  !< flag for initialization
REAL,POINTER,INTENT(IN)     :: Values_p(:,:,:,:,:)     !< Statevector
!REAL,POINTER,INTENT(INOUT)  :: Values_p_out(:,:,:,:,:) !< Statevector
REAL,ALLOCATABLE,INTENT(INOUT)  :: Values_out(:,:,:,:,:) !< Statevector
INTEGER,INTENT(IN)          :: NVisu,Nvisu_k,nElems,nVal
INTEGER,ALLOCATABLE,SAVE    :: range(:,:)
INTEGER:: test(3)
! LOCAL VARIABLES
!===================================================================================================================================

! Initialize in the first run
IF(doINIT) THEN
  ALLOCATE(Values_out(0:Nvisu,0:Nvisu,0:NVisu_k,nElems,2*nVal))
  Values_out=0.
END IF
Values_out(:,:,:,:,1:nVal)=Values_out(:,:,:,:,1:nVal)+Values_p                   ! mean
Values_out(:,:,:,:,1+nVal:2*nVal)=Values_out(:,:,:,:,1+nVal:2*nVal)+Values_p**2  ! mean square
END SUBROUTINE postistat


END MODULE MOD_Statistics


!===================================================================================================================================
!> Posti statistics tool to generate statistic evaluation of derived quantities. Read in parameter file, loop over all given State
!> files and call the visu3D routine to calculate instantaneous values, from which mean, variance etc are evaluated
!>
!> Usage: posti parameter_posti.ini [parameter_flexi.ini] State1.h5 State2.h5 ...
!> The optional parameter_flexi.ini is used for FLEXI parameters instead of the ones that are found in the userblock of the 
!> State file.
!===================================================================================================================================
PROGRAM Posti_statistics
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_Commandline_Arguments
USE MOD_Visu3D                ,ONLY: visu3D
USE MOD_Statistics      
USE MOD_ISO_VARYING_STRING
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_VTK                   ,ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet
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
REAL,ALLOCATABLE,TARGET        :: Values_stat(:,:,:,:,:)
REAL,POINTER                   :: Values_stat_p(:,:,:,:,:)
CHARACTER(LEN=255)             :: FileString_DG
CHARACTER(LEN=255)             :: FileString_FV
CHARACTER(LEN=255)             :: FileString_multiblock
#if !USE_MPI
INTEGER                        :: MPI_COMM_WORLD = 0
#endif
INTEGER                        :: NVisu_k,NVisu_k_FV
INTEGER                        :: nstates
LOGICAL                        :: firstState
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
ELSE
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti [posti-prm-file [flexi-prm-file]] statefile [statefiles]')
END IF

firstState=.TRUE.
DO iArg=1+skipArgs,nArgs
  statefile = TRIM(Args(iArg))
  SWRITE(*,*) "Processing state-file: ",TRIM(statefile)
  
  CALL visu3D(MPI_COMM_WORLD, prmfile, postifile, statefile, &
      coordsDG_out,valuesDG_out,nodeidsDG_out, &
      coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)
  nVal = varnames_out%len/255
  CALL C_F_POINTER(varnames_out%data, VarNames_p, [nVal])
  NVisu_k = NVisu * (VisuDimension-2)
  CALL C_F_POINTER(coordsDG_out%data, Coords_p, [3,NVisu+1,NVisu+1,NVisu_k+1,nElems_DG])
  CALL C_F_POINTER(valuesDG_out%data, Values_p, [  NVisu+1,NVisu+1,NVisu_k+1,nElems_DG,nVal])
  ! Call statistical evaluation / sum here
  CALL postistat(firstState,Values_p,Values_stat,NVisu,Nvisu_k,nElems_DG,nVal)
  firstState=.FALSE.
END DO

nstates=nArgs-skipArgs
Values_stat=Values_stat/nstates
Values_stat(:,:,:,:,1+nVal:2*nVal)=SQRT(Values_stat(:,:,:,:,1+nVal:2*nVal)-Values_stat(:,:,:,:,1:nVal)**2)

FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Mean',OutputTime))//'.vtu'
Values_stat_p=>Values_stat(:,:,:,:,1:nVal)
CALL WriteDataToVTK(nVal,NVisu   ,nElems_DG,VarNames_p,Coords_p,Values_stat_p,FileString_DG,&
    dim=VisuDimension,DGFV=0,nValAtLastDimension=.TRUE.)

FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_MeanSq',OutputTime))//'.vtu'
Values_stat_p=>Values_stat(:,:,:,:,1+nVal:2*nVal)
CALL WriteDataToVTK(nVal,NVisu   ,nElems_DG,VarNames_p,Coords_p,Values_stat_p,FileString_DG,&
    dim=VisuDimension,DGFV=0,nValAtLastDimension=.TRUE.)
DEALLOCATE(Values_stat)

END PROGRAM 

