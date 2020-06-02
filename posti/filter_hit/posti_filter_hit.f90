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
!> ??????????
!===================================================================================================================================
PROGRAM posti_filter_hit
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Filter_HIT
USE MOD_Filter_HIT_Vars
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Interpolation_Vars      ,ONLY: NodeType
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,               ONLY: nElems,Elem_xGP,nElems_IJK,Elem_IJK
USE MOD_Mesh_ReadIn,             ONLY: ReadIJKSorting
USE MOD_Output,                  ONLY: DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Output_Vars,             ONLY: ProjectName,NOut
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5,OpenDataFile
USE MOD_HDF5_Input
USE MOD_HDF5_Output,             ONLY: WriteState
USE MOD_Commandline_Arguments
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_ReadInTools
USE MOD_FFT,                     ONLY: InitFFT,FinalizeFFT
USE MOD_FFT,                     ONLY: Interpolate_DG2FFT, Interpolate_FFT2DG
USE MOD_ANALYZE
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE FFTW3
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,iExt,iElem              !some loop indexes
INTEGER                            :: i,j,k
INTEGER                            :: iVar
CHARACTER(LEN=255)                 :: InputStateFile          !dummy for state file for analysis
CHARACTER(LEN=255)                 :: StateName               !state file for analysis
INTEGER                            :: N_HDF5_old              !Polynominal degree
CHARACTER(LEN=255)                 :: NodeType_HDF5_old
CHARACTER(LEN=255)                 :: MeshFile_HDF5_old
LOGICAL                            :: changedMeshFile=.FALSE. !True if mesh between states changed
LOGICAL                            :: changedN       =.FALSE. !True if N between states changes
INTEGER                            :: HSize_proc(5)          !no idea
REAL,ALLOCATABLE                   :: U_local(:,:,:,:,:)     !solution per cell
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()

CALL ParseCommandlineArguments()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
    " ||==========================================||                                                           "
SWRITE(UNIT_stdOut,'(A)') &
    " ||         Spectral filter for HIT          ||                                                           "
SWRITE(UNIT_stdOut,'(A)') &
    " ||==========================================||                                                           "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Define Parameters
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
!=====================================
CALL prms%SetSection("filterHIT")
CALL prms%CreateStringOption(  "MeshFile"       , "Overwrite mesh file saved to state files' userblock.")
CALL prms%CreateIntOption(     "N_Filter"       , "Cutoff filter")
CALL prms%CreateIntOption(     "N_Visu"         , "Polynomial degree to perform DFFT on")

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: init_hit [prm-file]')
END IF
! Parse parameters
CALL prms%read_options(Args(1))
ParameterFile = Args(1)

! Readin Parameters
N_Filter = GETINT('N_Filter')
N_Visu   = GETINT('N_Visu')

! Initialize IO
CALL InitIOHDF5()

! Loop over all files specified on commandline
DO iArg=2,nArgs
  InputStateFile = Args(iArg)

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputStateFile),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')

  ! Get start index of file extension to check if it is a h5 file
  iExt=INDEX(InputStateFile,'.',BACK = .TRUE.)
  IF(InputStateFile(iExt+1:iExt+2) .NE. 'h5') THEN
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid file extension!')
  END IF

  ! Check if state file is a valid state
  validHDF5 = ISVALIDHDF5FILE(InputStateFile)
  IF(.NOT.validHDF5) &
    CALL CollectiveStop(__STAMP__,'ERROR - Restart file not a valid state file.')

  ! Read attributes and solution from state file
  CALL ReadOldStateFile(InputStateFile)

  ! Check if input attributes have changed since last state file
  changedMeshFile = .FALSE.
  changedN        = .FALSE.
  IF(TRIM(MeshFile_HDF5).NE.TRIM(MeshFile_HDF5_old)) changedMeshFile =.TRUE.
  IF(N_HDF5_old.NE.N_HDF5)                           changedN        =.TRUE.

  ! Re-initialize interpolation if N changed
  IF(changedN) THEN
    CALL FinalizeInterpolation()
    CALL InitInterpolation(N_HDF5)
  END IF

  ! Re-initialize mesh if it has changed
  IF(changedMeshFile) THEN
    SWRITE(UNIT_stdOUT,*) "INITIALIZING MESH FROM FILE """,TRIM(MeshFile_HDF5),""""
    CALL FinalizeMesh()
    CALL DefineParametersMesh()
    CALL InitMesh(MeshMode=0,MeshFile_IN=MeshFile_HDF5)
    CALL ReadIJKSorting() ! Read global xyz sorting of structured mesh

    ! Get new number of points for fourier analysis
    N_FFT=(N_Visu+1)*nElems_IJK(1)
    SDEALLOCATE(U_Global)
    ALLOCATE(U_Global(1:nVar_HDF5,1:N_FFT,1:N_FFT,1:N_FFT))
  END IF

  !! 2D currently not supported
  !IF (HSize(4).EQ.1) CALL CollectiveStop(__STAMP__,'ERROR - Tool does not support 2D data')
  !! FV currently not supported
  !IF (ANY(FV_Elems(:).EQ.1)) CALL CollectiveStop(__STAMP__, &
  !                                         'ERROR - Programm cannot handle FV Subcells!')
  !! Number of elements has to be equal in all three dimensions
  !IF(.NOT.((nElems_IJK(1).EQ.nElems_IJK(2)).AND.(nElems_IJK(1).EQ.nElems_IJK(3)))) THEN
  !  CALL ABORT(__STAMP__,'Mesh has not the same amount of elements in xyz!')
  !END IF

  IF(changedMeshFile .OR. changedN) THEN
    SWRITE(UNIT_stdOut,'(A)') 'FFT SETUP'
    CALL FinalizeFFT()
    CALL InitFFT()
    SWRITE(UNIT_stdOut,'(A)') 'FFT SETUP DONE'
    SWRITE(UNIT_StdOut,'(132("-"))')

    ! Allocate new solution array for FFT
    SDEALLOCATE(U_FFT)
    ALLOCATE(U_FFT(1:nVar_HDF5,1:Endw(1),1:Endw(2),1:Endw(3)))
  END IF

  ! Evaluate DG solution at equidistant points
  CALL Interpolate_DG2FFT(U_HDF5,U_Global)

  ! Apply Fourier-Transform on solution from state file
  DO iVar=1,nVar_HDF5
    CALL DFFTW_PLAN_DFT_R2C_3D(plan,N_FFT,N_FFT,N_FFT,U_Global(iVar,:,:,:),U_FFT(iVar,:,:,:),FFTW_ESTIMATE)
    CALL DFFTW_Execute(plan,U_Global(iVar,:,:,:),U_FFT(iVar,:,:,:))
  END DO
  CALL DFFTW_DESTROY_PLAN(plan)

  ! Normalize Data
  U_FFT=U_FFT/REAL(N_FFT**3)

  ! Apply Fourier Filter
  DO k=1,endw(3); DO j=1,endw(2); DO i=1,endw(1)
    IF(localk(4,i,j,k).GT.N_Filter) U_FFT(:,i,j,k) = 0.
  END DO; END DO; END DO

  ! Evaluate Fourier basis at DG interpolation points
  ALLOCATE(U(1:nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,nElems_HDF5))

  ! Apply inverse Fourier-Transform on solution from state file
  DO iVar=1,nVar_HDF5
    CALL DFFTW_PLAN_DFT_C2R_3D(plan,N_FFT,N_FFT,N_FFT,U_Global(iVar,:,:,:),U_FFT(iVar,:,:,:),FFTW_ESTIMATE)
    CALL DFFTW_Execute(plan,U_Global(iVar,:,:,:),U_FFT(iVar,:,:,:))
  END DO
  CALL DFFTW_DESTROY_PLAN(plan)

  ! Interpolate global solution at equidistant points back to DG solution
  CALL Interpolate_FFT2DG(U_Global,U)

  ! Some dirty stuff for output
  MeshFile = MeshFile_HDF5
  ProjectName = "Test23"
  NOut = N_HDF5

  ! Write State-File
  CALL WriteNewStateFile()

  DEALLOCATE(U)
  DEALLOCATE(U_HDF5)

END DO !iArg=1,nArgs

! Finalize everything
CALL FinalizeParameters()
CALL FinalizeInterpolation()
CALL FinalizeMesh()
CALL FinalizeFFT()
#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
CALL FinalizeMPI()
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' filterHIT FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM posti_filter_hit
