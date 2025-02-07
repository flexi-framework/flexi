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
!===================================================================================================================================
!> Standalone version of the Visu tool. Read in parameter file, loop over all given State files and call the visu routine for
!> all of them.
!>
!> Usage: posti parameter_posti.ini [parameter_flexi.ini] State1.h5 State2.h5 ...
!> The optional parameter_flexi.ini is used for FLEXI parameters instead of the ones that are found in the userblock of the
!> State file.
!===================================================================================================================================
PROGRAM Posti_Visu
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_Globals_Vars         ,ONLY: StartTime
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_ISO_VARYING_STRING
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_Output_Vars           ,ONLY: ProjectName,doPrintStatusLine
USE MOD_StringTools           ,ONLY: STRICMP,GetFileExtension
USE MOD_Visu                  ,ONLY: visu,FinalizeVisu
USE MOD_Visu_HDF5_Output      ,ONLY: visu_WriteHDF5
USE MOD_Visu_Vars
USE MOD_VTK                   ,ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet
#if USE_MPI
USE MOD_MPI                   ,ONLY: FinalizeMPI
USE MOD_VTK                   ,ONLY: WriteParallelVTK
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iArg,iVar,iExt
CHARACTER(LEN=255),TARGET      :: prmfile
CHARACTER(LEN=255),TARGET      :: postifile
CHARACTER(LEN=255),TARGET      :: statefile
INTEGER                        :: skipArgs
CHARACTER(LEN=255)             :: FileString_DG
CHARACTER(LEN=255)             :: FileString_SurfDG
REAL                           :: percent
INTEGER                        :: tmpLength
CHARACTER(LEN=255)             :: tmpString1,tmpString2
#if FV_ENABLED
CHARACTER(LEN=255)             :: FileString_FV
CHARACTER(LEN=255)             :: FileString_SurfFV
CHARACTER(LEN=255)             :: FileString_multiblock
#endif
#if USE_MPI
LOGICAL                        :: InitMPI_loc=.FALSE.
#endif /*USE_MPI*/
CHARACTER(LEN=255),ALLOCATABLE :: VarNames_loc(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesSurf_loc(:)
!==================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL ParseCommandlineArguments()

IF (doPrintHelp.GT.0) THEN
  prmfile = ''
  CALL visu(prmfile, prmfile, Args(1) &
#if USE_MPI
           ,MPI_COMM_WORLD &
#endif /*USE_MPI*/
           ) !pass first arg (section etc.) instead of statefile
END IF

IF (nArgs.LT.1) &
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti [posti-prm-file [flexi-prm-file]] statefile [statefiles]')

StartTime = FLEXITIME()
prmfile   = ''
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
  skipArgs = 0 ! do not skip arguments. first argument is a h5 file
  !create empty dummy prm file
  postifile = ".posti.ini"
  IF(MPIRoot)THEN
    IF(FILEEXISTS(postifile))THEN
      OPEN(UNIT=31, FILE=postifile, STATUS="old")
      CLOSE(31, STATUS="delete")
    END IF
    OPEN (UNIT=31, FILE=postifile, STATUS="new")
    CLOSE (UNIT=31)
  END IF
ELSE
  CALL CollectiveStop(__STAMP__,&
        'ERROR - Invalid syntax. Please use: posti [posti-prm-file [flexi-prm-file]] statefile [statefiles]')
END IF

DO iArg=1+skipArgs,nArgs
  statefile = TRIM(Args(iArg))

  WRITE(tmpString1,'(I0)') iArg-1-skipArgs
  WRITE(tmpString2,'(I0)') nArgs -skipArgs
  tmpLength = 96 - LEN(TRIM(tmpString1)) - LEN(TRIM(tmpString2))
  percent   = REAL(iArg-1-skipArgs)/REAL(nArgs-skipArgs)*100.
  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I0,A,I0,A)',ADVANCE='NO') ' Processing file ',iArg-skipArgs,' of ',nArgs-skipArgs,' |'
  SWRITE(UNIT_stdOut,'(A,A1,A)'     ,ADVANCE='NO')  REPEAT('=',MAX(CEILING(percent*(tmpLength+1)/100.),0)),'>',&
                                                    REPEAT(' ',(tmpLength+1)-MAX(CEILING(percent*(tmpLength+1)/100.),0))
  SWRITE(UNIT_stdOut,'(A3,F6.2,A3)' ,ADVANCE='YES') '| [',percent,'%] '
  SWRITE(UNIT_stdOut,'(132("-"))')

  ! Enable progress indicator
  doPrintStatusLine = .TRUE.

  CALL visu(prmfile, postifile, statefile &
#if USE_MPI
           ,MPI_COMM_WORLD &
#endif /*USE_MPI*/
           )

  IF (MeshFileMode) THEN
    ! Remove file extension
    iExt = INDEX(MeshFile,'.',BACK = .TRUE.) ! Position of file extension
    FileString_DG=MeshFile(:iExt-1)
    FileString_DG=TRIM(FileString_DG)//'_visu'
  ELSE
#if FV_ENABLED
    FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DG',OutputTime))
    IF (.NOT.MeshFileMode) FileString_FV=TRIM(TIMESTAMP(TRIM(ProjectName)//'_FV',OutputTime))
#else
    FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))
#endif
  END IF

  ALLOCATE(varnames_loc(nVarVisu))
  ALLOCATE(varnamesSurf_loc(nVarSurfVisuAll))
  DO iVar=1,nVarAll
    IF (mapAllVarsToVisuVars(iVar).GT.0) THEN
      VarNames_loc(mapAllVarsToVisuVars(iVar)) = VarnamesAll(iVar)
    END IF
    IF (mapAllVarsToSurfVisuVars(iVar).GT.0) THEN
      VarNamesSurf_loc(mapAllVarsToSurfVisuVars(iVar)) = VarnamesAll(iVar)
    END IF
  END DO

  SELECT CASE(OutputFormat)
    CASE(OUTPUTFORMAT_PARAVIEW)
      IF (Avg2D) THEN
        CALL WriteDataToVTK(nVarVisu,NVisu,nElemsAvg2D_DG,VarNames_loc,CoordsVisu_DG,UVisu_DG,FileString_DG,&
            dim=2,DGFV=0,nValAtLastDimension=.TRUE.,PostiParallel=.TRUE.,HighOrder=HighOrder,OutputDirectory=OutputDirectory)

#if FV_ENABLED
        CALL WriteDataToVTK(nVarVisu,NVisu_FV,nElemsAvg2D_FV,VarNames_loc,CoordsVisu_FV,UVisu_FV,FileString_FV,&
            dim=2,DGFV=1,nValAtLastDimension=.TRUE.,PostiParallel=.TRUE.,HighOrder=HighOrder,OutputDirectory=OutputDirectory)

        IF (MPIRoot) THEN
          ! write multiblock file
          FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))
          CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_DG,FileString_FV,OutputDirectory=OutputDirectory)
        ENDIF
#endif

      ! .NOT.Avg2D
      ELSE
        CALL WriteDataToVTK(nVarVisu,NVisu,nElems_DG,VarNames_loc,CoordsVisu_DG,UVisu_DG,FileString_DG,&
            dim=PP_dim,DGFV=0,nValAtLastDimension=.TRUE.,PostiParallel=.TRUE.,HighOrder=HighOrder,OutputDirectory=OutputDirectory)

#if FV_ENABLED
        IF (.NOT.MeshFileMode) THEN
          CALL WriteDataToVTK(nVarVisu,NVisu_FV,nElems_FV,VarNames_loc,CoordsVisu_FV,UVisu_FV,FileString_FV,&
              dim=PP_dim,DGFV=1,nValAtLastDimension=.TRUE.,PostiParallel=.TRUE.,HighOrder=HighOrder,OutputDirectory=OutputDirectory)

          IF (MPIRoot) THEN
            ! write multiblock file
            FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))
            CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_DG,FileString_FV,OutputDirectory=OutputDirectory)
          ENDIF
        END IF
#endif

        IF (doSurfVisu) THEN
          ! Surface data
#if FV_ENABLED
          FileString_SurfDG = TRIM(TIMESTAMP(TRIM(ProjectName)//'_SurfDG',OutputTime))
#else
          FileString_SurfDG = TRIM(TIMESTAMP(TRIM(ProjectName)//'_Surf',OutputTime))
#endif

          CALL WriteDataToVTK(nVarSurfVisuAll,NVisu,nBCSidesVisu_DG,VarNamesSurf_loc,CoordsSurfVisu_DG,USurfVisu_DG,&
              FileString_SurfDG,dim=PP_dim-1,DGFV=0,nValAtLastDimension=.TRUE.,HighOrder=HighOrder,OutputDirectory=OutputDirectory)
#if FV_ENABLED
          FileString_SurfFV = TRIM(TIMESTAMP(TRIM(ProjectName)//'_SurfFV',OutputTime))

          CALL WriteDataToVTK(nVarSurfVisuAll,NVisu_FV,nBCSidesVisu_FV,VarNamesSurf_loc,CoordsSurfVisu_FV,USurfVisu_FV,&
              FileString_SurfFV,dim=PP_dim-1,DGFV=1,nValAtLastDimension=.TRUE.,HighOrder=HighOrder,OutputDirectory=OutputDirectory)

          IF (MPIRoot) THEN
            ! write multiblock file
            FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_SurfSolution',OutputTime))
            CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_SurfDG,FileString_SurfFV,OutputDirectory=OutputDirectory)
          ENDIF
#endif
        END IF ! doSurfVisu
      END IF ! Avg2D

    CASE(OUTPUTFORMAT_HDF5)
#if USE_MPI
      ! Initialize communicator
      IF (MPI_COMM_NODE.EQ.MPI_COMM_NULL) THEN
        ! Create the worker MPI Comms
        CALL MPI_COMM_DUP(MPI_COMM_FLEXI,MPI_COMM_NODE,iError)
        CALL MPI_COMM_RANK(MPI_COMM_NODE,myLocalRank,iError)
        CALL MPI_COMM_SIZE(MPI_COMM_NODE,nLocalProcs,iError)
        MPILocalRoot=(myLocalRank .EQ. 0)

        ! Split global communicator into small group leaders and the others
        MPI_COMM_LEADERS = MPI_COMM_NULL
        MPI_COMM_WORKERS = MPI_COMM_NULL
        myLeaderRank     = -1
        myWorkerRank     = -1
        IF(myLocalRank.EQ.0)THEN
          CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,0,0,MPI_COMM_LEADERS,iError)
          CALL MPI_COMM_RANK( MPI_COMM_LEADERS,myLeaderRank,iError)
          CALL MPI_COMM_SIZE( MPI_COMM_LEADERS,nLeaderProcs,iError)
          nWorkerProcs = nProcessors-nLeaderProcs
        ELSE
          CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,1,0,MPI_COMM_WORKERS,iError)
          CALL MPI_COMM_RANK( MPI_COMM_WORKERS,myWorkerRank,iError)
          CALL MPI_COMM_SIZE( MPI_COMM_WORKERS,nWorkerProcs,iError)
          nLeaderProcs = nProcessors-nWorkerProcs
        END IF
      END IF

#endif /*USE_MPI*/
      FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.h5'
      CALL visu_WriteHDF5(nVarVisu     = nVarVisu      &
                         ,NVisu        = NVisu         &
                         ,nElems_loc   = nElems_DG     &
                         ,FileString   = FileString_DG &
                         ,MeshFileName = MeshFile      &
                         ,VarNames_loc = VarNames_loc  &
                         ,Coords_DG    = CoordsVisu_DG &
                         ,dim          = 3             &
                         ,UVisu_DG     = UVisu_DG      )
      IF (doSurfVisu) THEN
        FileString_SurfDG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.h5'
        CALL visu_WriteHDF5( nVarVisu     = nVarSurfVisuAll              &
                           , NVisu        = NVisu                        &
                           , nElems_loc   = nBCSidesVisu_DG              &
                           , FileString   = FileString_SurfDG            &
                           , MeshFileName = MeshFile                     &
                           , VarNames_loc = VarNamesSurf_loc             &
                           , Coords_DG2D  = CoordsSurfVisu_DG(:,:,:,0,:) &
                           , dim          = 2                            &
                           , UVisu_DG2D   = USurfVisu_DG(:,:,0,:,:))
      END IF

#if USE_MPI
      IF (InitMPI_loc) THEN
        ! Free MPI communicators
        IF(MPI_COMM_NODE   .NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_NODE   ,iError)
        IF(MPI_COMM_WORKERS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_WORKERS,iError)
        IF(MPI_COMM_LEADERS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_LEADERS,iError)
      END IF
#endif /*USE_MPI*/
  END SELECT ! OutputFormat

  DEALLOCATE(VarNames_loc)
  DEALLOCATE(VarNamesSurf_loc)
END DO

CALL FinalizeVisu()
#if USE_MPI
! For flexilib MPI init/finalize is controlled by main program
CALL FinalizeMPI()
CALL MPI_FINALIZE(iError)
IF (iError.NE.MPI_SUCCESS) STOP 'MPI finalize error'
#endif
END PROGRAM
