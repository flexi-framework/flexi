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
!> Module containing routines needed to average the solution in the zeta direction.
!> General idea of averaging with both FV and DG elements:
!> We calculate the amount of FV and DG elements in the average direction for each cell in the 2D output mesh. If more than 50%
!> of the cells are FV cells, then the averaged DG cells will be converted to the FV grid and the output of this cell will be
!> like for a pure FV cells. The other way round, if less than 50% are FV cells, then the FV averaged data will be converted
!> to the DG grid and the cell will be visualized as a pure DG cell.
!> The averaging process relies on the IJK sorting of the mesh and only works with meshes that include one. All elements with
!> the same k index willl be averaged into one cell.
!===================================================================================================================================
MODULE MOD_Visu_Avg2D
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitAverage2D
  MODULE PROCEDURE InitAverage2D
END INTERFACE

INTERFACE BuildVandermonds_Avg2D
  MODULE PROCEDURE BuildVandermonds_Avg2D
END INTERFACE

INTERFACE Average2D
  MODULE PROCEDURE Average2D
END INTERFACE

INTERFACE WriteAverageToHDF5
  MODULE PROCEDURE WriteAverageToHDF5
END INTERFACE

PUBLIC:: WriteAverageToHDF5
PUBLIC:: InitAverage2D
PUBLIC:: BuildVandermonds_Avg2D
PUBLIC:: Average2D

CONTAINS

!===================================================================================================================================
!> Initialization of the average routines. The IJK sorting of the elements (which is necessary for averaging) is read from the
!> mesh file and from that a mapping is built. This mapping will lead from the i and j indizes of the elements to the element
!> index in the averaged output arrays. Also amount of FV cells in the average direction for each 2D cell will be
!> calculated and stored.
!===================================================================================================================================
SUBROUTINE InitAverage2D()
USE MOD_PreProc
USE MOD_Globals
USE MOD_Visu_Vars
USE MOD_HDF5_Input         ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input         ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile,ReadArray,DatasetExists
USE MOD_Mesh_Vars          ,ONLY: nElems,nGlobalElems,offsetElem
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: exists
INTEGER                          :: ii,jj,iElem
!===================================================================================================================================
! Read the IJK sorting from the mesh file. This is required!!
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,'nElems_IJK',exists)
SDEALLOCATE(Elem_IJK)
SDEALLOCATE(Elem_IJK_glob)
ALLOCATE(Elem_IJK(3,nElems))
ALLOCATE(Elem_IJK_glob(3,nGlobalElems))
IF (exists) THEN
  CALL ReadArray('nElems_IJK',1,(/3/),0,1,IntArray=nElems_IJK)
  CALL ReadArray('Elem_IJK',2,(/3,nElems/),offsetElem,2,IntArray=Elem_IJK)
  CALL ReadArray('Elem_IJK',2,(/3,nGlobalElems/),0,2,IntArray=Elem_IJK_glob)
ELSE
  nElems_IJK(1) = 1
  nElems_IJK(2) = nElems
  nElems_IJK(3) = 1
  DO iElem=1,nElems
    Elem_IJK(1,iElem) = 1
    Elem_IJK(2,iElem) = iElem
    Elem_IJK(3,iElem) = 1
  END DO ! iElem
END IF
CALL CloseDataFile()

! For each cell in the 2D grid, store the amount of FV cells in the average direction
SDEALLOCATE(FVAmountAvg2D)
ALLOCATE(FVAmountAvg2D(nElems_IJK(1),nElems_IJK(2)))
FVAmountAvg2D = 0.
DO iElem=1,nElems
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  FVAmountAvg2D(ii,jj) = FVAmountAvg2D(ii,jj) + FV_Elems_loc(iElem)
END DO
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,FVAmountAvg2D,nElems_IJK(1)*nElems_IJK(2),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif
FVAmountAvg2D = FVAmountAvg2D / REAL(nElems_IJK(3))

! Create a mapping from the IJ indizes of each cell to the 2D visu elements, for DG and FV seperately
nElemsAvg2D_DG = 0
nElemsAvg2D_FV = 0
SDEALLOCATE(mapElemIJToDGElemAvg2D)
SDEALLOCATE(mapElemIJToFVElemAvg2D)
ALLOCATE(mapElemIJToDGElemAvg2D(nElems_IJK(1),nElems_IJK(2)))
ALLOCATE(mapElemIJToFVElemAvg2D(nElems_IJK(1),nElems_IJK(2)))
mapElemIJToDGElemAvg2D = 0
mapElemIJToDGElemAvg2D = 0
DO jj=1,nElems_IJK(2); DO ii=1,nElems_IJK(1)
  IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN
    nElemsAvg2D_DG = nElemsAvg2D_DG + 1
    mapElemIJToDGElemAvg2D(ii,jj) = nElemsAvg2D_DG
  ELSE
    nElemsAvg2D_FV = nElemsAvg2D_FV + 1
    mapElemIJToFVElemAvg2D(ii,jj) = nElemsAvg2D_FV
  END IF
END DO; END DO ! ii,jj=1,nElems_IJK
END SUBROUTINE InitAverage2D

!===================================================================================================================================
!> Builds Vandermonde matrizes thar are needed by the averaging routines. These matrizes are:
!> * Vdm_DGToFV: Conversion of DG solution to FV solution on the calc polynomial degree
!> * Vdm_FVToDG: Inverse of the above
!> * Vdm_DGToVisu: Conversion from calc to visu for DG
!> * Vdm_FVToVisu: Conversion from calc to visu for FV
!===================================================================================================================================
SUBROUTINE BuildVandermonds_Avg2D(NCalc_DG,NCalc_FV)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars          ,ONLY: Vdm_DGToFV,Vdm_FVToDG,Vdm_DGToVisu,Vdm_FVToVisu,NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu,NVisu_FV
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
#if FV_ENABLED
USE MOD_FV_Basis           ,ONLY: FV_GetVandermonde
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: NCalc_DG,NCalc_FV
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if FV_ENABLED
REAL,ALLOCATABLE   :: FV_Vdm(:,:),FV_sVdm(:,:)
REAL,ALLOCATABLE   :: FVdouble(:,:)
INTEGER            :: i
#endif
!===================================================================================================================================
SDEALLOCATE(Vdm_DGToFV)
SDEALLOCATE(Vdm_FVToDG)
SDEALLOCATE(Vdm_DGToVisu)
SDEALLOCATE(Vdm_FVToVisu)
#if FV_ENABLED
ALLOCATE(Vdm_DGToFV  (0:NCalc_FV,0:NCalc_DG))
ALLOCATE(Vdm_FVToDG  (0:NCalc_DG,0:NCalc_FV))
ALLOCATE(Vdm_DGToVisu(0:NVisu   ,0:NCalc_DG))
ALLOCATE(Vdm_FVToVisu(0:NVisu_FV,0:NCalc_FV))
ALLOCATE(FVdouble(0:NVisu_FV,0:PP_N))
FVdouble = 0.
DO i = 0, PP_N
  FVdouble(i*2  ,i) = 1.
  FVdouble(i*2+1,i) = 1.
END DO ! i = 0, PP_N

IF (NCalc_FV.EQ.NCalc_DG) THEN
  CALL FV_GetVandermonde(PP_N,NodeType, Vdm_DGToFV, Vdm_FVToDG)
  Vdm_FVToVisu = FVdouble
ELSE IF (NCalc_FV.EQ.(PP_N+1)*2-1) THEN
  ALLOCATE(FV_Vdm  (0:PP_N    ,0:PP_N))
  ALLOCATE(FV_sVdm (0:PP_N    ,0:PP_N))
  CALL FV_GetVandermonde(PP_N,NodeType, FV_Vdm, FV_sVdm)
  Vdm_DGToFV = MATMUL(FVdouble,FV_Vdm)
  Vdm_FVToDG = MATMUL(FV_sVdm,0.5*TRANSPOSE(FVdouble))
  Vdm_FVToVisu = 0.
  DO i = 0, NVisu_FV
    Vdm_FVToVisu(i,i) = 1.
  END DO
  DEALLOCATE(FV_Vdm)
  DEALLOCATE(FV_sVdm)
ELSE
  CALL CollectiveStop(__STAMP__,&
      "BuildVandermonds_Avg2D called with wrong NCalc_DG and NCalc_FV")
END IF
DEALLOCATE(FVdouble)
#else
ALLOCATE(Vdm_DGToFV  (0:0       ,0:NCalc_DG))
ALLOCATE(Vdm_FVToDG  (0:NCalc_DG,0:0       ))
ALLOCATE(Vdm_DGToVisu(0:NVisu   ,0:NCalc_DG))
ALLOCATE(Vdm_FVToVisu(0:NVisu_FV,0:0       ))
#endif
CALL GetVandermonde(NCalc_DG,NodeType,NVisu,NodeTypeVisuPosti,Vdm_DGToVisu,modal=.FALSE.)

END SUBROUTINE BuildVandermonds_Avg2D

!===================================================================================================================================
!> Main routine to do the averaging.
!===================================================================================================================================
SUBROUTINE Average2D(nVarCalc_DG,nVarCalc_FV,NCalc_DG,NCalc_FV,nElems_DG,nElems_FV,&
    NodeTypeCalc_DG,UCalc_DG,UCalc_FV,&
    Vdm_DGToFV,Vdm_FVToDG,Vdm_DGToVisu,Vdm_FVToVisu, &
    startIndexMapVarCalc,endIndexMapVarCalc,mapVarCalc, &
    UVisu_DG,UVisu_FV)
USE MOD_PreProc
USE MOD_Globals
USE MOD_Visu_Vars          ,ONLY: Elem_IJK,FVAmountAvg2D,Avg2DHDF5Output
USE MOD_Visu_Vars          ,ONLY: mapDGElemsToAllElems,mapFVElemsToAllElems
USE MOD_Visu_Vars          ,ONLY: mapElemIJToDGElemAvg2D,mapElemIJToFVElemAvg2D,mapAllVarsToVisuVars
USE MOD_Visu_Vars          ,ONLY: nVarVisu,NVisu,NVisu_FV,nElemsAvg2D_FV,nElemsAvg2D_DG
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVISUFVEqui,xGP
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nVarCalc_DG,nVarCalc_FV,NCalc_DG,NCalc_FV,nElems_DG,nElems_FV
CHARACTER(LEN=255),INTENT(IN) :: NodeTypeCalc_DG
REAL,INTENT(IN)               :: Vdm_DGToFV(0:NCalc_FV,0:NCalc_DG)
REAL,INTENT(IN)               :: Vdm_FVToDG(0:NCalc_DG,0:NCalc_FV)
REAL,INTENT(IN)               :: Vdm_DGToVisu(0:NVisu   ,0:NCalc_DG)
REAL,INTENT(IN)               :: Vdm_FVToVisu(0:NVisu_FV,0:NCalc_FV)
INTEGER,INTENT(IN)            :: startIndexMapVarCalc,endIndexMapVarCalc,mapVarCalc(startIndexMapVarCalc:endIndexMapVarCalc)
REAL,INTENT(IN)               :: UCalc_DG(0:NCalc_DG,0:NCalc_DG,0:NCalc_DG,1:nElems_DG     ,1:nVarCalc_DG)
REAL,INTENT(IN)               :: UCalc_FV(0:NCalc_FV,0:NCalc_FV,0:NCalc_FV,1:nElems_FV     ,1:nVarCalc_FV)
REAL,INTENT(INOUT)            :: UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,1:nElemsAvg2D_DG,1:nVarVisu)
REAL,INTENT(INOUT)            :: UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,1:nElemsAvg2D_FV,1:nVarVisu)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

REAL,ALLOCATABLE :: Utmp(:,:),Utmp2(:,:)
INTEGER          :: iElem, iElem_DG, iElem_FV, iElemAvg, ii,jj,k,iVar,iVarCalc,iVarVisu
REAL             :: xGP_loc(0:NCalc_DG),wGP_loc(0:NCalc_DG)
REAL,ALLOCATABLE :: UAvg_DG(:,:,:,:)
REAL,ALLOCATABLE :: UAvg_FV(:,:,:,:)
REAL             :: dx
REAL             :: dxSum_DG(nElemsAvg2D_DG),dxSum_FV(nElemsAvg2D_FV)
!===================================================================================================================================
ALLOCATE(UAvg_DG(0:NCalc_DG,0:NCalc_DG,nElemsAvg2D_DG,nVarVisu))
ALLOCATE(UAvg_FV(0:NCalc_FV,0:NCalc_FV,nElemsAvg2D_FV,nVarVisu))
UAvg_DG = 0.
UAvg_FV = 0.
dxSum_DG = 0.
dxSum_FV = 0.
CALL GetNodesAndWeights(NCalc_DG,NodeTypeCalc_DG,xGP_loc,wGP_loc)

! average all DG elements first
SWRITE(*,*) " [Avg2D] Average DG elements"
ALLOCATE(Utmp (0:NCalc_DG,0:NCalc_DG))
ALLOCATE(Utmp2(0:NCalc_FV,0:NCalc_FV))
DO iElem_DG = 1,nElems_DG                ! iterate over all DG elements
  iElem = mapDGElemsToAllElems(iElem_DG) ! get global element index
  IF (iElem.EQ.0) CYCLE
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  dx = (Elem_xGP(3,0,0,PP_N,iElem) - Elem_xGP(3,0,0,0,iElem))*2.0/(xGP(PP_N)-xGP(0))
  IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN  ! the averaged ii,jj-th element is rather a DG element and this
                                         ! element is a DG element => just average this element and add to UAvg_DG
    iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
    dxSum_DG(iElemAvg) = dxSum_DG(iElemAvg) + dx
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        DO k=0,NCalc_DG
          UAvg_DG(:,:,iElemAvg,iVarVisu) = UAvg_DG(:,:,iElemAvg,iVarVisu) + wGP_loc(k)/2.*dx * UCalc_DG(:,:,k,iElem_DG,iVarCalc)
        END DO
      END IF
    END DO
  ELSE ! the averaged ii,jj-th element is rather a FV element, but this element is a DG element
       ! => average this element as DG element and convert the average (Utmp) to FV and add to UAvg_FV
    iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
    dxSum_FV(iElemAvg) = dxSum_FV(iElemAvg) + dx
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        Utmp = 0.
        DO k=0,NCalc_DG
          Utmp(:,:) = Utmp(:,:) + wGP_loc(k)/2.*dx * UCalc_DG(:,:,k,iElem_DG,iVarCalc)
        END DO
        CALL ChangeBasis2D(NCalc_DG,NCalc_FV,Vdm_DGToFV,Utmp,Utmp2)
        UAvg_FV(:,:,iElemAvg,iVarVisu) = UAvg_FV(:,:,iElemAvg,iVarVisu) + Utmp2
      END IF
    END DO
  END IF
END DO

! all DG elements are averaged and we will now average all FV elements, which requires conversions
! from FV to DG, the other way round then before.
SWRITE(*,*) " [Avg2D] Average FV elements"
DEALLOCATE(Utmp)
DEALLOCATE(Utmp2)
ALLOCATE(Utmp( 0:NCalc_FV,0:NCalc_FV))
ALLOCATE(Utmp2(0:NCalc_DG,0:NCalc_DG))
DO iElem_FV = 1,nElems_FV                ! iterate over all FV elements
  iElem = mapFVElemsToAllElems(iElem_FV) ! get global element index
  IF (iElem.EQ.0) CYCLE
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  dx = (Elem_xGP(3,0,0,PP_N,iElem) - Elem_xGP(3,0,0,0,iElem))*2.0/(xGP(PP_N)-xGP(0))
  IF (FVAmountAvg2D(ii,jj).GT.0.5) THEN  ! the averaged ii,jj-th element is rather a FV element and this
                                         ! element is a FV element => just average this element and add to UAvg_FV
    iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
    dxSum_FV(iElemAvg) = dxSum_FV(iElemAvg) + dx
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        DO k=0,NCalc_FV
          UAvg_FV(:,:,iElemAvg,iVarVisu) = UAvg_FV(:,:,iElemAvg,iVarVisu) + 1./(NCalc_FV+1)*dx * UCalc_FV(:,:,k,iElem_FV,iVarCalc)
        END DO
      END IF
    END DO
  ELSE ! the averaged ii,jj-th element is rather a DG element, but this element is a FV element
       ! => average this element as FV element and convert the average (Utmp) to DG and add to UAvg_DG
    iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
    dxSum_DG(iElemAvg) = dxSum_DG(iElemAvg) + dx
    DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
      iVarVisu = mapAllVarsToVisuVars(iVar)
      IF (iVarVisu.GT.0) THEN
        iVarCalc = mapVarCalc(iVar)
        Utmp = 0.
        DO k=0,NCalc_FV
          Utmp(:,:) = Utmp(:,:) +  1./(NCalc_FV+1)*dx * UCalc_FV(:,:,k,iElem_FV,iVarCalc)
        END DO
        CALL ChangeBasis2D(NCalc_FV,NCalc_DG,Vdm_FVToDG,Utmp,Utmp2)
        UAvg_DG(:,:,iElemAvg,iVarVisu) = UAvg_DG(:,:,iElemAvg,iVarVisu) + Utmp2
      END IF
    END DO
  END IF
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,UAvg_DG,nElemsAvg2D_DG*(NCalc_DG+1)**2*nVarVisu,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,UAvg_FV,nElemsAvg2D_FV*(NCalc_FV+1)**2*nVarVisu,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,dxSum_DG,nElemsAvg2D_DG,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,dxSum_FV,nElemsAvg2D_FV,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(UAvg_DG      ,0     ,nElemsAvg2D_DG*(NCalc_DG+1)**2*nVarVisu,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(UAvg_FV      ,0     ,nElemsAvg2D_FV*(NCalc_FV+1)**2*nVarVisu,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(dxSum_DG     ,0     ,nElemsAvg2D_DG,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(dxSum_FV     ,0     ,nElemsAvg2D_FV,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif
IF (MPIRoot) THEN
  DO iElem_DG = 1,nElemsAvg2D_DG
    UAvg_DG(:,:,iElem_DG,:) = UAvg_DG(:,:,iElem_DG,:) / dxSum_DG(iElem_DG)
  END DO
  DO iElem_FV = 1,nElemsAvg2D_FV
    UAvg_FV(:,:,iElem_FV,:) = UAvg_FV(:,:,iElem_FV,:) / dxSum_FV(iElem_FV)
  END DO
END IF

#if USE_MPI
IF (Avg2DHDF5Output) THEN
  ! Distribute the averaged data back
  CALL MPI_BCAST(UAvg_DG,nElemsAvg2D_DG*(NCalc_DG+1)**2*nVarVisu,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
  CALL MPI_BCAST(UAvg_FV,nElemsAvg2D_FV*(NCalc_FV+1)**2*nVarVisu,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
END IF
#endif

IF ((MPIRoot).OR.(Avg2DHDF5Output)) THEN
  ! Convert the averaged data to the visu grid
  DO iVar=startIndexMapVarCalc,endIndexMapVarCalc
    iVarVisu = mapAllVarsToVisuVars(iVar)
    IF (iVarVisu.GT.0) THEN
      DO iElemAvg=1,nElemsAvg2D_DG
        CALL ChangeBasis2D(NCalc_DG,NVisu   ,Vdm_DGToVisu  ,UAvg_DG(:,:,iElemAvg,iVarVisu),UVisu_DG(:,:,0,iElemAvg,iVarVisu))
      END DO
      IF (NCalc_FV.EQ.NVisu_FV) THEN
        UVisu_FV(:,:,0,:,iVarVisu) = UAvg_FV(:,:,:,iVarVisu)
      ELSE
        DO iElemAvg=1,nElemsAvg2D_FV
          CALL ChangeBasis2D(NCalc_FV,NVisu_FV,Vdm_FVToVisu,UAvg_FV(:,:,iElemAvg,iVarVisu),UVisu_FV(:,:,0,iElemAvg,iVarVisu))
        END DO
      END IF
    END IF
  END DO
END IF


! clear local variables
DEALLOCATE(Utmp)
DEALLOCATE(Utmp2)
DEALLOCATE(UAvg_DG)
DEALLOCATE(UAvg_FV)
END SUBROUTINE Average2D

!===================================================================================================================================
!> HDF5 output routine for 2D averaged data. The very simple idea is to use the 3D mesh (since creating a new 2D mesh would be
!> very complicated) and copy the averaged solution to the third dimension.
!> The output must be done on PP_N and on the calculation NodeType, this is enforced in the InitFile routine.
!===================================================================================================================================
SUBROUTINE WriteAverageToHDF5(nVar,NVisu,NVisu_FV,NodeType,OutputTime,MeshFileName,UVisu_DG,UVisu_FV)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE HDF5
USE MOD_HDF5_Output,    ONLY: GatheredWriteArray,WriteAttribute,WriteHeader,WriteAdditionalElemData
USE MOD_IO_HDF5,        ONLY: OpenDataFile,CloseDataFile,AddToElemData
USE MOD_HDF5_Input,     ONLY: File_ID
USE MOD_Output_Vars,    ONLY: ProjectName
USE MOD_Visu_Vars,      ONLY: Elem_IJK,VarnamesAll,mapAllVarsToVisuVars,nVarAll
USE MOD_Visu_Vars,      ONLY: nElemsAvg2D_DG,mapElemIJToDGElemAvg2D,nElemsAvg2D_FV,mapElemIJToFVElemAvg2D
USE MOD_Mesh_Vars,      ONLY: nGlobalElems,offsetElem,nElems
#if FV_ENABLED
#if FV_RECONSTRUCT
USE MOD_Visu_Vars,      ONLY: NCalc_FV
#endif
USE MOD_ChangeBasis,    ONLY: ChangeBasis2D
USE MOD_FV_Vars,        ONLY: FV_Elems
USE MOD_IO_HDF5,        ONLY: ElementOut
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVar,NVisu,NVisu_FV
CHARACTER(LEN=255)             :: NodeType
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN)                :: UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,nElemsAvg2D_DG,nVar)
REAL,INTENT(IN)                :: UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV,nVar)
CHARACTER(LEN=255),INTENT(IN)  :: MeshFileName
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: UVisu3D(nVar,0:NVisu,0:NVisu,0:NVisu,nElems)
INTEGER             :: iElemAvg2D_DG,iElemAvg2D_FV,iElem,ii,jj,k,iVar,iVarVisu
INTEGER(HSIZE_T)    :: Dimsf(5)
INTEGER(HID_T)      :: DSet_ID,FileSpace,HDF5DataType
CHARACTER(LEN=255)  :: StrVarNames(nVar)
CHARACTER(LEN=255)  :: FileName
#if FV_ENABLED
#if FV_RECONSTRUCT
INTEGER             :: i
REAL                :: Vdm_FVRecon_PP_N(0:PP_N,0:NVisu_FV)
#endif
#endif
!===================================================================================================================================
#if FV_RECONSTRUCT
Vdm_FVRecon_PP_N = 0.
DO i=0,PP_N
  Vdm_FVRecon_PP_N(i,2*i) = 1.
END DO ! i=0,PP_N
#endif

!================== Prepare output data =============================!

! Create full three dimensional array to write to the 3D mesh file.
DO iElem=1,nElems
  ii = Elem_IJK(1,iElem)
  jj = Elem_IJK(2,iElem)
  iElemAvg2D_DG = mapElemIJToDGElemAvg2D(ii,jj)
  iElemAvg2D_FV = mapElemIJToFVElemAvg2D(ii,jj)
  IF (iElemAvg2D_DG.GT.0) THEN
    ! Averaged as a DG element
    DO k=0,NVisu
      DO iVar=1,nVar
        UVisu3D(iVar,:,:,k,iElem) = UVisu_DG(:,:,0,iElemAvg2D_DG,iVar)
      END DO ! iVar=1,nVar
    END DO ! k=0,NVisu
#if FV_ENABLED
    FV_Elems(iElem) = 0
  ELSE
    ! Averaged as a FV Element
    DO k=0,NVisu
      DO iVar=1,nVar
#if FV_RECONSTRUCT
        ! For reconstruction, FV elements will be on 2*(PP_N+1) points
        CALL ChangeBasis2D(NCalc_FV,PP_N,Vdm_FVRecon_PP_N,UVisu_FV(:,:,0,iElemAvg2D_FV,iVar),UVisu3D(iVar,:,:,k,iElem))
        !UVisu3D(iVar,:,:,k,iElem) = UVisu_FV(:,:,0,iElemAvg2D_FV,iVar)
#else
        UVisu3D(iVar,:,:,k,iElem) = UVisu_FV(:,:,0,iElemAvg2D_FV,iVar)
#endif
      END DO ! iVar=1,nVar
    END DO ! k=0,NVisu
    FV_Elems(iElem) = 1
#endif
  END IF
END DO ! iElem

! Build string array with the visualize variable names
DO iVar=1,nVarAll
  iVarVisu = mapAllVarsToVisuVars(iVar)
  IF (iVarVisu.GT.0) THEN
    StrVarNames(iVarVisu) = VarnamesAll(iVar)
  END IF
END DO ! iVar=1,nVarAll
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif


!================= Create and prepare HDF5 file =======================!
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_avg2D',OutputTime))//'.h5'
IF (MPIRoot) THEN
  ! Create file
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

  ! Write file header with file type 'Avg2D'
  CALL WriteHeader('Avg2D',File_ID)

  ! Preallocate the data space for the dataset.
  Dimsf=(/nVar,NVisu+1,NVisu+1,NVisu+1,nGlobalElems/)

  CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
  ! Create the dataset with default properties.
  HDF5DataType=H5T_NATIVE_DOUBLE
  CALL H5DCREATE_F(File_ID,'DG_Solution', HDF5DataType, FileSpace, DSet_ID, iError)
  ! Close the filespace and the dataset
  CALL H5DCLOSE_F(Dset_id, iError)
  CALL H5SCLOSE_F(FileSpace, iError)

  ! Write dataset properties "N","Time","MeshFile","NodeType","VarNames","NComputation"
  CALL WriteAttribute(File_ID,'N',1,IntScalar=NVisu)
  CALL WriteAttribute(File_ID,'Time',1,RealScalar=OutputTime)
  CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=(/MeshFileName/))
  CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttribute(File_ID,'VarNames',nVar,StrArray=StrVarNames)
  CALL WriteAttribute(File_ID,'NComputation',1,IntScalar=PP_N)

  CALL CloseDataFile()
END IF

!================= Actual data output =======================!
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif
CALL GatheredWriteArray(TRIM(FileName),create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/nVar,NVisu+1,NVisu+1,NVisu+1,nGlobalElems/),&
                        nVal=      (/nVar,NVisu+1,NVisu+1,NVisu+1,nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=UVisu3D)

#if FV_ENABLED
!=========== FV/DG element distribution output ===============!
CALL WriteAdditionalElemData(FileName,ElementOut)
#endif

END SUBROUTINE WriteAverageToHDF5

END MODULE MOD_Visu_Avg2D
