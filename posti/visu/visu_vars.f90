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
!==================================================================================================================================
!> Contains global variables provided by the posti routines 
!==================================================================================================================================
MODULE MOD_Visu_Vars
USE ISO_C_BINDING
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
CHARACTER(LEN=255)                :: fileType = ""               ! possible values: 
                                                                 ! * 'State' for FLEXI state files matching the compiled EOS
                                                                 ! * 'Generic'  
                                                                 ! * 'Mesh' 
CHARACTER(LEN=255)                :: prmfile_old = ""            ! saves the filename of the previous FLEXI parameter file 
CHARACTER(LEN=255)                :: statefile_old = ""          ! saves the filename of the previous state (*.h5)
CHARACTER(LEN=255)                :: MeshFile = ""               ! acutal filename of the mesh used for visualization
CHARACTER(LEN=255)                :: MeshFile_state = ""         ! filename of the mesh given in the state file
CHARACTER(LEN=255)                :: MeshFile_old = ""           ! saves  previous MeshFile
CHARACTER(LEN=255)                :: NodeTypeVisuPosti = "VISU"  ! NodeType used for visualization output
CHARACTER(LEN=255)                :: NodeTypeVisuPosti_old = ""  ! saves previous NodeType
INTEGER                           :: NVisu                       ! polynomial degree of the visualization
INTEGER                           :: NVisu_old = -1              ! saves previous NVisu
INTEGER                           :: NVisu_FV                    ! number of output points for FV elements (always == 2*(PP_N+1))
INTEGER                           :: nVar_State                  ! number of variables in the state file 
INTEGER                           :: nVar_State_old = -1         ! saves previous nVar_State_old
INTEGER                           :: nElems_DG                   ! number of DG elements in state
INTEGER                           :: nElems_FV                   ! number of FV elements in state
LOGICAL                           :: withDGOperator              ! flag indicating if call of 'DGTimeDerivative' is required
LOGICAL                           :: withDGOperator_old = .FALSE.! saves previous withDGOperator
REAL                              :: OutputTime                  ! simulation time of actual state file 
LOGICAL                           :: hasFV_Elems = .FALSE.       ! flag indicating if state file contains any FV elements
LOGICAL                           :: DGonly = .FALSE.            ! flag to force visualization of FV elements as DG elements
LOGICAL                           :: DGonly_old = .TRUE.         ! saves previous DGonly
INTEGER,ALLOCATABLE               :: mapDGElemsToAllElems(:)     ! maps element index of DG elements to all elements
INTEGER,ALLOCATABLE               :: mapFVElemsToAllElems(:)
INTEGER,ALLOCATABLE               :: FV_Elems_loc(:)             ! current distribution of FV/DG elems
INTEGER,ALLOCATABLE               :: FV_Elems_old(:)             ! saves previous FV_Elems, which holds DG/FV elements distribution
INTEGER                           :: VisuDimension               ! TODO: Avg2D
INTEGER                           :: meshMode_old=0              ! Used to check if InitMesh must be called again with different
                                                                 ! mesh mode
LOGICAL                           :: doSurfVisu                  ! Flag indicating if any surfaces need to be visualized

                
! The following flags indicate if during successive visualizations of (different) state files the respective properties
! changed. For example the mesh file of different state files in a timeseries is the same ...
LOGICAL                           :: changedStateFile            ! 
LOGICAL                           :: changedMeshFile             ! 
LOGICAL                           :: changedNVisu                ! 
LOGICAL                           :: changedVarNames             ! variables selected for visualization changed (ParaView plugin) 
LOGICAL                           :: changedFV_Elems             ! different distribution of DG and FV elements 
LOGICAL                           :: changedWithDGOperator       ! 
LOGICAL                           :: changedDGonly               ! 
LOGICAL                           :: changedBCnames              ! BCnames selected for visualization changed (ParaView plugin)

CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarNamesHDF5(:)         ! varnames in state file (DG_Solution, not including generic 
                                                                 ! element- or pointwise)
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarnamesAll(:)          ! all available varnames (state file + dependent vars + generic)
INTEGER                               :: nVarAll                 ! number of all available visu variables
INTEGER                               :: nVarDep                 ! number of dependent variables, that EOS can calculate
INTEGER                               :: nVarVisu                ! number of variables selected for visualization
INTEGER,ALLOCATABLE                   :: mapAllVarsToVisuVars(:) ! maps all available variable index to visualization variable index
INTEGER,ALLOCATABLE                   :: DepTable(:,:)           ! table holding the EOS dependencies required to calculate 
                                                                 ! variables, that depend on other variables (e.g. primitive ...)
                                                                 ! The i-th line of this table holds the dependency informations of
                                                                 ! the i-th quantity on the previous quantities. The j-th column
                                                                 ! is 0 if the i-th quantity does NOT depends on the j-th quantity
                                                                 ! is 1 if the i-th quantity DOES depends on the j-th quantity

REAL,ALLOCATABLE                      :: UCalc_DG(:,:,:,:,:)     ! dependet variables require the computation of intermediate
                                                                 ! variables, that may not be visualized. Therefore the whole
                                                                 ! computation process takes place on this array and is afterwards
                                                                 ! converted to the visualization array (UVisu_DG)
REAL,ALLOCATABLE                      :: UCalc_FV(:,:,:,:,:)
INTEGER                               :: nVarCalc                ! number of (intermediate) variables that must be calculated
INTEGER,ALLOCATABLE                   :: mapDepToCalc(:)         ! maps all dependend variable index to calc variable index
#if FV_ENABLED && FV_RECONSTRUCT
INTEGER                               :: nVarCalc_FV             ! since FV reconstruction is done in primitive quantities, the 
INTEGER,ALLOCATABLE                   :: mapDepToCalc_FV(:)      ! dependencies are different to the DG case, where everything is
                                                                 ! based on conservative quantities
#endif

REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: UVisu_DG(:,:,:,:,:)     ! solution that is written to VTK or send to ParaView
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: UVisu_FV(:,:,:,:,:)     ! 
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsVisu_DG(:,:,:,:,:)! coordinates of UVisu_DG
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsVisu_FV(:,:,:,:,:)! 
INTEGER,ALLOCATABLE,TARGET            :: nodeids_DG(:)           ! nodeids for CoordsVisu_DG
INTEGER,ALLOCATABLE,TARGET            :: nodeids_FV(:)           ! 



! ==============================================================================================================================
! Surface visualization
! ==============================================================================================================================
INTEGER,ALLOCATABLE                   :: DepSurfaceOnly(:)       ! same but for quantities that are exclusively available on BCs

INTEGER                               :: nBCNamesAll                  ! number of all BC names in mesh file
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: BCNamesAll(:)                ! all BC names in mesh file 
INTEGER                               :: nBCNamesVisu                 ! number of BC names selected for visualization
INTEGER,ALLOCATABLE                   :: mapAllBCNamesToVisuBCNames(:)! maps global BCName index to visu BCName index
INTEGER,ALLOCATABLE                   :: mapAllBCNamesToVisuBCNames_old(:) 

INTEGER                               :: nVarSurfVisuAll                 ! number of all avail. vars that are visualized on surf.
INTEGER,ALLOCATABLE                   :: mapAllVarsToSurfVisuVars(:)     ! maps all avail. var index to surf. visu. var index
INTEGER,ALLOCATABLE                   :: mapAllVarsToSurfVisuVars_old(:) ! saves previous mapAllVarsToSurfVisuVars
INTEGER                               :: nBCSidesVisu_DG                 ! number of DG BCsides selected for visualization
INTEGER                               :: nBCSidesVisu_FV                 ! number of FV BCsides selected for visualization
INTEGER,ALLOCATABLE                   :: mapAllBCSidesToDGVisuBCSides(:) ! map global BC side index to DG visu BC sides
INTEGER,ALLOCATABLE                   :: mapAllBCSidesToFVVisuBCSides(:) ! map global BC side index to FV visu BC sides
INTEGER,ALLOCATABLE                   :: nSidesPerBCNameVisu_DG(:)       ! holds number of DG BCsides for each BCName
INTEGER,ALLOCATABLE                   :: nSidesPerBCNameVisu_FV(:)       ! holds number of FV BCsides for each BCName

REAL,ALLOCATABLE                      :: USurfCalc_DG(:,:,:,:)        ! array on which dependent quantities are calculated
REAL,ALLOCATABLE                      :: USurfCalc_FV(:,:,:,:)        ! 

REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: USurfVisu_DG(     :,:,:,:,:) ! surf. solution that is written to VTK or send to ParaView
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: USurfVisu_FV(     :,:,:,:,:) ! 
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsSurfVisu_DG(:,:,:,:,:) ! coordinates of surface solution 
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsSurfVisu_FV(:,:,:,:,:) ! 
INTEGER,ALLOCATABLE,TARGET            :: nodeidsSurf_DG(:)            ! nodeids for surface coordinates
INTEGER,ALLOCATABLE,TARGET            :: nodeidsSurf_FV(:)            ! 


! ==============================================================================================================================
! TODO: Avg2D
! ==============================================================================================================================
INTEGER,ALLOCATABLE,TARGET        :: nodeids_DG_2D(:)           ! visu nodeids
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: CoordsVisu_DG_2D(:,:,:,:,:)! visu coordinates
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: UVisu_DG_2D(:,:,:,:,:)     ! state at visu points
INTEGER,ALLOCATABLE,TARGET        :: nodeids_FV_2D(:)           ! visu nodeids
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: CoordsVisu_FV_2D(:,:,:,:,:)! visu coordinates
REAL(C_DOUBLE),ALLOCATABLE,TARGET :: UVisu_FV_2D(:,:,:,:,:)     ! state at visu points



END MODULE MOD_Visu_Vars
