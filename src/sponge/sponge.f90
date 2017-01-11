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

!==================================================================================================================================
!> Subroutines for applying a sponge layer to the solution, to reduce reflections at the boundaries.
!> This is done by forcing the solution towards a previously known state, the baseflow state.
!> The forcing is applied within the sponge region, where the forcing strength is increased till the end of the region.
!==================================================================================================================================
MODULE MOD_Sponge
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTEGER,PARAMETER :: SPONGESHAPE_RAMP        = 1
INTEGER,PARAMETER :: SPONGESHAPE_CYLINDRICAL = 2

INTEGER,PARAMETER :: SPONGEBASEFLOW_CONSTANT  = 1
INTEGER,PARAMETER :: SPONGEBASEFLOW_EXACTFUNC = 2
INTEGER,PARAMETER :: SPONGEBASEFLOW_FILE      = 3
INTEGER,PARAMETER :: SPONGEBASEFLOW_PRUETT    = 4

INTERFACE InitSponge
  MODULE PROCEDURE InitSponge
END INTERFACE

INTERFACE Sponge
  MODULE PROCEDURE Sponge
END INTERFACE

INTERFACE FinalizeSponge
  MODULE PROCEDURE FinalizeSponge
END INTERFACE


PUBLIC::InitSponge,Sponge,FinalizeSponge
!==================================================================================================================================

PUBLIC::DefineParametersSponge
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersSponge()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Sponge")
CALL prms%CreateLogicalOption('SpongeLayer',    "Turn on to use sponge regions for reducing reflections at boundaries.",'.FALSE.')
CALL prms%CreateRealOption(   'damping',        "Damping factor of sponge (0..1).", '0.01')
CALL prms%CreateIntFromStringOption( 'SpongeShape',    "Set shape of sponge: (1) ramp : cartesian / vector-aligned, (2) "//&
                                                       " cylindrical", multiple=.TRUE.)
CALL addStrListEntry('SpongeShape','ramp',       SPONGESHAPE_RAMP)
CALL addStrListEntry('SpongeShape','cylindrical',SPONGESHAPE_CYLINDRICAL)
!CALL prms%CreateIntOption( 'SpongeShape',    "Set shape of sponge: (1) ramp : cartesian / vector-aligned, (2) "//&
!                                                       " cylindrical", multiple=.TRUE.)
CALL prms%CreateRealOption(   'SpongeDistance', "Length of sponge ramp. The sponge will have maximum strength at the end "//&
                                                "of the ramp and after that point.", multiple=.TRUE.)
CALL prms%CreateRealArrayOption('xStart',       "Coordinates of start postion of sponge ramp (SpongeShape=ramp) "//&
                                                "or center (SpongeShape=cylindrical).", multiple=.TRUE.)
CALL prms%CreateRealArrayOption('SpongeDir',    "Direction vector of the sponge ramp (SpongeShape=ramp)", multiple=.TRUE.)
CALL prms%CreateRealOption(   'SpongeRadius',   "Radius of the sponge zone (SpongeShape=cylindrical)", multiple=.TRUE.)
#if (PP_dim==3)
CALL prms%CreateRealArrayOption('SpongeAxis',   "Axis vector of cylindrical sponge (SpongeShape=cylindrical)", multiple=.TRUE.)
#endif
CALL prms%CreateLogicalOption('SpongeViz',      "Turn on to write a visualization file of sponge region and strength.",'.FALSE.')
CALL prms%CreateIntFromStringOption( 'SpongeBaseFlow', "Type of baseflow to be used for sponge. (1) constant: fixed state,"//&
                                                "(2) exactfunction: exact function, (3) file: read baseflow file, (4) pruett: "//&
                                                "temporally varying, solution adaptive Pruett baseflow",'1')
CALL addStrListEntry('SpongeBaseFlow','constant',     SPONGEBASEFLOW_CONSTANT)
CALL addStrListEntry('SpongeBaseFlow','exactfunction',SPONGEBASEFLOW_EXACTFUNC)
CALL addStrListEntry('SpongeBaseFlow','file',         SPONGEBASEFLOW_FILE)
CALL addStrListEntry('SpongeBaseFlow','pruett',       SPONGEBASEFLOW_PRUETT)
CALL prms%CreateIntOption(    'SpongeRefState', "Index of refstate in ini-file (SpongeBaseFlow=cartesian)")
CALL prms%CreateIntOption(    'SpongeExactFunc',"Index of exactfunction (SpongeBaseFlow=exactfunction)")
CALL prms%CreateStringOption( 'SpongeBaseFlowFile',"FLEXI solution (e.g. TimeAvg) file from which baseflow is read.")
CALL prms%CreateRealOption(   'tempFilterWidth',"Temporal filter width used to advance Pruett baseflow in time.)")
END SUBROUTINE DefineParametersSponge

!==================================================================================================================================
!> \brief Initialize sponge region (get parameters, allocate arrays).
!>
!> Important parameters are:
!>  - damping: strength of the sponge
!>  - SpongeShape: Either a ramp (1) or a cylinder (2).
!>  - SpBaseFlowType: Constant (1), evaluation of exact function (2), from a .h5 file (3) or a pruett type base flow (4).
!>
!> Depending on the chosen sponge types, the necessary parameters are read in and the strength of the sponge is calculated.
!> After this, the base flow of the sponge is initialized. The base flow is stored in the whole computational domain,
!> not just in the sponge region. This allows e.g. for easier change of the sponge region.
!==================================================================================================================================
SUBROUTINE InitSponge
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_Sponge_Vars
USE MOD_Exactfunc,    ONLY:ExactFunc
USE MOD_Equation_Vars,ONLY:RefStateCons
USE MOD_Mesh_Vars,    ONLY:Elem_xGP,nElems
USE MOD_Output_Vars,  ONLY:ProjectName
USE MOD_PruettDamping,ONLY:InitPruettDamping
USE MOD_Restart_Vars, ONLY:DoRestart,RestartTime,RestartFile
USE MOD_Equation_Vars,ONLY:IniExactFunc
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k
INTEGER             :: SpongeExactFunc,SpongeRefState,SpBaseFlowType
CHARACTER(LEN=300)  :: BaseFlowFile
LOGICAL             :: validBaseFlowFile
!==================================================================================================================================
doSponge=GETLOGICAL('SpongeLayer','.FALSE.')
IF(.NOT.doSponge) RETURN
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SPONGE...'

damping   = GETREAL('damping','0.01')
IF(damping .LE. 0.)  THEN
  doSponge = .FALSE.
  RETURN
END IF

SpongeViz=GETLOGICAL('SpongeViz','.FALSE.')

CALL CalcSpongeRamp()

CalcPruettDamping=.FALSE.
! Readin of Baseflow parameters
SpBaseFlowType = GETINTFROMSTR('SpongeBaseFlow')
SELECT CASE(SpBaseflowType)
CASE(SPONGEBASEFLOW_CONSTANT) ! constant baseflow from refstate
  spongeRefState  = GETINT('SpongeRefState')
CASE(SPONGEBASEFLOW_EXACTFUNC) ! Exactfunction
  spongeExactFunc = GETINT('SpongeExactFunc')
CASE(SPONGEBASEFLOW_FILE) ! Base Flow from .h5 File
  BaseFlowFile = GETSTR('SpongeBaseFlowFile')
CASE(SPONGEBASEFLOW_PRUETT) ! Pruett
  CalcPruettDamping=.TRUE.
  CALL InitPruettDamping()
  IF(DoRestart)THEN
    BaseFlowFile    = GETSTR('SpongeBaseFlowFile','')
    IF (TRIM(BaseFlowFile) .EQ. '') THEN
      ! If no base flow file has been specified, assume a standard name for the base flow file
      BaseFlowFile=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BaseFlow',RestartTime))//'.h5'
      ! Check if this file exists
      INQUIRE(FILE=TRIM(BaseFlowFile),EXIST=validBaseFlowFile)
      IF (.NOT.validBaseFlowFile) THEN
        ! If the assumed base flow file does not exist, use the restart state to initialize the sponge base flow
        BaseFlowFile = RestartFile
        SWRITE(UNIT_StdOut,'(A)') 'WARNING: No baseflow file found! Using the restart state to initialize sponge base flow.'
      END IF
    ELSE
      ! check if baseflow exists
      INQUIRE(FILE=TRIM(BaseFlowFile),EXIST=validBaseFlowFile)
      IF (.NOT.validBaseFlowFile) THEN
        CALL CollectiveStop(__STAMP__,&
          'ERROR: Sponge base flow file '//TRIM(BaseFlowFile)//' does not exist.')
      END IF
    END IF
  ELSE
    spongeExactFunc = GETINT('SpongeExactFunc','-1')
  END IF
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    "Undefined SpongeBaseFlow!")   
END SELECT

! Preparation of the baseflow on each Gauss Point
SWRITE(UNIT_StdOut,'(A)') '  Initialize Sponge Base Flow...'
ALLOCATE(SpBaseFlow(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
SELECT CASE(SpBaseflowType)
CASE(SPONGEBASEFLOW_CONSTANT) ! constant baseflow from refstate
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      SpBaseFlow(:,i,j,k,iElem)=RefStateCons(spongeRefState,:)
    END DO; END DO; END DO
  END DO
CASE(SPONGEBASEFLOW_EXACTFUNC) ! Exactfunction
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      !Save exactFunc state for later use
      CALL ExactFunc(SpongeExactFunc,0.,Elem_xGP(:,i,j,k,iElem),SpBaseFlow(:,i,j,k,iElem))
    END DO; END DO; END DO
  END DO
CASE(SPONGEBASEFLOW_FILE) ! Base Flow from .h5 File
  CALL ReadBaseFlow(BaseFlowfile)
!readin of the hdf5 base flow solution
CASE(SPONGEBASEFLOW_PRUETT) ! Pruett: RefState for computation from scratch, Base flow file for restart
  IF(DoRestart)THEN
    CALL ReadBaseFlow(BaseFlowfile)
  ELSE
    IF (SpongeExactFunc.LT.0) THEN
      SWRITE(UNIT_StdOut,'(A)') 'WARNING: No sponge exact func given! Use ini exact func instead.'
      SpongeExactFunc = IniExactFunc
    END IF
    DO iElem=1,nElems
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        !Save exactFunc state for later use
        CALL ExactFunc(SpongeExactFunc,0.,Elem_xGP(:,i,j,k,iElem),SpBaseFlow(:,i,j,k,iElem))
      END DO; END DO; END DO
    END DO
  END IF
END SELECT

SWRITE(UNIT_StdOut,'(A)')' INIT SPONGE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitSponge



!==================================================================================================================================
!> \brief Compute sponge shape and strength at each solution point. Visualize sponge.
!>
!> First, depending on the shape (linear or cylindrical), the strength  of the shape without the damping factor (x_star) is
!> calculated on the solution points.
!> From this, a mapping is build which contains only the elements with x_star > 0 somewhere, which is used to later apply
!> the sponge only to regions where it is needed.
!> In this sponge region, the final strength of the sponge is build by limiting x_star to [0,1] and mutiply it by the damping.
!> If set in the parameter file, a visualization of the sponge strength is written as .vtu files.
!> At the end, the sponge is pre-multiplied by the Jacobian since we need to do this anyway when the sponge is applied.
!==================================================================================================================================
SUBROUTINE CalcSpongeRamp()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_Sponge_Vars
USE MOD_Output_Vars       ,ONLY:ProjectName
USE MOD_Mesh_Vars         ,ONLY:Elem_xGP
USE MOD_Interpolation_Vars,ONLY:NodeTypeCL,NodeType
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Output_Vars       ,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_ChangeBasisByDim  ,ONLY:ChangeBasisVolume
USE MOD_Mesh_Vars         ,ONLY:sJ,nElems
USE MOD_VTK               ,ONLY:WriteDataToVTK
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                                 :: applySponge(nElems)
INTEGER                                 :: iElem,iSpongeElem,i,j,k,iRamp
CHARACTER(LEN=255)                      :: FileString,VarNameSponge(1)
REAL,DIMENSION(  0:PP_N,0:PP_N,0:PP_NZ) :: sigma, x_star
REAL                                    :: r_vec(PP_dim)
REAL,ALLOCATABLE,TARGET                 :: SpDummy(:,:,:,:)
REAL,ALLOCATABLE,TARGET                 :: SpongeMat_NVisu(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET                 :: Coords_NVisu(:,:,:,:,:)
REAL,POINTER                            :: SpongeMat_NVisu_p(:,:,:,:,:)
REAL,POINTER                            :: Coords_NVisu_p(:,:,:,:,:)
INTEGER                                 :: nSpongeRamps
INTEGER,ALLOCATABLE                     :: SpongeShape(:)  
REAL,ALLOCATABLE                        :: xStart(:,:)                 ! Starting Point for Sponge Ramp
REAL,ALLOCATABLE                        :: SpVec(:,:)                  ! Vector defining the ramp direction
REAL,ALLOCATABLE                        :: SpDistance(:)               ! Distance of the sponge layer
REAL,ALLOCATABLE                        :: SpRadius(:)                 ! Radius of the cylindrical (3D) / radial (2D) sponge layer
#if(PP_dim==3)
REAL,ALLOCATABLE                        :: SpAxis(:,:)                 ! Axis of the cylindrical sponge layer (only 3D)
#endif
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(A)') '  Initialize Sponge Ramping Function...'

! Precalculation of the sponge strength on the whole domain to determine actual sponge region

nSpongeRamps= CountOption('SpongeShape')
ALLOCATE(SpongeShape(nSpongeRamps))
ALLOCATE(SpDistance(nSpongeRamps))
ALLOCATE(xStart(3,nSpongeRamps))
ALLOCATE(SpVec(3,nSpongeRamps))
ALLOCATE(SpRadius(nSpongeRamps))
#if(PP_dim==3)
ALLOCATE(SpAxis(3,nSpongeRamps))
#endif

DO iRamp=1,nSpongeRamps
  ! readin geometrical parameters of the sponge ramp
  SpongeShape(iRamp)=GETINTFROMSTR('SpongeShape')
!  SpongeShape(iRamp)=GETINT('SpongeShape')
  ! Readin of the sponge Ramp thickness  
  SpDistance(iRamp) = GETREAL('SpongeDistance')
  ! start Sponge Ramp at xStart 
  xStart(:,iRamp)= GETREALARRAY('xStart',3,'(/0.,0.,0./)')
#if PP_dim==2
    IF(xStart(3,iRamp).NE.0) THEN
      CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set xStart(3) = 0!') 
    END IF
#endif
  ! Readin of geometrical parameters for different sponge shapes
  SELECT CASE(SpongeShape(iRamp))
  CASE(SPONGESHAPE_RAMP) ! ramp aligned with a vector
    SpVec(:,iRamp)= GETREALARRAY('SpongeDir',3,'(/1.,0.,0./)')
#if PP_dim==2
    IF(SpVec(3,iRamp).NE.0) THEN
      CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set SpVec(3) = 0!') 
    END IF
#endif
    SpVec(:,iRamp)=SpVec(:,iRamp)/NORM2(SpVec(:,iRamp)) ! Normalize SpVec
  CASE(SPONGESHAPE_CYLINDRICAL) ! circular sponge 
    SpRadius(iRamp)=GETREAL('SpongeRadius')
!    SpRadius(iRamp)=0
#if PP_dim==3
    SpAxis(:,iRamp)=GETREALARRAY('SpongeAxis',3,'(/0.,0.,1./)')
#endif
  END SELECT 
END DO!iRamp

applySponge=.FALSE.
DO iElem=1,nElems
  x_star=0.
  DO iRamp=1,nSpongeRamps
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    SELECT CASE(SpongeShape(iRamp))
      CASE(SPONGESHAPE_RAMP) ! ramp aligned with a vector
      x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      CASE(SPONGESHAPE_CYLINDRICAL) ! cylindrical sponge
      r_vec(:) = Elem_xGP(:,i,j,k,iElem)-xStart(1:PP_dim,iRamp)
#if(PP_dim==3)
      r_vec = r_vec  -SUM((Elem_xGP(:,i,j,k,iElem)-xStart(:,iRamp))*SpAxis(:,iRamp))*SpAxis(:,iRamp)
#endif
      x_star(i,j,k) = (SQRT(SUM(r_vec*r_vec))-SpRadius(iRamp))/SpDistance(iRamp)
    END SELECT
  END DO; END DO; END DO
  IF(ANY(x_star.GT.0.)) THEN
    applySponge(iElem)=.TRUE.
     CYCLE
  END IF
  END DO !iRamp
END DO !iElem=1,nElems

! Get sponge count and build sponge mappings
nSpongeElems=COUNT(applySponge)
ALLOCATE(SpongeMat(0:PP_N,0:PP_N,0:PP_NZ,nSpongeElems))
ALLOCATE(SpongeMap(nSpongeElems))
iSpongeElem=0
DO iElem=1,nElems
  IF(applySponge(iElem))THEN
    iSpongeElem = iSpongeElem+1
    spongeMap(iSpongeElem) = iElem
  END IF
END DO

! Calculate the final sponge strength in the sponge region
SpongeMat=0.
DO iSpongeElem=1,nSpongeElems
  iElem=spongeMap(iSpongeElem)
  sigma=0.
  DO iRamp=1,nSpongeRamps
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    SELECT CASE(SpongeShape(iRamp))
      CASE(SPONGESHAPE_RAMP) ! ramp aligned with a vector
      x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      CASE(SPONGESHAPE_CYLINDRICAL) ! cylindrical sponge
      r_vec(:) = Elem_xGP(:,i,j,k,iElem)-xStart(1:PP_dim,iRamp)
#if(PP_dim==3)
      r_vec = r_vec  -SUM((Elem_xGP(:,i,j,k,iElem)-xStart(:,iRamp))*SpAxis(:,iRamp))*SpAxis(:,iRamp)
#endif
      x_star(i,j,k) = (SQRT(SUM(r_vec*r_vec))-SpRadius(iRamp))/SpDistance(iRamp)
    END SELECT
  END DO; END DO; END DO
  ! Limit to [0,1]
  x_star = MAX(0.,x_star)
  x_star = MIN(1.,x_star)
  ! Sponge Ramping Function ala Babucke
  sigma  = MIN(1.,sigma+6.*x_star**5. - 15.*x_star**4. + 10.*x_star**3.)
  END DO !iRamp
  ! Apply damping factor
  SpongeMat(:,:,:,iSpongeElem) = damping*sigma(:,:,:)
END DO !iSpongeElem=1,nSpongeElems

DEALLOCATE(SpongeShape)
DEALLOCATE(SpDistance)
DEALLOCATE(xStart)
DEALLOCATE(SpVec)
DEALLOCATE(SpRadius)

! Visualize the Sponge Ramp - until now only 3D visualization!
IF(SpongeViz) THEN
  FileString=TRIM(INTSTAMP(TRIM(ProjectName),myRank))//'_SpongeRamp.vtu'
  ALLOCATE(Coords_NVisu(1:3, 0:NVisu,0:NVisu,0:PP_NVisuZ,nElems))
  ALLOCATE(SpongeMat_NVisu(1,0:NVisu,0:NVisu,0:PP_NVisuZ,nElems))
  ALLOCATE(SpDummy(1,0:PP_N,0:PP_N,0:PP_NZ))
  ! Create coordinates of visualization points
  DO iElem=1,nElems
    CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
  END DO
  ! Interpolate solution onto visu grid
  SpongeMat_NVisu=0.
  DO iSpongeElem=1,nSpongeElems
    iElem=spongeMap(iSpongeElem)
    SpDummy(1,:,:,:)=SpongeMat(:,:,:,iSpongeElem)
    CALL ChangeBasisVolume(1,PP_N,NVisu,Vdm_GaussN_NVisu,SpDummy(1:1,:,:,:),SpongeMat_NVisu(1:1,:,:,:,iElem))
  END DO !SpongeElem=1,nSpongeElems
  VarNameSponge(1)='dSponge'
  Coords_NVisu_p => Coords_NVisu
  SpongeMat_NVisu_p => SpongeMat_NVisu
  CALL WriteDataToVTK(1,NVisu,nElems,VarNameSponge,Coords_NVisu_p,SpongeMat_NVisu_p,TRIM(FileString),dim=PP_dim,DGFV=0)
  WRITE(*,*) '********************************************************************'
  DEALLOCATE(Coords_NVisu)
  DEALLOCATE(SpongeMat_NVisu)
  DEALLOCATE(SpDummy)
END IF !SpongeViz

! Finally add the contribution of the Jacobian to SpongeMat (JU_src = (U-UBase)*SpMat)
DO iSpongeElem=1,nSpongeElems
  iElem=spongeMap(iSpongeElem)
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    SpongeMat(i,j,k,iSpongeElem) = SpongeMat(i,j,k,iSpongeElem)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
END DO

END SUBROUTINE CalcSpongeRamp


!==================================================================================================================================
!> \brief Read sponge baseflow from HDF5 file
!> 
!> This routine reads the base flow for the sponge from a .h5 file. This is used for the pruett damping or for the base
!> flow type three. It is checked if the base flow is using the same polynomial degree and node type as the current solution.
!> If not, a interpolation is performed first.
!> Since the base flow is stored on the whole domain, there are no problems if we e.g. change the shape of the sponge region.
!==================================================================================================================================
SUBROUTINE ReadBaseFlow(FileName)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Sponge_Vars
USE MOD_Mesh_Vars  ,       ONLY: offsetElem,nGlobalElems,nElems
USE MOD_HDF5_input ,       ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataProps
USE MOD_ChangeBasisByDim,  ONLY: ChangeBasisVolume
USE MOD_Interpolation,     ONLY: GetVandermonde
USE MOD_Interpolation_Vars,ONLY: NodeType
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=300),INTENT(IN) :: FileName                 !< HDF5 filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem
INTEGER            :: N_Base,nVar_Base,nElems_Base
CHARACTER(LEN=255) :: NodeType_Base
REAL,ALLOCATABLE   :: UTmp(:,:,:,:,:),Vdm_NBase_N(:,:)
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(A,A)')'  Read Sponge Base Flow from file "',FileName
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataProps(nVar_Base,N_Base,nElems_Base,NodeType_Base)

IF(nElems_Base.NE.nGlobalElems)THEN
  CALL abort(__STAMP__,&
             'Baseflow file does not match solution. Elements,nVar',nElems_Base,REAL(nVar_Base))
END IF

! Read in state
IF((N_Base.EQ.PP_N).AND.(TRIM(NodeType_Base).EQ.TRIM(NodeType)))THEN
  ! No interpolation needed, read solution directly from file
  CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_NZ+1,nElems/),OffsetElem,5,RealArray=SpBaseFlow)
ELSE
  ! We need to interpolate the solution to the new computational grid
  SWRITE(UNIT_stdOut,*)'Interpolating base flow from file with N_Base=',N_Base,' to N=',PP_N
  ALLOCATE(UTmp(PP_nVar,0:N_Base,0:N_Base,0:N_Base,nElems))
  ALLOCATE(Vdm_NBase_N(0:N_Base,0:PP_N))
  CALL GetVandermonde(N_Base,NodeType_Base,PP_N,NodeType,Vdm_NBase_N,modal=.TRUE.)
  CALL ReadArray('DG_Solution',5,(/PP_nVar,N_Base+1,N_Base+1,N_Base+1,nElems/),OffsetElem,5,RealArray=UTmp)
  DO iElem=1,nElems
    CALL ChangeBasisVolume(PP_nVar,N_Base,PP_N,Vdm_NBase_N,UTmp(:,:,:,:,iElem),SpBaseFlow(:,:,:,:,iElem))
  END DO
  DEALLOCATE(UTmp,Vdm_NBase_N)
END IF
CALL CloseDataFile()
SWRITE(UNIT_stdOut,*)'DONE READING BASE FLOW!'

END SUBROUTINE ReadBaseFlow


!==================================================================================================================================
!> \brief Apply the sponge to the solution vector (compute contribution to Ut).
!>
!> This routine adds the contribution of the sponge to the time derivatve Ut. 
!> \f$ U_t = U_t - \sigma(x)*(U-U_B) \f$, where \f$ \sigma(x) \f$ is the sponge strength and \f$ U_B \f$ is the base flow.
!> The operation will be performed in the sponge region only using the sponge mapping. The sponge is already pre-multiplied
!> by the Jacobian since we are working in the reference space at this point (at the end of DGTimeDerivative_weakForm).
!==================================================================================================================================
SUBROUTINE Sponge(Ut)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Sponge_Vars ,ONLY: SpongeMap,SpongeMat,SpBaseFlow,nSpongeElems
USE MOD_DG_Vars     ,ONLY: U
USE MOD_Mesh_Vars   ,ONLY: nElems
#if FV_ENABLED
USE MOD_ChangeBasis ,ONLY: ChangeBasis3D
USE MOD_FV_Vars     ,ONLY: FV_Vdm,FV_Elems
USE MOD_Mesh_Vars   ,ONLY: sJ
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iSpongeElem,i,j,k
#if FV_ENABLED    
REAL                :: SpongeMatTmp(1,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: SpongeMat_FV(1,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: SpBaseFlow_FV(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================
DO iSpongeElem=1,nSpongeElems
  iElem=spongeMap(iSpongeElem)
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) THEN ! FV elem     
    ! Remove DG Jacobi from SpongeMat
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      SpongeMatTmp(1,i,j,k) = sJ(i,j,k,iElem,0)*SpongeMat(i,j,k,iSpongeElem)
    END DO; END DO; END DO ! i,j,k
    ! Change Basis of SpongeMat and SpongeBaseFlow to FV grid
    CALL ChangeBasis3D(1,PP_N,PP_N,FV_Vdm,SpongeMatTmp(:,:,:,:),SpongeMat_FV(:,:,:,:))
    CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_Vdm,SpBaseFlow(:,:,:,:,iElem),SpBaseFlow_FV(:,:,:,:))
    ! Calc and add source, take the FV Jacobian into account
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) - SpongeMat_Fv(1,i,j,k)/sJ(i,j,k,iElem,1) * &
                          (U(:,i,j,k,iElem) - SpBaseFlow_FV(:,i,j,k))
    END DO; END DO; END DO ! i,j,k
  ELSE
#endif
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) - SpongeMat(   i,j,k,iSpongeElem) * &
                          (U(:,i,j,k,iElem) - SpBaseFlow(:,i,j,k,iElem))
    END DO; END DO; END DO
#if FV_ENABLED
  END IF
#endif

END DO
END SUBROUTINE Sponge


!==================================================================================================================================
!> Deallocate sponge arrays
!==================================================================================================================================
SUBROUTINE FinalizeSponge()
! MODULES
USE MOD_Sponge_Vars  ,ONLY:SpongeMat,SpongeMap,SpBaseFlow
USE MOD_PruettDamping,ONLY:FinalizePruettDamping
IMPLICIT NONE
!==================================================================================================================================
CALL FinalizePruettDamping()
SDEALLOCATE(SpongeMap)
SDEALLOCATE(SpongeMat)
SDEALLOCATE(SpBaseFlow)
END SUBROUTINE FinalizeSponge

END MODULE MOD_Sponge
