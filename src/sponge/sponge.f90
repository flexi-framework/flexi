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

INTEGER,PARAMETER :: SPONGEBASEFLOW_CONSTANT  = 1
INTEGER,PARAMETER :: SPONGEBASEFLOW_EXACTFUNC = 2
INTEGER,PARAMETER :: SPONGEBASEFLOW_FILE      = 3
INTEGER,PARAMETER :: SPONGEBASEFLOW_PRUETT    = 4

PUBLIC:: DefineParametersSponge
PUBLIC:: InitSponge
PUBLIC:: Sponge
PUBLIC:: FinalizeSponge
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersSponge()
! MODULES
USE MOD_Areas_Vars
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================

CALL prms%SetSection(                'Sponge')
CALL prms%CreateLogicalOption(       'SpongeLayer'     , "Turn on to use sponge regions for reducing reflections at boundaries.",'.FALSE.')
CALL prms%CreateRealOption(          'damping'         , "Damping factor of sponge. U_t=U_t-damping*(U-U_base) in fully damped " //&
                                                         "regions."                                                                &
                                                       , multiple=.TRUE.)
CALL prms%CreateIntFromStringOption( 'SpongeShape'     , "Set shape of sponge: (1) ramp : cartesian / vector-aligned, (2) "      //&
                                                         " cylindrical"                                                            &
                                                       , multiple=.TRUE.)
CALL addStrListEntry('SpongeShape',  'ramp'            , SHAPE_REGION)
CALL addStrListEntry('SpongeShape',  'cuboid'          , SHAPE_CUBOID_CARTESIAN)
CALL addStrListEntry('SpongeShape',  'cylindrical'     , SHAPE_CYLINDRICAL_OUTER)
CALL addStrListEntry('SpongeShape',  'polygon'         , SHAPE_POLYGON)
CALL prms%CreateIntOption(           'nSpongeVertices' , "Define number of vertices per Polygon sponge Zone defining the Polygon"  &
                                                       , multiple=.TRUE.)
CALL prms%CreateRealArrayOption(     'SpongeVertex'    , "Sponge Vertex that defines polygon"                                      &
                                                       , multiple=.TRUE.)
CALL prms%CreateRealOption(          'SpongeDistance'  , "Length of sponge ramp. The sponge will have maximum strength at the "  //&
                                                         "end of the ramp and after that point."                                   &
                                                       , multiple=.TRUE.)
CALL prms%CreateRealArrayOption(     'SpongeXStart'    , "Coordinates of start position of sponge ramp (SpongeShape=ramp) "      //&
                                                          "or center (SpongeShape=cylindrical)."                                   &
                                                       , multiple=.TRUE.)
CALL prms%CreateRealArrayOption(     'SpongeXEnd'      , "Coordinates of second point to define cartesian aligned cube.", multiple=.TRUE.)
CALL prms%CreateRealArrayOption(     'SpongeDir'       , "Direction vector of the sponge ramp (SpongeShape=ramp)", multiple=.TRUE.)
CALL prms%CreateRealOption(          'SpongeRadius'    , "Radius of the sponge zone (SpongeShape=cylindrical)", multiple=.TRUE.)
#if (PP_dim==3)
CALL prms%CreateRealArrayOption(     'SpongeAxis'      , "Axis vector of cylindrical sponge (SpongeShape=cylindrical)", multiple=.TRUE.)
#else
CALL prms%CreateRealArrayOption(     'SpongeXCenter'   , "Center coordinates of cylindrical sponge (SpongeShape=cylindrical)", multiple=.TRUE.)
#endif
CALL prms%CreateLogicalOption(       'SpongeViz'            ,"Turn on to write a visualization file of sponge region and strength."   &
                                                            ,'.FALSE.')
CALL prms%CreateLogicalOption(       'WriteSponge'          ,"Turn on to write the sponge region and strength to the first state "  //&
                                                             "file."     &
                                                            ,'.FALSE.')
CALL prms%CreateIntFromStringOption( 'SpongeBaseFlow'       ,"Type of baseflow to be used for sponge. (1) constant: fixed state,"   //&
                                                             "(2) exactfunction: exact function, (3) file: read baseflow file, "    //&
                                                             "(4) pruett: temporally varying, solution adaptive Pruett baseflow"      &
                                                            ,'1')
CALL addStrListEntry(                'SpongeBaseFlow'       ,'constant',     SPONGEBASEFLOW_CONSTANT)
CALL addStrListEntry(                'SpongeBaseFlow'       ,'exactfunction',SPONGEBASEFLOW_EXACTFUNC)
CALL addStrListEntry(                'SpongeBaseFlow'       ,'file',         SPONGEBASEFLOW_FILE)
CALL addStrListEntry(                'SpongeBaseFlow'       ,'pruett',       SPONGEBASEFLOW_PRUETT)
CALL prms%CreateIntOption(           'SpongeRefState'       ,"Index of refstate in ini-file (SpongeBaseFlow=constant)")
CALL prms%CreateIntOption(           'SpongeExactFunc'      ,"Index of exactfunction (SpongeBaseFlow=exactfunction)")
CALL prms%CreateStringOption(        'SpongeRefFile'        ,"FLEXI solution (e.g. TimeAvg) file from which sponge is read.")
CALL prms%CreateRealOption(          'tempFilterWidthSponge',"Temporal filter width used to advance Pruett baseflow in time.)")
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
USE MOD_BaseFlow          ,ONLY: InitBaseFlow
USE MOD_BaseFlow_Vars     ,ONLY: BaseFlow,doBaseFlow
USE MOD_Equation_Vars     ,ONLY: RefStateCons
USE MOD_Exactfunc         ,ONLY: ExactFunc
USE MOD_Mesh_Vars         ,ONLY: Elem_xGP,nElems
USE MOD_Output_Vars       ,ONLY: ProjectName
USE MOD_ReadInTools       ,ONLY: GETLOGICAL,GETINT,GETINTFROMSTR,GETREAL,GETSTR
USE MOD_Restart_Vars      ,ONLY: DoRestart,RestartTime,RestartFile
USE MOD_Sponge_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k
INTEGER             :: SpongeExactFunc,SpongeRefState
LOGICAL             :: validBaseFlowFile
!==================================================================================================================================

doSponge=GETLOGICAL('SpongeLayer')
IF(.NOT.doSponge) RETURN

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SPONGE...'

SpongeViz      = GETLOGICAL('SpongeViz')
! Readin of BaseFlow parameters
SpBaseFlowType = GETINTFROMSTR('SpongeBaseFlow')

SELECT CASE(SpBaseFlowType)
  CASE(SPONGEBASEFLOW_CONSTANT)  ! Constant baseflow from refstate
    spongeRefState  = GETINT('SpongeRefState')
  CASE(SPONGEBASEFLOW_EXACTFUNC) ! Exact function
    spongeExactFunc = GETINT('SpongeExactFunc')
  CASE(SPONGEBASEFLOW_FILE)      ! Base Flow from .h5 File
    IF (DoRestart) THEN
      SpRefFile  = GETSTR('SpongeRefFile')
      IF (TRIM(SpRefFile).EQ.'none') THEN
        ! If no base flow file has been specified, assume a standard name for the base flow file
        SpRefFile = TRIM(TIMESTAMP(TRIM(ProjectName)//'_BaseFlow',RestartTime))//'.h5'
        ! Check if this file exists
        validBaseFlowFile = FILEEXISTS(SpRefFile)
        IF (.NOT.validBaseFlowFile) THEN
          ! If the assumed base flow file does not exist, use the restart state to initialize the sponge base flow
          SpRefFile = RestartFile
          SWRITE(UNIT_stdOut,'(A)') 'WARNING: No baseflow file found! Using the restart state to initialize sponge base flow.'
        END IF
      ELSE
        ! check if baseflow exists
        validBaseFlowFile = FILEEXISTS(SpRefFile)
        IF (.NOT.validBaseFlowFile) CALL CollectiveStop(__STAMP__,'ERROR: Sponge base flow file '//TRIM(SpRefFile)//' does not exist.')
      END IF
    END IF
  CASE(SPONGEBASEFLOW_PRUETT)    ! Pruett sponge
    IF (.NOT.doBaseFlow) THEN
      CALL PrintWarning('Trying to use Pruett Sponge without enabling BaseFlow!\n'//&
                        'To avoid crash, BaseFlow functionality is now enabled and reinitialized!')
      ! Enable BaseFlow and initialize
      doBaseFlow = .TRUE.
      CALL InitBaseFlow()
    END IF
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,"Undefined SpongeBaseFlow!")
END SELECT

CALL CalcSpongeRamp()

! Preparation of the baseflow on each Gauss Point
SWRITE(UNIT_stdOut,'(A)') '  Initialize Sponge Base Flow...'

SELECT CASE(SpBaseFlowType)
  CASE(SPONGEBASEFLOW_CONSTANT)  ! Constant baseflow from refstate
    ALLOCATE(SpRefState(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

    DO iElem=1,nElems
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        SpRefState(:,i,j,k,iElem) = RefStateCons(:,spongeRefState)
      END DO; END DO; END DO
    END DO
    SpBaseFlow_p(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) => SpRefState

  CASE(SPONGEBASEFLOW_EXACTFUNC) ! Exact function
    ALLOCATE(SpRefState(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

    DO iElem=1,nElems
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        !Save exactFunc state for later use
        CALL ExactFunc(SpongeExactFunc,0.,Elem_xGP(:,i,j,k,iElem),SpRefState(:,i,j,k,iElem))
      END DO; END DO; END DO
    END DO
    SpBaseFlow_p(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) => SpRefState

  CASE(SPONGEBASEFLOW_FILE)      ! Base Flow from .h5 File
    ALLOCATE(SpRefState(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

    CALL ReadBaseFlowSp(SpRefFile)
    SpBaseFlow_p(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) => SpRefState
    ! readin of the hdf5 base flow solution

  CASE(SPONGEBASEFLOW_PRUETT)    ! Pruett: RefState for computation from scratch, Base flow file for restart
    SpBaseFlow_p(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) => BaseFlow
END SELECT

SWRITE(UNIT_stdOut,'(A)')' INIT SPONGE DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitSponge


!==================================================================================================================================
!> \brief Compute sponge shape and strength at each solution point. Visualize sponge.
!>
!> First, depending on the shape (linear or cylindrical), the strength  of the shape without the damping factor (x_star) is
!> calculated on the solution points.
!> From this, a mapping is built which contains only the elements with x_star > 0 somewhere, which is used to later apply
!> the sponge only to regions where it is needed.
!> In this sponge region, the final strength of the sponge is built by limiting x_star to [0,1] and mutiply it by the damping.
!> If set in the parameter file, a visualization of the sponge strength is written as .vtu files.
!> At the end, the sponge is pre-multiplied by the Jacobian since we need to do this anyway when the sponge is applied.
!==================================================================================================================================
SUBROUTINE CalcSpongeRamp()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Areas             ,ONLY:InitArea,PointInPoly
USE MOD_Areas_Vars
USE MOD_BaseFlow_Vars     ,ONLY:TimeFilterWidthBaseFlow
USE MOD_ChangeBasisByDim  ,ONLY:ChangeBasisVolume
USE MOD_Interpolation_Vars,ONLY:NodeTypeCL,NodeType
USE MOD_Interpolation     ,ONLY:GetVandermonde
! USE MOD_IO_HDF5           ,ONLY:AddToFieldData,FieldOut
USE MOD_Mesh_Vars         ,ONLY:sJ,nElems
USE MOD_Mesh_Vars         ,ONLY:Elem_xGP
USE MOD_Output_Vars       ,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_Output_Vars       ,ONLY:ProjectName
USE MOD_ReadInTools
USE MOD_Sponge_Vars
USE MOD_VTK               ,ONLY:WriteDataToVTK
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tArea),POINTER                     :: locSponge
TYPE(tShape),POINTER                    :: locShape
INTEGER                                 :: iElem,iSpongeElem,i,j,k,iRamp,nTimeFilter,nDamping
CHARACTER(LEN=255)                      :: FileString,VarNameSponge(3)
REAL,DIMENSION(  0:PP_N,0:PP_N,0:PP_NZ) :: sigma, x_star
REAL                                    :: SpongeMat_Temp(0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL                                    :: SpongeCount_Temp(1:nElems)
REAL,ALLOCATABLE                        :: SpongeMat_loc(:,:,:,:,:)
REAL                                    :: r_vec(3),tmp1,tmp2
REAL,ALLOCATABLE,TARGET                 :: SpDummy(:,:,:,:)
REAL,ALLOCATABLE,TARGET                 :: SpongeMat_NVisu(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET                 :: Coords_NVisu(:,:,:,:,:)
REAL,POINTER                            :: SpongeMat_NVisu_p(:,:,:,:,:)
REAL,POINTER                            :: Coords_NVisu_p(:,:,:,:,:)
INTEGER,ALLOCATABLE                     :: SpongeShape(:)
REAL,ALLOCATABLE                        :: PruettTimeFilterWidth(:)    ! Filter width of each sponge zone
REAL,ALLOCATABLE                        :: dampingFac(:,:)
INTEGER                                 :: iVertex
REAL                                    :: a(2),b(2),distance
LOGICAL                                 :: applyPolygonSponge
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' | Initialize Sponge Ramping Function...'

! Precalculation of the sponge strength on the whole domain to determine actual sponge region
nSpongeRamps = CountOption('SpongeShape')
ALLOCATE(SpongeShape(   nSpongeRamps))
ALLOCATE(Sponges(       nSpongeRamps))
ALLOCATE(damping(       nSpongeRamps))
ALLOCATE(SpongeDistance(nSpongeRamps))
ALLOCATE(dampingFac(    nSpongeRamps,nElems))

IF (SpBaseFlowType.EQ.SPONGEBASEFLOW_PRUETT) THEN
  ALLOCATE(tempFilterWidthSp(nElems))
  ALLOCATE(PruettTimeFilterWidth(nSpongeRamps))

  nTimeFilter = CountOption('tempFilterWidthSponge')
  IF (nTimeFilter.EQ.1 ) THEN
    PruettTimeFilterWidth = 1./GETREAL("tempFilterWidthSponge")
  ELSE IF (nTimeFilter .EQ. nSpongeRamps ) THEN
    DO iRamp = 1,nSpongeRamps
      PruettTimeFilterWidth(iRamp) = 1./GETREAL("tempFilterWidthSponge")
    END DO
  ELSE
    CALL CollectiveStop(__STAMP__,'Number of Pruett time filter width given does not match number of sponge ramps')
  END IF
  ! Set initial value
  tempFilterWidthSp = 0. !PruettTimeFilterWidth(1)
END IF

nDamping = CountOption('damping')
IF (nDamping .EQ. 1 ) THEN
  damping = GETREAL("damping")
ELSE IF (nDamping .EQ. nSpongeRamps ) THEN
  DO iRamp=1,nSpongeRamps
    damping(iRamp) = GETREAL("damping")
  END DO
ELSE
  CALL CollectiveStop(__STAMP__,'Number of damping factor given does not match number of sponge ramps')
END IF

! Create sponge ramps
DO iRamp=1,nSpongeRamps
  ! readin geometrical parameters of the sponge ramp
  SpongeShape(iRamp) = GETINTFROMSTR('SpongeShape')

  ! Read in Sponge Distance
  SELECT CASE(SpongeShape(iRamp))
  CASE(SHAPE_REGION,SHAPE_CYLINDRICAL_OUTER,SHAPE_CUBOID_CARTESIAN)
    SpongeDistance(iRamp) = GETREAL("SpongeDistance")
  END SELECT

  ! Initialize sponge areas
  CALL InitArea('Sponge',Sponges(iRamp),SpongeShape(iRamp))

  ! Assign damping and time filters
  locSponge => Sponges(iRamp)
  DO iSpongeElem=1,locSponge%nAreaElems
    iElem = locSponge%AreaMap(iSpongeElem)
    IF (SpBaseFlowType.EQ.SPONGEBASEFLOW_PRUETT) THEN
      ! Warning: This is defined per element. Gets overwritten for overlapping sponges!!!
      tempFilterWidthSp(iElem)       = PruettTimeFilterWidth(iRamp)
      TimeFilterWidthBaseFlow(iElem) = PruettTimeFilterWidth(iRamp)
    END IF
    dampingFac(iRamp,iElem) = damping(iRamp)
    CYCLE
  END DO
END DO ! iRamp

! Calculate the sponge strength in every sponge region
ALLOCATE(SpongeMat_loc(nSpongeRamps,0:PP_N,0:PP_N,0:PP_NZ,MAXVAL(Sponges(:)%nAreaElems)))
DO iRamp = 1,nSpongeRamps
  locSponge => Sponges(iRamp)
  locShape  => locSponge%Shape

  DO iSpongeElem = 1,locSponge%nAreaElems
    iElem  = locSponge%AreaMap(iSpongeElem)
    sigma  = 0.
    x_star = 0.
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      SELECT CASE(locSponge%AreaShape)
        CASE(SHAPE_REGION) ! ramp aligned with a vector
          ! Region between xStart and xEnd
          IF (SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-locShape%xStart  (1:PP_dim))   *locShape%Vec(1:PP_dim)).GE.0 .AND. &
              SUM((locShape%xEnd(1:PP_dim)       -Elem_xGP(1:PP_dim,i,j,k,iElem))*locShape%Vec(1:PP_dim)).GE.0) THEN
            x_star(i,j,k) = SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-locShape%xStart(1:PP_dim))*locShape%Vec(1:PP_dim))/SpongeDistance(iRamp)
          END IF

        CASE(SHAPE_CUBOID_CARTESIAN) ! cuboid cartesian aligned defined by two points
          IF(ABS(Elem_xGP(1,i,j,k,iElem)-locShape%xCenter(1)).LT.ABS(locShape%xStart(1)-locShape%xCenter(1))) THEN
            IF(ABS(Elem_xGP(2,i,j,k,iElem)-locShape%xCenter(2)).LT.ABS(locShape%xStart(2)-locShape%xCenter(2))) THEN
              IF(ABS(Elem_xGP(3,i,j,k,iElem)-locShape%xCenter(3)).LT.ABS(locShape%xStart(3)-locShape%xCenter(3))) THEN
                tmp1 = MINVAL((ABS(Elem_xGP(:,i,j,k,iElem)-locShape%xStart(:)))/SpongeDistance(iRamp))
                tmp2 = MINVAL((ABS(Elem_xGP(:,i,j,k,iElem)-locShape%xEnd  (:)))/SpongeDistance(iRamp))
                x_star(i,j,k) = MIN(tmp1,tmp2)
                x_star(i,j,k) = MIN(1.,x_star(i,j,k))
              ELSE
                x_star(i,j,k) =  0.
              END IF
            END IF
          END IF

      CASE(SHAPE_CYLINDRICAL_OUTER) ! cylindrical sponge
      r_vec(:) = Elem_xGP(:,i,j,k,iElem)-locShape%xStart(:)
#if(PP_dim==3)
      r_vec    = r_vec - SUM((Elem_xGP(:,i,j,k,iElem)-locShape%xStart(:))*locShape%Axis(:))*locShape%Axis(:)
#endif
      x_star(i,j,k) = (SQRT(SUM(r_vec*r_vec))-locShape%Radius)/SpongeDistance(iRamp)

      CASE(SHAPE_POLYGON)
        CALL PointInPoly(Elem_xGP(1,i,j,k,iElem),Elem_xGP(2,i,j,k,iElem),locShape%AreaVertex(1:locShape%nAreaVertices,1), &
                         locShape%AreaVertex(1:locShape%nAreaVertices,2),locShape%nAreaVertices,applyPolygonSponge)
        IF(applyPolygonSponge) THEN
          x_star(i,j,k) = 1.0
          DO iVertex=1,locShape%nAreaVertices
            a = Elem_xGP(1:2,i,j,k,iElem)-locShape%AreaVertex(iVertex,1:2)
            b = locShape%AreaVertex(MODULO(iVertex,locShape%nAreaVertices)+1,1:2)-locShape%AreaVertex(iVertex,1:2)
            distance      = ABS((a(1)*b(2)-a(2)*b(1)))/NORM2(b)/locShape%AreaVertex(iVertex,3)
            x_star(i,j,k) = MIN(x_star(i,j,k),distance)
          END DO
        END IF
    END SELECT
  END DO; END DO; END DO

  ! Limit to [0,1]
  x_star = MAX(0.,x_star)
  x_star = MIN(1.,x_star)
  ! Sponge Ramping Function ala Babucke
  sigma  = MIN(1.,sigma+6.*x_star**5. - 15.*x_star**4. + 10.*x_star**3.)
  ! Apply damping factor
  SpongeMat_loc(iRamp,:,:,:,iSpongeElem) = dampingFac(iRamp,iElem)*sigma(:,:,:)
  END DO ! iSpongeElem=1,nSpongeElems
END DO ! iRamp

! Build non reduced Mapping
SpongeCount_Temp = 0
SpongeMat_Temp = 0.
DO iRamp = 1,nSpongeRamps
  locSponge => Sponges(iRamp)
  DO iSpongeElem=1,locSponge%nAreaElems
    iElem = locSponge%AreaMap(iSpongeElem)
    SpongeCount_Temp(iElem)     = SpongeCount_Temp(iElem) + 1
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      SpongeMat_Temp(i,j,k,iElem) = MAX(SpongeMat_Temp(i,j,k,iElem),SpongeMat_loc(iRamp,i,j,k,iSpongeElem))
    END DO; END DO; END DO
  END DO
END DO

! Build up global mapping without duplicate elements
nSpongeElems = 0
DO iElem = 1,nElems
  IF (SpongeCount_Temp(iElem).NE.0) THEN
    nSpongeElems = nSpongeElems + 1
  END IF
END DO

! Allocate global scaling Matrix
ALLOCATE(SpongeMat(0:PP_N,0:PP_N,0:PP_NZ,1:nSpongeElems))
ALLOCATE(SpongeMap(1:nSpongeElems))
iSpongeElem = 1
DO iElem = 1,nElems
  IF (SpongeCount_Temp(iElem).NE.0) THEN
    SpongeMap(iSpongeElem)       = iElem
    SpongeMat(:,:,:,iSpongeElem) = SpongeMat_Temp(:,:,:,iElem)
    iSpongeElem = iSpongeElem +1
  END IF
END DO

SDEALLOCATE(SpongeShape)
SDEALLOCATE(dampingFac)
SDEALLOCATE(PruettTimeFilterWidth)
SDEALLOCATE(SpongeMat_loc)

! Visualize the Sponge Ramp - until now only 3D visualization!
IF(SpongeViz) THEN
  ! Create visu dir, where all vtu files are placed
#if USE_MPI
  IF(nProcessors.GT.1) CALL SYSTEM('mkdir -p visu')
#endif

  FileString=TRIM(ProjectName)//'_SpongeRamp'
  ALLOCATE(Coords_NVisu(1:3, 0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems))
  ALLOCATE(SpongeMat_NVisu(3,0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems))
  ALLOCATE(SpDummy(3,0:PP_N,0:PP_N,0:PP_NZ))
  SpDummy(1,:,:,:) = 0.
  SpDummy(2,:,:,:) = 0.
  SpDummy(3,:,:,:) = 0.

  ! Create coordinates of visualization points
  DO iElem=1,nElems
    CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
  END DO

  ! Interpolate solution onto visu grid
  IF (SpBaseFlowType.EQ.SPONGEBASEFLOW_PRUETT) THEN
    DO iElem=1,nElems
      SpongeMat_NVisu(2,:,:,:,iElem) = tempFilterWidthSp(iElem)
    END DO
  ELSE
    SpongeMat_NVisu(2,:,:,:,:) = 0.
  END IF

  SpongeMat_NVisu(1,:,:,:,:) = 0.
  DO iSpongeElem = 1,nSpongeElems
    iElem = spongeMap(iSpongeElem)
    SpDummy(1,:,:,:) = SpongeMat(:,:,:,iSpongeElem)
    IF (SpBaseFlowType.EQ.SPONGEBASEFLOW_PRUETT) THEN
      SpDummy(2,:,:,:) = tempFilterWidthSp(iElem)
    END IF
    SpDummy(3,:,:,:) = 1.
    CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_GaussN_NVisu,SpDummy(1:3,:,:,:),SpongeMat_NVisu(1:3,:,:,:,iElem))
  END DO ! SpongeElem=1,nSpongeElems

  VarNameSponge(1) = 'dSponge'
  VarNameSponge(2) = 'PruettFilterWidth'
  VarNameSponge(3) = 'SpongeElems'
  Coords_NVisu_p    => Coords_NVisu
  SpongeMat_NVisu_p => SpongeMat_NVisu
  CALL WriteDataToVTK(3,NVisu,nElems,VarNameSponge,Coords_NVisu_p,SpongeMat_NVisu_p,TRIM(FileString),dim=PP_dim)
  DEALLOCATE(Coords_NVisu)
  DEALLOCATE(SpongeMat_NVisu)
  DEALLOCATE(SpDummy)
END IF !SpongeViz

! Write the SpongeMat into the output file
! IF (WriteSponge) THEN
!   ALLOCATE(SpongeMat_Out(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
!   SpongeMat_Out = 0.
!   DO iSpongeElem = 1,nSpongeElems
!     iElem = spongeMap(iSpongeElem)
!     SpongeMat_Out(1,:,:,:,iElem) = SpongeMat(:,:,:,iSpongeElem)
!   END DO
!   CALL AddToFieldData(FieldOut,(/1,PP_N+1,PP_N+1,PP_NZ+1,nElems/),'SpongeMat',(/'SpongeMat'/),RealArray=SpongeMat_Out)
! END IF

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
SUBROUTINE ReadBaseFlowSp(FileName)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasisByDim,  ONLY: ChangeBasisVolume
USE MOD_HDF5_Input,        ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataProps,DatasetExists
USE MOD_Interpolation,     ONLY: GetVandermonde
USE MOD_Interpolation_Vars,ONLY: NodeType
USE MOD_IO_HDF5,           ONLY: File_ID
USE MOD_Mesh_Vars,         ONLY: offsetElem,nGlobalElems,nElems
USE MOD_Sponge_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileName                 !< HDF5 filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: WriteSuccessful
INTEGER            :: iElem
INTEGER            :: N_Base,nVar_Base,nElems_Base
CHARACTER(LEN=255) :: NodeType_Base
REAL,ALLOCATABLE   :: UTmp(:,:,:,:,:),Vdm_NBase_N(:,:)
! Timers
REAL               :: StartT,EndT
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A,A)')' |> Reading sponge base flow from file "',TRIM(FileName)
GETTIME(StartT)

CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Check if restart file was written successfully
CALL DatasetExists(File_ID,'TIME',WriteSuccessful,attrib=.TRUE.)
IF (.NOT.WriteSuccessful) &
  CALL Abort(__STAMP__,'BaseFlow file missing WriteSuccessful marker. Aborting...')

CALL GetDataProps(nVar_Base,N_Base,nElems_Base,NodeType_Base)

IF(nElems_Base.NE.nGlobalElems)THEN
  CALL Abort(__STAMP__,&
             'BaseFlow file does not match solution. Elements,nVar',nElems_Base,REAL(nVar_Base))
END IF

! Read in state
IF((N_Base.EQ.PP_N).AND.(TRIM(NodeType_Base).EQ.TRIM(NodeType)))THEN
  ! No interpolation needed, read solution directly from file
  CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_NZ+1,nElems/),OffsetElem,5,RealArray=SpRefState)
ELSE
  ! We need to interpolate the solution to the new computational grid
  SWRITE(UNIT_stdOut,*)'Interpolating base flow from file with N_Base=',N_Base,' to N=',PP_N
  ALLOCATE(UTmp(PP_nVar,0:N_Base,0:N_Base,0:N_Base,nElems))
  ALLOCATE(Vdm_NBase_N(0:N_Base,0:PP_N))
  CALL GetVandermonde(N_Base,NodeType_Base,PP_N,NodeType,Vdm_NBase_N,modal=.TRUE.)
  CALL ReadArray('DG_Solution',5,(/PP_nVar,N_Base+1,N_Base+1,N_Base+1,nElems/),OffsetElem,5,RealArray=UTmp)
  DO iElem=1,nElems
    CALL ChangeBasisVolume(PP_nVar,N_Base,PP_N,Vdm_NBase_N,UTmp(:,:,:,:,iElem),SpRefState(:,:,:,:,iElem))
  END DO
  DEALLOCATE(UTmp,Vdm_NBase_N)
END IF
CALL CloseDataFile()

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, '|> Reading sponge base flow from file DONE!', DisplayLine=.TRUE.)

END SUBROUTINE ReadBaseFlowSp


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
USE MOD_DG_Vars,           ONLY: U
USE MOD_IO_HDF5,           ONLY: RemoveFromFieldData,FieldOut
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_Sponge_Vars,       ONLY: SpongeMap,SpongeMat,SpBaseFlow_p,nSpongeElems
USE MOD_Sponge_Vars,       ONLY: WriteSponge,SpongeMat_Out
USE MOD_TimeDisc_Vars,     ONLY: iter
#if FV_ENABLED == 1
USE MOD_ChangeBasisByDim,  ONLY: ChangeBasisVolume
USE MOD_FV_Vars,           ONLY: FV_Vdm,FV_Elems
USE MOD_Mesh_Vars,         ONLY: sJ
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iSpongeElem,i,j,k
#if FV_ENABLED == 1
REAL                :: SpongeMatTmp(1,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: SpongeMat_FV(1,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: SpBaseFlow_FV(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================

! Remove the SpongeMat output from all but the first state file
IF (iter.GT.0 .and. WriteSponge) THEN
  CALL RemoveFromFieldData(FieldOut,'SpongeMat')
  SDEALLOCATE(SpongeMat_Out)
END IF

DO iSpongeElem=1,nSpongeElems
  iElem=spongeMap(iSpongeElem)
#if FV_ENABLED == 1
  IF (FV_Elems(iElem).GT.0) THEN ! FV elem
    ! Remove DG Jacobi from SpongeMat
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      SpongeMatTmp(1,i,j,k) = sJ(i,j,k,iElem,0)*SpongeMat(i,j,k,iSpongeElem)
    END DO; END DO; END DO ! i,j,k

    ! Change Basis of SpongeMat and SpongeBaseFlow to FV grid
    CALL ChangeBasisVolume(1      ,PP_N,PP_N,FV_Vdm,SpongeMatTmp(:,:,:,:)      ,SpongeMat_FV( :,:,:,:))
    CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,SpBaseFlow_p(:,:,:,:,iElem),SpBaseFlow_FV(:,:,:,:))

    ! Calc and add source, take the FV Jacobian into account
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) - SpongeMat_Fv(1,i,j,k)/sJ(i,j,k,iElem,1) * &
                          (U(:,i,j,k,iElem) - SpBaseFlow_FV(:,i,j,k))
    END DO; END DO; END DO ! i,j,k

  ELSE
#endif
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) - SpongeMat(     i,j,k,iSpongeElem) * &
                          (U(:,i,j,k,iElem) - SpBaseFlow_p(:,i,j,k,iElem))
    END DO; END DO; END DO
#if FV_ENABLED == 1
  END IF
#endif
END DO

END SUBROUTINE Sponge


!==================================================================================================================================
!> Deallocate sponge arrays
!==================================================================================================================================
SUBROUTINE FinalizeSponge()
! MODULES
USE MOD_Areas
USE MOD_Sponge_Vars      ,ONLY:Sponges,nSpongeRamps
USE MOD_Sponge_Vars      ,ONLY:SpongeMat,SpongeMap,SpRefState
USE MOD_Sponge_Vars      ,ONLY:SpongeMat_Out
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iRamp
!==================================================================================================================================
DO iRamp=1,nSpongeRamps
  CALL FinalizeArea(Sponges(iRamp))
END DO
SDEALLOCATE(Sponges)
SDEALLOCATE(SpongeMap)
SDEALLOCATE(SpongeMat)
SDEALLOCATE(SpongeMat_Out)
SDEALLOCATE(SpRefState)

END SUBROUTINE FinalizeSponge

END MODULE MOD_Sponge
