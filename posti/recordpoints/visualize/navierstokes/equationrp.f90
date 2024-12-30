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
!> Module containing the routines that are used needed to calculate derived quantities from the variables in the RP file.
!===================================================================================================================================
MODULE MOD_EquationRP
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitEquationRP
  MODULE PROCEDURE InitEquationRP
END INTERFACE

INTERFACE CalcEquationRP
  MODULE PROCEDURE CalcEquationRP
END INTERFACE

INTERFACE Plane_BLProps
  MODULE PROCEDURE Plane_BLProps
END INTERFACE

INTERFACE Box_BLProps
  MODULE PROCEDURE Box_BLProps
END INTERFACE

INTERFACE Line_TransformVel
  MODULE PROCEDURE Line_TransformVel
END INTERFACE

INTERFACE Plane_TransformVel
  MODULE PROCEDURE Plane_TransformVel
END INTERFACE

INTERFACE Box_TransformVel
  MODULE PROCEDURE Box_TransformVel
END INTERFACE

INTERFACE FinalizeEquationRP
  MODULE PROCEDURE FinalizeEquationRP
END INTERFACE

PUBLIC::InitEquationRP,CalcEquationRP,Plane_BLProps,Box_BLProps,Line_TransformVel,Plane_TransformVel,Box_TransformVel,FinalizeEquationRP
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the visualization and map the variable names to classify these in conservative and derived quantities.
!===================================================================================================================================
SUBROUTINE InitEquationRP()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars       ,ONLY:VarNames_HDF5,nVar_HDF5
USE MOD_ParametersVisu
USE MOD_EquationRP_Vars
USE MOD_EOS               ,ONLY:InitEOS
USE MOD_EOS_Posti_Vars    ,ONLY:nVarDepEOS,DepTableEOS,DepNames
USE MOD_Readintools       ,ONLY:CountOption,GETSTR
USE MOD_StringTools       ,ONLY:STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: iVar,iVar2,strlen,countCons
INTEGER                               :: mapCand(3)
CHARACTER(LEN=255)                    :: tmp255,tmp255_2
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT EquationRP ...'

justVisualizeState=.FALSE.

! In case no output variables specified: take instead the variable names out of the HDF5 file (used for timeavg-files)
nVarVisu=CountOption("VarName")
IF(nVarVisu .LT. 1) THEN
  justVisualizeState=.TRUE.
  WRITE(*,*) 'No output variables specified, using existing variables from state file.'
  nVarVisu   = nVar_HDF5
  nVarDep=nVar_HDF5
  ALLOCATE(VarNameVisu(nVar_HDF5))
  ALLOCATE(VarNamesAll(nVar_HDF5))
  VarNameVisu = VarNames_HDF5
  VarNamesAll = VarNames_HDF5
!  ALLOCATE(DepTable(nVarVisu,0:nVarVisu))
!  DepTable=0
!  DO iVar=1,nVarVisu
!    DepTable(iVar,iVar)=1
!  END DO
ELSE
  ALLOCATE(VarNameVisu(nVarVisu))
  DO iVar=1,nVarVisu
    VarNameVisu(iVar) = GETSTR("VarName")
  END DO
  !check if Varnames in HDF5 file are conservatives
  ! TODO: also generate mapping in case conservatives are there but not in the correct order
  countCons=0
  IF(nVar_HDF5.GE.PP_nVar)THEN
    DO iVar=1,PP_nVar
      IF (STRICMP(VarNames_HDF5(iVar), DepNames(iVar))) THEN
        countCons=countCons+1
      END IF
    END DO
  END IF
  IF(countCons.NE.PP_nVar) THEN
    CALL CollectiveStop(__STAMP__,'Not all necessary variables are present in HDF5 files')
  END IF
  nVarDep=nVarDepEOS
  ALLOCATE(VarNamesAll(nVarDep))
  VarNamesAll=DepNames
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  DepTable=DepTableEOS

  ! generate mappings
  CALL Build_mapCalc_mapVisu()

  ! initialize EOS
  CALL InitEOS()
END IF


! For local transforms, we need a mapping to all physical vector quantities
! with 3 available components to transform them
ALLOCATE(TransMap(3,INT(nVarVisu/2)))
ALLOCATE(is2D(INT(nVarVisu/2)))
is2D=.FALSE.
TransMap=-1
iVelocity=-1
nVecTrans=0
DO iVar=1,nVarVisu
  tmp255=TRIM(VarNameVisu(iVar))
  strlen=LEN(TRIM(ADJUSTL(tmp255)))
  IF(tmp255(strlen:strlen).EQ.TRIM('X')) THEN
    WRITE(tmp255_2,'(A)')TRIM(tmp255(1:strlen-1))//'Y'
    mapCand(2)=GETMAPBYNAME(TRIM(tmp255_2),VarNameVisu,nVarVisu)
    WRITE(tmp255_2,'(A)')TRIM(tmp255(1:strlen-1))//'Z'
    mapCand(3)=GETMAPBYNAME(TRIM(tmp255_2),VarNameVisu,nVarVisu)
    mapCand(1)=iVar
    IF(.NOT.ANY(mapCand.LE.0)) THEN
      nVecTrans=nVecTrans+1
      TransMap(:,nVecTrans)=mapCand
    ELSEIF((.NOT.ANY(mapCand(1:2).LE.0))) THEN ! 2D Case
      nVecTrans=nVecTrans+1
      TransMap(:,nVecTrans)=mapCand
      is2D(nVecTrans)=.TRUE.
    END IF
  END IF
END DO
WRITE(UNIT_stdOut,'(A)')' Quantities to transform to local coordinate system:'
DO iVar=1,nVecTrans
  DO iVar2=1,3
    WRITE(UNIT_stdOut,'(A,A)')'  ',TRIM(VarNameVisu(TransMap(iVar2,iVar)))
    IF("VelocityX".EQ.TRIM(VarNameVisu(TransMap(iVar2,iVar))))THEN
      iVelocity = iVar
    END IF
  END DO
END DO

IF(Plane_doBLProps.OR.Box_doBLProps) THEN
  IF(iVelocity.EQ.-1) CALL CollectiveStop(__STAMP__,'To calculate BL properties, VelocityX needs to be provided as output variable')
  iDensity=-1
  iDensity=GETMAPBYNAME('Density',VarNameVisu,nVarVisu)
  IF(iDensity.EQ.-1) CALL CollectiveStop(__STAMP__,'To calculate BL properties, Density needs to be provided as output variable')
  iPressure=-1
  iPressure=GETMAPBYNAME('Pressure',VarNameVisu,nVarVisu)
  IF(iPressure.GT.0) THEN
    nBLProps=12
  ELSE
    nBLProps=11
    WRITE(UNIT_stdOut,'(A)') '!!! No Pressure data provided, cp is not computed !!!'
  END IF
  ALLOCATE(VarNames_BLProps(nBLProps))
  VarNames_BLProps(1)='delta99'
  VarNames_BLProps(2)='u_delta'
  VarNames_BLProps(3)='delta1'
  VarNames_BLProps(4)='theta'
  VarNames_BLProps(5)='H12'
  VarNames_BLProps(6)='Re_delta'
  VarNames_BLProps(7)='Re_theta'
  VarNames_BLProps(8)='u_reverse'
  VarNames_BLProps(9)='tau_w'
  VarNames_BLProps(10)='Re_tau'
  VarNames_BLProps(11)='u_tau'
  IF(iPressure.GT.0) VarNames_BLProps(12)='c_p'
END IF

WRITE(UNIT_stdOut,'(A)')' INIT EquationRP DONE!'
WRITE(UNIT_stdOut,'(132("-"))')
EquationRPInitIsDone=.TRUE.
END SUBROUTINE InitEquationRP



!===================================================================================================================================
!> This routine computes the state on the visualization grid
!===================================================================================================================================
SUBROUTINE CalcEquationRP()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars        ,ONLY: RPData
USE MOD_RPData_Vars        ,ONLY: nVar_HDF5,VarNames_HDF5
USE MOD_RPSetVisuVisu_Vars ,ONLY: nRP_global
USE MOD_OutputRPVisu_Vars  ,ONLY: nSamples_out
USE MOD_ParametersVisu     ,ONLY: nVarDep,nVarCalc,mapCalc,mapVisu,VarNamesAll,justVisualizeState
USE MOD_ParametersVisu     ,ONLY: Line_LocalVel,Plane_LocalVel
USE MOD_OutputRPVisu_Vars  ,ONLY: RPData_out
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_StringTools        ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nVal(2)
INTEGER            :: iVarOut,iVarIn,iVar,iVarCalc,iVarVisu
REAL,ALLOCATABLE   :: UCalc(:,:,:)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')" CONVERT DERIVED QUANTITIES..."
! CALCULATE DERIVED QUATITIES -----------------------------------------------------------------------------------------------------!
IF(justVisualizeState)THEN
  RPData_out=RPData
ELSE
  maskCalc=1
  nVal=(/nRP_global,nSamples_out/)
  ! Copy existing variables from solution array
  ! Attention: nVarCalc must be last dimension (needed for CalcQuantities from flexilib!)
  ALLOCATE(UCalc(nRP_global,nSamples_out,nVarCalc))

  DO iVarOut=1,nVarDep ! iterate over all out variables
    IF (mapCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
    DO iVarIn=1,nVar_HDF5 ! iterate over all in variables
      IF( STRICMP(VarNamesAll(iVarOut),VarNames_HDF5(iVarIn))) THEN
        UCalc(:,:,mapCalc(iVarOut))=RPData(iVarIn,:,:)
        maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
      END IF
    END DO
  END DO

  ! calculate all quantities
  CALL CalcQuantities(nVarCalc,nVal,(/1/),mapCalc,UCalc,maskCalc)

  ! fill output array
  DO iVar=1,nVarDep
    IF (mapVisu(iVar).GT.0) THEN
      iVarCalc = mapCalc(iVar)
      iVarVisu = mapVisu(iVar)
      RPData_out(iVarVisu,:,:)=UCalc(:,:,iVarCalc)
    END IF
  END DO
  DEALLOCATE(UCalc)
END IF

! Coordinate Transform
IF(Line_LocalVel) &
  CALL Line_TransformVel(RPData_out,nSamples_out)
IF(Plane_LocalVel) &
  CALL Plane_TransformVel(RPData_out,nSamples_out)

WRITE(UNIT_stdOut,'(A)')" CONVERT DERIVED QUANTITIES DONE!"
END SUBROUTINE CalcEquationRP



!===================================================================================================================================
!> This routine transformes velocities to the line-local coordinate system
!===================================================================================================================================
SUBROUTINE Line_TransformVel(RPData,nSamples)
! MODULES
USE MOD_Globals
USE MOD_ParametersVisu     ,ONLY: nVarVisu
USE MOD_RPSetVisuVisu_Vars ,ONLY: nLines,Lines,tLine,nRP_global
USE MOD_EquationRP_Vars    ,ONLY: nVecTrans,TransMap
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,   INTENT(INOUT)      :: RPData(nVarVisu,nRP_global,nSamples)
INTEGER,INTENT(IN)         :: nSamples
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iLine,iRP,iSample,iVec
TYPE(tLine),POINTER :: aLine
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')" coordinate transform of velocity along lines.."
DO iLine=1,nLines
  aLine=>Lines(iLine)
  DO iRP=1,aLine%nRP
    DO iSample=1,nSamples
      DO iVec=1,nVecTrans
        RPData(TransMap(1:3,iVec),aLine%IDlist(iRP),iSample)=&
           MATMUL(aLine%Tmat,RPData(TransMap(1:3,iVec),aLine%IDlist(iRP),iSample))
      END DO
    END DO !iSample
  END DO !iRP
END DO !iLine
WRITE(UNIT_stdOut,'(A)')" done!"
END SUBROUTINE Line_TransformVel


!===================================================================================================================================
!> This routine transformes velocities to the plane-local coordinate system
!===================================================================================================================================
SUBROUTINE Plane_TransformVel(RPData,nSamples)
! MODULES
USE MOD_Globals
USE MOD_Mathtools          ,ONLY:CROSS
USE MOD_ParametersVisu     ,ONLY:nVarVisu
USE MOD_RPSetVisuVisu_Vars ,ONLY:nPlanes,Planes,tPlane,nRP_global
USE MOD_EquationRP_Vars    ,ONLY:nVecTrans,TransMap,is2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,   INTENT(INOUT)      :: RPData(nVarVisu,nRP_global,nSamples)
INTEGER,INTENT(IN)         :: nSamples
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iPlane,i,j,iSample,iVec
REAL                 :: Tmat(3,3)
TYPE(tPlane),POINTER :: Plane
!===================================================================================================================================

WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Coordinate transform of velocity along Plane coordinates..."
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF(Plane%Type.EQ.2) THEN ! BLPlane
    DO i=1,Plane%nRP(1)
      ! T1= tangential vector
      Tmat(1,1:3) = Plane%TangVec(:,i)
      ! T2= normalized Line Vector
      Tmat(2,1:3) = Plane%NormVec(:,i)
      ! T3= crossprod to guarantee right hand system
      Tmat(3,1:3) = CROSS(Tmat(1,1:3),Tmat(2,1:3))
      Tmat(3,1:3) = Tmat(3,:)/NORM2(Tmat(3,:))
      DO iVec=1,nVecTrans
        IF(.NOT.is2D(iVec))THEN
          DO j=1,Plane%nRP(2)
            DO iSample=1,nSamples
              RPData(TransMap(1:3,iVec),Plane%IDlist(i,j),iSample)=&
                 MATMUL(Tmat,RPData(TransMap(1:3,iVec),Plane%IDlist(i,j),iSample))
            END DO !iSample
          END DO !j
        ELSE ! 2D case
          DO j=1,Plane%nRP(2)
            DO iSample=1,nSamples
              RPData(TransMap(1:2,iVec),Plane%IDlist(i,j),iSample)=&
                 MATMUL(Tmat(1:2,1:2),RPData(TransMap(1:2,iVec),Plane%IDlist(i,j),iSample))
            END DO !iSample
          END DO !j
        END IF !.NOT.is2D
      END DO ! iVec
    END DO !i
  END IF!(Plane%Type.EQ.2) THEN ! BLPlane
END DO !iPlane
WRITE(UNIT_stdOut,'(A)')" done!"

END SUBROUTINE Plane_TransformVel

!===================================================================================================================================
!> This routine transformes velocities to the Box-local coordinate system
!===================================================================================================================================
SUBROUTINE Box_TransformVel(RPData,nSamples)
! MODULES
USE MOD_Globals
USE MOD_Mathtools          ,ONLY:CROSS
USE MOD_ParametersVisu     ,ONLY:nVarVisu
USE MOD_RPSetVisuVisu_Vars ,ONLY:nBoxes,Boxes,tBox,nRP_global
USE MOD_EquationRP_Vars    ,ONLY:nVecTrans,TransMap,is2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,   INTENT(INOUT)      :: RPData(nVarVisu,nRP_global,nSamples)
INTEGER,INTENT(IN)         :: nSamples
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iBox,i,j,k,iSample,iVec
REAL                 :: Tmat(3,3)
TYPE(tBox),POINTER   :: Box
!===================================================================================================================================

WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Coordinate transform of velocity along Box coordinates..."
DO iBox=1,nBoxes
  Box=>Boxes(iBox)
  IF(Box%Type.EQ.1) THEN ! BLBox
    DO i=1,Box%nRP(1)
      DO k=1,Box%nRP(3)
        ! T1= tangential vector
        Tmat(1,1:3) = Box%TangVec(:,i,k)
        ! T2= normalized Line Vector
        Tmat(2,1:3) = Box%NormVec(:,i,k)
        ! T3= crossprod to guarantee right hand system
        Tmat(3,1:3) = CROSS(Tmat(1,1:3),Tmat(2,1:3))
        Tmat(3,1:3) = Tmat(3,:)/NORM2(Tmat(3,:))
        DO iVec=1,nVecTrans
          IF(.NOT.is2D(iVec))THEN
            DO j=1,Box%nRP(2)
              DO iSample=1,nSamples
                RPData(TransMap(1:3,iVec),Box%IDlist(i,j,k),iSample)=&
                      MATMUL(Tmat,RPData(TransMap(1:3,iVec),Box%IDlist(i,j,k),iSample))
              END DO !iSample
            END DO !j
          ELSE ! 2D case
            DO j=1,Box%nRP(2)
              DO iSample=1,nSamples
                RPData(TransMap(1:2,iVec),Box%IDlist(i,j,k),iSample)=&
                    MATMUL(Tmat(1:2,1:2),RPData(TransMap(1:2,iVec),Box%IDlist(i,j,k),iSample))
              END DO !iSample
            END DO !j
          END IF !.NOT.is2D
        END DO ! iVec
      END DO !k
    END DO !i
  END IF!(Box%Type.EQ.2) THEN ! BLBox
END DO !iBox
WRITE(UNIT_stdOut,'(A)')" done!"

END SUBROUTINE Box_TransformVel

!===================================================================================================================================
!> This routine calculates the boundary layer specific quantities on all boundary layer planes.
!===================================================================================================================================
SUBROUTINE Plane_BLProps()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars ,ONLY: nPlanes,Planes,tPlane
USE MOD_EquationRP_Vars    ,ONLY: nBLProps
USE MOD_ParametersVisu     ,ONLY: Plane_BLvelScaling
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iPlane,i
TYPE(tPlane),POINTER :: Plane
!===================================================================================================================================

WRITE(UNIT_stdOut,'(A)')" Calculate integral boundary layer properties for planes"
SELECT CASE(Plane_BLvelScaling)
  CASE(0) ! do nothing
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Velocity profiles are not scaled."
  CASE(1) ! laminar scaling
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Velocity profiles are scaled with boundary layer edge velocity, PlaneY with BL thickness."
  CASE(2) ! turbulent scaling
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Velocity profiles and PlaneY are scaled with friction velocity (= u+ over y+)."
END SELECT

DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF(Plane%Type.EQ.2) THEN ! BLPlane
    ALLOCATE(Plane%BLProps(nBLProps,1:Plane%nRP(1)))
    DO i=1,Plane%nRP(1)
      CALL Calc_BLProps(Plane%LocalCoord(:,i,:),Plane%IDlist(i,:),Plane%nRP(2),Plane_BLvelScaling,Plane%BLProps(:,i))
    END DO
  END IF!(Plane%Type.EQ.2) THEN ! BLPlane
END DO !iPlane
WRITE(UNIT_stdOut,'(A)')" done!"

END SUBROUTINE Plane_BLProps

!===================================================================================================================================
!> This routine calculates the boundary layer specific quantities on all boundary layer boxes.
!===================================================================================================================================
SUBROUTINE Box_BLProps()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars ,ONLY: nBoxes,Boxes,tBox
USE MOD_EquationRP_Vars    ,ONLY: nBLProps
USE MOD_ParametersVisu     ,ONLY: Box_BLvelScaling
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iBox,i,k
TYPE(tBox),POINTER   :: Box
!===================================================================================================================================

WRITE(UNIT_stdOut,'(A)')" Calculate integral boundary layer properties for boxes"
SELECT CASE(Box_BLvelScaling)
  CASE(0) ! do nothing
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Velocity profiles are not scaled."
  CASE(1) ! laminar scaling
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Velocity profiles are scaled with boundary layer edge velocity, BoxY with BL thickness."
  CASE(2) ! turbulent scaling
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" Velocity profiles and BoxY are scaled with friction velocity (= u+ over y+)."
END SELECT

DO iBox=1,nBoxes
  Box=>Boxes(iBox)
  IF(Box%Type.EQ.1) THEN ! BLBox
    ALLOCATE(Box%BLProps(nBLProps,1:Box%nRP(1),1:Box%nRP(3)))
    DO k=1,Box%nRP(3)
      DO i=1,Box%nRP(1)
        CALL Calc_BLProps(Box%LocalCoord(1:2,i,:,k),Box%IDlist(i,:,k),Box%nRP(2),Box_BLvelScaling,Box%BLProps(:,i,k))
      END DO
    END DO
  END IF!(Box%Type.EQ.2) THEN ! BLBox
END DO !iBox
WRITE(UNIT_stdOut,'(A)')" done!"

END SUBROUTINE Box_BLProps

SUBROUTINE Calc_BLProps(LocalCoord,IDlist,nRP,Scaling,BLProps)
! MODULES
USE MOD_Globals
USE MOD_OutputRPVisu_Vars  ,ONLY: RPDataTimeAvg_out
USE MOD_EquationRP_Vars    ,ONLY: TransMap,pInf,uInf,rhoInf,iPressure,iVelocity,is2D,nBLProps
USE MOD_ParametersVisu     ,ONLY: Mu0
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)        :: LocalCoord(2,nRP)
INTEGER,INTENT(IN)        :: IDlist(nRP)
INTEGER                   :: nRP
INTEGER                   :: Scaling
REAL,INTENT(OUT)          :: BLProps(nBLProps)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: j,j_max,ndim
REAL                 :: u_max,y_max,diffu1,diffu2,y1,y2,rho1,rho2,rho_delta,dudy,dudy1
REAL                 :: u1,u2,u3,dy,dy2
REAL                 :: u_loc_bledge
REAL                 :: u_delta,delta99,delta1,theta,u_r,tau_W,u_tau
REAL                 :: u_loc(nRP),y_loc(nRP),u_star(nRP)
!==================================================================================================================================

! Initialize
BLProps = 0.

! get wall friction tau_W = mu*(du/dy+dv/dx)
dy =LocalCoord(2,2)-LocalCoord(2,1)
dy2=LocalCoord(2,3)-LocalCoord(2,2)
u1=RPDataTimeAvg_out(TransMap(1,iVelocity),IDlist(1))
u2=RPDataTimeAvg_out(TransMap(1,iVelocity),IDlist(2))
u3=RPDataTimeAvg_out(TransMap(1,iVelocity),IDlist(3))
dudy=(u2-u1)/dy- (dy*(u3-u2)-dy2*(u2-u1))/((dy+dy2)*dy2)
tau_W = mu0*dudy
! get delta99
! Kloker method: integrate vorticity to get u_star.
u_star(1)=RPDataTimeAvg_out(TransMap(1,iVelocity),IDlist(1)) ! Set to wall velocity <- slip wall
DO j=2,nRP-1
  !calculate vorticity using central FD O2
  dy2=LocalCoord(2,j+1)-LocalCoord(2,j-1)
  u2=RPDataTimeAvg_out(TransMap(1,1),IDlist(j+1))
  u1=RPDataTimeAvg_out(TransMap(1,1),IDlist(j-1))
  dudy1=(u2-u1)/dy2
  dy=LocalCoord(2,j)-LocalCoord(2,j-1)
  ! integration with trapezoidal rule O2
  u_star(j)=u_star(j-1)+0.5*(dudy+dudy1)*dy
  dudy=dudy1
END DO
! find point with u_star=0.99*u_starmax
! define this as delta99 and u_delta
u_max=u_star(nRP-1)
DO j=2,nRP-1 ! loop starting from wall (j=1)
  diffu1=(u_max-u_star(j))/u_max
  IF(diffu1.LT.0.01) THEN
    y1=LocalCoord(2,j)
    y2=LocalCoord(2,j-1)
    u1=u_star(j)
    u2=u_star(j-1)
    rho1=RPDataTimeAvg_out(1,IDlist(j))
    rho2=RPDataTimeAvg_out(1,IDlist(j-1))
    diffu2=(u_max-u_star(j-1))/u_max
    ! interpolate linearly between two closest points around u/u_int=0.99 to get
    ! delta99 and u_delta
    delta99=y1+(y2-y1)/(diffu2-diffu1)*(0.01-diffu1)   ! delta99
    u_delta=u_max
    rho_delta=rho1+(rho2-rho1)/(diffu2-diffu1)*(0.01-diffu1)   ! rho_delta
    j_max=j
    y_max=y1
    EXIT
  END IF
END DO !j=1,RP
! scale the velocity with the local u_delta,wall distance with delta99
DO j=1,nRP
  u_loc(j)= &
          RPDataTimeAvg_out(TransMap(1,1),IDlist(j))/u_delta
  y_loc(j)=LocalCoord(2,j)/delta99
END DO! j=1,nRP
! calculate integral BL properties (trapezoidal rule)
delta1=0.
theta=0.
DO j=2,j_max
  IF (j.LT.j_max) THEN
    dy=y_loc(j)-y_loc(j-1)
    ! displacement thickness
    delta1=delta1+ 0.5*dy*(1.-u_loc(j) +1.-u_loc(j-1))
    ! momentum thickness
    theta=theta+ 0.5*dy*(u_loc(j)*(1.-u_loc(j)) +u_loc(j-1)*(1.-u_loc(j-1)))
  ELSE !Interpolate y_loc and u_loc at at the boundarylayer edge to be consitents
    !y_loc already scaled with delta99: y_loc at BL edge = 1
    dy=1-y_loc(j-1)
    u_loc_bledge=(u_loc(j)-u_loc(j-1))/(y_loc(j)-y_loc(j-1))*(1-y_loc(j-1))+u_loc(j-1)
    ! displacement thickness
    delta1=delta1+ 0.5*dy*(1.-u_loc_bledge +1.-u_loc(j-1))
    ! momentum thickness
    theta=theta+ 0.5*dy*(u_loc_bledge*(1.-u_loc_bledge) +u_loc(j-1)*(1.-u_loc(j-1)))
  END IF !j.LT.j_max
END DO
! get max. reverse flow if any
u_r=MIN(MINVAL(u_loc(:)),0.)
delta1=delta1*delta99 ! integration has been performed in non-dimensional units
theta=theta*delta99   !
! write to solution array
BLProps(1)=delta99
BLProps(2)=u_delta
BLProps(3)=delta1
BLProps(4)=theta
BLProps(5)=delta1/theta  !H12
BLProps(6)=rho_delta*u_delta*delta1/Mu0                  !Re_delta
BLProps(7)=rho_delta*u_delta*theta/Mu0                   !Re_theta
BLProps(8)=ABS(u_r)                                      !u_reverse
BLProps(9)=tau_W                                         !tau_W
BLProps(10)=rho_delta*SQRT(ABS(tau_W)/rho_delta)*delta99/Mu0  !Re_tau
BLProps(11)=SQRT(ABS(tau_W)/rho_delta)                   !u_tau
IF(iPressure.GT.0) BLProps(12)=(RPDataTimeAvg_out(iPressure,IDlist(1))-pInf)/(0.5*rhoInf*uInf**2) !c_p

! perform scaling if required
ndim = MERGE(3,2,.NOT.is2D(1))
SELECT CASE(Scaling)
CASE(0) ! do nothing
CASE(1) ! laminar scaling
  DO j=1,nRP
    RPDataTimeAvg_out(TransMap(1:ndim,1),IDlist(j)) = &
      RPDataTimeAvg_out(TransMap(1:ndim,1),IDlist(j))/u_delta
    LocalCoord(2,j)=LocalCoord(2,j)/delta99
  END DO
CASE(2) ! turbulent scaling
  u_tau=SQRT(ABS(tau_w)/rho_delta)
  DO j=1,nRP
    RPDataTimeAvg_out(TransMap(1:ndim,1),IDlist(j)) = &
      RPDataTimeAvg_out(TransMap(1:ndim,1),IDlist(j))/u_tau
    LocalCoord(2,j)=LocalCoord(2,j)*u_tau*rho_delta/mu0
  END DO
CASE(3) ! testing
  DO j=1,nRP
    RPDataTimeAvg_out(TransMap(1:ndim,1),IDlist(j)) = &
      u_star(j)/u_max
  END DO
END SELECT

END SUBROUTINE Calc_BLProps


!===================================================================================================================================
!> This function returns the index of the variable "VarName" in the array of variables "VarNameList"
!===================================================================================================================================
FUNCTION GETMAPBYNAME(VarName,VarNameList,nVarList)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: VarName,VarNameList(nVarList)
INTEGER,INTENT(IN)             :: nVarList
INTEGER                        :: GETMAPBYNAME
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i
!===================================================================================================================================
GETMAPBYNAME=-1
DO i=1,nVarList
  IF(TRIM(VarName).EQ.TRIM(VarNameList(i)))THEN
    GETMAPBYNAME=i
    RETURN
  END IF
END DO
RETURN
end function GETMAPBYNAME


!===================================================================================================================================
!> Deallocate the global variables
!===================================================================================================================================
SUBROUTINE FinalizeEquationRP()
! MODULES
USE MOD_Globals
USE MOD_EquationRP_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(TransMap)
SDEALLOCATE(is2D)
WRITE(UNIT_stdOut,'(A)') ' EquationRP FINALIZED'

END SUBROUTINE FinalizeEquationRP

END MODULE MOD_EquationRP
