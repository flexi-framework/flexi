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

INTERFACE Line_TransformVel
  MODULE PROCEDURE Line_TransformVel
END INTERFACE

INTERFACE Plane_TransformVel
  MODULE PROCEDURE Plane_TransformVel
END INTERFACE

INTERFACE FinalizeEquationRP
  MODULE PROCEDURE FinalizeEquationRP
END INTERFACE

PUBLIC::InitEquationRP,CalcEquationRP,Plane_BLProps,Line_TransformVel,Plane_TransformVel,FinalizeEquationRP
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
WRITE(UNIT_StdOut,'(132("-"))')
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
WRITE(UNIT_StdOut,'(A)')' Quantities to transform to local coordinate system:'
DO iVar=1,nVecTrans
  DO iVar2=1,3
    WRITE(UNIT_StdOut,'(A,A)')'  ',TRIM(VarNameVisu(TransMap(iVar2,iVar)))
  END DO
END DO

IF(Plane_doBLProps) THEN
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
WRITE(UNIT_StdOut,'(132("-"))')
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
WRITE(UNIT_StdOut,'(132("-"))')
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
WRITE(UNIT_stdOut,'(A)')" coordinate transform of velocity along Plane coordinates.."
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
!> This routine calculates the boundary layer specific quantities on all boundary layer planes.
!===================================================================================================================================
SUBROUTINE Plane_BLProps()
! MODULES
USE MOD_Globals
USE MOD_OutputRPVisu_Vars  ,ONLY: RPDataTimeAvg_out
USE MOD_RPSetVisuVisu_Vars ,ONLY: nPlanes,Planes,tPlane
USE MOD_EquationRP_Vars    ,ONLY: is2D,TransMap,nBLProps,pInf,uInf,rhoInf,iPressure
USE MOD_ParametersVisu     ,ONLY: Plane_BLvelScaling
USE MOD_ParametersVisu     ,ONLY: Mu0
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iPlane,i,j,jj,j_max,ndim
REAL                 :: u_max,y_max,diffu1,diffu2,y1,y2,u1,u2,u3,rho1,rho2,rho_delta,dudy,dudy1
REAL                 :: dy,dy2,u_loc_bledge
REAL                 :: u_delta,delta99,delta1,theta,u_r,tau_W,u_tau
REAL,ALLOCATABLE     :: u_loc(:),y_loc(:),u_star(:)
TYPE(tPlane),POINTER :: Plane
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')" Calculate integral boundary layer properties"
SELECT CASE(Plane_BLvelScaling)
CASE(0) ! do nothing
WRITE(UNIT_stdOut,'(A)')" Velocity profiles are not scaled."
CASE(1) ! laminar scaling
WRITE(UNIT_stdOut,'(A)')" Velocity profiles are scaled with boundary layer edge velocity, PlaneY with BL thickness."
CASE(2) ! turbulent scaling
WRITE(UNIT_stdOut,'(A)')" Velocity profiles and PlaneY are scaled with friction velocity (= u+ over y+)."
END SELECT
ndim=3
IF(is2D(1)) ndim=2
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF(Plane%Type.EQ.2) THEN ! BLPlane
    ALLOCATE(Plane%BLProps(nBLProps,1:Plane%nRP(1)))
    ALLOCATE(u_loc(Plane%nRP(2)),y_loc(Plane%nRP(2)),u_star(Plane%nRP(2)))
    DO i=1,Plane%nRP(1)
      ! get wall friction tau_W = mu*(du/dy+dv/dx)
      dy=Plane%LocalCoord(2,i,2)-Plane%LocalCoord(2,i,1)
      dy2=Plane%LocalCoord(2,i,3)-Plane%LocalCoord(2,i,2)
      u1=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,1))
      u2=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,2))
      u3=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,3))
      dudy=(u2-u1)/dy- (dy*(u3-u2)-dy2*(u2-u1))/((dy+dy2)*dy2)
      tau_W = mu0*dudy
      ! get delta99
      ! Kloker method: integrate vorticity to get u_star.
      jj=Plane%nRP(2)
      u_star(1)=0.
      DO j=2,jj-1
        !calculate vorticity using central FD O2
        dy2=Plane%LocalCoord(2,i,j+1)-Plane%LocalCoord(2,i,j-1)
        u2=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j+1))
        u1=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j-1))
        dudy1=(u2-u1)/dy2
        dy=Plane%LocalCoord(2,i,j)-Plane%LocalCoord(2,i,j-1)
        ! integration with trapezoidal rule O2
        u_star(j)=u_star(j-1)+0.5*(dudy+dudy1)*dy
        dudy=dudy1
      END DO
      ! find point with u_star=0.99*u_starmax
      ! define this as delta99 and u_delta
      u_max=u_star(jj-1)
      DO j=2,jj-1 ! loop starting from wall (j=1)
        diffu1=(u_max-u_star(j))/u_max
        IF(diffu1.LT.0.01) THEN
          y1=Plane%LocalCoord(2,i,j)
          y2=Plane%LocalCoord(2,i,j-1)
          u1=u_star(j)
          u2=u_star(j-1)
          rho1=RPDataTimeAvg_out(1,Plane%IDlist(i,j))
          rho2=RPDataTimeAvg_out(1,Plane%IDlist(i,j-1))
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
      END DO !j=1,Plane%nRP(2)
!      jj=Plane%nRP(2)
!      uu=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,jj))
!      yy=Plane%LocalCoord(2,i,jj)
!      ! find local maximum in VelocityX
!      u_max=0.
!      DO j=1,jj
!       IF(RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j)).GT.u_max) THEN
!          u_max=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j))
!          y_max=Plane%LocalCoord(2,i,j)
!          j_max=j
!        END IF
!      END DO !j=1,Plane%nRP(2)
!      ! find point with u=0.995*u_max ala Wuerz
!      ! define this as delta99 and u_delta
!      DO j=1,jj ! loop starting from wall (j=1)
!        IF(y_max.LT.yy) THEN
!          ! interpolate linearly between uppermost point and maximum
!          u_int=uu + (u_max-uu)/(y_max-yy)*(Plane%LocalCoord(2,i,j)-yy)
!        ELSE
!          ! if the maximum is at the edge of the discretization, use a vertical line
!          u_int=u_max
!        END IF
!        diffu1=(u_int-RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j)))/u_max
!        IF(diffu1.LT.0.01) THEN
!          y1=Plane%LocalCoord(2,i,j)
!          y2=Plane%LocalCoord(2,i,j-1)
!          u1=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j))
!          u2=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j-1))
!          rho1=RPDataTimeAvg_out(1,Plane%IDlist(i,j))
!          rho2=RPDataTimeAvg_out(1,Plane%IDlist(i,j-1))
!          diffu2=(u_int-RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j-1)))/u_max
!          ! interpolate linearly between two closest points around u/u_int=0.99 to get
!          ! delta99 and u_delta
!          delta99=y1+(y2-y1)/(diffu2-diffu1)*(0.01-diffu1)   ! delta99
!!          u_delta=u1+(u2-u1)/(diffu2-diffu1)*(0.01-diffu1)   ! u_delta
!          u_delta=u_max
!          rho_delta=rho1+(rho2-rho1)/(diffu2-diffu1)*(0.01-diffu1)   ! rho_delta
!!          u_delta(i)=RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j))
!!          delta99(i)=Plane%LocalCoord(2,i,j)
!          EXIT
!        END IF
!      END DO !j=1,Plane%nRP(2)
      ! scale the velocity with the local u_delta,wall distance with delta99
      DO j=1,jj
        u_loc(j)= &
                RPDataTimeAvg_out(TransMap(1,1),Plane%IDlist(i,j))/u_delta
        y_loc(j)=Plane%LocalCoord(2,i,j)/delta99
      END DO! j=1,jj
      ! calculate integral BL properties (trapezoidal rule)
      delta1=0.
      theta=0.
!      j_max=jj
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
      Plane%BLProps(1,i)=delta99
      Plane%BLProps(2,i)=u_delta
      Plane%BLProps(3,i)=delta1
      Plane%BLProps(4,i)=theta
      Plane%BLProps(5,i)=delta1/theta  !H12
      Plane%BLProps(6,i)=rho_delta*u_delta*delta1/Mu0                  !Re_delta
      Plane%BLProps(7,i)=rho_delta*u_delta*theta/Mu0                   !Re_theta
      Plane%BLProps(8,i)=ABS(u_r)                                      !u_reverse
      Plane%BLProps(9,i)=tau_W                                         !tau_W
      Plane%BLProps(10,i)=rho_delta*SQRT(ABS(tau_W)/rho_delta)*delta99/Mu0  !Re_tau
      Plane%BLProps(11,i)=SQRT(ABS(tau_W)/rho_delta)                   !u_tau
      IF(iPressure.GT.0)  Plane%BLProps(12,i)=(RPDataTimeAvg_out(iPressure,Plane%IDlist(i,1))-pInf)/(0.5*rhoInf*uInf**2) !c_p
      ! perform scaling if required
      SELECT CASE(Plane_BLvelScaling)
      CASE(0) ! do nothing
      CASE(1) ! laminar scaling
        DO j=1,jj
          RPDataTimeAvg_out(TransMap(1:ndim,1),Plane%IDlist(i,j)) = &
            RPDataTimeAvg_out(TransMap(1:ndim,1),Plane%IDlist(i,j))/u_delta
          Plane%LocalCoord(2,i,j)=Plane%LocalCoord(2,i,j)/delta99
        END DO
      CASE(2) ! turbulent scaling
        u_tau=SQRT(ABS(tau_w)/rho_delta)
        DO j=1,jj
          RPDataTimeAvg_out(TransMap(1:ndim,1),Plane%IDlist(i,j)) = &
            RPDataTimeAvg_out(TransMap(1:ndim,1),Plane%IDlist(i,j))/u_tau
          Plane%LocalCoord(2,i,j)=Plane%LocalCoord(2,i,j)*u_tau*rho_delta/mu0
        END DO
      CASE(3) ! testing
        DO j=1,jj
          RPDataTimeAvg_out(TransMap(1:ndim,1),Plane%IDlist(i,j)) = &
            u_star(j)/u_max
        END DO
      END SELECT
    END DO !i
    DEALLOCATE(u_loc,y_loc,u_star)
  END IF!(Plane%Type.EQ.2) THEN ! BLPlane
END DO !iPlane
WRITE(UNIT_stdOut,'(A)')" done!"

END SUBROUTINE Plane_BLProps


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
END FUNCTION


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
WRITE(UNIT_stdOut,'(A)') '  EquationRP FINALIZED'
END SUBROUTINE FinalizeEquationRP

END MODULE MOD_EquationRP

