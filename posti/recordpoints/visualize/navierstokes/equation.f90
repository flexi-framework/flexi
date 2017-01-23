#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_Equation
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE

INTERFACE CalcEquation
  MODULE PROCEDURE CalcEquation
END INTERFACE

INTERFACE Plane_BLProps
  MODULE PROCEDURE Plane_BLProps
END INTERFACE


INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC::InitEquation,CalcEquation,Plane_BLProps,FinalizeEquation
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the visualization and map the variable names to classify these in conservative and derived quantities.
!===================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_EOS               ,ONLY:InitEOS
USE MOD_RPData_Vars       ,ONLY:VarNames_HDF5,nVar_HDF5
USE MOD_Parameters        ,ONLY:nVar_visu,VarNameVisu,usePrims
USE MOD_Parameters        ,ONLY:Plane_doBLProps
USE MOD_Equation_Vars
USE MOD_VarNameMappingsRP
USE MOD_VarNameMappingsRP_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: iVar1,iVar2,strlen
INTEGER                               :: mapCand(3)
CHARACTER(LEN=255)                    :: VarName(5),tmp255,tmp255_2
LOGICAL                               :: justVisualizeState=.FALSE.
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT EQUATION ...'

! Check which conservative variables are to be visualized
! ReadFromFile detects if any derived quantity already exists in the state file
! (i.e. with time avg state files)  
max_nVar_visu= 1

! In case no output variables specified: take instead the variable names out of the hDF5 file (used for timeavg-files)
IF(nVar_visu .LT. 1) THEN
  WRITE(*,*) 'No output variables specified, using existing variables from state file.'
  nVar_visu   = nVar_HDF5
  SDEALLOCATE(VarNameVisu)
  ALLOCATE(VarNameVisu(nVar_visu))
  VarNameVisu = VarNames_HDF5
END IF


justVisualizeState=.FALSE.
IF(usePrims) THEN
  ! USE THE PRIMITIVE STATE VECTOR FOR ALL DERIVED QUANTITIES 
  VarName(1) ='Density'
  VarName(2) ='VelocityX'
  VarName(3) ='VelocityY'
  VarName(4) ='VelocityZ'
  VarName(5) ='Pressure'
  DO iVar1=1,5
    ! check if primitive variables are in the state file and create mapping
    iVar2 = GETMAPBYNAME(VarName(iVar1),VarNames_HDF5,nVar_HDF5)
    IF(iVar2.NE.-1)THEN
      PrimMap(iVar1)=iVar2
    ELSE 
      justVisualizeState=.TRUE.
      WRITE(UNIT_StdOut,*) 'WARNING: Not all Primitive Variables available in State File,'
      WRITE(UNIT_StdOut,*) '         cannot calculate any derived quantities.'
      WRITE(UNIT_StdOut,*) '         Visualizing State File Variables instead!'
      EXIT     
    END IF
  END DO
ELSE
  ! USE THE CONSERVATIVE STATE VECTOR FOR ALL DERIVED QUANTITIES 
  !  in this case, the variables in the state file need to be in the order below
  VarName(1) ='Density'
  VarName(2) ='MomentumX'
  VarName(3) ='MomentumY'
  VarName(4) ='MomentumZ'
  VarName(5) ='EnergyStagnationDensity'
  DO iVar1=1,5
    ! check if main conservative variables  are in the state file
    iVar2 = GETMAPBYNAME(VarName(iVar1),VarNames_HDF5,nVar_HDF5)
      PrimMap(iVar1)=iVar2
    IF(iVar2.EQ.-1)THEN
      justVisualizeState=.TRUE.
      WRITE(UNIT_StdOut,*) 'WARNING: Not all Conservative Variables available in State File,'
      WRITE(UNIT_StdOut,*) '         cannot calculate any derived quantities.'
      WRITE(UNIT_StdOut,*) '         Visualizing State File Variables instead!'
      EXIT
    END IF
  END DO
  ! check order of the Conservatives
  DO iVar1=1,5
    IF(iVar1.NE.PrimMap(iVar1)) THEN
      justVisualizeState=.TRUE.
      WRITE(UNIT_StdOut,*) 'WARNING: Conservative Variables are not in the right order,'
      WRITE(UNIT_StdOut,*) '         cannot calculate any derived quantities.'
      WRITE(UNIT_StdOut,*) '         Visualizing State File Variables instead!'
      EXIT
    END IF
  END DO
END IF

IF(justVisualizeState) THEN
  ! visualize all state file variables
  nVar_visu   = nVar_HDF5
  SDEALLOCATE(VarNameVisu)
  ALLOCATE(VarNameVisu(nVar_visu))
  VarNameVisu = VarNames_HDF5
  CALL CreateStateMappings(nVar_HDF5,VarNames_HDF5,Cons)
ELSE
! Create Mappings to the variables to be visualized:
  ! Variables directly available in State file, usually conservative variables
  WRITE(UNIT_StdOut,*) 'Preparing variables from state file...'
  CALL CreateStateMappings(nVar_HDF5,VarNames_HDF5,Cons)
  
  ! Primitive variables
  WRITE(UNIT_StdOut,*) 'Preparing primitive variables...'
  ! Initialize the primitive Variables calculation
  CALL InitEOS()
END IF

! For local transforms, we need a mapping to all physical vector quantities
! with 3 available components to transform them
ALLOCATE(TransMap(3,INT(nVar_visu/2)))
ALLOCATE(is2D(INT(nVar_visu/2)))
is2D=.FALSE.
TransMap=-1
nVecTrans=0
DO iVar1=1,nVar_visu
  tmp255=TRIM(VarNameVisu(iVar1))
  strlen=LEN(TRIM(ADJUSTL(tmp255)))
  IF(tmp255(strlen:strlen).EQ.TRIM('X')) THEN
    WRITE(tmp255_2,'(A)')TRIM(tmp255(1:strlen-1))//'Y'
    mapCand(2)=GETMAPBYNAME(TRIM(tmp255_2),VarNameVisu,nVar_visu) 
    WRITE(tmp255_2,'(A)')TRIM(tmp255(1:strlen-1))//'Z'
    mapCand(3)=GETMAPBYNAME(TRIM(tmp255_2),VarNameVisu,nVar_visu) 
    mapCand(1)=iVar1
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
DO iVar1=1,nVecTrans
  DO iVar2=1,3
    WRITE(UNIT_StdOut,'(A,A)')'  ',TRIM(VarNameVisu(TransMap(iVar2,iVar1)))
  END DO
END DO

IF(Plane_doBLProps) THEN
  nBLProps=10
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
END IF

WRITE(UNIT_stdOut,'(A)')' INIT EQUATION DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
EquationInitIsDone=.TRUE.
END SUBROUTINE InitEquation



!===================================================================================================================================
!> This routine computes the state on the visualization grid 
!===================================================================================================================================
SUBROUTINE CalcEquation()
! MODULES
USE MOD_Globals
USE MOD_VarNameMappingsRP_Vars
USE MOD_EOS
USE MOD_EOS_vars
USE MOD_RPData_Vars       ,ONLY:RPData       
USE MOD_RPData_Vars       ,ONLY:nVar_HDF5              
USE MOD_RPSet_Vars        ,ONLY:nRP_global        
USE MOD_OutputRPVisu_Vars       ,ONLY:nSamples_out
USE MOD_Parameters        ,ONLY:usePrims
USE MOD_Parameters        ,ONLY:nVar_visu
USE MOD_Parameters        ,ONLY:Line_LocalVel,Plane_LocalVel,Plane_doBLProps
USE MOD_OutputRPVisu_Vars       ,ONLY:RPData_out 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iSample
REAL,ALLOCATABLE            :: U_RP_out(:,:)
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')" CONVERT DERIVED QUANTITIES..."
! CALCULATE DERIVED QUATITIES -----------------------------------------------------------------------------------------------------!
ALLOCATE(U_RP_out(1:nVar_visu,1:nRP_global))

! Conservative variables
IF(Cons%nVar_Visu .GT. 0)THEN
  IF(usePrims)THEN
    RPData_out(Cons%IndGlobal(Cons%Ind),:,1:nSamples_out)=RPData(Cons%Ind,:,1:nSamples_out)
  ELSE
    RPData_out(Cons%IndGlobal(Cons%Ind),:,1:nSamples_out)=RPData(Cons%Ind,:,1:nSamples_out)
  END IF
END IF !(nConsVisu .GT. 0)

! Primitive Variables
IF(Prim%nVar_visu .GT. 0)THEN
  DO iSample=1,nSamples_out
    U_RP_out=0.
    CALL CalcPrims(nVar_HDF5,nRP_global,RPData(:,:,iSample),U_RP_out(1:Prim%nVar_visu,:))
    RPData_out(Prim%IndGlobal(Prim%Ind),:,iSample)=U_RP_out(1:Prim%nVar_visu,:)
  END DO
END IF  !(nPrimVisu .GT. 0)

! Coordinate Transform
IF(Line_LocalVel) &
  CALL Line_TransformVel()
IF(Plane_LocalVel) &
  CALL Plane_TransformVel()

WRITE(UNIT_stdOut,'(A)')" CONVERT DERIVED QUANTITIES DONE!"
END SUBROUTINE CalcEquation



!===================================================================================================================================
!> This routine computes the state on the visualization grid 
!===================================================================================================================================
SUBROUTINE Line_TransformVel()
! MODULES
USE MOD_Globals
USE MOD_OutputRPVisu_Vars           ,ONLY:RPData_out,nSamples_out
USE MOD_RPSet_Vars            ,ONLY:nLines,Lines,tLine
USE MOD_Equation_Vars         ,ONLY:nVecTrans,TransMap
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iLine,iRP,iSample,iVec
TYPE(tLine),POINTER :: aLine
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')" coordinate transform of velocity along lines.."
DO iLine=1,nLines
  aLine=>Lines(iLine)
  DO iRP=1,aLine%nRP
    DO iSample=1,nSamples_out
      DO iVec=1,nVecTrans
        RPData_out(TransMap(1:3,iVec),aLine%IDlist(iRP),iSample)=&
           MATMUL(aLine%Tmat,RPData_out(TransMap(1:3,iVec),aLine%IDlist(iRP),iSample))
      END DO
    END DO !iSample
  END DO !iRP
END DO !iLine
WRITE(UNIT_stdOut,'(A)')" done!"
END SUBROUTINE Line_TransformVel


!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE Plane_TransformVel()
! MODULES
USE MOD_Globals
USE MOD_OutputRPVisu_Vars           ,ONLY:nSamples_out,RPData_out
USE MOD_RPSet_Vars            ,ONLY:nPlanes,Planes,tPlane
USE MOD_Equation_Vars         ,ONLY:nVecTrans,TransMap,is2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
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
            DO iSample=1,nSamples_out
              RPData_out(TransMap(1:3,iVec),Plane%IDlist(i,j),iSample)=&
                 MATMUL(Tmat,RPData_out(TransMap(1:3,iVec),Plane%IDlist(i,j),iSample)) 
            END DO !iSample
          END DO !j
        ELSE ! 2D case
          DO j=1,Plane%nRP(2)
            DO iSample=1,nSamples_out
              RPData_out(TransMap(1:2,iVec),Plane%IDlist(i,j),iSample)=&
                 MATMUL(Tmat(1:2,1:2),RPData_out(TransMap(1:2,iVec),Plane%IDlist(i,j),iSample)) 
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
!> 
!===================================================================================================================================
SUBROUTINE Plane_BLProps()
! MODULES
USE MOD_Globals
USE MOD_OutputRPVisu_Vars           ,ONLY:nSamples_out,RPDataTimeAvg_out
USE MOD_RPSet_Vars            ,ONLY:nPlanes,Planes,tPlane,xF_RP
USE MOD_Equation_Vars         ,ONLY:is2D,TransMap,nBLProps
USE MOD_Parameters        ,ONLY:Plane_BLvelScaling
USE MOD_Parameters        ,ONLY:Mu0
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iPlane,i,j,jj,j_max,ndim
REAL                 :: uu,yy,u_max,y_max,u_int,diffu1,diffu2,y1,y2,u1,u2,u3,rho1,rho2,rho_delta,dudy,dudy1
REAL                 :: dy,dy2
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
        dy=y_loc(j)-y_loc(j-1)
        ! displacement thickness
        delta1=delta1+ 0.5*dy*(1.-u_loc(j) +1.-u_loc(j-1))
        ! momentum thickness
        theta=theta+ 0.5*dy*(u_loc(j)*(1.-u_loc(j)) +u_loc(j-1)*(1.-u_loc(j-1)))
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
!> 
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
!>
!===================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Globals
USE MOD_EOS  ,ONLY:FinalizeEOS
USE MOD_Equation_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DEALLOCATE(TransMap,is2D) 
CALL FinalizeEOS()
WRITE(UNIT_stdOut,'(A)') '  EQUATION FINALIZED'
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation

