!===================================================================================================================================
!> Contains global variables provided by the output routines
!===================================================================================================================================
MODULE MOD_Parameters
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(len=255)                  :: ProjectName
CHARACTER(len=255)                  :: RP_DefFile 
CHARACTER(len=255),ALLOCATABLE      :: GroupNames_visu(:)
INTEGER                             :: nGroups_visu
LOGICAL                             :: OutputTimeData
LOGICAL                             :: doFluctuations
LOGICAL                             :: equiTimeSpacing
LOGICAL                             :: OutputTimeAverage
LOGICAL                             :: OutputLines,OutputPlanes,OutputPoints
LOGICAL                             :: Line_GlobalCoords
LOGICAL                             :: Line_LocalCoords
LOGICAL                             :: Line_LocalVel
LOGICAL                             :: Plane_LocalCoords
LOGICAL                             :: Plane_LocalVel
LOGICAL                             :: Plane_doBLProps
INTEGER                             :: Plane_BLvelScaling
LOGICAL                             :: usePrims
LOGICAL                             :: RP_SET_defined
!--------------------------------------------------
! Filter
LOGICAL                             :: doFilter
INTEGER                             :: FilterMode
REAL                                :: FilterWidth
!--------------------------------------------------
! Spectral Analysis, FFT, PSD
LOGICAL                             :: doSpec
LOGICAL                             :: doPSD
LOGICAL                             :: doFFT
INTEGER                             :: nBlocks 
INTEGER                             :: BlockSize
REAL                                :: cutoffFreq,samplingFreq
LOGICAL                             :: doHanning
LOGICAL                             :: fourthDeriv,ThirdOct
REAL                                :: u_infPhys,chordPhys
REAL                                :: Line_LocalVel_vec(3)
REAL                                :: Mu0
!--------------------------------------------------
! Turbulence
LOGICAL                             :: doTurb
INTEGER                             :: nVarVisu  
CHARACTER(len=255),ALLOCATABLE      :: VarNamevisu(:)
INTEGER                             :: OutputFormat
INTEGER                             :: Skip ! nur jeder skipte RP sample wird eingelesen



INTEGER,ALLOCATABLE               :: DepTable(:,:)
CHARACTER(LEN=255),ALLOCATABLE    :: DepNames(:)
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarNamesAll(:)
INTEGER                           :: nVarDep                 ! 
INTEGER                           :: nVarCalc
INTEGER                           :: nVarVisuTotal
INTEGER,ALLOCATABLE               :: mapCalc(:)
INTEGER,ALLOCATABLE               :: mapVisu(:)
!===================================================================================================================================
END MODULE MOD_Parameters
