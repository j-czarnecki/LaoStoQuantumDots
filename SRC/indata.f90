
MODULE indata
  USE constants
  IMPLICIT  NONE
  SAVE

  !parameters for numerical calculations
  INTEGER :: Nx             !number of sites in the x direction
  INTEGER :: Ny             !number of sites in the  direction
  REAL*8  :: dx             !size of grid
  INTEGER :: norbs          !number of orbitals including spin
  INTEGER :: nstate         !number of single electron states taken for calculations

  !physical paramters for LAO/STO
  REAL :: tl  
  REAL :: th
  REAL :: td
  REAL :: dso
  REAL :: drso
  REAL :: dE
  REAL :: g

  !external parameters
  REAL :: Bx  
  REAL :: By
  REAL :: Bz

  REAL :: a !lattice constant

  
  NAMELIST /calculation_parameters/             &
       &  Nx,                                   &
       &  Ny,                                   &
       &  dx,                                   &
       &  norbs,                                &
       &  nstate
  
  
  NAMELIST /physical_parameters/               &
       &  tl,                                  &
       &  th,                                  &
       &  td,			                             &
       &  dso,			                           &
       &  drso,                                &
       &  dE,                                  &
       &  g
  
  ! calculations flags
       NAMELIST /external_parameters/          &
       &  Bx,                                  &
       &  By,                                  &
       &  Bz
       
CONTAINS


!===================================================================
SUBROUTINE INDATA_GET(nmlfile)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: nmlfile

  INTEGER :: i
  
  OPEN(33, FILE=TRIM(nmlfile), FORM="FORMATTED", ACTION="READ",  &
       &   STATUS="OLD")

  !default values
  Nx=0
  Ny=0

  tl=0.0
  th=0.0
  td=0.0
  dso=0.0
  drso=0.0
  dx=0.0
  dE=0.0
  g=0
  norbs=1
  nstate=10

  Bx=0.0
  By=0.0
  Bz=0.0
  
! read namelist
  READ(33,NML=calculation_parameters)  
  
  IF(Nx.le.0) then
    print *, "Nx has to be a positive interger. STOP"
    STOP
  ENDIF
  
  IF(Ny.le.0) then
    print *, "Ny has to be a positive interger. STOP"
    STOP
  ENDIF

  IF(dx.le.0) then
    print *, "dx has to be a positive interger. STOP"
    STOP
  ENDIF

  ! read namelist
  READ(33,NML=physical_parameters)  

  a=0.39*nm2au
  th=th*eV2au
  tl=tl*eV2au
  td=td*eV2au
  dso=dso*eV2au
  drso=drso*eV2au
  dx=dx*nm2au
  dE=dE*eV2au

  ! read namelist
  READ(33,NML=external_parameters) 
  Bx=Bx*T2au
  By=By*T2au
  Bz=Bz*T2au


END SUBROUTINE INDATA_GET

END MODULE indata
