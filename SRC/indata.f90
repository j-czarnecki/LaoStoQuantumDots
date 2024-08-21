
MODULE indata
  USE constants
  IMPLICIT NONE
  SAVE

  !parameters for numerical calculations
  INTEGER :: Nx             !number of sites in the x direction
  INTEGER :: Ny             !number of sites in the  direction
  REAL*8  :: dx             !size of grid
  INTEGER :: norbs          !number of orbitals including spin
  INTEGER :: nstate         !number of single electron states taken for calculations

  !physical paramters for LAO/STO
  REAL*8 :: tl
  REAL*8 :: th
  REAL*8 :: td
  REAL*8 :: dso
  REAL*8 :: drso
  REAL*8 :: dE
  REAL*8 :: g

  !external parameters
  REAL*8 :: Bx
  REAL*8 :: By
  REAL*8 :: Bz

  !Derived parameters
  REAL*8 :: a !lattice constant

  NAMELIST /calculation_parameters/             &
       &  Nx,                                   &
       &  Ny,                                   &
       &  dx,                                   &
       &  norbs,                                &
       &  nstate

  NAMELIST /physical_parameters/               &
       &  tl,                                  &
       &  th,                                  &
       &  td,                                                     &
       &  dso,                                                   &
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

    OPEN (33, FILE=TRIM(nmlfile), FORM="FORMATTED", ACTION="READ",  &
         &   STATUS="OLD")

    !default values
    Nx = 0
    Ny = 0

    tl = 0.0
    th = 0.0
    td = 0.0
    dso = 0.0
    drso = 0.0
    dx = 0.0
    dE = 0.0
    g = 0
    norbs = 1
    nstate = 10

    Bx = 0.0
    By = 0.0
    Bz = 0.0

! read namelist
    READ (33, NML=calculation_parameters)

    IF (Nx .LE. 0) THEN
      PRINT *, "Nx has to be a positive interger. STOP"
      STOP
    END IF

    IF (Ny .LE. 0) THEN
      PRINT *, "Ny has to be a positive interger. STOP"
      STOP
    END IF

    IF (dx .LE. 0) THEN
      PRINT *, "dx has to be a positive interger. STOP"
      STOP
    END IF

    ! read namelist
    READ (33, NML=physical_parameters)

    a = 0.39 * nm2au !Julian: is this meant to be a fixed value?
    th = th * eV2au
    tl = tl * eV2au
    td = td * eV2au
    dso = dso * eV2au
    drso = drso * eV2au
    dx = dx * nm2au
    dE = dE * eV2au

    ! read namelist
    READ (33, NML=external_parameters)
    Bx = Bx * T2au
    By = By * T2au
    Bz = Bz * T2au

  END SUBROUTINE INDATA_GET

END MODULE indata
