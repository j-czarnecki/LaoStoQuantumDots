#include "macros_def.f90"
MODULE indata
  USE constants
  USE logger
  IMPLICIT NONE
  SAVE

  !parameters for numerical calculations
  INTEGER :: Nx             !number of sites in the x direction
  INTEGER :: Ny             !number of sites in the  direction
  REAL*8  :: dx             !size of grid
  INTEGER :: norbs          !number of orbitals including spin
  INTEGER :: nstate_1       !number of single electron states taken for calculations
  INTEGER :: nstate_2       !number of many electron states taken for calculations
  INTEGER :: k_electrons    !number of electrons
  REAL*8 :: dt              !Time step [ns]
  REAL*8 :: t_max           !Time span [ns]

  !physical paramters for LAO/STO
  REAL*8 :: tl
  REAL*8 :: th
  REAL*8 :: td
  REAL*8 :: dso
  REAL*8 :: drso
  REAL*8 :: dE
  REAL*8 :: g
  REAL*8 :: eps_r

  !external parameters
  REAL*8 :: omega
  REAL*8 :: Bx
  REAL*8 :: By
  REAL*8 :: Bz
  REAL*8 :: domega_ac
  REAL*8 :: omega_ac_max
  REAL*8 :: f_ac
  REAL*8 :: Vb
  REAL*8 :: V0

  !Derived parameters
  REAL*8 :: a !lattice constant
  INTEGER*4 :: N_t_steps    !Number of time steps
  INTEGER*4 :: N_omega_ac_steps !number of steps of omega_ac
  NAMELIST /calculation_parameters/             &
       &  Nx,                                   &
       &  Ny,                                   &
       &  dx,                                   &
       &  norbs,                                &
       &  nstate_1,                             &
       &  nstate_2,                             &
       &  k_electrons,                          &
       &  dt,                                   &
       &  t_max

  NAMELIST /physical_parameters/               &
       &  tl,                                  &
       &  th,                                  &
       &  td,                                                     &
       &  dso,                                                   &
       &  drso,                                &
       &  dE,                                  &
       &  g,                                   &
       &  eps_r

  ! calculations flags
  NAMELIST /external_parameters/          &
  &  omega,                               &
  &  Bx,                                  &
  &  By,                                  &
  &  Bz,                                  &
  &  domega_ac,                           &
  &  omega_ac_max,                        &
  &  f_ac,                                &
  &  Vb,                                  &
  &  V0

CONTAINS

!===================================================================
  SUBROUTINE INDATA_GET(nmlfile)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: nmlfile

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
    nstate_1 = 6
    nstate_2 = 2

    omega = 0.0
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

    a = 0.39 * nm2au !Julian: is this meant to be a fixed value?
    dx = dx * nm2au
    dt = dt * ns2au
    t_max = t_max * ns2au
    N_t_steps = INT(t_max / dt)


    ! read namelist
    READ (33, NML=physical_parameters)

    th = th * eV2au
    tl = tl * eV2au
    td = td * eV2au
    dso = dso * eV2au
    drso = drso * eV2au
    dE = dE * eV2au

    IF (eps_r == 0.0d0) STOP "eps_r cannot be 0.0"

    ! read namelist
    READ (33, NML=external_parameters)
    omega = omega * eV2au
    Bx = Bx * T2au
    By = By * T2au
    Bz = Bz * T2au
    domega_ac = domega_ac * eV2au
    omega_ac_max = omega_ac_max * eV2au
    f_ac = f_ac * F2au
    Vb = Vb * eV2au
    V0 = V0 * eV2au

    N_omega_ac_steps = INT(omega_ac_max / domega_ac)


  END SUBROUTINE INDATA_GET


  SUBROUTINE READ_SINGLE_ELECTRON_WAVEFUNCTIONS(path, Psi, psi_size, nstates, norbs)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: path
    INTEGER*4, INTENT(IN) :: psi_size, nstates, norbs
    COMPLEX*16, INTENT(OUT) :: Psi(psi_size, nstates)
    CHARACTER(LEN=200) :: filename_state
    INTEGER*4 :: n, i, j, m
    REAL*8 :: x, y
    REAL*8, ALLOCATABLE :: Line_as_array(:)

    ALLOCATE(Line_as_array(norbs*2)) ! *2 due to separation of real and imaginary part

    DO n = 1, nstates
      WRITE(log_string,*) "Reading Psi_1_n", n
      LOG_INFO(log_string)

      WRITE(filename_state, '(2A, I0, A)') path, '/OutputData/Psi_1_n', n, '.dat'
      OPEN(10, FILE = filename_state, ACTION = 'READ', FORM = 'FORMATTED')
      READ(10, '(A)') !Skip first line, since it is a comment
      DO i = 1, psi_size, norbs
        READ(10,'(14E20.8)') x, y, Line_as_array
        DO j = 0, norbs*2 - 1, 2 !Shifted to 0, because that offset with respect to current subvector position
          Psi(i + j, n) = DCMPLX(Line_as_array(j + 1), Line_as_array(j + 2))
        END DO
      END DO
      CLOSE(10)
    END DO

    DEALLOCATE(Line_as_array)
  END SUBROUTINE READ_SINGLE_ELECTRON_WAVEFUNCTIONS

END MODULE indata
