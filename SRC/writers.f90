MODULE writers
USE constants
USE indata
USE utility

IMPLICIT NONE
CONTAINS

SUBROUTINE WRITE_STATE_MAP(Psi, psi_size, nstates, norbs, Nx, Ny, dx, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi(psi_size, nstates)
  INTEGER*4, INTENT(IN) :: Nx, Ny, psi_size, nstates, norbs
  REAL*8, INTENT(IN) :: dx
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=200) :: filename_state

  INTEGER*4 :: nn, ix, iy, is, i, iorb

  DO is = 1, nstates
    WRITE(filename_state, '(2A, I0, A)') filename, '_n', is, '.dat'
    OPEN (1, FILE=filename_state)
    WRITE(1,*) '#x[nm] y[nm] |Psi(orb1, s = +)|**2 |Psi(orb1, s = -)|**2 |Psi(orb2, s = +)|**2 ...'

    nn = 1
    DO iy = -Ny, Ny
    DO ix = -Nx, Nx
      WRITE(1,'(2E20.8)', ADVANCE='NO') ix*dx/nm2au, iy*dx/nm2au

      DO iorb = 1, norbs / 2
        DO i = 1, 2
          WRITE(1, '(E20.8)', ADVANCE='NO') ABS(Psi(nn,is))**2
          nn = nn + 1
        END DO
      END DO
      WRITE(1,*)
    END DO
    END DO
    IF (nn /= psi_size + 1) PRINT*, "Error in printing state map, nn /= psi_size + 1"

    CLOSE(1)
  END DO


END SUBROUTINE WRITE_STATE_MAP

SUBROUTINE WRITE_ENERGIES(Energies, nstates, filename)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Energies(nstates)
  INTEGER*4, INTENT(IN) :: nstates
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: i


  OPEN (1, FILE=filename)
  WRITE(1,*) '#No. sate [-]   Energy [meV]'
  DO i = 1, nstates
    WRITE(1,*) i, Energies(i) / eV2au * 1e3
  END DO
  CLOSE(1)

END SUBROUTINE WRITE_ENERGIES

SUBROUTINE WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1, ham_1_size, nstates, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstates
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: i, j, n
  CHARACTER(LEN=200) :: format_string



END SUBROUTINE WRITE_SINGLE_ELECTRON_EXPECTATIONS

END MODULE writers