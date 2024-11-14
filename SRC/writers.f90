MODULE writers
USE constants
USE indata
USE utility
USE many_body

IMPLICIT NONE
CONTAINS

!########### Universal writers
SUBROUTINE WRITE_ENERGIES(Energies, nstates, filename)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Energies(nstates)
  INTEGER*4, INTENT(IN) :: nstates
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: i


  OPEN(unit = 1, FILE= filename, FORM = "FORMATTED", ACTION = "WRITE")
  WRITE(1,*) '#No. sate [-]   Energy [meV]'
  DO i = 1, nstates
    WRITE(1,*) i, Energies(i) / eV2au * 1e3
  END DO
  CLOSE(1)

END SUBROUTINE WRITE_ENERGIES


!########### Single-electron writers
SUBROUTINE WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi, psi_size, nstates, norbs, Nx, Ny, dx, filename)
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
    WRITE(1,*) '#x[nm] y[nm] Re(Psi(orb1, s = +)) Im(Psi(orb1, s = +)) Re(Psi(orb1, s = -)) ...'

    nn = 1
    DO iy = -Ny, Ny
    DO ix = -Nx, Nx
      WRITE(1,'(2E20.8)', ADVANCE='NO') ix*dx/nm2au, iy*dx/nm2au

      DO iorb = 1, norbs / 2
        DO i = 1, 2
          WRITE(1, '(2E20.8)', ADVANCE='NO') REAL(Psi(nn,is)), AIMAG(Psi(nn,is))
          nn = nn + 1
        END DO
      END DO
      WRITE(1,*)
    END DO
    END DO

    CLOSE(1)
  END DO


END SUBROUTINE WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS

SUBROUTINE WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1, ham_1_size, nstates, norbitals, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstates, norbitals
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: n
  CHARACTER(LEN=200) :: format_string

  !Writing expectations for a given state
  format_string = '(I10, 9E20.8)'
  OPEN(unit = 9, FILE= filename, FORM = "FORMATTED", ACTION = "WRITE")
  WRITE(9,*) "#No. state [-]  <s_x>   <s_y>   <s_z>   <d_xy>   <d_xz>   <d_yz>   <parity>   <x>"
  DO n = 1, nstates
    WRITE(9, format_string) n,&
      & REAL(sigma_x_expected_value(Psi_1(:,n), Psi_1(:,n), ham_1_size)),&
      & REAL(sigma_y_expected_value(Psi_1(:,n), Psi_1(:,n), ham_1_size)),&
      & REAL(sigma_z_expected_value(Psi_1(:,n), Psi_1(:,n), ham_1_size)), &
      & d_xy_share(Psi_1(:,n), ham_1_size, norbitals), d_xz_share(Psi_1(:,n), ham_1_size, norbitals), d_yz_share(Psi_1(:,n), ham_1_size, norbitals), &
      & REAL(single_electron_parity(Psi_1(:,n),Psi_1(:,n), ham_1_size, norbitals, Nx, Ny)), &
      & REAL(single_electron_x_expected_value(Psi_1(:,n), Psi_1(:,n), norbitals, Nx, dx, ham_1_size)), &
      & REAL(single_electron_y_expected_value(Psi_1(:,n), Psi_1(:,n), norbitals, Nx, Ny, dx, ham_1_size))
  END DO
  CLOSE(9)

END SUBROUTINE WRITE_SINGLE_ELECTRON_EXPECTATIONS


!########### Multi-electron writers
SUBROUTINE WRITE_SLATER_COEFFICIENTS(C_slater, ham_2_size, nstates, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates)
  INTEGER*4, INTENT(IN) :: ham_2_size, nstates
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: i, j
  CHARACTER(LEN=200) :: format_string

  !Writing expectations for a given state
  format_string = '(2I10, 2E20.8)'
  OPEN(unit = 9, FILE= filename, FORM = "FORMATTED", ACTION = "WRITE")
  WRITE(9,*) "#No. multielectron state [-]    No. combination [-]    Re(c_slater)    Im(c_slater)"
  DO i = 1, nstates
    DO j = 1, ham_2_size
      WRITE(9, format_string) i, j, REAL(C_slater(j,i)), AIMAG(C_slater(j,i))
    END DO
    WRITE(9,*)
    WRITE(9,*)
  END DO
  CLOSE(9)

END SUBROUTINE WRITE_SLATER_COEFFICIENTS


SUBROUTINE WRITE_MULTI_ELECTRON_EXPECTATIONS(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, dx, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstate_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: nstate_1, nstate_2, norbs
  INTEGER*4, INTENT(IN) :: Nx, Ny
  REAL*8, INTENT(IN) :: dx
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: n
  CHARACTER(LEN=200) :: format_string

  !Writing expectations for a given state
  format_string = '(I10, 5E20.8)'
  OPEN(unit = 9, FILE= filename, FORM = "FORMATTED", ACTION = "WRITE")
  WRITE(9,*) "#No. state [-]  <x>    <S_x>   <S_y>   <S_z>   <parity>"
  DO n = 1, nstate_2
    WRITE(9, format_string) n,&
    & REAL(many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n, Nx, dx, norbs)) / nm2au / k_electrons,&
    & REAL(many_body_sigma_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n)),&
    & REAL(many_body_sigma_y_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n)),&
    & REAL(many_body_sigma_z_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n)),&
    & REAL(many_body_parity_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n))
  END DO
  CLOSE(9)

END SUBROUTINE WRITE_MULTI_ELECTRON_EXPECTATIONS

!########### Time-dependent calculations writers




END MODULE writers