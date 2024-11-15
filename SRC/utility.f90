MODULE utility
USE constants
IMPLICIT NONE
CONTAINS

PURE RECURSIVE REAL*8 FUNCTION v_ac(f_ac, omega_ac, t)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: f_ac, t, omega_ac
  v_ac = -f_ac*SIN(omega_ac*t)
END FUNCTION v_ac

PURE COMPLEX*16 FUNCTION energy_phase_offset(En, Em, t)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: En, Em, t
  energy_phase_offset = EXP(imag*t*(En - Em))
END FUNCTION energy_phase_offset

PURE RECURSIVE COMPLEX*16 FUNCTION sigma_x_expected_value(Psi1, Psi2, psi_size)
  !! Calculates expectation value of sigma_x Pauli matrix.
  !! Assumes that order of psi is (up, down, up, down, ...)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size
  COMPLEX*16, INTENT(IN) :: Psi1(psi_size), Psi2(psi_size)
  INTEGER*4 :: i
  sigma_x_expected_value = 0.0d0
  DO i = 1, psi_size, 2
    sigma_x_expected_value = sigma_x_expected_value + CONJG(Psi1(i))*Psi2(i + 1) + CONJG(Psi1(i + 1))*Psi2(i)
  END DO
  RETURN
END FUNCTION sigma_x_expected_value

PURE RECURSIVE COMPLEX*16 FUNCTION sigma_y_expected_value(Psi1, Psi2, psi_size)
  !! Calculates expectation value of sigma_x Pauli matrix.
  !! Assumes that order of psi is (up, down, up, down, ...)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size
  COMPLEX*16, INTENT(IN) :: Psi1(psi_size), Psi2(psi_size)
  INTEGER*4 :: i
  sigma_y_expected_value = 0.0d0
  DO i = 1, psi_size, 2
    sigma_y_expected_value = sigma_y_expected_value - imag*CONJG(Psi1(i))*Psi2(i + 1) + imag*CONJG(Psi1(i + 1))*Psi2(i)
  END DO
  RETURN
END FUNCTION sigma_y_expected_value

PURE RECURSIVE COMPLEX*16 FUNCTION sigma_z_expected_value(Psi1, Psi2, psi_size)
  !! Calculates expectation value of sigma_z Pauli matrix.
  !! Assumes that order of psi is (up, down, up, down, ...)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size
  COMPLEX*16, INTENT(IN) :: Psi1(psi_size), Psi2(psi_size)
  INTEGER*4 :: i
  sigma_z_expected_value = 0.0d0
  DO i = 1, psi_size
    sigma_z_expected_value = sigma_z_expected_value + (-1)**(i+1)*CONJG(Psi1(i))*Psi2(i)
  END DO
  RETURN
END FUNCTION sigma_z_expected_value

PURE RECURSIVE REAL*8 FUNCTION d_xy_share(Psi, psi_size, norbs)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size, norbs
  COMPLEX*16, INTENT(IN) :: Psi(psi_size)
  INTEGER*4 :: i
  d_xy_share = 0.0d0
  DO i = 1, psi_size, norbs
    d_xy_share = d_xy_share + REAL(CONJG(Psi(i))*Psi(i)) !Spin up
    d_xy_share = d_xy_share + REAL(CONJG(Psi(i + 1))*Psi(i + 1)) !Spin down
  END DO
  RETURN
END FUNCTION d_xy_share


PURE RECURSIVE REAL*8 FUNCTION d_xz_share(Psi, psi_size, norbs)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size, norbs
  COMPLEX*16, INTENT(IN) :: Psi(psi_size)
  INTEGER*4 :: i, offset
  offset = 2
  d_xz_share = 0.0d0
  DO i = 1, psi_size, norbs
    d_xz_share = d_xz_share + REAL(CONJG(Psi(i + offset))*Psi(i + offset)) !Spin up
    d_xz_share = d_xz_share + REAL(CONJG(Psi(i + 1 + offset))*Psi(i + 1 + offset)) !Spin down
  END DO
  RETURN
END FUNCTION d_xz_share

PURE RECURSIVE REAL*8 FUNCTION d_yz_share(Psi, psi_size, norbs)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size, norbs
  COMPLEX*16, INTENT(IN) :: Psi(psi_size)
  INTEGER*4 :: i, offset
  offset = 4
  d_yz_share = 0.0d0
  DO i = 1, psi_size, norbs
    d_yz_share = d_yz_share + REAL(CONJG(Psi(i + offset))*Psi(i + offset)) !Spin up
    d_yz_share = d_yz_share + REAL(CONJG(Psi(i + 1 + offset))*Psi(i + 1 + offset)) !Spin down
  END DO
  RETURN
END FUNCTION d_yz_share

PURE RECURSIVE COMPLEX*16 FUNCTION single_electron_parity(Psi_1, Psi_2, psi_size, norbs, Nx, Ny)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size, norbs, Nx, Ny
  COMPLEX*16, INTENT(IN) :: Psi_1(psi_size), Psi_2(psi_size)
  INTEGER*4 :: i, offset
  single_electron_parity = DCMPLX(0.0d0, 0.0d0)
  DO i = 1, psi_size, norbs
    offset = 0
    single_electron_parity = single_electron_parity  + CONJG(Psi_1(i + offset))*Psi_2(get_opposite_r_index(Nx, Ny, norbs, i + offset))
    offset = 1
    single_electron_parity = single_electron_parity  - CONJG(Psi_1(i + offset))*Psi_2(get_opposite_r_index(Nx, Ny, norbs, i + offset))
    offset = 2
    single_electron_parity = single_electron_parity  - CONJG(Psi_1(i + offset))*Psi_2(get_opposite_r_index(Nx, Ny, norbs, i + offset))
    offset = 3
    single_electron_parity = single_electron_parity  + CONJG(Psi_1(i + offset))*Psi_2(get_opposite_r_index(Nx, Ny, norbs, i + offset))
    offset = 4
    single_electron_parity = single_electron_parity  - CONJG(Psi_1(i + offset))*Psi_2(get_opposite_r_index(Nx, Ny, norbs, i + offset))
    offset = 5
    single_electron_parity = single_electron_parity  + CONJG(Psi_1(i + offset))*Psi_2(get_opposite_r_index(Nx, Ny, norbs, i + offset))
  END DO
  RETURN
END FUNCTION

PURE RECURSIVE COMPLEX*16 FUNCTION single_electron_x_expected_value(Psi1, Psi2, norbs, Nx, dx, ham_1_size)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi1(ham_1_size), Psi2(ham_1_size)
  INTEGER*4, INTENT(IN) :: ham_1_size
  INTEGER*4, INTENT(IN) :: norbs, Nx
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: i
  single_electron_x_expected_value = (0.0d0, 0.0d0)
  DO i = 1, ham_1_size
    single_electron_x_expected_value = single_electron_x_expected_value + CONJG(Psi1(i))*Psi2(i)*get_x_from_psi_index(i, norbs, Nx, dx)
  END DO
  RETURN
END FUNCTION single_electron_x_expected_value

PURE RECURSIVE COMPLEX*16 FUNCTION single_electron_y_expected_value(Psi1, Psi2, norbs, Nx, Ny, dx, ham_1_size)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi1(ham_1_size), Psi2(ham_1_size)
  INTEGER*4, INTENT(IN) :: ham_1_size
  INTEGER*4, INTENT(IN) :: norbs, Nx, Ny
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: i
  single_electron_y_expected_value = (0.0d0, 0.0d0)
  DO i = 1, ham_1_size
    single_electron_y_expected_value = single_electron_y_expected_value + CONJG(Psi1(i))*Psi2(i)*get_y_from_psi_index(i, norbs, Nx, Ny, dx)
  END DO
  RETURN
END FUNCTION single_electron_y_expected_value

PURE RECURSIVE REAL*8 FUNCTION get_y_from_psi_index(i, norbs, Nx, Ny, dx)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: i, norbs, Nx, Ny
  REAL*8, INTENT(IN) :: dx
  get_y_from_psi_index = (((i - 1)/norbs)/(2*Nx + 1) - Ny) * dx !-1 Due to indexing from 1
END FUNCTION get_y_from_psi_index

PURE RECURSIVE REAL*8 FUNCTION get_x_from_psi_index(i, norbs, Nx, dx)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: i, norbs, Nx
  REAL*8, INTENT(IN) :: dx
  get_x_from_psi_index = (MOD((i - 1)/norbs, 2*Nx + 1) - Nx) * dx !-1 Due to indexing from 1
END FUNCTION get_x_from_psi_index

PURE RECURSIVE INTEGER*4 FUNCTION get_opposite_r_index(Nx, Ny, nstates, current_index)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: Nx, Ny, nstates, current_index
  INTEGER*4 :: so_index !Spin-orbital index from current index
  INTEGER*4 :: current_nx, current_ny, current_nx_opposite, current_ny_opposite
  so_index = MOD(current_index - 1, nstates)
  current_nx = MOD((current_index - 1)/nstates, 2*Nx + 1) - Nx
  current_ny = (((current_index - 1)/nstates)/(2*Nx + 1) - Ny)

  current_nx_opposite = -current_nx
  current_ny_opposite = -current_ny
  get_opposite_r_index = 1 + (current_nx_opposite + Nx)*nstates + (current_ny_opposite + Ny) * (2*Nx + 1) * nstates + so_index
  RETURN
END FUNCTION

PURE RECURSIVE INTEGER*4 FUNCTION get_upper_hermitian_index(i,j, size)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: i, j, size
  IF (j < i) THEN
    get_upper_hermitian_index = -1
    RETURN
  END IF
  get_upper_hermitian_index = ((i - 1) * (2*size - i)) / 2 + j !Mind the order!!!
  RETURN
END FUNCTION get_upper_hermitian_index

PURE INTEGER*4 FUNCTION factorial(n)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: n
  INTEGER*4 :: i
  factorial = 1
  DO i = 1, n
    factorial = factorial * i
  END DO
END FUNCTION

RECURSIVE SUBROUTINE GET_SLICE_FROM_HERMITIAN_MATRIX(Slice, Matrix, slice_size, original_dim, matrix_size,  i, j)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: slice_size, original_dim, matrix_size, i, j
  COMPLEX*16, INTENT(IN) :: Matrix(matrix_size, slice_size)
  COMPLEX*16, INTENT(OUT) :: Slice(slice_size)
  IF (j >= i) THEN
    Slice(:) = Matrix(get_upper_hermitian_index(i, j, original_dim), :)
  ELSE
    Slice(:) = CONJG(Matrix(get_upper_hermitian_index(j, i, original_dim), :))
  END IF
END SUBROUTINE

END MODULE utility