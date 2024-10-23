MODULE utility
USE constants
IMPLICIT NONE
CONTAINS

PURE REAL*8 FUNCTION v_ac(f_ac, omega_ac, t)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: f_ac, t, omega_ac
  v_ac = -f_ac*SIN(omega_ac*t)
END FUNCTION v_ac

PURE COMPLEX*16 FUNCTION energy_phase_offset(En, Em, t)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: En, Em, t
  energy_phase_offset = EXP(imag*t*(En - Em))
END FUNCTION energy_phase_offset

PURE COMPLEX*16 FUNCTION sigma_x_expected_value(Psi1, Psi2, psi_size)
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

PURE COMPLEX*16 FUNCTION sigma_y_expected_value(Psi1, Psi2, psi_size)
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

PURE COMPLEX*16 FUNCTION sigma_z_expected_value(Psi1, Psi2, psi_size)
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

PURE REAL*8 FUNCTION d_xy_share(Psi, psi_size, norbs)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size, norbs
  COMPLEX*16, INTENT(IN) :: Psi(psi_size)
  INTEGER*4 :: i
  d_xy_share = 0.0d0
  DO i = 1, psi_size, norbs
    d_xy_share = d_xy_share + CONJG(Psi(i))*Psi(i) !Spin up
    d_xy_share = d_xy_share + CONJG(Psi(i + 1))*Psi(i + 1) !Spin down
  END DO
  RETURN
END FUNCTION d_xy_share


PURE REAL*8 FUNCTION d_xz_share(Psi, psi_size, norbs)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size, norbs
  COMPLEX*16, INTENT(IN) :: Psi(psi_size)
  INTEGER*4 :: i, offset
  offset = 2
  d_xz_share = 0.0d0
  DO i = 1, psi_size, norbs
    d_xz_share = d_xz_share + CONJG(Psi(i + offset))*Psi(i + offset) !Spin up
    d_xz_share = d_xz_share + CONJG(Psi(i + 1 + offset))*Psi(i + 1 + offset) !Spin down
  END DO
  RETURN
END FUNCTION d_xz_share

PURE REAL*8 FUNCTION d_yz_share(Psi, psi_size, norbs)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: psi_size, norbs
  COMPLEX*16, INTENT(IN) :: Psi(psi_size)
  INTEGER*4 :: i, offset
  offset = 4
  d_yz_share = 0.0d0
  DO i = 1, psi_size, norbs
    d_yz_share = d_yz_share + CONJG(Psi(i + offset))*Psi(i + offset) !Spin up
    d_yz_share = d_yz_share + CONJG(Psi(i + 1 + offset))*Psi(i + 1 + offset) !Spin down
  END DO
  RETURN
END FUNCTION d_yz_share

COMPLEX*16 FUNCTION single_electron_x_expected_value(Psi1, Psi2, norbs, Nx, Ny, dx, ham_1_size)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi1(ham_1_size), Psi2(ham_1_size)
  INTEGER*4, INTENT(IN) :: ham_1_size
  INTEGER*4, INTENT(IN) :: norbs, Nx, Ny
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: i,j
  single_electron_x_expected_value = (0.0d0, 0.0d0)
  DO i = 1, ham_1_size
    single_electron_x_expected_value = single_electron_x_expected_value + CONJG(Psi1(i))*Psi2(i)*get_x_from_psi_index(i, norbs, Nx, dx)
  END DO
  RETURN
END FUNCTION single_electron_x_expected_value

PURE REAL*8 FUNCTION get_y_from_psi_index(i, norbs, Nx, Ny, dx)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: i, norbs, Nx, Ny
  REAL*8, INTENT(IN) :: dx
  get_y_from_psi_index = ((i/norbs)/(2*Nx + 1) - Ny) * dx
END FUNCTION get_y_from_psi_index

PURE REAL*8 FUNCTION get_x_from_psi_index(i, norbs, Nx, dx)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: i, norbs, Nx
  REAL*8, INTENT(IN) :: dx
  get_x_from_psi_index = (MOD(i/norbs, 2*Nx + 1) - Nx) * dx
END FUNCTION get_x_from_psi_index

END MODULE utility