#include "macros_def.f90"
MODULE writers
USE constants
USE indata
USE utility
USE many_body
USE logger

IMPLICIT NONE
CONTAINS

!########### Universal writers
SUBROUTINE WRITE_ENERGIES(Energies, nstates, filename)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Energies(nstates)
  INTEGER*4, INTENT(IN) :: nstates
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: i

  WRITE (log_string, *) "Saving energies to file "//TRIM(filename)
  LOG_INFO(log_string)

  OPEN (unit=1, FILE=filename, FORM="FORMATTED", ACTION="WRITE")
  WRITE (1, *) '#No. sate [-]   Energy [meV]'
  DO i = 1, nstates
    WRITE (1, *) i, Energies(i) / eV2au * 1e3
  END DO
  CLOSE (1)

END SUBROUTINE WRITE_ENERGIES

SUBROUTINE WRITE_POTENTIAL(potential, Nx, Ny, filename)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: potential(-Nx:Nx, -Ny:Ny)
  INTEGER*4, INTENT(IN) :: Nx, Ny
  INTEGER*4 :: iy, ix
  CHARACTER(LEN=*), INTENT(IN) :: filename

  WRITE (log_string, *) "Saving potential to file "//TRIM(filename)
  LOG_INFO(log_string)

  OPEN (unit=1, FILE=filename, FORM="FORMATTED", ACTION="WRITE")
  WRITE (1, *) '#x[nm] y[nm] potential[meV]'
  DO iy = -Ny, Ny
  DO ix = -Nx, Nx
    WRITE (1, '(3E20.8)') ix * dx / nm2au, iy * dx / nm2au, potential(ix, iy) / eV2au
  END DO
  END DO
  CLOSE (1)
END SUBROUTINE

!########### Single-electron writers
SUBROUTINE WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi, psi_size, nstates, norbs, Nx, Ny, dx, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi(psi_size, nstates)
  INTEGER*4, INTENT(IN) :: Nx, Ny, psi_size, nstates, norbs
  REAL*8, INTENT(IN) :: dx
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=200) :: filename_state

  INTEGER*4 :: nn, ix, iy, is, i, iorb

  WRITE (log_string, *) "Saving wavefunctions to file "//TRIM(filename)
  LOG_INFO(log_string)

  DO is = 1, nstates
    WRITE (filename_state, '(2A, I0, A)') filename, '_n', is, '.dat'
    OPEN (1, FILE=filename_state)
    WRITE (1, *) '#x[nm] y[nm] Re(Psi(orb1, s = +)) Im(Psi(orb1, s = +)) Re(Psi(orb1, s = -)) ...'

    nn = 1
    DO iy = -Ny, Ny
    DO ix = -Nx, Nx
      WRITE (1, '(2E20.8)', ADVANCE='NO') ix * dx / nm2au, iy * dx / nm2au

      DO iorb = 1, norbs / 2
        DO i = 1, 2
          WRITE (1, '(2E20.8)', ADVANCE='NO') REAL(Psi(nn, is)), AIMAG(Psi(nn, is))
          nn = nn + 1
        END DO
      END DO
      WRITE (1, *)
    END DO
    END DO

    CLOSE (1)
  END DO

END SUBROUTINE WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS

SUBROUTINE WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1, ham_1_size, nstates, norbitals, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstates, norbitals
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: n
  CHARACTER(LEN=200) :: format_string

  WRITE (log_string, *) "Saving 1e expectations to file "//TRIM(filename)
  LOG_INFO(log_string)

  !Writing expectations for a given state
  format_string = '(I10, 14E20.8)'
  OPEN (unit=9, FILE=filename, FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) "#No. state [-]  <s_x>   <s_y>   <s_z>   <d_xy_up>   <d_xy_down>   <d_xz_up>   <d_xz_down>   <d_yz_up>   <d_yz_down>   <parity>   <x>   <y>"
  DO n = 1, nstates
    WRITE (9, format_string) n,&
      & REAL(sigma_x_expected_value(Psi_1(:, n), Psi_1(:, n), ham_1_size)),&
      & REAL(sigma_y_expected_value(Psi_1(:, n), Psi_1(:, n), ham_1_size)),&
      & REAL(sigma_z_expected_value(Psi_1(:, n), Psi_1(:, n), ham_1_size)), &
      & REAL(d_xy_up_share(Psi_1(:, n), ham_1_size, norbitals)), REAL(d_xy_down_share(Psi_1(:, n), ham_1_size, norbitals)),&
      & REAL(d_xz_up_share(Psi_1(:, n), ham_1_size, norbitals)), REAL(d_xz_down_share(Psi_1(:, n), ham_1_size, norbitals)),&
      & REAL(d_yz_up_share(Psi_1(:, n), ham_1_size, norbitals)), REAL(d_yz_down_share(Psi_1(:, n), ham_1_size, norbitals)),&
      & REAL(single_electron_parity(Psi_1(:, n), Psi_1(:, n), ham_1_size, norbitals, Nx, Ny)), &
      & REAL(single_electron_x_expected_value(Psi_1(:, n), Psi_1(:, n), norbitals, Nx, dx, ham_1_size)), &
      & REAL(single_electron_y_expected_value(Psi_1(:, n), Psi_1(:, n), norbitals, Nx, Ny, dx, ham_1_size)), &
      & REAL(sigma_z_expected_value_L(Psi_1(:, n), Psi_1(:, n), ham_1_size, Nx, Ny, 6)), &
      & REAL(sigma_z_expected_value_R(Psi_1(:, n), Psi_1(:, n), ham_1_size, Nx, Ny, 6))
  END DO
  CLOSE (9)

END SUBROUTINE WRITE_SINGLE_ELECTRON_EXPECTATIONS

!########### Multi-electron writers
SUBROUTINE WRITE_SLATER_COEFFICIENTS(C_slater, ham_2_size, nstates, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates)
  INTEGER*4, INTENT(IN) :: ham_2_size, nstates
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: i, j
  CHARACTER(LEN=200) :: format_string

  WRITE (log_string, *) "Saving Slater coeffs to file "//TRIM(filename)
  LOG_INFO(log_string)

  !Writing expectations for a given state
  format_string = '(2I10, 2E20.8)'
  OPEN (unit=9, FILE=filename, FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) "#No. multielectron state [-]    No. combination [-]    Re(c_slater)    Im(c_slater)"
  DO i = 1, nstates
    DO j = 1, ham_2_size
      WRITE (9, format_string) i, j, REAL(C_slater(j, i)), AIMAG(C_slater(j, i))
    END DO
    WRITE (9, *)
    WRITE (9, *)
  END DO
  CLOSE (9)

END SUBROUTINE WRITE_SLATER_COEFFICIENTS

SUBROUTINE WRITE_MULTI_ELECTRON_EXPECTATIONS(Ham_1_crs, col_crs, row_crs, nonzero_ham_1, V_image, V_confinement, Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, dx, filename)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nonzero_ham_1
  COMPLEX*16, INTENT(IN) :: Ham_1_crs(nonzero_ham_1)
  INTEGER*4, INTENT(IN) :: col_crs(nonzero_ham_1)
  INTEGER*4, INTENT(IN) :: row_crs(ham_1_size + 1)
  REAL*8, INTENT(IN) :: V_image(-Nx:Nx, -Ny:Ny)
  REAL*8, INTENT(IN) :: V_confinement(-Nx:Nx, -Ny:Ny)
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstate_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: nstate_1, nstate_2, norbs
  INTEGER*4, INTENT(IN) :: Nx, Ny
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: tilde_upper_triangle_elems
  COMPLEX*16, ALLOCATABLE :: R_tilde_upper(:, :)
  COMPLEX*16, ALLOCATABLE :: V_tilde_upper(:, :)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: n
  CHARACTER(LEN=200) :: format_string
  REAL*8, PARAMETER :: imag_tol = 1.0D-6  ! tolerance for imaginary part check

  ! Complex buffers for each expectation value (imaginary part must be ~0 for physical observables)
  COMPLEX*16 :: x_expected        ! mean position <x> [a.u.]
  COMPLEX*16 :: spin_x_expected   ! total spin x-component <S_x>
  COMPLEX*16 :: spin_y_expected   ! total spin y-component <S_y>
  COMPLEX*16 :: spin_z_expected   ! total spin z-component <S_z>
  COMPLEX*16 :: d_xy_up_expected   ! d_xy spin-up orbital population
  COMPLEX*16 :: d_xy_down_expected ! d_xy spin-down orbital population
  COMPLEX*16 :: d_xz_up_expected   ! d_xz spin-up orbital population
  COMPLEX*16 :: d_xz_down_expected ! d_xz spin-down orbital population
  COMPLEX*16 :: d_yz_up_expected   ! d_yz spin-up orbital population
  COMPLEX*16 :: d_yz_down_expected ! d_yz spin-down orbital population
  COMPLEX*16 :: parity_expected    ! spatial parity
  COMPLEX*16 :: spin_z_L_expected  ! left-dot spin z-component <S_z^L>
  COMPLEX*16 :: spin_z_R_expected  ! right-dot spin z-component <S_z^R>
  COMPLEX*16 :: v_image_expected   ! mean interaction energy with image <V_{image}>
  COMPLEX*16 :: v_confinement_expected ! mean interaction energy with confinement potential  <V_{conf}>
  COMPLEX*16 :: coulomb_interaction_expected ! mean value of Coulomb interaction
  COMPLEX*16 :: relative_distance_expected ! mean relative distance between particles
  COMPLEX*16 :: h1_no_potential_expected_value ! Mean value of the single-electron Hamiltonian without potentials
  tilde_upper_triangle_elems = (nstate_1 * (nstate_1 + 1)) / 2 !Number of elements in upper triangle of hermitian matrix V_tilde

  ALLOCATE (R_tilde_upper(tilde_upper_triangle_elems, ham_1_size))
  ALLOCATE (V_tilde_upper(tilde_upper_triangle_elems, ham_1_size))

  WRITE (log_string, *) "Calculating R_tilde for expectation values"
  LOG_INFO(log_string)
  CALL CALCULATE_R_TILDE(Psi_1, ham_1_size, nstate_1, R_tilde_upper, tilde_upper_triangle_elems, norbs, Nx, Ny, dx)

  WRITE (log_string, *) "Calculating V_tilde for expectation values"
  LOG_INFO(log_string)
  CALL CALCULATE_V_TILDE(Psi_1, ham_1_size, nstate_1, V_tilde_upper, tilde_upper_triangle_elems, norbs, Nx, Ny, dx)

  WRITE (log_string, *) "Saving 2e expectations to file "//TRIM(filename)
  LOG_INFO(log_string)

  !Writing expectations for a given state
  format_string = '(I10, 20E20.8)'
  OPEN (unit=9, FILE=filename, FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) "#No. state [-]  <x>[nm]    <S_x>   <S_y>   <S_z>    <d_xy_up>   <d_xy_down>   <d_xz_up>   <d_xz_down>   <d_yz_up>   <d_yz_down>   <parity>   <S_z^L>   <S_z^R>   <V_image>[eV]   <V_confinement>[eV]   <V_Coulomb>[eV]   <r_12>   <H_0(V=0)>[eV]   <H_tot>[eV]"
  DO n = 1, nstate_2

    ! Compute all expectation values for state n
    x_expected = many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n, Nx, dx, norbs)
    spin_x_expected = many_body_sigma_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n)
    spin_y_expected = many_body_sigma_y_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n)
    spin_z_expected = many_body_sigma_z_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n)
    d_xy_up_expected = many_body_d_xy_up_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n)
    d_xy_down_expected = many_body_d_xy_down_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n)
    d_xz_up_expected = many_body_d_xz_up_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n)
    d_xz_down_expected = many_body_d_xz_down_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n)
    d_yz_up_expected = many_body_d_yz_up_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n)
    d_yz_down_expected = many_body_d_yz_down_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n)
    parity_expected = many_body_parity_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, n, n)
    spin_z_L_expected = many_body_sigma_z_expected_value_L(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n, Nx, Ny, norbs)
    spin_z_R_expected = many_body_sigma_z_expected_value_R(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n, Nx, Ny, norbs)
    v_image_expected = many_body_potential_expected_value(V_image, Psi_1, C_slater, Combinations, N_Changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n, Nx, Ny, dx, norbs)
    v_confinement_expected = many_body_potential_expected_value(V_confinement, Psi_1, C_slater, Combinations, N_Changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n, Nx, Ny, dx, norbs)
    coulomb_interaction_expected = many_body_coulomb_expected_value(V_tilde_upper, Psi_1, C_slater, N_Changed_indeces, Changed_indeces, Combinations, tilde_upper_triangle_elems, ham_1_size, ham_2_size, nstate_1, nstate_2, n, n, k_electrons, norbs, Nx, Ny, dx, eps_r)
    relative_distance_expected = many_body_relative_distance_expected_value(R_tilde_upper, Psi_1, C_slater, N_Changed_indeces, Changed_indeces, Combinations, tilde_upper_triangle_elems, ham_1_size, ham_2_size, nstate_1, nstate_2, n, n, k_electrons, norbs, Nx, Ny, dx)
    h1_no_potential_expected_value = many_body_hamiltonian_expected_value(Ham_1_crs, col_crs, row_crs, nonzero_ham_1, Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, n)
    ! Verify that imaginary parts are negligible (physical observables must be real)
    IF (ABS(AIMAG(x_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(x_expected) for state n =', n, ':', AIMAG(x_expected), 'Re(x_expected): ', REAL(x_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(spin_x_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(spin_x_expected) for state n =', n, ':', AIMAG(spin_x_expected), 'Re(spin_x_expected): ', REAL(spin_x_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(spin_y_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(spin_y_expected) for state n =', n, ':', AIMAG(spin_y_expected), 'Re(spin_y_expected): ', REAL(spin_y_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(spin_z_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(spin_z_expected) for state n =', n, ':', AIMAG(spin_z_expected), 'Re(spin_z_expected): ', REAL(spin_z_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(d_xy_up_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(d_xy_up_expected) for state n =', n, ':', AIMAG(d_xy_up_expected), 'Re(d_xy_up_expected): ', REAL(d_xy_up_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(d_xy_down_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(d_xy_down_expected) for state n =', n, ':', AIMAG(d_xy_down_expected), 'Re(d_xy_down_expected): ', REAL(d_xy_down_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(d_xz_up_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(d_xz_up_expected) for state n =', n, ':', AIMAG(d_xz_up_expected), 'Re(d_xz_up_expected): ', REAL(d_xz_up_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(d_xz_down_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(d_xz_down_expected) for state n =', n, ':', AIMAG(d_xz_down_expected), 'Re(d_xz_down_expected): ', REAL(d_xz_down_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(d_yz_up_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(d_yz_up_expected) for state n =', n, ':', AIMAG(d_yz_up_expected), 'Re(d_yz_up_expected): ', REAL(d_yz_up_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(d_yz_down_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(d_yz_down_expected) for state n =', n, ':', AIMAG(d_yz_down_expected), 'Re(d_yz_down_expected): ', REAL(d_yz_down_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(parity_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(parity_expected) for state n =', n, ':', AIMAG(parity_expected), 'Re(parity_expected): ', REAL(parity_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(spin_z_L_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(spin_z_L_expected) for state n =', n, ':', AIMAG(spin_z_L_expected), 'Re(spin_z_L_expected): ', REAL(spin_z_L_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(spin_z_R_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(spin_z_R_expected) for state n =', n, ':', AIMAG(spin_z_R_expected), 'Re(spin_z_R_expected): ', REAL(spin_z_R_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(v_image_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(v_image_expected) for state n =', n, ':', AIMAG(v_image_expected), 'Re(v_image_expected): ', REAL(v_image_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(v_confinement_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(v_confinement_expected) for state n =', n, ':', AIMAG(v_confinement_expected), 'Re(v_confinement_expected): ', REAL(v_confinement_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(coulomb_interaction_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(coulomb_interaction_expected) for state n =', n, ':', AIMAG(coulomb_interaction_expected), 'Re(coulomb_interaction_expected): ', REAL(coulomb_interaction_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(relative_distance_expected)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(relative_distance_expected) for state n =', n, ':', AIMAG(relative_distance_expected), 'Re(relative_distance_expected): ', REAL(relative_distance_expected)
      LOG_ABNORMAL(log_string)
    END IF
    IF (ABS(AIMAG(h1_no_potential_expected_value)) > imag_tol) THEN
      WRITE (log_string, *) 'Non-zero Im(h1_no_potential_expected_value) for state n =', n, ':', AIMAG(h1_no_potential_expected_value), 'Re(h1_no_potential_expected_value): ', REAL(h1_no_potential_expected_value)
      LOG_ABNORMAL(log_string)
    END IF
    ! Write all expectation values for state n in a single statement;
    ! x_expected is converted from a.u. to nm and divided by k_electrons (per-electron average)
    WRITE (9, format_string) n,&
    & REAL(x_expected) / nm2au / k_electrons,&
    & REAL(spin_x_expected),&
    & REAL(spin_y_expected),&
    & REAL(spin_z_expected),&
    & REAL(d_xy_up_expected),&
    & REAL(d_xy_down_expected),&
    & REAL(d_xz_up_expected),&
    & REAL(d_xz_down_expected),&
    & REAL(d_yz_up_expected),&
    & REAL(d_yz_down_expected),&
    & REAL(parity_expected),&
    & REAL(spin_z_L_expected),&
    & REAL(spin_z_R_expected),&
    & REAL(v_image_expected) / eV2au,&
    & REAL(v_confinement_expected) / eV2au,&
    & REAL(coulomb_interaction_expected) / eV2au,&
    & REAL(relative_distance_expected) / nm2au, &
    & REAL(h1_no_potential_expected_value) / eV2au, &
    & REAL(h1_no_potential_expected_value + coulomb_interaction_expected) / eV2au

  END DO
  CLOSE (9)

  DEALLOCATE (R_tilde_upper)
  DEALLOCATE (V_tilde_upper)

END SUBROUTINE WRITE_MULTI_ELECTRON_EXPECTATIONS

!########### Time-dependent calculations writers

SUBROUTINE WRITE_TIME_EVOLUTION(Spin_t, t_max_int, nmax, filename)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Spin_t(t_max_int, nmax)
  INTEGER*4, INTENT(IN) ::t_max_int, nmax
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=200) :: format_string
  INTEGER*4 :: ti, ni
  REAL*8 :: t

  WRITE (log_string, *) "Saving time evolution to file "//TRIM(filename)
  LOG_INFO(log_string)

  format_string = '(3E20.8)'
  OPEN (unit=9, FILE=filename, FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) "#Time [ns] s_x(L) s_x(R) s_y(L) s_y(R) s_z(L) s_z(R)"
  DO ti = 1, t_max_int
    t = ti * dt
    WRITE (9, '(E20.8)', ADVANCE="NO") t / ns2au
    DO ni = 1, nmax
      WRITE (9, '(E20.8)', ADVANCE="NO") REAL(Spin_t(ti, ni))
    END DO
    WRITE (9, *)   ! nowa linia
  END DO

  CLOSE (9)

END SUBROUTINE WRITE_TIME_EVOLUTION

END MODULE writers
