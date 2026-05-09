#include "macros_def.f90"
MODULE many_body
USE combinatory
USE utility
USe logger
USE omp_lib
IMPLICIT NONE
CONTAINS

SUBROUTINE CREATE_MANY_BODY_HAMILTONIAN_CRS(Hamiltonian_2_crs, column_2_crs, row_2_crs, N_changed_indeces, Changed_indeces, Combinations,&
  & N_ham_2_elems_in_prev_rows, Psi_1, Energies_1, ham_1_size, nstate_1, nonzero_ham_2, ham_2_size, k_electrons, norbs, Nx, Ny, dx, eps_r)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: ham_2_size, norbs, Nx, Ny, k_electrons, nonzero_ham_2
  REAL*8, INTENT(IN) :: dx, eps_r
  COMPLEX*16, INTENT(OUT) :: Hamiltonian_2_crs(nonzero_ham_2)
  INTEGER*4, INTENT(OUT) :: column_2_crs(nonzero_ham_2)
  INTEGER*4, INTENT(OUT) :: row_2_crs(ham_2_size + 1)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: N_ham_2_elems_in_prev_rows(ham_2_size)
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1)
  REAL*8, INTENT(IN) :: Energies_1(nstate_1)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstate_1

  COMPLEX*16, ALLOCATABLE :: V_tilde_upper(:, :) !This matrix contains matrix elements of potential
  !\tilde{V}(r_2) = < i(r_1) | V(r_1, r_2) | j(r_1) >.
  !Since this matrix is hermitian, only upper triangle is stored.
  COMPLEX*16, ALLOCATABLE :: V_tilde_slice(:) !This is an element that should be passed to calculate interaction element
  INTEGER*4 :: v_tilde_elems
  INTEGER*4 :: i, j, k, l, nn
  COMPLEX*16 :: interaction_element
  INTEGER*4 :: phase

  WRITE (log_string, *) "Creating many-body Hamiltonian"
  LOG_INFO(log_string)

  v_tilde_elems = (nstate_1 * (nstate_1 + 1)) / 2 !Number of elements in upper triangle of hermitian matrix V_tilde
  ALLOCATE (V_tilde_upper(v_tilde_elems, ham_1_size))
  ALLOCATE (V_tilde_slice(ham_1_size))
  CALL CALCULATE_V_TILDE(Psi_1, ham_1_size, nstate_1, V_tilde_upper, v_tilde_elems, norbs, Nx, Ny, dx)

  Hamiltonian_2_crs = DCMPLX(0.0d0, 0.0d0)
  column_2_crs = 0
  row_2_crs = 0

  OPEN (10, FILE='./OutputData/Coulomb_integrals.dat', ACTION='WRITE', FORM='FORMATTED')
  !$omp parallel private(nn, interaction_element, V_tilde_slice, phase, log_string)
  !$omp do schedule(dynamic, 1)
  DO i = 1, ham_2_size
    WRITE (log_string, *) "Ham_2 i = ", i
    LOG_INFO(log_string)

    nn = N_ham_2_elems_in_prev_rows(i)
    row_2_crs(i) = nn
    DO j = i, ham_2_size
        !!PRINT*, i, j
      !If statements' order determined by frequency of given elements
      IF (N_changed_indeces(i, j) == 2) THEN
        phase = get_parity_phase(Combinations(i, :), Combinations(j, :), N_changed_indeces(i, j), Changed_indeces(i, j, 1, 1),&
                                &Changed_indeces(i, j, 1, 2), Changed_indeces(i, j, 2, 1), Changed_indeces(i, j, 2, 2), k_electrons)

        !Interaction elements
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 1, 2))
        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Changed_indeces(i, j, 2, 1)),&
                                              & Psi_1(:, Changed_indeces(i, j, 1, 2)), Psi_1(:, Changed_indeces(i, j, 2, 2)),&
                                              & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
        Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element * phase

        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 2, 2))
        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Changed_indeces(i, j, 2, 1)),&
                                              & Psi_1(:, Changed_indeces(i, j, 2, 2)), Psi_1(:, Changed_indeces(i, j, 1, 2)),&
                                              & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
        Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) - interaction_element * phase
        WRITE (10, *) i, j, REAL(Hamiltonian_2_crs(nn)), AIMAG(Hamiltonian_2_crs(nn))
        column_2_crs(nn) = j
        nn = nn + 1

      ELSE IF (N_changed_indeces(i, j) == 1) THEN
        !Interaction elements
        DO k = 1, k_electrons
          ! !PRINT*, Combinations(i, k)
          phase = get_parity_phase(Combinations(i, :), Combinations(j, :), N_changed_indeces(i, j), Changed_indeces(i, j, 1, 1),&
                                  &Changed_indeces(i, j, 1, 2), Changed_indeces(i, j, 2, 1), Changed_indeces(i, j, 2, 2), k_electrons)

          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 1, 2))
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Combinations(i, k)),&
                                                & Psi_1(:, Changed_indeces(i, j, 1, 2)), Psi_1(:, Combinations(i, k)),&
                                                & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
          Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element * phase

          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Combinations(i, k))
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Combinations(i, k)),&
                                                & Psi_1(:, Combinations(i, k)), Psi_1(:, Changed_indeces(i, j, 1, 2)),&
                                                & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
          Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) - interaction_element * phase
        END DO
        WRITE (10, *) i, j, REAL(Hamiltonian_2_crs(nn)), AIMAG(Hamiltonian_2_crs(nn))
        column_2_crs(nn) = j
        nn = nn + 1
        !Diagonal elements
      ELSE IF (N_changed_indeces(i, j) == 0) THEN
        !Interaction elements
        DO k = 1, k_electrons - 1
          DO l = k + 1, k_electrons !Check whether sum could be reduced to sum_{k, l>k}. Then I will get rid of 0.5*
            IF (Combinations(i, k) /= Combinations(i, l)) THEN
              ! !PRINT*, Combinations(i,k), Combinations(i,l)
              CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Combinations(i, k), Combinations(i, k))
              CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i, k)), Psi_1(:, Combinations(i, l)),&
                                                    & Psi_1(:, Combinations(i, k)), Psi_1(:, Combinations(i, l)),&
                                                    & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
              Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element!Check whether it should be  * 0.5.

              CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Combinations(i, k), Combinations(i, l))
              CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i, k)), Psi_1(:, Combinations(i, l)),&
                                                    & Psi_1(:, Combinations(i, l)), Psi_1(:, Combinations(i, k)),&
                                                    & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
              Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) - interaction_element!Check whether it should be  * 0.5

            END IF
          END DO
        END DO
        WRITE (10, *) i, j, REAL(Hamiltonian_2_crs(nn)), AIMAG(Hamiltonian_2_crs(nn))

        DO k = 1, k_electrons
          Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + Energies_1(Combinations(i, k))
        END DO

        column_2_crs(nn) = j
        nn = nn + 1
      ELSE IF (N_changed_indeces(i, j) /= 3) THEN
        STOP
      END IF
    END DO
  END DO
  !$omp end do
  !$omp end parallel
  CLOSE (10)
  !PRINT*, "Nonzero elements traversed in ham 2", nn
  row_2_crs(ham_2_size + 1) = nonzero_ham_2 + 1 !Sanity check already done in INIT_PREV_ELEMS()

    !!!!!!!! ONLY FOR TESTING !!!!!!!!!!!!
  OPEN (unit=1, FILE='./OutputData/V_integrals.dat', FORM="FORMATTED", ACTION="WRITE")
  WRITE (1, *) '#No. state    V_{iiii}'
  DO i = 1, nstate_1
    WRITE (log_string, *) "V_{iiii} integral for i = ", i
    LOG_INFO(log_string)

    CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, i, i)
    CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, i), Psi_1(:, i),&
                                          & Psi_1(:, i), Psi_1(:, i),&
                                          & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
    WRITE (1, *) i, interaction_element
  END DO
  CLOSE (1)

  DEALLOCATE (V_tilde_upper)
  DEALLOCATE (V_tilde_slice)
END SUBROUTINE CREATE_MANY_BODY_HAMILTONIAN_CRS

SUBROUTINE CALCULATE_INTERACTION_ELEMENTS(Psi_1, Psi_2, Psi_3, Psi_4, V_tilde, ham_1_size, matrix_element, norbs, eps_r)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: ham_1_size, norbs
  REAL*8, INTENT(IN) :: eps_r
  COMPLEX*16, INTENT(OUT) :: matrix_element
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size), Psi_2(ham_1_size), Psi_3(ham_1_size), Psi_4(ham_1_size)
  COMPLEX*16, INTENT(IN) :: V_tilde(ham_1_size)

  INTEGER*4 :: r2, o2, oi, oj, ok, ol, i_so, j_so, k_so, l_so, si, sj, s2, so2, i_so_2, j_so_2
  REAL*8, PARAMETER :: epsilon(3) = [0.336, 0.306, 0.015] !eps(i,i,i,i), eps(i,j,i,j), eps(i,j,j,i)

  matrix_element = DCMPLX(0.0d0, 0.0d0)

  DO r2 = 1, ham_1_size, norbs
    DO s2 = 0, 1
      DO o2 = 0, norbs / 2 - 1
        so2 = s2 + 2 * o2
        matrix_element = matrix_element +&
        & CONJG(Psi_2(r2 + so2)) * Psi_4(r2 + so2) * V_tilde(r2 + so2)
      END DO
    END DO

    !Calculating onsite integrals
    !Calculating all the same orbitals
    DO oi = 0, norbs / 2 - 1
      DO si = 0, 1
        DO sj = 0, 1
          i_so = si + 2 * oi
          j_so = sj + 2 * oi
          matrix_element = matrix_element + epsilon(1) * CONJG(Psi_1(r2 + i_so) * Psi_2(r2 + j_so)) * Psi_3(r2 + i_so) * Psi_4(r2 + j_so)
        END DO
      END DO
    END DO
    !Calculating two different orbitals
    DO oi = 0, norbs / 2 - 1
      DO oj = 0, norbs / 2 - 1
        IF (oi /= oj) THEN
          DO si = 0, 1
            DO sj = 0, 1
              i_so = si + 2 * oi
              j_so = sj + 2 * oj
              matrix_element = matrix_element + epsilon(2) * CONJG(Psi_1(r2 + i_so) * Psi_2(r2 + j_so)) * Psi_3(r2 + i_so) * Psi_4(r2 + j_so)

              !To avoid interaction elements < psi_1 || psi_3 > with different spins
              i_so = si + 2 * oi
              j_so = si + 2 * oj
              i_so_2 = sj + 2 * oi
              j_so_2 = sj + 2 * oj
              matrix_element = matrix_element + epsilon(3) * CONJG(Psi_1(r2 + i_so) * Psi_2(r2 + j_so_2)) * Psi_3(r2 + j_so) * Psi_4(r2 + i_so_2)
            END DO
          END DO
        END IF
      END DO
    END DO

  END DO
  matrix_element = matrix_element / eps_r
END SUBROUTINE CALCULATE_INTERACTION_ELEMENTS

SUBROUTINE CALCULATE_V_TILDE(Psi_1, ham_1_size, nstate_1, V_tilde_upper, v_tilde_elems, norbs, Nx, Ny, dx)
  !! This subroutine integrates potential V(r_1, r_2) over r_1, treating r_2 as a parameter.
  !! \tilde{V}(r_2) = < i(r_1) | V(r_1, r_2) | j(r_1) >
  !! It stores the results in 1-D array, since resulting matrix is hermitian.
  IMPLICIT NONE
  COMPLEX*16, INTENT(OUT) :: V_tilde_upper(v_tilde_elems, ham_1_size)
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstate_1, v_tilde_elems
  INTEGER*4, INTENT(IN) :: norbs, Nx, Ny
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: r1, r2, i, j, s1, s2, o1, o2, so1, so2
  REAL*8 :: r12, x1, y1, x2, y2
  CHARACTER(LEN=500) :: filename
  CHARACTER(LEN=1000) :: filename_state
  INTEGER*4 :: nn, ix, iy, iorb, ii, is

  V_tilde_upper = (0.0d0, 0.0d0)
  !Loops over state
  !$omp parallel private (so1, so2, x1, y1, x2, y2, r12, log_string)
  !$omp do
  DO i = 1, nstate_1
    WRITE (log_string, *) "V_tilde i = ", i
    LOG_INFO(log_string)

    DO j = i, nstate_1
      !Loop over position r_2 of second electron
      DO r2 = 1, ham_1_size, norbs
        x2 = get_x_from_psi_index(r2, norbs, Nx, dx)
        y2 = get_y_from_psi_index(r2, norbs, Nx, Ny, dx)

        DO s2 = 0, 1
          DO o2 = 0, norbs / 2 - 1
            so2 = s2 + 2 * o2
            !Loop over position r_1 of first electron
            DO r1 = 1, ham_1_size, norbs
              x1 = get_x_from_psi_index(r1, norbs, Nx, dx)
              y1 = get_y_from_psi_index(r1, norbs, Nx, Ny, dx)

              r12 = SQRT((x1 - x2)**2 + (y1 - y2)**2)
              !r12 = 1.0d0 !JULIAN: CHANGED FOR TESTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (r1 /= r2) THEN
                DO s1 = 0, 1
                  DO o1 = 0, norbs / 2 - 1
                    so1 = s1 + 2 * o1
                    V_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), r2 + so2) = V_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), r2 + so2) + &
                    & 1./r12 * CONJG(Psi_1(r1 + so1, i)) * Psi_1(r1 + so1, j) !Not dividing by eps_r, since it will be handled in matrix elements
                  END DO
                END DO
              END IF
              !Not calculating onsite integrals, since they can be determined in integration over r_2.
              !Then we assume that each r_2 has had an r_1 hat was equal to it and calculate contribution to matrix element.
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  !$omp end do
  !$omp end parallel

  ! OPEN (10, FILE='./OutputData/V_tilde_integrated.dat', ACTION='WRITE', FORM='FORMATTED')
  ! filename = './OutputData/V_tilde_position_dependent'
  ! DO i = 1, nstate_1
  !   DO j = i, nstate_1
  !     WRITE (10, *) i, j, SUM(V_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), :))

  !     WRITE (filename_state, '(A, A, I0, A, I0, A)') TRIM(filename), '_i', i, '_j', j, '.dat'
  !     OPEN (1, FILE=filename_state)
  !     WRITE (1, *) '#x[nm] y[nm] Re(V_tilde(orb1, s = +)) Im(V_tilde(orb1, s = +)) Re(V_tilde(orb1, s = -)) ...'

  !     nn = 1
  !     DO iy = -Ny, Ny
  !     DO ix = -Nx, Nx
  !       WRITE (1, '(2E20.8)', ADVANCE='NO') ix * dx / nm2au, iy * dx / nm2au

  !       DO iorb = 1, norbs / 2
  !         DO ii = 1, 2
  !           WRITE (1, '(2E20.8)', ADVANCE='NO') REAL(V_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), nn)), AIMAG(V_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), nn))
  !           nn = nn + 1
  !         END DO
  !       END DO
  !       WRITE (1, *)
  !     END DO
  !     END DO

  !     CLOSE (1)

  !   END DO
  ! END DO
  ! CLOSE (10)

END SUBROUTINE CALCULATE_V_TILDE

SUBROUTINE CALCULATE_R_TILDE(Psi_1, ham_1_size, nstate_1, R_tilde_upper, r_tilde_elems, norbs, Nx, Ny, dx)
  !! This subroutine integrates relative position |r_1 - r_2| over r_1, treating r_2 as a parameter.
  !! \tilde{R}(r_2) = < i(r_1) | |r_1 - r_2| | j(r_1) >
  !! It stores the results in 1-D array, since resulting matrix is hermitian.
  IMPLICIT NONE
  COMPLEX*16, INTENT(OUT) :: R_tilde_upper(r_tilde_elems, ham_1_size)
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstate_1, r_tilde_elems
  INTEGER*4, INTENT(IN) :: norbs, Nx, Ny
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: r1, r2, i, j, s1, s2, o1, o2, so1, so2
  REAL*8 :: r12, x1, y1, x2, y2
  CHARACTER(LEN=500) :: filename
  CHARACTER(LEN=1000) :: filename_state
  INTEGER*4 :: nn, ix, iy, iorb, ii, is

  R_tilde_upper = (0.0d0, 0.0d0)
  !Loops over state
  !$omp parallel private (so1, so2, x1, y1, x2, y2, r12, log_string)
  !$omp do
  DO i = 1, nstate_1
    WRITE (log_string, *) "R_tilde i = ", i
    LOG_INFO(log_string)

    DO j = i, nstate_1
      !Loop over position r_2 of second electron
      DO r2 = 1, ham_1_size, norbs
        x2 = get_x_from_psi_index(r2, norbs, Nx, dx)
        y2 = get_y_from_psi_index(r2, norbs, Nx, Ny, dx)

        DO s2 = 0, 1
          DO o2 = 0, norbs / 2 - 1
            so2 = s2 + 2 * o2
            !Loop over position r_1 of first electron
            DO r1 = 1, ham_1_size, norbs
              x1 = get_x_from_psi_index(r1, norbs, Nx, dx)
              y1 = get_y_from_psi_index(r1, norbs, Nx, Ny, dx)

              r12 = SQRT((x1 - x2)**2 + (y1 - y2)**2)
              DO s1 = 0, 1
                DO o1 = 0, norbs / 2 - 1
                  so1 = s1 + 2 * o1
                  R_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), r2 + so2) = R_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), r2 + so2) + &
                  & r12 * CONJG(Psi_1(r1 + so1, i)) * Psi_1(r1 + so1, j)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  !$omp end do
  !$omp end parallel

END SUBROUTINE CALCULATE_R_TILDE

RECURSIVE COMPLEX * 16 FUNCTION many_body_relative_distance_expected_value(R_tilde_upper, Psi_1, C_slater, N_changed_indeces, Changed_indeces, Combinations, r_tilde_elems, ham_1_size, ham_2_size, nstate_1, nstate_2, n, m, k_electrons, norbs, Nx, Ny, dx)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, Nx, Ny, norbs, r_tilde_elems
  INTEGER*4, INTENT(IN) :: n, m
  COMPLEX*16, INTENT(IN) :: R_tilde_upper(r_tilde_elems, ham_1_size) !! This matrix contains matrix elements of relative positions
                                                                      !! $$
                                                                      !!\tilde{R}(r_2) = \int \varphi^*_n(r_1, r_2) |r_1 - r_2| \varphi_m(r_1,r_2) dr_1
                                                                      !!$$
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1) !! Single-electron wavefunctions $$\Psi_n(x,y)$$
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstate_2) !! Coefficient of given Slater determinant
  ! (as ordered in many-body Hamiltonian)
  ! for a given many-body state
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size) ! Number of single-particle wavefunction indices changed
  ! between  i-th and j-th Slater determinant
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2) !Specifies which index has been changed to which for example 1 -> 2 and 3 -> 4.
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: k_electrons, nstate_1, nstate_2
  REAL*8, INTENT(IN) :: dx

  COMPLEX*16 :: R_tilde_slice(ham_1_size) !This is an element that should be passed to calculate interaction element
  REAL*8 :: phase

  INTEGER*4 :: i, j, k, l, nstate

  many_body_relative_distance_expected_value = CMPLX(0.0d0, 0.0d0)

  DO i = 1, ham_2_size
    DO j = 1, ham_2_size
      IF (N_changed_indeces(i, j) == 2) THEN
        phase = get_parity_phase(Combinations(i, :), Combinations(j, :), N_changed_indeces(i, j), Changed_indeces(i, j, 1, 1),&
          & Changed_indeces(i, j, 1, 2), Changed_indeces(i, j, 2, 1), Changed_indeces(i, j, 2, 2), k_electrons)

        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(R_tilde_slice, R_tilde_upper, ham_1_size, nstate_1, r_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 1, 2))

        many_body_relative_distance_expected_value = many_body_relative_distance_expected_value + &
          & phase * SUM(CONJG(Psi_1(:, Changed_indeces(i, j, 2, 1))) * Psi_1(:, Changed_indeces(i, j, 2, 2)) * R_tilde_slice) * &
          & CONJG(C_slater(i, n)) * C_slater(j, m)

        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(R_tilde_slice, R_tilde_upper, ham_1_size, nstate_1, r_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 2, 2))

        many_body_relative_distance_expected_value = many_body_relative_distance_expected_value - &
          & phase * SUM(CONJG(Psi_1(:, Changed_indeces(i, j, 2, 1))) * Psi_1(:, Changed_indeces(i, j, 1, 2)) * R_tilde_slice) * &
          & CONJG(C_slater(i, n)) * C_slater(j, m)

      ELSE IF (N_Changed_indeces(i, j) == 1) THEN
        DO k = 1, k_electrons
          phase = get_parity_phase(Combinations(i, :), Combinations(j, :), N_changed_indeces(i, j), Changed_indeces(i, j, 1, 1),&
            & Changed_indeces(i, j, 1, 2), Changed_indeces(i, j, 2, 1), Changed_indeces(i, j, 2, 2), k_electrons)

          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(R_tilde_slice, R_tilde_upper, ham_1_size, nstate_1, r_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 1, 2))

          many_body_relative_distance_expected_value = many_body_relative_distance_expected_value + &
            & phase * SUM(CONJG(Psi_1(:, Combinations(i, k))) * Psi_1(:, Combinations(i, k)) * R_tilde_slice) * &
            & CONJG(C_slater(i, n)) * C_slater(j, m)

          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(R_tilde_slice, R_tilde_upper, ham_1_size, nstate_1, r_tilde_elems, Changed_indeces(i, j, 1, 1), Combinations(i, k))

          many_body_relative_distance_expected_value = many_body_relative_distance_expected_value - &
            & phase * SUM(CONJG(Psi_1(:, Combinations(i, k))) * Psi_1(:, Changed_indeces(i, j, 1, 2)) * R_tilde_slice) * &
            & CONJG(C_slater(i, n)) * C_slater(j, m)

        END DO
      ELSE IF (N_changed_indeces(i, j) == 0) THEN
        DO k = 1, k_electrons
          DO l = k + 1, k_electrons
            IF (Combinations(i, k) /= Combinations(i, l)) THEN
              phase = 1
              CALL GET_SLICE_FROM_HERMITIAN_MATRIX(R_tilde_slice, R_tilde_upper, ham_1_size, nstate_1, r_tilde_elems, Combinations(i, k), Combinations(i, k))

              many_body_relative_distance_expected_value = many_body_relative_distance_expected_value + &
                & phase * SUM(CONJG(Psi_1(:, Combinations(i, l))) * Psi_1(:, Combinations(i, l)) * R_tilde_slice) * &
                & CONJG(C_slater(i, n)) * C_slater(j, m)

              CALL GET_SLICE_FROM_HERMITIAN_MATRIX(R_tilde_slice, R_tilde_upper, ham_1_size, nstate_1, r_tilde_elems, Combinations(i, k), Combinations(i, l))

              many_body_relative_distance_expected_value = many_body_relative_distance_expected_value - &
                & phase * SUM(CONJG(Psi_1(:, Combinations(i, l))) * Psi_1(:, Combinations(i, k)) * R_tilde_slice) * &
                & CONJG(C_slater(i, n)) * C_slater(j, m)
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO

END FUNCTION many_body_relative_distance_expected_value

RECURSIVE COMPLEX * 16 FUNCTION many_body_coulomb_expected_value(V_tilde_upper, Psi_1, C_slater, N_changed_indeces, Changed_indeces, Combinations, v_tilde_elems, ham_1_size, ham_2_size, nstate_1, nstate_2, n, m, k_electrons, norbs, Nx, Ny, dx, eps_r)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, Nx, Ny, norbs, v_tilde_elems
  INTEGER*4, INTENT(IN) :: n, m
  REAL*8, INTENT(IN) :: eps_r
  COMPLEX*16, INTENT(IN) :: V_tilde_upper(v_tilde_elems, ham_1_size)
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1) ! Single-electron wavefunctions \Psi_n(x,y)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstate_2) ! Coefficient of given Slater determinant
  ! (as ordered in many-body Hamiltonian)
  ! for a given many-body state
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size) ! Number of single-particle wavefunction indices changed
  ! between  i-th and j-th Slater determinant
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2) !Specifies which index has been changed to which for example 1 -> 2 and 3 -> 4.
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: k_electrons, nstate_1, nstate_2
  REAL*8, INTENT(IN) :: dx

  COMPLEX*16 :: V_tilde_slice(ham_1_size)
  INTEGER*4 :: i, j, k, l
  COMPLEX*16 :: interaction_element
  INTEGER*4 :: phase

  many_body_coulomb_expected_value = 0.0d0

  DO i = 1, ham_2_size
    DO j = i, ham_2_size
      !If statements' order determined by frequency of given elements
      IF (N_changed_indeces(i, j) == 2) THEN
        phase = get_parity_phase(Combinations(i, :), Combinations(j, :), N_changed_indeces(i, j), Changed_indeces(i, j, 1, 1),&
                                &Changed_indeces(i, j, 1, 2), Changed_indeces(i, j, 2, 1), Changed_indeces(i, j, 2, 2), k_electrons)

        !Interaction elements
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 1, 2))
        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Changed_indeces(i, j, 2, 1)),&
                                              & Psi_1(:, Changed_indeces(i, j, 1, 2)), Psi_1(:, Changed_indeces(i, j, 2, 2)),&
                                              & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
        many_body_coulomb_expected_value = many_body_coulomb_expected_value + interaction_element * phase * CONJG(C_slater(i, n)) * C_slater(j, m)

        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 2, 2))
        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Changed_indeces(i, j, 2, 1)),&
                                              & Psi_1(:, Changed_indeces(i, j, 2, 2)), Psi_1(:, Changed_indeces(i, j, 1, 2)),&
                                              & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
        many_body_coulomb_expected_value = many_body_coulomb_expected_value - interaction_element * phase * CONJG(C_slater(i, n)) * C_slater(j, m)

      ELSE IF (N_changed_indeces(i, j) == 1) THEN
        !Interaction elements
        DO k = 1, k_electrons
          phase = get_parity_phase(Combinations(i, :), Combinations(j, :), N_changed_indeces(i, j), Changed_indeces(i, j, 1, 1),&
                                  &Changed_indeces(i, j, 1, 2), Changed_indeces(i, j, 2, 1), Changed_indeces(i, j, 2, 2), k_electrons)

          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Changed_indeces(i, j, 1, 2))
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Combinations(i, k)),&
                                                & Psi_1(:, Changed_indeces(i, j, 1, 2)), Psi_1(:, Combinations(i, k)),&
                                                & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
          many_body_coulomb_expected_value = many_body_coulomb_expected_value + interaction_element * phase * CONJG(C_slater(i, n)) * C_slater(j, m)

          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i, j, 1, 1), Combinations(i, k))
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i, j, 1, 1)), Psi_1(:, Combinations(i, k)),&
                                                & Psi_1(:, Combinations(i, k)), Psi_1(:, Changed_indeces(i, j, 1, 2)),&
                                                & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
          many_body_coulomb_expected_value = many_body_coulomb_expected_value - interaction_element * phase * CONJG(C_slater(i, n)) * C_slater(j, m)
        END DO
        !Diagonal elements
      ELSE IF (N_changed_indeces(i, j) == 0) THEN
        !COulomb interaction elements
        DO k = 1, k_electrons - 1
          DO l = k + 1, k_electrons !Check whether sum could be reduced to sum_{k, l>k}. Then I will get rid of 0.5*
            IF (Combinations(i, k) /= Combinations(i, l)) THEN
              CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Combinations(i, k), Combinations(i, k))
              CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i, k)), Psi_1(:, Combinations(i, l)),&
                                                    & Psi_1(:, Combinations(i, k)), Psi_1(:, Combinations(i, l)),&
                                                    & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
              many_body_coulomb_expected_value = many_body_coulomb_expected_value + interaction_element * CONJG(C_slater(i, n)) * C_slater(j, m)!Check whether it should be  * 0.5.

              CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Combinations(i, k), Combinations(i, l))
              CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i, k)), Psi_1(:, Combinations(i, l)),&
                                                    & Psi_1(:, Combinations(i, l)), Psi_1(:, Combinations(i, k)),&
                                                    & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
              many_body_coulomb_expected_value = many_body_coulomb_expected_value - interaction_element * CONJG(C_slater(i, n)) * C_slater(j, m)!Check whether it should be  * 0.5
            END IF
          END DO
        END DO
      ELSE IF (N_changed_indeces(i, j) /= 3) THEN
        STOP
      END IF
    END DO
  END DO

END FUNCTION many_body_coulomb_expected_value

SUBROUTINE CALCULATE_PARTICLE_DENSITY(Particle_density, Psi_1, C_slater, N_changed_indeces, Changed_indeces, Combinations, ham_1_size, ham_2_size, nstate_1, nstate_2, k_electrons, dx)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size
  REAL*8, INTENT(OUT) :: Particle_density(ham_1_size, nstate_2) ! Particle density \rho(x,y) for each many-body state
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1) ! Single-electron wavefunctions \Psi_n(x,y)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstate_2) ! Coefficient of given Slater determinant
  ! (as ordered in many-body Hamiltonian)
  ! for a given many-body state
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size) ! Number of single-particle wavefunction indices changed
  ! between  i-th and j-th Slater determinant
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2) !Specifies which index has been changed to which for example 1 -> 2 and 3 -> 4.
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: k_electrons, nstate_1, nstate_2
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: i, j, n, m, k
  REAL*8 :: phase
  REAL*8, PARAMETER :: tolerance = 1e-6
  COMPLEX*16 :: Density_matrix(ham_1_size, nstate_2, nstate_2)

  Particle_density = 0.0d0
  Density_matrix = CMPLX(0.0d0, 0.0d0)

  DO n = 1, nstate_2
    DO m = 1, nstate_2
      DO i = 1, ham_2_size
        DO j = 1, ham_2_size
          IF (N_changed_indeces(i, j) == 0) THEN
            DO k = 1, k_electrons
              Density_matrix(:, n, m) = Density_matrix(:, n, m) + CONJG(C_slater(i, n)) * C_slater(j, m) *&
                & CONJG(Psi_1(:, Combinations(i, k))) * Psi_1(:, Combinations(j, k))
            END DO
          ELSE IF (N_changed_indeces(i, j) == 1) THEN
            phase = get_parity_phase(Combinations(i, :), Combinations(j, :), N_changed_indeces(i, j), Changed_indeces(i, j, 1, 1),&
              & Changed_indeces(i, j, 1, 2), Changed_indeces(i, j, 2, 1), Changed_indeces(i, j, 2, 2), k_electrons)

            Density_matrix(:, n, m) = Density_matrix(:, n, m) + phase * CONJG(C_slater(i, n)) * C_slater(j, m) *&
              & CONJG(Psi_1(:, Changed_indeces(i, j, 1, 1))) * Psi_1(:, Changed_indeces(i, j, 1, 2))
          END IF
        END DO
      END DO
      WRITE (log_string, *) 'Density matrix for ', n, m, SUM(Density_matrix(:, n, m))
      LOG_INFO(log_string)
      IF (n .eq. m) THEN
        Particle_density(:, n) = REAL(Density_matrix(:, n, n))
        IF (MAXVAL(ABS(AIMAG(Density_matrix(:, n, n)))) > tolerance) THEN
          WRITE (log_string, *) "Imaginary part of particle density too big: ", MAXVAL(ABS(AIMAG(Density_matrix(:, n, n))))
          LOG_ABNORMAL(log_string)
        END IF
      END IF
    END DO
  END DO

END SUBROUTINE CALCULATE_PARTICLE_DENSITY

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_potential_expected_value(Potential, Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, Nx, Ny, dx, norbs)
  !! Calculates matrix element $$\langle n|V (r)|m \rangle$$, where n and m denote multi-body wavefunctions and
  !! $$ V(r)$$ is the potential

  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, Nx, Ny, norbs, n, m
  INTEGER*4, INTENT(IN) :: k_electrons, nstates_1, nstates_2
  REAL*8, INTENT(IN) :: Potential(-Nx:Nx, -Ny:Ny) !! Potential $$V(r)$$
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1) !! Single-electron wavefunctions $$\Psi_n(x,y)$$
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2) !! Coefficient of given Slater determinant
                                                            !! (as ordered in many-body Hamiltonian)
                                                            !! for a given many-body state
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size) !! Number of single-particle wavefunction indices changed
                                                                     !! between  i-th and j-th Slater determinant
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2) !! Specifies which index has been changed to which for example 1 -> 2 and 3 -> 4.
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: a, b, k
  REAL*8 :: phase

  many_body_potential_expected_value = CMPLX(0.0d0, 0.0d0)

  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_potential_expected_value = many_body_potential_expected_value + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & potential_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), Potential, ham_1_size, Nx, Ny, norbs)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
          & Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)
        many_body_potential_expected_value = many_body_potential_expected_value + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & potential_expected_value(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), Potential, ham_1_size, Nx, Ny, norbs)
      END IF
    END DO
  END DO

END FUNCTION many_body_potential_expected_value

RECURSIVE COMPLEX * 16 FUNCTION many_body_hamiltonian_expected_value(Ham_1_crs, col_crs, row_crs, nonzero_ham_1, Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: nonzero_ham_1
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: Ham_1_crs(nonzero_ham_1)
  INTEGER*4, INTENT(IN) :: col_crs(nonzero_ham_1)
  INTEGER*4, INTENT(IN) :: row_crs(ham_1_size + 1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2) !! Coefficient of given Slater determinant
                                                            !! (as ordered in many-body Hamiltonian)
                                                            !! for a given many-body state
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size) !! Number of single-particle wavefunction indices changed
                                                                     !! between  i-th and j-th Slater determinant
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2) !! Specifies which index has been changed to which for example 1 -> 2 and 3 -> 4.
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4 :: a, b, k
  REAL*8 :: phase

  many_body_hamiltonian_expected_value = CMPLX(0.0d0, 0.0d0)

  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_hamiltonian_expected_value = many_body_hamiltonian_expected_value + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & hamiltonian_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), Ham_1_crs, col_crs, row_crs, nonzero_ham_1, ham_1_size)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
          & Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)
        many_body_hamiltonian_expected_value = many_body_hamiltonian_expected_value + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & hamiltonian_expected_value(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), Ham_1_crs, col_crs, row_crs, nonzero_ham_1, ham_1_size)
      END IF
    END DO
  END DO

END FUNCTION many_body_hamiltonian_expected_value

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, Nx, dx, norbs)
  !! Calculates matrix element of <n|X|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! X = \sum_i^{k_electrons} x_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4, INTENT(IN) :: Nx, norbs
  REAL*8, INTENT(IN) :: dx
  INTEGER*4 :: phase
  INTEGER*4 :: a, b, k

  many_body_x_expected_value = (0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_x_expected_value = many_body_x_expected_value + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & single_electron_x_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), norbs, Nx, dx, ham_1_size)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_x_expected_value = many_body_x_expected_value + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & single_electron_x_expected_value(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), norbs, Nx, dx, ham_1_size)
      END IF
    END DO
  END DO

END FUNCTION many_body_x_expected_value

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m)
  !! Calculates matrix element of <n|s_x|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_x = \sum_i^{k_electrons} s_x_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase
  INTEGER*4 :: a, b, k

  many_body_sigma_x_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_x_expected_value = many_body_sigma_x_expected_value + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_x_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_x_expected_value = many_body_sigma_x_expected_value + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_x_expected_value(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_x_expected_value

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_y_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m)
  !! Calculates matrix element of <n|s_y|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_y = \sum_i^{k_electrons} s_y_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase
  INTEGER*4 :: a, b, k

  many_body_sigma_y_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_y_expected_value = many_body_sigma_y_expected_value + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_y_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_y_expected_value = many_body_sigma_y_expected_value + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_y_expected_value(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_y_expected_value

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_z_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m)
  !! Calculates matrix element of <n|s_z|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_z = \sum_i^{k_electrons} s_z_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase

  INTEGER*4 :: a, b, k

  many_body_sigma_z_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_z_expected_value = many_body_sigma_z_expected_value + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_z_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_z_expected_value = many_body_sigma_z_expected_value + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_z_expected_value(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_z_expected_value

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_parity_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, norbs, Nx, Ny, n, m)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, Nx, Ny
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: a, k
  COMPLEX*16 :: current_combination_parity
  many_body_parity_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    current_combination_parity = 1
    DO k = 1, k_electrons
      current_combination_parity = current_combination_parity * single_electron_parity(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(a, k)), ham_1_size, norbs, Nx, Ny)
    END DO
    many_body_parity_expected_value = many_body_parity_expected_value + CONJG(C_slater(a, n)) * C_slater(a, m) * current_combination_parity
  END DO

END FUNCTION

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_d_xy_up_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, norbs, Nx, Ny, n, m)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, Nx, Ny
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: a, k
  COMPLEX*16 :: share
  many_body_d_xy_up_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    share = 0
    DO k = 1, k_electrons
      share = share + d_xy_up_share(Psi_1(:, Combinations(a, k)), ham_1_size, norbs)
    END DO
    many_body_d_xy_up_expected_value = many_body_d_xy_up_expected_value + CONJG(C_slater(a, n)) * C_slater(a, m) * share
  END DO
END FUNCTION

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_d_xy_down_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, norbs, Nx, Ny, n, m)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, Nx, Ny
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: a, k
  COMPLEX*16 :: share
  many_body_d_xy_down_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    share = 0
    DO k = 1, k_electrons
      share = share + d_xy_down_share(Psi_1(:, Combinations(a, k)), ham_1_size, norbs)
    END DO
    many_body_d_xy_down_expected_value = many_body_d_xy_down_expected_value + CONJG(C_slater(a, n)) * C_slater(a, m) * share
  END DO
END FUNCTION

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_d_xz_up_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, norbs, Nx, Ny, n, m)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, Nx, Ny
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: a, k
  COMPLEX*16 :: share
  many_body_d_xz_up_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    share = 0
    DO k = 1, k_electrons
      share = share + d_xz_up_share(Psi_1(:, Combinations(a, k)), ham_1_size, norbs)
    END DO
    many_body_d_xz_up_expected_value = many_body_d_xz_up_expected_value + CONJG(C_slater(a, n)) * C_slater(a, m) * share
  END DO
END FUNCTION

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_d_xz_down_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, norbs, Nx, Ny, n, m)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, Nx, Ny
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: a, k
  COMPLEX*16 :: share
  many_body_d_xz_down_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    share = 0
    DO k = 1, k_electrons
      share = share + d_xz_down_share(Psi_1(:, Combinations(a, k)), ham_1_size, norbs)
    END DO
    many_body_d_xz_down_expected_value = many_body_d_xz_down_expected_value + CONJG(C_slater(a, n)) * C_slater(a, m) * share
  END DO
END FUNCTION

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_d_yz_up_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, norbs, Nx, Ny, n, m)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, Nx, Ny
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: a, k
  COMPLEX*16 :: share
  many_body_d_yz_up_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    share = 0
    DO k = 1, k_electrons
      share = share + d_yz_up_share(Psi_1(:, Combinations(a, k)), ham_1_size, norbs)
    END DO
    many_body_d_yz_up_expected_value = many_body_d_yz_up_expected_value + CONJG(C_slater(a, n)) * C_slater(a, m) * share
  END DO
END FUNCTION

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_d_yz_down_expected_value(Psi_1, C_slater, Combinations, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, norbs, Nx, Ny, n, m)
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons, Nx, Ny
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: a, k
  COMPLEX*16 :: share
  many_body_d_yz_down_expected_value = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    share = 0
    DO k = 1, k_electrons
      share = share + d_yz_down_share(Psi_1(:, Combinations(a, k)), ham_1_size, norbs)
    END DO
    many_body_d_yz_down_expected_value = many_body_d_yz_down_expected_value + CONJG(C_slater(a, n)) * C_slater(a, m) * share
  END DO
END FUNCTION

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_x_expected_value_R(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs)
  !! Calculates matrix element of <n|s_z|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_z = \sum_i^{k_electrons} s_z_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, nx, ny, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase

  INTEGER*4 :: a, b, k

  many_body_sigma_x_expected_value_R = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_x_expected_value_R = many_body_sigma_x_expected_value_R + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_x_expected_value_R(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size, nx, ny, norbs)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_x_expected_value_R = many_body_sigma_x_expected_value_R + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_x_expected_value_R(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size, nx, ny, norbs)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_x_expected_value_R

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_x_expected_value_L(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs)
  !! Calculates matrix element of <n|s_z|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_z = \sum_i^{k_electrons} s_z_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, nx, ny, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase

  INTEGER*4 :: a, b, k

  many_body_sigma_x_expected_value_L = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_x_expected_value_L = many_body_sigma_x_expected_value_L + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_x_expected_value_L(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size, nx, ny, norbs)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_x_expected_value_L = many_body_sigma_x_expected_value_L + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_x_expected_value_L(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size, nx, ny, norbs)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_x_expected_value_L

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_y_expected_value_R(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs)
  !! Calculates matrix element of <n|s_z|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_z = \sum_i^{k_electrons} s_z_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, nx, ny, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase

  INTEGER*4 :: a, b, k

  many_body_sigma_y_expected_value_R = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_y_expected_value_R = many_body_sigma_y_expected_value_R + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_y_expected_value_R(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size, nx, ny, norbs)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_y_expected_value_R = many_body_sigma_y_expected_value_R + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_y_expected_value_R(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size, nx, ny, norbs)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_y_expected_value_R

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_y_expected_value_L(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs)
  !! Calculates matrix element of <n|s_z|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_z = \sum_i^{k_electrons} s_z_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, nx, ny, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase

  INTEGER*4 :: a, b, k

  many_body_sigma_y_expected_value_L = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_y_expected_value_L = many_body_sigma_y_expected_value_L + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_y_expected_value_L(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size, nx, ny, norbs)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_y_expected_value_L = many_body_sigma_y_expected_value_L + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_y_expected_value_L(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size, nx, ny, norbs)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_y_expected_value_L

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_z_expected_value_R(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs)
  !! Calculates matrix element of <n|s_z|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_z = \sum_i^{k_electrons} s_z_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, nx, ny, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase

  INTEGER*4 :: a, b, k

  many_body_sigma_z_expected_value_R = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_z_expected_value_R = many_body_sigma_z_expected_value_R + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_z_expected_value_R(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size, nx, ny, norbs)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_z_expected_value_R = many_body_sigma_z_expected_value_R + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_z_expected_value_R(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size, nx, ny, norbs)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_z_expected_value_R

PURE RECURSIVE COMPLEX * 16 FUNCTION many_body_sigma_z_expected_value_L(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs)
  !! Calculates matrix element of <n|s_z|m>, where n and m denote multi-body wavefunctions and position operator x is defined as
  !! S_z = \sum_i^{k_electrons} s_z_i.
  IMPLICIT NONE
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
  INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, nx, ny, norbs
  INTEGER*4, INTENT(IN) :: n, m !Two many-body states which expected value should be calculated <n|x|m>
  INTEGER*4 :: phase

  INTEGER*4 :: a, b, k

  many_body_sigma_z_expected_value_L = DCMPLX(0.0d0, 0.0d0)
  DO a = 1, ham_2_size
    DO b = 1, ham_2_size !Probably I can sum over upper triangle
      IF (N_changed_indeces(a, b) == 0) THEN
        DO k = 1, k_electrons
          many_body_sigma_z_expected_value_L = many_body_sigma_z_expected_value_L + CONJG(C_slater(a, n)) * C_slater(b, m) *&
            & sigma_z_expected_value_L(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size, nx, ny, norbs)
        END DO
      ELSE IF (N_changed_indeces(a, b) == 1) THEN
        phase = get_parity_phase(Combinations(a, :), Combinations(b, :), N_changed_indeces(a, b), Changed_indeces(a, b, 1, 1),&
                                &Changed_indeces(a, b, 1, 2), Changed_indeces(a, b, 2, 1), Changed_indeces(a, b, 2, 2), k_electrons)

        many_body_sigma_z_expected_value_L = many_body_sigma_z_expected_value_L + phase * CONJG(C_slater(a, n)) * C_slater(b, m) *&
          & sigma_z_expected_value_L(Psi_1(:, Changed_indeces(a, b, 1, 1)), Psi_1(:, Changed_indeces(a, b, 1, 2)), ham_1_size, nx, ny, norbs)
      END IF
    END DO
  END DO

END FUNCTION many_body_sigma_z_expected_value_L
END MODULE many_body
