MODULE many_body
  USE combinatory
  USE utility
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE CREATE_MANY_BODY_HAMILTONIAN_CRS(Hamiltonian_2_crs, column_2_crs, row_2_crs, N_changed_indeces, Changed_indeces, Combinations,&
    & Psi_1, Energies_1, ham_1_size, nstate_1, nonzero_ham_2, ham_2_size, k_electrons, norbs, Nx, Ny, dx, eps_r)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: ham_2_size, norbs, Nx, Ny, k_electrons, nonzero_ham_2
    REAL*8, INTENT(IN) :: dx, eps_r
    COMPLEX*16, INTENT(OUT) :: Hamiltonian_2_crs(nonzero_ham_2)
    INTEGER*4, INTENT(OUT) :: column_2_crs(nonzero_ham_2)
    INTEGER*4, INTENT(OUT) :: row_2_crs(ham_2_size + 1)
    INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
    INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
    INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
    COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1)
    REAL*8, INTENT(IN) :: Energies_1(nstate_1)
    INTEGER*4, INTENT(IN) :: ham_1_size, nstate_1

    INTEGER*4 :: i, j, k, l, nn
    COMPLEX*16 :: interaction_element
    Hamiltonian_2_crs = DCMPLX(0.0d0, 0.0d0)
    column_2_crs = 0
    row_2_crs = 0
    nn = 1
    DO i = 1, ham_2_size
      !PRINT*, "i = ", i
      WRITE(66,*) "nn / nonzero ", nn, nonzero_ham_2
      row_2_crs(i) = nn
      DO j = i, ham_2_size
        !!PRINT*, i, j
        !If statements' order determined by frequency of given elements
        IF (N_changed_indeces(i,j) == 2) THEN
          !Interaction elements
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)),&
                                            & Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:, Changed_indeces(i,j, 2, 2)),&
                                            & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
          Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element


          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)),&
                                            & Psi_1(:, Changed_indeces(i,j, 2, 2)), Psi_1(:, Changed_indeces(i,j, 1, 2)),&
                                            & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
          Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) - interaction_element
          column_2_crs(nn) = j
          nn = nn + 1
        ELSE IF (N_changed_indeces(i,j) == 1) THEN
          !Interaction elements
          DO k = 1, k_electrons
            ! !PRINT*, Combinations(i, k)
            CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)),&
                                              & Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:, Combinations(i,k)),&
                                              & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
            Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element

            CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)),&
                                              & Psi_1(:, Combinations(i,k)), Psi_1(:, Changed_indeces(i,j, 1, 2)),&
                                              & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
            Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) - interaction_element
          END DO
          column_2_crs(nn) = j
          nn = nn + 1
        !Diagonal elements
        ELSE IF (N_changed_indeces(i,j) == 0) THEN
          DO k = 1, k_electrons
            Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + Energies_1(Combinations(i,k))
          END DO

          !Interaction elements
          DO k = 1, k_electrons
            DO l = k + 1, k_electrons !Check whether sum could be reduced to sum_{k, l>k}. Then I will get rid of 0.5*
              IF (Combinations(i,k) /= Combinations(i,l)) THEN
                ! !PRINT*, Combinations(i,k), Combinations(i,l)
                CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                                  & Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                                  &ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
                Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element!Check whether it should be  * 0.5.
                CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                                  & Psi_1(:, Combinations(i,l)), Psi_1(:, Combinations(i,k)),&
                                                  & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
                Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) - interaction_element!Check whether it should be  * 0.5
              END IF
            END DO
          END DO
          column_2_crs(nn) = j
          nn = nn + 1
        ELSE IF (N_changed_indeces(i,j) /= 3) THEN
          STOP
        END IF
      END DO
    END DO
    !PRINT*, "Nonzero elements traversed in ham 2", nn
    row_2_crs(ham_2_size + 1) = nn

  END SUBROUTINE CREATE_MANY_BODY_HAMILTONIAN_CRS

  SUBROUTINE CALCULATE_INTERACTION_ELEMENTS(Psi_1, Psi_2, Psi_3, Psi_4, size, matrix_element, norbs, Nx, Ny, dx, eps_r)
    !! Computes integral of the form <psi_1 psi_2 | 1/r_{ij} | psi_3 psi_4>
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: size, norbs, Nx, Ny
    REAL*8, INTENT(IN) :: dx, eps_r
    COMPLEX*16, INTENT(OUT) :: matrix_element
    COMPLEX*16, INTENT(IN) :: Psi_1(size), Psi_2(size), Psi_3(size), Psi_4(size)

    !REAL*8, EXTERNAL :: get_x_from_psi_index, get_y_from_psi_index

    INTEGER*4 :: ri,rj,k,l,m,n, si, sj, oi, oj, ok, ol
    INTEGER*4 :: i_so, j_so, k_so, l_so !index of spin-orbitals
    REAL*8 :: xi, yi, xj, yj, r_ij
    REAL*8, PARAMETER :: epsilon(3) = [0.336, 0.306, 0.015] !eps(i,i,i,i), eps(i,j,i,j), eps(i,j,j,i)

    matrix_element = (0.0d0, 0.0d0)

    !Only two loops determining position due to two center approach.
    !delta(r_j1, r_j3) and delta(r_j2, r_j4)
    DO ri = 1, size, norbs
      DO rj = 1, size, norbs
        ! xi = get_x_from_psi_index(ri,norbs, Nx, Ny, dx)
        ! yi = get_y_from_psi_index(ri,norbs, Ny, dx)

        ! xj = get_x_from_psi_index(rj,norbs, Nx, Ny, dx)
        ! yj = get_y_from_psi_index(rj,norbs, Ny, dx)
        xi = get_x_from_psi_index(ri,norbs, Nx, dx)
        yi = get_y_from_psi_index(ri,norbs, Nx, Ny, dx)

        xj = get_x_from_psi_index(rj,norbs, Nx, dx)
        yj = get_y_from_psi_index(rj,norbs, Nx, Ny, dx)
        r_ij = SQRT((xi - xj)**2 + (yi - yj)**2)
        !Only iteratinng over two spins, because we assume orthogonality  of spin states
        !delta(s_j1, s_j3) and delta(s_j2, s_j4)
        DO si = 0, 1
          DO sj = 0, 1

            IF (r_ij /= 0.0d0) THEN

              !delta(o_j1, o_j3) and delta(o_j2, o_j4)
              DO oi = 0, norbs/2 - 1
                DO oj = 0, norbs/2 - 1
                  i_so = si + 2*oi
                  j_so = sj + 2*oj
                  !!PRINT*, si, sj, oi, oj, i_so, j_so
                  matrix_element = matrix_element + 1/r_ij*CONJG(Psi_1(ri + i_so)*Psi_2(rj + j_so))*Psi_3(ri + i_so)*Psi_4(rj + j_so)
                END DO
              END DO
            ELSE

              !This loop can be unrolled to and ifs can be removed
              DO oi = 0, norbs/2 - 1
                DO oj = 0, norbs/2 - 1
                  DO ok = 0, norbs/2 - 1
                    DO ol = 0, norbs/2 - 1
                      i_so = si + 2*oi
                      j_so = sj + 2*oj
                      k_so = si + 2*ok
                      l_so = sj + 2*ol

                      !Those ifs can be unrolled
                      IF (oi == ok .AND. oj == ol .AND. oi /= oj) THEN
                        matrix_element = matrix_element + epsilon(2) * CONJG(Psi_1(ri + i_so)*Psi_2(rj + j_so))*Psi_3(ri + k_so)*Psi_4(rj + l_so)
                      ELSE IF (oi == ol .AND. oj == ok .AND. oi /= oj) THEN
                        matrix_element = matrix_element + epsilon(3) * CONJG(Psi_1(ri + i_so)*Psi_2(rj + j_so))*Psi_3(ri + k_so)*Psi_4(rj + l_so)
                      ELSE IF (oi == oj .AND.oi == ok .AND. oi == ol) THEN
                        matrix_element = matrix_element + epsilon(1) * CONJG(Psi_1(ri + i_so)*Psi_2(rj + j_so))*Psi_3(ri + k_so)*Psi_4(rj + l_so)
                      END IF

                    END DO
                  END DO
                END DO
              END DO

            END IF

          END DO
        END DO

      END DO
    END DO

    matrix_element = matrix_element/eps_r
    !!PRINT*, matrix_element

  END SUBROUTINE CALCULATE_INTERACTION_ELEMENTS

  SUBROUTINE CALCULATE_PARTICLE_DENSITY(Particle_density, Psi_1, C_slater, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, Nx, Ny, n_states_1, n_states_2, dx, k_electrons)
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: Particle_density(ham_1_size, n_states_2)
    COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, n_states_1), C_slater(ham_2_size, n_states_2)
    INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
    INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
    INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size
    INTEGER*4, INTENT(IN) :: Nx, Ny, k_electrons, n_states_1, n_states_2
    REAL*8, INTENT(IN) :: dx
    INTEGER*4 :: i, j, n, nstate
    REAL*8 :: norm

    Particle_density = 0.0d0
    norm = GAMMA(REAL(k_electrons)) !Since all sums are multiplied by (N - 1)!

    DO nstate = 1, n_states_2
      DO i = 1, ham_2_size
        DO n = 1, k_electrons
          Particle_density(:,nstate) = Particle_density(:,nstate) + norm*CONJG(C_slater(i,nstate))*C_slater(i,nstate)*CONJG(Psi_1(:,n))*Psi_1(:,n)
        END DO
      END DO

      DO i = 1, ham_2_size
        DO j = i + 1, ham_2_size
          IF (N_changed_indeces(i,j) == 1) THEN
            Particle_density(:,nstate) = Particle_density(:,nstate) + REAL(norm*CONJG(C_slater(i,nstate))*C_slater(j,nstate)*CONJG(C_slater(j, Changed_indeces(i,j,1,1)))*C_slater(j, Changed_indeces(i,j,1,2)))
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE CALCULATE_PARTICLE_DENSITY

  COMPLEX*16 FUNCTION many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, Nx, Ny, dx, norbs)
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
    INTEGER*4, INTENT(IN) :: Nx, Ny, norbs
    REAL*8, INTENT(IN) :: dx
    INTEGER*4 :: i, j, a, b, k

    many_body_x_expected_value = (0.0d0, 0.0d0)
    DO a = 1, ham_2_size
      DO b = 1, ham_2_size !Probably I can sum over upper triangle
        IF (N_changed_indeces(a,b) == 0) THEN
          DO k = 1, k_electrons
            many_body_x_expected_value = many_body_x_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
              & single_electron_x_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), norbs, Nx, Ny, dx, ham_1_size)
          END DO
        ELSE IF (N_changed_indeces(a,b) == 1) THEN
          many_body_x_expected_value = many_body_x_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
            & single_electron_x_expected_value(Psi_1(:, Changed_indeces(a,b,1,1)), Psi_1(:, Changed_indeces(a,b,1,2)), norbs, Nx, Ny, dx, ham_1_size)
        END IF
      END DO
    END DO

  END FUNCTION many_body_x_expected_value


  COMPLEX*16 FUNCTION many_body_sigma_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m)
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
    INTEGER*4 :: i, j, a, b, k

    many_body_sigma_x_expected_value = DCMPLX(0.0d0, 0.0d0)
    DO a = 1, ham_2_size
      DO b = 1, ham_2_size !Probably I can sum over upper triangle
        IF (N_changed_indeces(a,b) == 0) THEN
          DO k = 1, k_electrons
            many_body_sigma_x_expected_value = many_body_sigma_x_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
              & sigma_x_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size)
          END DO
        ELSE IF (N_changed_indeces(a,b) == 1) THEN
          many_body_sigma_x_expected_value = many_body_sigma_x_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
            & sigma_x_expected_value(Psi_1(:, Changed_indeces(a,b,1,1)), Psi_1(:, Changed_indeces(a,b,1,2)), ham_1_size)
        END IF
      END DO
    END DO

  END FUNCTION many_body_sigma_x_expected_value

  COMPLEX*16 FUNCTION many_body_sigma_y_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m)
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
    INTEGER*4 :: i, j, a, b, k

    many_body_sigma_y_expected_value = DCMPLX(0.0d0, 0.0d0)
    DO a = 1, ham_2_size
      DO b = 1, ham_2_size !Probably I can sum over upper triangle
        IF (N_changed_indeces(a,b) == 0) THEN
          DO k = 1, k_electrons
            many_body_sigma_y_expected_value = many_body_sigma_y_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
              & sigma_y_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size)
          END DO
        ELSE IF (N_changed_indeces(a,b) == 1) THEN
          many_body_sigma_y_expected_value = many_body_sigma_y_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
            & sigma_y_expected_value(Psi_1(:, Changed_indeces(a,b,1,1)), Psi_1(:, Changed_indeces(a,b,1,2)), ham_1_size)
        END IF
      END DO
    END DO

  END FUNCTION many_body_sigma_y_expected_value

  COMPLEX*16 FUNCTION many_body_sigma_z_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m)
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

    INTEGER*4 :: i, j, a, b, k

    many_body_sigma_z_expected_value = DCMPLX(0.0d0, 0.0d0)
    DO a = 1, ham_2_size
      DO b = 1, ham_2_size !Probably I can sum over upper triangle
        IF (N_changed_indeces(a,b) == 0) THEN
          DO k = 1, k_electrons
            many_body_sigma_z_expected_value = many_body_sigma_z_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
              & sigma_z_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), ham_1_size)
          END DO
        ELSE IF (N_changed_indeces(a,b) == 1) THEN
          many_body_sigma_z_expected_value = many_body_sigma_z_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
            & sigma_z_expected_value(Psi_1(:, Changed_indeces(a,b,1,1)), Psi_1(:, Changed_indeces(a,b,1,2)), ham_1_size)
        END IF
      END DO
    END DO

  END FUNCTION many_body_sigma_z_expected_value

END MODULE many_body