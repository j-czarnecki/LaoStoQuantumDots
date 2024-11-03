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


    COMPLEX*16, ALLOCATABLE :: V_tilde_upper(:,:) !This matrix contains matrix elements of potential
                                                  !\tilde{V}(r_2) = < i(r_1) | V(r_1, r_2) | j(r_1) >.
                                                  !Since this matrix is hermitian, only upper triangle is stored.
    COMPLEX*16, ALLOCATABLE :: V_tilde_slice(:) !This is an element that should be passed to calculate interaction element
    INTEGER*4 :: v_tilde_elems
    INTEGER*4 :: i, j, k, l, nn
    COMPLEX*16 :: interaction_element


    v_tilde_elems = nstate_1*(nstate_1 + 1)/2 !Number of elements in upper triangle of hermitian matrix V_tilde
    ALLOCATE(V_tilde_upper(v_tilde_elems, ham_1_size))
    ALLOCATE(V_tilde_slice(ham_1_size))
    CALL CALCULATE_V_TILDE(Psi_1, ham_1_size, nstate_1, V_tilde_upper, v_tilde_elems, norbs, Nx, Ny, dx)

    Hamiltonian_2_crs = DCMPLX(0.0d0, 0.0d0)
    column_2_crs = 0
    row_2_crs = 0
    nn = 1
    DO i = 1, ham_2_size
      !PRINT*, "i = ", i
      WRITE(66,*) "nn / nonzero ", nn, nonzero_ham_2
      FLUSH(66)
      row_2_crs(i) = nn
      DO j = i, ham_2_size
        !!PRINT*, i, j
        !If statements' order determined by frequency of given elements
        IF (N_changed_indeces(i,j) == 2) THEN
          !Interaction elements
          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i,j,1,1),  Changed_indeces(i,j, 1, 2))
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)),&
                                                & Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:, Changed_indeces(i,j, 2, 2)),&
                                                & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
          Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element


          CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper,  ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i,j,1,1),  Changed_indeces(i,j, 2, 2))
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)),&
                                                & Psi_1(:, Changed_indeces(i,j, 2, 2)), Psi_1(:, Changed_indeces(i,j, 1, 2)),&
                                                & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
          Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) - interaction_element
          column_2_crs(nn) = j
          nn = nn + 1

        ELSE IF (N_changed_indeces(i,j) == 1) THEN
          !Interaction elements
          DO k = 1, k_electrons
            ! !PRINT*, Combinations(i, k)
            CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i,j,1,1),  Changed_indeces(i,j, 1, 2))
            CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)),&
                                                  & Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:,  Combinations(i,k)),&
                                                  & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
            Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element

            CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Changed_indeces(i,j,1,1),  Combinations(i,k))
            CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)),&
                                                  & Psi_1(:, Combinations(i,k)), Psi_1(:,  Changed_indeces(i,j, 1, 2)),&
                                                  & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
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
          DO k = 1, k_electrons - 1
            DO l = k + 1, k_electrons !Check whether sum could be reduced to sum_{k, l>k}. Then I will get rid of 0.5*
              IF (Combinations(i,k) /= Combinations(i,l)) THEN
                ! !PRINT*, Combinations(i,k), Combinations(i,l)
                CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper,  ham_1_size, nstate_1, v_tilde_elems, Combinations(i,k), Combinations(i,k))
                CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                                      & Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                                      & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
                Hamiltonian_2_crs(nn) = Hamiltonian_2_crs(nn) + interaction_element!Check whether it should be  * 0.5.

                CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, Combinations(i,k), Combinations(i,l))
                CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                                      & Psi_1(:, Combinations(i,l)), Psi_1(:, Combinations(i,k)),&
                                                      & V_tilde_slice(:), ham_1_size, interaction_element, norbs, eps_r)
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

    DEALLOCATE(V_tilde_upper)
    DEALLOCATE(V_tilde_slice)
  END SUBROUTINE CREATE_MANY_BODY_HAMILTONIAN_CRS


  SUBROUTINE CALCULATE_INTERACTION_ELEMENTS(Psi_1, Psi_2, Psi_3, Psi_4, V_tilde, ham_1_size, matrix_element, norbs, eps_r)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: ham_1_size, norbs
    REAL*8, INTENT(IN) :: eps_r
    COMPLEX*16, INTENT(OUT) :: matrix_element
    COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size), Psi_2(ham_1_size), Psi_3(ham_1_size), Psi_4(ham_1_size)
    COMPLEX*16, INTENT(IN) :: V_tilde(ham_1_size)

    INTEGER*4 :: r2, o2, oi, oj, ok, ol, i_so, j_so, k_so, l_so, si, sj, s2, so2
    REAL*8, PARAMETER :: epsilon(3) = [0.336, 0.306, 0.015] !eps(i,i,i,i), eps(i,j,i,j), eps(i,j,j,i)

    matrix_element = DCMPLX(0.0d0, 0.0d0)

    DO r2 = 1, ham_1_size, norbs
      DO s2 = 0, 1
        DO o2 = 0, norbs/2 - 1
          so2 = s2 + 2*o2
          matrix_element = matrix_element +&
          & CONJG(Psi_2(r2 + so2))*Psi_4(r2 + so2)*V_tilde(r2 + so2)
        END DO
      END DO

      !Calculating onsite integrals
      DO si = 0, 1
        DO sj = 0, 1

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
                    matrix_element = matrix_element + epsilon(2) * CONJG(Psi_1(r2 + i_so)*Psi_2(r2 + j_so))*Psi_3(r2 + k_so)*Psi_4(r2 + l_so)
                  ELSE IF (oi == ol .AND. oj == ok .AND. oi /= oj) THEN
                    matrix_element = matrix_element + epsilon(3) * CONJG(Psi_1(r2 + i_so)*Psi_2(r2 + j_so))*Psi_3(r2 + k_so)*Psi_4(r2 + l_so)
                  ELSE IF (oi == oj .AND.oi == ok .AND. oi == ol) THEN
                    matrix_element = matrix_element + epsilon(1) * CONJG(Psi_1(r2 + i_so)*Psi_2(r2 + j_so))*Psi_3(r2 + k_so)*Psi_4(r2 + l_so)
                  END IF
                  
                END DO
              END DO
            END DO
          END DO

        END DO
      END DO

    END DO
    matrix_element = matrix_element/eps_r
    !WRITE(66,*) 'Matrix element', matrix_element
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

    V_tilde_upper = (0.0d0, 0.0d0)
    !Loops over state
    !!$omp parallel do private(j, r2, s2, o2, so2, r1, x1, y1, x2, y2, r12, s1, o1, so1)
    DO i = 1, nstate_1
      WRITE(66,*) "V_tilde i = ", i
      FLUSH(66)
      DO j = i, nstate_1
        !Loop over position r_2 of second electron
        DO r2 = 1, ham_1_size, norbs
          x2 = get_x_from_psi_index(r2,norbs, Nx, dx)
          y2 = get_y_from_psi_index(r2,norbs, Nx, Ny, dx)

          DO s2 = 0, 1
            DO o2 = 0, norbs/2 - 1 
              so2 = s2 + 2*o2
              !Loop over position r_1 of first electron
              DO r1 = 1, ham_1_size, norbs
                x1 = get_x_from_psi_index(r1,norbs, Nx, dx)
                y1 = get_y_from_psi_index(r1,norbs, Nx, Ny, dx)

                r12 = SQRT((x1 - x2)**2 + (y1 - y2)**2)
                IF (r1 /= r2) THEN
                  DO s1 = 0, 1
                    DO o1 = 0, norbs/2 - 1
                      so1 = s1 + 2*o1
                      V_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), r2 + so2) =  V_tilde_upper(get_upper_hermitian_index(i, j, nstate_1), r2 + so2) + &
                      & 1. / r12 * CONJG(Psi_1(r1 + so1, i)) * Psi_1(r1 + so1, j) !Not dividing by eps_r, since it will be handled in matrix elements
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
    !!$omp end parallel do
  END SUBROUTINE CALCULATE_V_TILDE

  SUBROUTINE CALCULATE_PARTICLE_DENSITY(Particle_density, Psi_1, C_slater, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, n_states_1, n_states_2, k_electrons)
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: Particle_density(ham_1_size, n_states_2)
    COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, n_states_1), C_slater(ham_2_size, n_states_2)
    INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
    INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
    INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size
    INTEGER*4, INTENT(IN) :: k_electrons, n_states_1, n_states_2
    INTEGER*4 :: i, j, n, nstate
    REAL*8 :: norm

    Particle_density = 0.0d0
    norm = GAMMA(REAL(k_electrons)) !Since all sums are multiplied by (N - 1)!

    DO nstate = 1, n_states_2
      DO i = 1, ham_2_size
        DO n = 1, k_electrons
          Particle_density(:,nstate) = Particle_density(:,nstate) + norm*REAL(CONJG(C_slater(i,nstate))*C_slater(i,nstate)*CONJG(Psi_1(:,n))*Psi_1(:,n))
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

  COMPLEX*16 FUNCTION many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, Nx, dx, norbs)
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
    INTEGER*4 :: a, b, k

    many_body_x_expected_value = (0.0d0, 0.0d0)
    DO a = 1, ham_2_size
      DO b = 1, ham_2_size !Probably I can sum over upper triangle
        IF (N_changed_indeces(a,b) == 0) THEN
          DO k = 1, k_electrons
            many_body_x_expected_value = many_body_x_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
              & single_electron_x_expected_value(Psi_1(:, Combinations(a, k)), Psi_1(:, Combinations(b, k)), norbs, Nx, dx, ham_1_size)
          END DO
        ELSE IF (N_changed_indeces(a,b) == 1) THEN
          many_body_x_expected_value = many_body_x_expected_value + CONJG(C_slater(a,n))*C_slater(b,m)*&
            & single_electron_x_expected_value(Psi_1(:, Changed_indeces(a,b,1,1)), Psi_1(:, Changed_indeces(a,b,1,2)), norbs, Nx, dx, ham_1_size)
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
    INTEGER*4 :: a, b, k

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
    INTEGER*4 :: a, b, k

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

    INTEGER*4 :: a, b, k

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