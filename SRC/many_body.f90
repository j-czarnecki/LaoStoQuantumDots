MODULE many_body
  USE combinatory
  IMPLICIT NONE
  CONTAINS


  SUBROUTINE CALCULATE_INTERACTION_ELEMENTS(Psi_1, Psi_2, Psi_3, Psi_4, size, matrix_element, norbs, Nx, Ny, dx, eps_r)
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
        xi = get_x_from_psi_index(ri,norbs, Nx, Ny, dx)
        yi = get_y_from_psi_index(ri,norbs, Ny, dx)


        xj = get_x_from_psi_index(rj,norbs, Nx, Ny, dx)
        yj = get_y_from_psi_index(rj,norbs, Ny, dx)

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
                  !PRINT*, si, sj, oi, oj, i_so, j_so
                  matrix_element = matrix_element + 1/r_ij*CONJG(Psi_1(ri + i_so)*Psi_2(rj + j_so))*Psi_3(ri + i_so)*Psi_4(rj + j_so)
                END DO
              END DO
            ELSE

              DO oi = 0, norbs/2 - 1
                DO oj = 0, norbs/2 - 1
                  DO ok = 0, norbs/2 - 1
                    DO ol = 0, norbs/2 - 1
                      i_so = si + 2*oi
                      j_so = sj + 2*oj
                      k_so = si + 2*ok
                      l_so = sj + 2*ol
          
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

    matrix_element = matrix_element*dx**2/eps_r
    !PRINT*, matrix_element

  END SUBROUTINE CALCULATE_INTERACTION_ELEMENTS

  SUBROUTINE CALCULATE_PARTICLE_DENSITY(Particle_density, Psi_1, Psi_2, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, Nx, Ny, n_states_1, n_states_2, dx, k_electrons)
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: Particle_density(ham_1_size, n_states_2)
    COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, n_states_1), Psi_2(ham_2_size, n_states_2)
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
          Particle_density(:,nstate) = Particle_density(:,nstate) + norm*CONJG(Psi_2(i,nstate))*Psi_2(i,nstate)*CONJG(Psi_1(:,n))*Psi_1(:,n)
        END DO
      END DO

      DO i = 1, ham_2_size
        DO j = i + 1, ham_2_size
          IF (N_changed_indeces(i,j) == 1) THEN
            Particle_density(:,nstate) = Particle_density(:,nstate) + REAL(norm*CONJG(Psi_2(i,nstate))*Psi_2(j,nstate)*CONJG(Psi_2(j, Changed_indeces(i,j,1,1)))*Psi_2(j, Changed_indeces(i,j,1,2)))
          END IF
        END DO
      END DO
    END DO
    
  END SUBROUTINE CALCULATE_PARTICLE_DENSITY

  REAL*8 FUNCTION get_x_from_psi_index(i, norbs, Nx, Ny, dx)
    IMPLICIT NONE
    INTEGER*4 :: i, norbs, Nx, Ny
    REAL*8 :: dx
    get_x_from_psi_index = ((i/norbs)/(2*Ny + 1) - Nx) * dx
  END FUNCTION get_x_from_psi_index

  REAL*8 FUNCTION get_y_from_psi_index(i, norbs, Ny, dx)
    IMPLICIT NONE
    INTEGER*4 :: i, norbs, Ny
    REAL*8 :: dx
    get_y_from_psi_index = (MOD(i/norbs, 2*Ny + 1) - Ny) * dx
  END FUNCTION get_y_from_psi_index
END MODULE many_body