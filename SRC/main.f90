
PROGRAM MAIN
  USE indata
  USE hamiltonian
  USE diagonalize

  IMPLICIT NONE

  COMPLEX*16, ALLOCATABLE :: Hamiltonian_1(:,:) !One electron hamiltonian of LAO-STO
  COMPLEX*16, ALLOCATABLE :: Hamiltonian_2(:,:) !Many-body hamiltonian
  REAL*8, ALLOCATABLE :: potential(:, :)
  INTEGER :: state_to_write
  !COMPLEX*16, ALLOCATABLE :: psi_lapack(:, :, :, :, :)
  !REAL*8, ALLOCATABLE :: ev_lapack(:)
  COMPLEX*16, ALLOCATABLE :: Psi_1(:, :), Psi_2(:,:)
  REAL*8, ALLOCATABLE :: Energies_1(:), Energies_2(:)
  INTEGER*4, ALLOCATABLE :: Combination_current(:)
  INTEGER*4,  ALLOCATABLE :: Combinations(:,:)
  INTEGER*1, ALLOCATABLE :: N_changed_indeces(:,:) !Number of changed indeces in given element of multi-body hamiltonian.
                                                 !Needed for Slater-Condon rules implementation.
                                                 !0 - no indeces chaned,
                                                 !1 - one index changed,
                                                 !2 - two indeces changed,
                                                 !3 - three or more indeces chnged.
  INTEGER*4, ALLOCATABLE :: Changed_indeces(:,:,:,:)
  LOGICAL, ALLOCATABLE :: Index_found(:)
  INTEGER*4 :: k_electrons, total
  INTEGER :: ix, iy, iorb
  INTEGER*4 :: i,j,k,l, n
  REAL*8 :: x, y, omega
  INTEGER*4 :: ham_1_size
  INTEGER*4 :: ham_2_size
  INTEGER*4 :: n_changed
  LOGICAL :: same_index
  REAL*8 :: slater_norm
  COMPLEX*16 :: interaction_element

  !read a file
  CALL INDATA_GET("../RUNS/quantum_dot.nml")
  WRITE (*, *) "Physical parameters:"
  WRITE (*, *) "th=", th / eV2au
  WRITE (*, *) "tl=", tl / eV2au
  WRITE (*, *) "td=", td / eV2au
  WRITE (*, *) "dso=", dso / eV2au
  WRITE (*, *) "drso=", drso / eV2au
  WRITE (*, *) "dE=", dE / eV2au
  WRITE (*, *) "g=", g
  WRITE (*, *) "Nx=", Nx
  WRITE (*, *) "Ny=", Ny
  WRITE (*, *) "dx=", dx / nm2au
  WRITE (*, *) "norbs=", norbs
  WRITE (*, *) "nstate=", nstate

  WRITE (*, *)
  WRITE (*, *) "External parameters:"
  WRITE (*, *) "Bx=", Bx / T2au
  WRITE (*, *) "By=", By / T2au
  WRITE (*, *) "Bz=", Bz / T2au

  k_electrons = 2
  ham_1_size = (2*Nx + 1)*(2*Ny + 1)*norbs
  !Size of two-electron hamiltonian is equal to number of combinations of one-electron wavefunctions
  !That we put into Slater determinants.
  ham_2_size = GAMMA(nstate + 1.0d0)/(GAMMA(k_electrons + 1.0d0) * GAMMA(nstate - k_electrons + 1.0d0))
  slater_norm = SQRT(GAMMA(k_electrons + 1.0d0))
  PRINT*, "Number of combinations: ", ham_2_size
  
  
  ALLOCATE(Combination_current(k_electrons))
  ALLOCATE (potential(-Nx:Nx, -Ny:Ny))
  !ALLOCATE (psi_lapack(nstate, -Nx:Nx, -Ny:Ny, norbs / 2, 2))
  !ALLOCATE (ev_lapack(nstate))
  !ALLOCATE (psi_arpack(nstate, -Nx:Nx, -Ny:Ny, norbs / 2, 2))
  ALLOCATE (Psi_1(ham_1_size, nstate))
  ALLOCATE (Energies_1(nstate))
  ALLOCATE (Psi_2(ham_2_size, nstate))
  ALLOCATE (Energies_2(nstate))
  ALLOCATE(Hamiltonian_1(ham_1_size, ham_1_size))
  ALLOCATE(Hamiltonian_2(ham_2_size, ham_2_size))
  ALLOCATE (N_changed_indeces(ham_2_size, ham_2_size))
  ALLOCATE (Combinations(ham_2_size, k_electrons))
  ALLOCATE (Changed_indeces(ham_2_size, ham_2_size, 2, 2))
  ALLOCATE( Index_found(k_electrons))


  Hamiltonian_1(:,:) = DCMPLX(0.0d0, 0.0d0)
  Hamiltonian_2(:,:) = DCMPLX(0.0d0, 0.0d0)
  N_changed_indeces(:,:) = 0
  Changed_indeces(:,:,:,:) = 0


  !#################### ONE-ELECTRON PROBLEM ####################
  omega = 3.0E-3 * eV2au !Julian: is this meant to be a fixed value? Maybe should be in input .nml file?
  potential = 0.0
  DO ix = -Nx, Nx
    DO iy = -Ny, Ny
      x = ix * dx
      y = iy * dx
      potential(ix, iy) = 0.5 * omega**2 * (x**2 + y**2)
    END DO
  END DO

  CALL HAMILTONIAN_CREATE(Hamiltonian_1(:,:), ham_1_size, Nx, Ny, norbs, potential)

  !Julian: I will split it a bit differently - hamiltonian should be initialized here in main.f90 and passed to DIAGONALIZE_XXX.
  !Then DIAGONALIZE_ARPACK will be more universal. It could reshape 1D eigenvector from ARPACK to a map like psi_arpack in main.f90 in the end of the subroutine.
  !CALL DIAGONALIZE_LAPACK(psi_lapack, ev_lapack, Nx, Ny, norbs, nstate, potential)
  !CALL DIAGONALIZE_ARPACK(psi_arpack, ev_arpack, Nx, Ny, norbs, nstate, potential)
  CALL DIAGONALIZE_ARPACK(Hamiltonian_1, ham_1_size, Psi_1, Energies_1, Nx, Ny, norbs, nstate)
  !###############################################################


  !#################### TWO-ELECTRON PROBLEM #####################
  !Many problem is solved using Galerkin method, where basis
  !is chosen as combinnation of Slater determinants.


  !Initialize combination_row to generate next sets
  CALL INIT_COMBINATION(Combination_current, k_electrons)
  DO i = 1, ham_2_size
    CALL GET_COMBINATION(Combination_current, nstate, k_electrons)
    Combinations(i, :) = Combination_current(:)
  END DO

  DO i = 1, ham_2_size
    DO j = 1, ham_2_size
      n_changed = 0

      !Check how many indeces are changed between combinations specifying row and column
      !For each element of row-combination check heter it exists in column-combination.
      !If not, increment N_changed_indexes by 1 and write index without match to Changed_indeces(i,j,n_chnanged,1)
      DO k = 1, k_electrons
        same_index = .FALSE.
        DO l = 1, k_electrons
          IF(Combinations(i,k) == Combinations(j,l)) THEN
            same_index = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. same_index) THEN
          n_changed = n_changed + 1
          IF (n_changed < 3) THEN
            Changed_indeces(i,j,n_changed,1) = Combinations(i,k)
          END IF
        END IF
      END DO
      N_changed_indeces(i,j) = MIN(n_changed, 3)

      n_changed = 0
      !For each element of column-combination check wheter it exists in row-combination.
      !TODO: Consider whther it could be done in a more efficient way
      DO k = 1, k_electrons
        same_index = .FALSE.
        DO l = 1, k_electrons
          IF(Combinations(j,k) == Combinations(i,l)) THEN
            same_index = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. same_index) THEN
          n_changed = n_changed + 1
          IF (n_changed < 3) THEN
            Changed_indeces(i,j,n_changed,2) = Combinations(j,k)
          END IF
        END IF
      END DO

      !PRINT*, N_changed_indeces(i,j)
      !PRINT*
      !PRINT*
      !WRITE(*,'(I2)', ADVANCE = 'NO') N_changed_indeces(i,j)
      !WRITE(*, *) 'i = ', i, ' j = ', j, ' N_changed  = ', N_changed_indeces(i,j), ' Unpaired row: ', (Changed_indeces(i,j,k,1), k = 1, 2), ' Unparied column: ', (Changed_indeces(i,j,k,2), k = 1, 2)
    END DO
    !WRITE(*,*)
  END DO

  PRINT*, "Constructing mult-body hamiltonian..."
  !CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:,1), Psi_1(:,2), Psi_1(:,3), Psi_1(:,4), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
  !Computing diagonal elements of Hamiltonian_2
  DO i = 1, ham_2_size
    !Single electron energies
    DO n = 1, k_electrons
      Hamiltonian_2(i,i) = Hamiltonian_2(i,i) + Energies_1(Combinations(i,n))
    END DO

    !Interaction elements
    DO k = 1, k_electrons
      DO l = 1, k_electrons !Check whether sum could be reduced to sum_{k, l>k}. Then I will get rid of 0.5*
        IF (Combinations(i,k) /= Combinations(i,l)) THEN
          !PRINT*, Combinations(i,k), Combinations(i,l)
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)), Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
          Hamiltonian_2(i,i) = Hamiltonian_2(i,i) + 0.5*interaction_element !Check whether it should be  * 0.5.
          !PRINT*, interaction_element

          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)), Psi_1(:, Combinations(i,l)), Psi_1(:, Combinations(i,k)), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
          Hamiltonian_2(i,i) = Hamiltonian_2(i,i) - 0.5*interaction_element !Check whether it should be  * 0.5
          !PRINT*, interaction_element
          !PRINT*
        END IF
      END DO
    END DO


  END DO


  DO i = 1, ham_2_size
    DO j = i + 1, ham_2_size
      PRINT*, i, j
      !If statements' order determined by frequency of given elements      
      IF (N_changed_indeces(i,j) == 2) THEN
        !Interaction elements
        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)), Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:, Changed_indeces(i,j, 2, 2)), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
        Hamiltonian_2(i,j) = Hamiltonian_2(i,j) + interaction_element

        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)), Psi_1(:, Changed_indeces(i,j, 2, 2)), Psi_1(:, Changed_indeces(i,j, 1, 2)), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
        Hamiltonian_2(i,j) = Hamiltonian_2(i,j) - interaction_element

      ELSE IF (N_changed_indeces(i,j) == 1) THEN
        !Single electron energies
        !Computes <psi_i|H_1|psi_p> where psi is one-electron wavefunction and index i was changed to index p.
        Hamiltonian_2(i,j) = Hamiltonian_2(i,j) + dx**2*DOT_PRODUCT(Psi_1(:, Changed_indeces(i,j,1,1)), MATMUL(Hamiltonian_1, Psi_1(:, Changed_indeces(i,j,1,2))))

        !Interaction elements
        DO k = 1, k_electrons
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)), Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:, Combinations(i,k)), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
          Hamiltonian_2(i,j) = Hamiltonian_2(i,j) + interaction_element

          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,k)), Psi_1(:, Changed_indeces(i,j, 1, 2)), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
          Hamiltonian_2(i,j) = Hamiltonian_2(i,j) - interaction_element

        END DO
      ELSE IF (N_changed_indeces(i,j) /= 3) THEN
        STOP
        !ERROR
      END IF


    END DO
  END DO


  ! DO i = 1, ham_2_size
  !   WRITE(*,*) (REAL(Hamiltonian_2(i,j)), j = 1, ham_2_size)
  ! END DO

  !Does not work for small size of ham_2_size
  CALL DIAGONALIZE_ARPACK(Hamiltonian_2, ham_2_size, Psi_2, Energies_2, Nx, Ny, norbs, nstate)




  !##########################################################




  !Julian: This is only data output. I will wrap it in some subroutines to make the main function cleaner.
  WRITE (*, *) "Porownanie energii LAPACK i ARPACK"
  DO i = 1, nstate
    WRITE (*, '(200e20.12)') Energies_1(i) / eV2au! ,ev_lapack(i) / eV2au
  END DO

  WRITE (*, *) "Multi-body energies"
  DO i = 1, nstate
    WRITE (*, '(200e20.12)') Energies_2(i) / eV2au! ,ev_lapack(i) / eV2au
  END DO

  ! state_to_write = 3
  ! OPEN (1, FILE="psi_lapack.dat")
  ! DO ix = -Nx, Nx
  !   DO iy = -Ny, Ny
  !     x = ix * dx
  !     y = iy * dx
  !     WRITE (1, *) x, y, (ABS(psi_lapack(state_to_write, ix, iy, iorb, 1))**2 + ABS(psi_lapack(state_to_write, ix, iy, iorb, 2))**2, iorb=1, 3)
  !   END DO
  ! END DO
  ! CLOSE (1)

  ! OPEN (1, FILE="psi_arpack.dat")
  ! DO ix = -Nx, Nx
  !   DO iy = -Ny, Ny
  !     x = ix * dx
  !     y = iy * dx
  !     WRITE (1, *) x, y, (ABS(psi_arpack(state_to_write, ix, iy, iorb, 1))**2 + ABS(psi_arpack(state_to_write, ix, iy, iorb, 2))**2, iorb=1, 3)
  !   END DO
  ! END DO
  ! CLOSE (1)

  DEALLOCATE (potential)
  !DEALLOCATE (psi_lapack)
  !DEALLOCATE (ev_lapack)
  DEALLOCATE (Psi_1)
  DEALLOCATE (Energies_1)
  DEALLOCATE (Psi_2)
  DEALLOCATE (Energies_2)
  DEALLOCATE (Hamiltonian_1)
  DEALLOCATE (Hamiltonian_2)
  DEALLOCATE (N_changed_indeces)
  DEALLOCATE (Combinations)
  DEALLOCATE (Combination_current)
END PROGRAM MAIN


SUBROUTINE GET_COMBINATION(Combination, N, k)
  !! Returns set of indices that form a k-element combination from set {1, 2, ..., N}.
  !! Should be initialized with an array conataining first combination with last element diminished by 1.
  !! Example: k=3,then in first iteration Combination = [1, 2, 2] should be passed and will return [1,2,3]
  !! It Could be viewed as traversing an upper triangle of N-dimensional tensor and writing all sets of indices.
  IMPLICIT NONE
  INTEGER*4, INTENT(OUT) :: Combination(k) !! Array containg set of indices forming a k-element combination
  INTEGER*4, INTENT(IN) :: N !! Number of elements in the set
  INTEGER*4, INTENT(IN) :: k !! Number of elements in combination
  
  INTEGER*4 :: i, j, m

  Combination(k) = Combination(k) + 1

  j = k
  i = 0
  DO  WHILE (j .GT. 1 .AND. Combination(j) .GT. N - i)
    Combination(j - 1) = Combination(j - 1) + 1
    !TODO: This probably could be done only once to reduce some spare operations.
    DO m = j, k
      Combination(m) = Combination(m - 1) + 1
    END DO
    j = j - 1
    i = i + 1
  END DO
  
END SUBROUTINE

SUBROUTINE INIT_COMBINATION(Combination, k)
  !! Initializes COmbination array so that it could be passed to GET_COBINATION subroutine
  IMPLICIT NONE
  INTEGER*4, INTENT(OUT) :: Combination(k) !! Array containg set of indices forming a k-element combination
                                           !! with last eleent shifted by -1
  INTEGER*4, INTENT(IN) :: k !! Number of elements in combination
  INTEGER*4 :: i

  DO i = 1, k - 1
    Combination(i) = i
  END DO
  Combination(k) = Combination(k - 1)

END SUBROUTINE

SUBROUTINE CALCULATE_INTERACTION_ELEMENTS(Psi_1, Psi_2, Psi_3, Psi_4, size, matrix_element, norbs, Nx, Ny, dx)
  USE constants
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: size, norbs, Nx, Ny
  REAL*8, INTENT(IN) :: dx
  COMPLEX*16, INTENT(OUT) :: matrix_element
  COMPLEX*16, INTENT(IN) :: Psi_1(size), Psi_2(size), Psi_3(size), Psi_4(size)

  REAL*8, EXTERNAL :: get_x_from_psi_index, get_y_from_psi_index

  INTEGER*4 :: i,j,k,l,m,n, si, sj, oi, oj, ok, ol
  INTEGER*4 :: i_so, j_so, k_so, l_so !index of spin-orbitals
  REAL*8 :: xi, yi, ri, xj, yj, rj, r_ij
  REAL*8, PARAMETER :: epsilon(3) = [0.336, 0.306, 0.015] !eps(i,i,i,i), eps(i,j,i,j), eps(i,j,j,i) 
  
  matrix_element = (0.0d0, 0.0d0)

  !Only two loops determining position due to two center approach.
  !delta(r_j1, r_j3) and delta(r_j2, r_j4)
  DO i = 1, size, norbs
    DO j = 1, size, norbs
      xi = get_x_from_psi_index(i,norbs, Nx, Ny, dx)
      yi = get_y_from_psi_index(i,norbs, Nx, Ny, dx)


      xj = get_x_from_psi_index(j,norbs, Nx, Ny, dx)
      yj = get_y_from_psi_index(j,norbs, Nx, Ny, dx)

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
                matrix_element = matrix_element + 1/r_ij*CONJG(Psi_1(i + i_so)*Psi_2(j + j_so))*Psi_3(i + i_so)*Psi_4(j + j_so)
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
        
                    IF (oi == ok .AND. oj == ol) THEN
                      matrix_element = matrix_element + epsilon(2) * CONJG(Psi_1(i + i_so)*Psi_2(j + j_so))*Psi_3(i + k_so)*Psi_4(j + l_so)
                    ELSE IF (oi == ol .AND. oj == ok) THEN
                      matrix_element = matrix_element + epsilon(3) * CONJG(Psi_1(i + i_so)*Psi_2(j + j_so))*Psi_3(i + k_so)*Psi_4(j + l_so)
                    ELSE IF (oi == oj .AND.oi == ok .AND. oi == ol) THEN
                      matrix_element = matrix_element + epsilon(1) * CONJG(Psi_1(i + i_so)*Psi_2(j + j_so))*Psi_3(i + k_so)*Psi_4(j + l_so)
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

  !PRINT*, matrix_element

END SUBROUTINE CALCULATE_INTERACTION_ELEMENTS

REAL*8 FUNCTION get_x_from_psi_index(i, norbs, Nx, Ny, dx)
  IMPLICIT NONE
  INTEGER*4 :: i, norbs, Nx, Ny
  REAL*8 :: dx
  get_x_from_psi_index = ((i/norbs)/(2*Ny + 1) - Nx) * dx
END FUNCTION get_x_from_psi_index

REAL*8 FUNCTION get_y_from_psi_index(i, norbs, Nx, Ny, dx)
  IMPLICIT NONE
  INTEGER*4 :: i, norbs, Nx, Ny
  REAL*8 :: dx
  get_y_from_psi_index = (MOD(i/norbs, 2*Ny + 1) - Ny) * dx
END FUNCTION get_y_from_psi_index