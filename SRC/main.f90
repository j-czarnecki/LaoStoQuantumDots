
PROGRAM MAIN
  USE indata
  USE hamiltonian
  USE diagonalize
  USE combinatory
  USE many_body
  USE writers
  USE constants
  USE utility


  IMPLICIT NONE

  !COMPLEX*16, ALLOCATABLE :: Hamiltonian_1(:,:) !One electron hamiltonian of LAO-STO
  COMPLEX*16, ALLOCATABLE :: Hamiltonian_2(:,:) !Many-body hamiltonian
  COMPLEX*16, ALLOCATABLE :: Hamiltonian_1_crs(:) !One electron hamiltonian in compressed row storage
  INTEGER*4, ALLOCATABLE :: column_1_crs(:)
  INTEGER*4, ALLOCATABLE :: row_1_crs(:)

  COMPLEX*16, ALLOCATABLE :: Hamiltonian_2_crs(:) !Many-body hamiltonian in compressed row storage
  INTEGER*4, ALLOCATABLE :: column_2_crs(:)
  INTEGER*4, ALLOCATABLE :: row_2_crs(:)  
  REAL*8, ALLOCATABLE :: potential(:, :)
  INTEGER :: state_to_write
  !COMPLEX*16, ALLOCATABLE :: psi_lapack(:, :, :, :, :)
  !REAL*8, ALLOCATABLE :: ev_lapack(:)
  COMPLEX*16, ALLOCATABLE :: Psi_1(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
  COMPLEX*16, ALLOCATABLE :: C_slater(:,:) !Coefficients of consecutive Slater determinants (\varphi_i) for many-body hamiltonian such that many-body wavefunction
                                           !is given by \Psi_2 = \sum_{i=1}^{n} [c_i \varphi_i(r_1, r_2, ..., r_n)]
  REAL*8, ALLOCATABLE :: Particle_density(:,:)
  REAL*8, ALLOCATABLE :: Energies_1(:) !Single electron energies
  REAL*8, ALLOCATABLE :: Energies_2(:) !Many-body energies in the same order as multi-body wavefunctions in \Psi_2
  INTEGER*4, ALLOCATABLE :: Combination_current(:)
  INTEGER*4,  ALLOCATABLE :: Combinations(:,:)
  INTEGER*1, ALLOCATABLE :: N_changed_indeces(:,:) !Number of changed indeces in given element of multi-body hamiltonian.
                                                 !Needed for Slater-Condon rules implementation.
                                                 !0 - no indeces chaned,
                                                 !1 - one index changed,
                                                 !2 - two indeces changed,
                                                 !3 - three or more indeces chnged.
  INTEGER*4, ALLOCATABLE :: Changed_indeces(:,:,:,:) !Specifies which index has been changed to whic for example 1 -> 2 and 3 -> 4.
  LOGICAL, ALLOCATABLE :: Index_found(:)
  INTEGER*4 :: k_electrons, total
  INTEGER :: ix, iy, iorb, is, nn
  INTEGER*4 :: i,j,k,l, n
  REAL*8 :: x, y
  INTEGER*4 :: ham_1_size, nonzero_ham_1
  INTEGER*4 :: ham_2_size, nonzero_ham_2
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
  nonzero_ham_1 = (2*Nx - 1)*2*Ny*N_INTERIOR_ELEMENTS + (2*Ny)*N_RIGHT_FACET_ELEMENTS +&
                & (2*Nx)*N_TOP_FACET_ELEMENTS + (2*Ny)*N_LEFT_FACET_ELEMENTS + N_RIGHT_TOP_CORNER_ELEMENTS
  !Size of two-electron hamiltonian is equal to number of combinations of one-electron wavefunctions
  !That we put into Slater determinants.
  ham_2_size = GAMMA(nstate + 1.0d0)/(GAMMA(k_electrons + 1.0d0) * GAMMA(nstate - k_electrons + 1.0d0))
  slater_norm = SQRT(GAMMA(k_electrons + 1.0d0))
  PRINT*, "Number of combinations: ", ham_2_size
  PRINT*, "Size of one-body hamiltonian: ", ham_1_size
  PRINT*, 'Nonzero elements in hamiltonian 1: ', nonzero_ham_1
  
  
  ALLOCATE(Combination_current(k_electrons))
  ALLOCATE (potential(-Nx:Nx, -Ny:Ny))
  !ALLOCATE (psi_lapack(nstate, -Nx:Nx, -Ny:Ny, norbs / 2, 2))
  !ALLOCATE (ev_lapack(nstate))
  !ALLOCATE (psi_arpack(nstate, -Nx:Nx, -Ny:Ny, norbs / 2, 2))
  ALLOCATE(Hamiltonian_1_crs(nonzero_ham_1))
  ALLOCATE(column_1_crs(nonzero_ham_1))
  ALLOCATE(row_1_crs(ham_1_size + 1))
  ALLOCATE (Particle_density(ham_1_size, nstate))
  ALLOCATE (Psi_1(ham_1_size, nstate))
  ALLOCATE (Energies_1(nstate))
  ALLOCATE (C_slater(ham_2_size, nstate))
  ALLOCATE (Energies_2(nstate))
  !ALLOCATE (Hamiltonian_1(ham_1_size, ham_1_size))
  ALLOCATE (Hamiltonian_2(ham_2_size, ham_2_size))
  ALLOCATE (N_changed_indeces(ham_2_size, ham_2_size))
  ALLOCATE (Combinations(ham_2_size, k_electrons))
  ALLOCATE (Changed_indeces(ham_2_size, ham_2_size, 2, 2))
  ALLOCATE (Index_found(k_electrons))
  PRINT*, "Memory check..."


  !Hamiltonian_1(:,:) = DCMPLX(0.0d0, 0.0d0)
  Hamiltonian_2(:,:) = DCMPLX(0.0d0, 0.0d0)
  N_changed_indeces(:,:) = 0
  Changed_indeces(:,:,:,:) = 0

  PRINT*, "Memory check 2..."

  !#################### ONE-ELECTRON PROBLEM ####################
  !omega = 37.378e-3 * eV2au !Julian: is this meant to be a fixed value? Maybe should be in input .nml file?
  potential = 0.0d0
  DO ix = -Nx, Nx
    DO iy = -Ny, Ny
      x = ix * dx
      y = iy * dx
      potential(ix, iy) = 0.5 * m_eff * omega**2 * (x**2 + y**2)
    END DO
  END DO

  PRINT*, "Creating one electron hamiltonian..."
  CALL CREATE_ONE_ELECTRON_HAMILTONIAN_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Nx, Ny, norbs, potential)
  !CALL HAMILTONIAN_CREATE(Hamiltonian_1(:,:), ham_1_size, Nx, Ny, norbs, potential)

  PRINT*, "Diagonalizing one electron hamiltonian..."
  CALL DIAGONALIZE_ARPACK_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Psi_1, Energies_1, Nx, Ny, norbs, nstate)
  !CALL DIAGONALIZE_ARPACK(Hamiltonian_1, ham_1_size, Psi_1, Energies_1, Nx, Ny, norbs, nstate)
  !###############################################################

  CALL WRITE_STATE_MAP(Psi_1, ham_1_size, nstate, norbs, Nx, Ny, dx, '../RUNS/Psi_1')

  !Julian: This is only data output. I will wrap it in some subroutines to make the main function cleaner.
  WRITE (*, *) "Single electron energies:"
  DO i = 1, nstate
    WRITE (*, *) Energies_1(i) / eV2au * 1e3, sigma_x_expected_value(Psi_1(:,i), ham_1_size),&
      sigma_y_expected_value(Psi_1(:,i), ham_1_size), sigma_z_expected_value(Psi_1(:,i), ham_1_size)
      WRITE(*,*) d_xy_share(Psi_1(:,i), ham_1_size, norbs), d_xz_share(Psi_1(:,i), ham_1_size, norbs), d_yz_share(Psi_1(:,i), ham_1_size, norbs)
      WRITE(*,*)
  END DO
  CALL WRITE_ENERGIES(Energies_1, nstate, '../RUNS/Energies1.dat')


  !#################### TWO-ELECTRON PROBLEM #####################
  !Many problem is solved using Galerkin method, where basis
  !is chosen as combinnation of Slater determinants.


  !Initialize combination_row to generate next sets
  CALL INIT_COMBINATION(Combination_current, k_electrons)
  DO i = 1, ham_2_size
    CALL GET_COMBINATION(Combination_current, nstate, k_electrons)
    Combinations(i, :) = Combination_current(:)
  END DO

  CALL GET_CHANGED_INDECES(Changed_indeces, Combinations, N_changed_indeces, ham_2_size, k_electrons)

  PRINT*, "Constructing multi-body hamiltonian..."
  !CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:,1), Psi_1(:,2), Psi_1(:,3), Psi_1(:,4), ham_1_size, interaction_element, norbs, Nx, Ny, dx)
  !Computing diagonal elements of Hamiltonian_2
  DO i = 1, ham_2_size
    !Single electron energies
    DO n = 1, k_electrons
      Hamiltonian_2(i,i) = Hamiltonian_2(i,i) + Energies_1(Combinations(i,n))
    END DO

    !Interaction elements
    DO k = 1, k_electrons
      DO l = k + 1, k_electrons !Check whether sum could be reduced to sum_{k, l>k}. Then I will get rid of 0.5*
        IF (Combinations(i,k) /= Combinations(i,l)) THEN
          !PRINT*, Combinations(i,k), Combinations(i,l)
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                            & Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                            &ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
          Hamiltonian_2(i,i) = Hamiltonian_2(i,i) + interaction_element !Check whether it should be  * 0.5.

          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Combinations(i,k)), Psi_1(:, Combinations(i,l)),&
                                            & Psi_1(:, Combinations(i,l)), Psi_1(:, Combinations(i,k)),&
                                            & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
          Hamiltonian_2(i,i) = Hamiltonian_2(i,i) - interaction_element !Check whether it should be  * 0.5
        END IF
      END DO
    END DO


  END DO


  DO i = 1, ham_2_size
    PRINT*, "i = ", i
    DO j = i + 1, ham_2_size
      !PRINT*, i, j
      !If statements' order determined by frequency of given elements
      IF (N_changed_indeces(i,j) == 2) THEN
        !Interaction elements
        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)),&
                                          & Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:, Changed_indeces(i,j, 2, 2)),&
                                          & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
        Hamiltonian_2(i,j) = Hamiltonian_2(i,j) + interaction_element

        CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Changed_indeces(i,j,2,1)),&
                                          & Psi_1(:, Changed_indeces(i,j, 2, 2)), Psi_1(:, Changed_indeces(i,j, 1, 2)),&
                                          & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
        Hamiltonian_2(i,j) = Hamiltonian_2(i,j) - interaction_element

      ELSE IF (N_changed_indeces(i,j) == 1) THEN
        !Interaction elements
        DO k = 1, k_electrons
          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)),&
                                            & Psi_1(:, Changed_indeces(i,j, 1, 2)), Psi_1(:, Combinations(i,k)),&
                                            & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
          Hamiltonian_2(i,j) = Hamiltonian_2(i,j) + interaction_element

          CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, Changed_indeces(i,j,1,1)), Psi_1(:, Combinations(i,k)),&
                                            & Psi_1(:, Combinations(i,k)), Psi_1(:, Changed_indeces(i,j, 1, 2)),&
                                            & ham_1_size, interaction_element, norbs, Nx, Ny, dx, eps_r)
          Hamiltonian_2(i,j) = Hamiltonian_2(i,j) - interaction_element

        END DO
      ELSE IF (N_changed_indeces(i,j) /= 3) THEN
        STOP
      END IF
    END DO
  END DO


  ! DO i = 1, ham_2_size
  !   WRITE(*,*) (REAL(Hamiltonian_2(i,j)), j = 1, ham_2_size)
  ! END DO

  !Does not work for small size of ham_2_size
  CALL DIAGONALIZE_ARPACK(Hamiltonian_2, ham_2_size, C_slater, Energies_2, Nx, Ny, norbs, 2)


  CALL CALCULATE_PARTICLE_DENSITY(Particle_density, Psi_1, C_slater, N_changed_indeces,&
  &Changed_indeces, ham_1_size, ham_2_size, Nx, Ny, nstate, nstate, dx, k_electrons)

  !##########################################################
  CALL WRITE_STATE_MAP(DCMPLX(Particle_density, 0.0d0), ham_1_size, nstate, norbs, Nx, Ny, dx, '../RUNS/Density')


  WRITE (*, *) "Multi-body energies"
  DO i = 1, nstate
    !WRITE (*, '(200e20.12)') Energies_2(i) / eV2au * 1e3! ,ev_lapack(i) / eV2au
    WRITE (*, *) Energies_2(i) / eV2au * 1e3! ,ev_lapack(i) / eV2au
  END DO
  CALL WRITE_ENERGIES(Energies_2, 2, '../RUNS/Energies2.dat')


  DEALLOCATE (potential)
  DEALLOCATE(Hamiltonian_1_crs)
  DEALLOCATE(row_1_crs)
  DEALLOCATE(column_1_crs)
  !DEALLOCATE (psi_lapack)
  !DEALLOCATE (ev_lapack)
  DEALLOCATE (Psi_1)
  DEALLOCATE (Energies_1)
  DEALLOCATE (C_slater)
  DEALLOCATE (Energies_2)
  !DEALLOCATE (Hamiltonian_1)
  DEALLOCATE (Hamiltonian_2)
  DEALLOCATE (N_changed_indeces)
  DEALLOCATE (Combinations)
  DEALLOCATE (Combination_current)
END PROGRAM MAIN

! SUBROUTINE WRITE_STATE_MAP(Psi, psi_size, nstates, norbs, Nx, Ny, dx, filename)
!   USE constants
!   IMPLICIT NONE
!   COMPLEX*16, INTENT(IN) :: Psi(psi_size, nstates)
!   INTEGER*4, INTENT(IN) :: Nx, Ny, psi_size, nstates, norbs
!   REAL*8, INTENT(IN) :: dx
!   CHARACTER(LEN=*), INTENT(IN) :: filename
!   CHARACTER(LEN=200) :: filename_state

!   INTEGER*4 :: nn, ix, iy, is, i, iorb

!   DO is = 1, nstates
!     WRITE(filename_state, '(2A, I0, A)') filename, '_n', is, '.dat'
!     OPEN (1, FILE=filename_state)
!     WRITE(1,*) '#x[nm] y[nm] |Psi(orb1, s = +)|**2 |Psi(orb1, s = -)|**2 |Psi(orb2, s = +)|**2 ...'

!     nn = 1
!     DO iy = -Ny, Ny
!     DO ix = -Nx, Nx
!       WRITE(1,'(2E20.8)', ADVANCE='NO') ix*dx/nm2au, iy*dx/nm2au

!       DO iorb = 1, norbs / 2
!         DO i = 1, 2
!           WRITE(1, '(E20.8)', ADVANCE='NO') ABS(Psi(nn,is))**2
!           nn = nn + 1
!         END DO
!       END DO
!       WRITE(1,*)
!     END DO
!     END DO
!     IF (nn /= psi_size + 1) PRINT*, "Error in printing state map, nn /= psi_size + 1"

!     CLOSE(1)
!   END DO


! END SUBROUTINE WRITE_STATE_MAP

! SUBROUTINE WRITE_ENERGIES(Energies, nstates, filename)
!   USE constants
!   IMPLICIT NONE
!   REAL*8, INTENT(IN) :: Energies(nstates)
!   INTEGER*4, INTENT(IN) :: nstates
!   CHARACTER(LEN=*), INTENT(IN) :: filename
!   INTEGER*4 :: i


!   OPEN (1, FILE=filename)
!   WRITE(1,*) '#No. sate [-]   Energy [meV]'
!   DO i = 1, nstates
!     WRITE(1,*) i, Energies(i) / eV2au * 1e3
!   END DO
!   CLOSE(1)

! END SUBROUTINE WRITE_ENERGIES