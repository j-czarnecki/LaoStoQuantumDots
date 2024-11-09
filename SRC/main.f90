#include "macros_def.f90"
PROGRAM MAIN
  USE indata
  USE hamiltonian
  USE diagonalize
  USE combinatory
  USE many_body
  USE writers
  USE constants
  USE utility
  USE logger

  USE omp_lib

  IMPLICIT NONE

  COMPLEX*16, ALLOCATABLE :: Hamiltonian_1_crs(:) !One electron hamiltonian in compressed row storage
  INTEGER*4, ALLOCATABLE :: column_1_crs(:)
  INTEGER*4, ALLOCATABLE :: row_1_crs(:)

  COMPLEX*16, ALLOCATABLE :: Hamiltonian_2_crs(:) !Many-body hamiltonian in compressed row storage
  INTEGER*4, ALLOCATABLE :: column_2_crs(:)
  INTEGER*4, ALLOCATABLE :: row_2_crs(:)
  REAL*8, ALLOCATABLE :: potential(:, :)
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
  INTEGER*4, ALLOCATABLE :: N_ham_2_elems_in_prev_rows(:) !Contains number of non-zero elements in previous rows of multi-body hamiltonian.
                                                          !Needed for multithreading, so that at each row we can start filling up new elements independently.
  COMPLEX*16, ALLOCATABLE :: Nxm_elems(:,:) !matrix elements of x-position operator <n|X|m> wheren n and m are multi-body eigenstates
  COMPLEX*16, ALLOCATABLE :: C_time(:) !Those are the time-dependent coefficients of the many-body wavefunction, such that
                                         !\Psi(t) = \sum_{m} [c_m(t) * exp(-i*E_n*t / \hbar) |m> ].
                                         ! Where m is a multi-body eigenstate.
                                         ! first index denotes time, second is the number of state. Time is stored such that t(1) > t(2) > ... > t(n)
  COMPLEX*16, ALLOCATABLE :: C_max_time(:)
  COMPLEX*16, ALLOCATABLE :: A_crank_nicolson(:,:)
  COMPLEX*16, ALLOCATABLE :: B_crank_nicolson(:,:)
  INTEGER*4, ALLOCATABLE :: IPIV(:)
  INTEGER*4 :: INFO
  INTEGER*4 :: it, iomega
  REAL*8 :: omega_ac
  INTEGER :: ix, iy
  INTEGER*4 :: i,n, m
  REAL*8 :: x, y
  INTEGER*4 :: ham_1_size, nonzero_ham_1
  INTEGER*4 :: ham_2_size, nonzero_ham_2
  REAL*8 :: slater_norm

  !OMP specific
  INTEGER*4 :: max_num_threads, used_threads
  max_num_threads = 0
  used_threads = 0
  CALL INIT_LOGGER()

  !read a file
  CALL INDATA_GET("./OutputData/quantum_dot.nml")
  WRITE(log_string, *) 'Calculation parameters: ',&
    & 'Nx = ', Nx,&
    & 'Ny = ', Ny,&
    & 'dx = ', dx / nm2au,&
    & 'norbs = ', norbs,&
    & 'nstate_1 = ', nstate_1,&
    & 'nstate_2 = ', nstate_2,&
    & 'k_electrons = ', k_electrons
  LOG_INFO(log_string)
  WRITE(log_string, *) 'Physical parameters: ',&
    & 'th = ', th / eV2au,&
    & 'tl = ', tl / eV2au,&
    & 'td = ', td / eV2au,&
    & 'dso = ', dso / eV2au,&
    & 'drso = ', drso / eV2au,&
    & 'dE = ', dE / eV2au,&
    & 'g = ', g
  LOG_INFO(log_string)
  WRITE(log_string, *) 'External parameters: ',&
    & 'omega = ', omega / eV2au,&
    & 'Bx = ', Bx / T2au,&
    & 'By = ', By / T2au,&
    & 'Bz = ', Bz / T2au,&
    & 'domega_ac = ', domega_ac / eV2au,&
    & 'omega_ac_max = ', omega_ac_max / eV2au,&
    & 'F = ', f_ac/F2au
  LOG_INFO(log_string)

  max_num_threads = omp_get_max_threads()
  CALL omp_set_num_threads(max_num_threads)
  WRITE(log_string, *) 'Number of threads: ', max_num_threads
  LOG_INFO(log_string)

  ! !This is to test if parallelization works
  ! !$omp parallel
  ! WRITE(66,*) "Thread", omp_get_thread_num()
  ! !$omp critical
  ! used_threads = omp_get_num_threads()
  ! !$omp end critical
  ! !$omp end parallel
  ! WRITE(66,*) "Used threads", used_threads

  ham_1_size = (2*Nx + 1)*(2*Ny + 1)*norbs
  nonzero_ham_1 = (2*Nx - 1)*2*Ny*N_INTERIOR_ELEMENTS + (2*Ny)*N_RIGHT_FACET_ELEMENTS +&
                & (2*Nx)*N_TOP_FACET_ELEMENTS + (2*Ny)*N_LEFT_FACET_ELEMENTS + N_RIGHT_TOP_CORNER_ELEMENTS
  !Size of two-electron hamiltonian is equal to number of combinations of one-electron wavefunctions
  !That we put into Slater determinants.
  ham_2_size = INT(GAMMA(nstate_1 + 1.0d0)/(GAMMA(k_electrons + 1.0d0) * GAMMA(nstate_1 - k_electrons + 1.0d0)))
  !Number of diagonal elements is equal to number of k_electrons-element combinations from nstates
  !Number of elements with one changed index: k*(n-k)
  !Number of elements with two changed index: (k    2) * (n-k   2) <- Newton symbols
  !First lets take off-diagonal terms - swapped one or two indeces. I only store upper triangle, therefore will only have half of them.
  !After all add ham_2_size which is the number of diagonal elements
  nonzero_ham_2 = ham_2_size*(k_electrons*(nstate_1 - k_electrons)&
  & + INT(GAMMA(k_electrons + 1.0d0)*GAMMA(nstate_1 - k_electrons + 1.0d0) / (4*GAMMA(k_electrons - 1.0d0)*GAMMA(nstate_1 - k_electrons - 1.0d0))))/2&
  & + ham_2_size

  slater_norm = SQRT(GAMMA(k_electrons + 1.0d0))
  ! PRINT*, "Number of combinations: ", ham_2_size
  ! PRINT*, "Size of one-body hamiltonian: ", ham_1_size
  ! PRINT*, 'Nonzero elements in hamiltonian 1: ', nonzero_ham_1
  ! PRINT*, 'Nonzero elements in hamiltonian 2: ', nonzero_ham_2

  ALLOCATE (Combination_current(k_electrons))
  ALLOCATE (potential(-Nx:Nx, -Ny:Ny))
  ALLOCATE (Hamiltonian_1_crs(nonzero_ham_1))
  ALLOCATE (column_1_crs(nonzero_ham_1))
  ALLOCATE (row_1_crs(ham_1_size + 1))
  ALLOCATE (Particle_density(ham_1_size, nstate_2))
  ALLOCATE (Psi_1(ham_1_size, nstate_1))
  ALLOCATE (Energies_1(nstate_1))
  ALLOCATE (C_slater(ham_2_size, nstate_2))
  ALLOCATE (Energies_2(nstate_2))
  ALLOCATE (Hamiltonian_2_crs(nonzero_ham_2))
  ALLOCATE (column_2_crs(nonzero_ham_2))
  ALLOCATE (row_2_crs(ham_2_size + 1))
  ALLOCATE (N_changed_indeces(ham_2_size, ham_2_size))
  ALLOCATE (Combinations(ham_2_size, k_electrons))
  ALLOCATE (Changed_indeces(ham_2_size, ham_2_size, 2, 2))
  ALLOCATE (N_ham_2_elems_in_prev_rows(ham_2_size))
  ALLOCATE (Nxm_elems(nstate_2, nstate_2))
  ALLOCATE (C_time(nstate_2))
  ALLOCATE (C_max_time(nstate_2))
  ALLOCATE (A_crank_nicolson(nstate_2, nstate_2))
  ALLOCATE (B_crank_nicolson(nstate_2, 1))
  ALLOCATE (IPIV(nstate_2))
  !PRINT*, "Memory check..."


  N_changed_indeces(:,:) = 0
  Changed_indeces(:,:,:,:) = 0

  !PRINT*, "Memory check 2..."
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

  CALL CREATE_ONE_ELECTRON_HAMILTONIAN_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Nx, Ny, norbs, potential)

  CALL DIAGONALIZE_ARPACK_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Psi_1, Energies_1, nstate_1)
  !###############################################################

  !Writing single-electron problem data to a file
  CALL WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi_1, ham_1_size, nstate_1, norbs, Nx, Ny, dx, './OutputData/Psi_1')
  CALL WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1, ham_1_size, nstate_1, norbs, './OutputData/Expectations_1.dat')
  CALL WRITE_ENERGIES(Energies_1, nstate_1, './OutputData/Energies1.dat')


  !#################### TWO-ELECTRON PROBLEM #####################
  !Many problem is solved using Galerkin method, where basis
  !is chosen as combinnation of Slater determinants.

  !Initialize combination_row to generate next sets
  CALL INIT_COMBINATION(Combination_current, k_electrons)
  DO i = 1, ham_2_size
    CALL GET_COMBINATION(Combination_current, nstate_1, k_electrons)
    Combinations(i, :) = Combination_current(:)
  END DO
  CALL GET_CHANGED_INDECES(Changed_indeces, Combinations, N_changed_indeces, ham_2_size, k_electrons)
  CALL INIT_PREV_ELEMS(N_ham_2_elems_in_prev_rows, N_changed_indeces, ham_2_size, nonzero_ham_2)

  CALL CREATE_MANY_BODY_HAMILTONIAN_CRS(Hamiltonian_2_crs, column_2_crs, row_2_crs, N_changed_indeces, Changed_indeces, Combinations,&
  & N_ham_2_elems_in_prev_rows, Psi_1, Energies_1, ham_1_size, nstate_1, nonzero_ham_2, ham_2_size, k_electrons, norbs, Nx, Ny, dx, eps_r)

  Energies_2(:) = 0.0d0
  CALL DIAGONALIZE_ARPACK_CRS(Hamiltonian_2_crs, column_2_crs, row_2_crs, nonzero_ham_2, ham_2_size, C_slater, Energies_2, nstate_2)

  ! CALL CALCULATE_PARTICLE_DENSITY(Particle_density, Psi_1, C_slater, N_changed_indeces,&
  ! &Changed_indeces, ham_1_size, ham_2_size, nstate_1, nstate_2, k_electrons)

  !##########################################################
  !CALL WRITE_STATE_MAP(DCMPLX(Particle_density, 0.0d0), ham_1_size, nstate_2, norbs, Nx, Ny, dx, './OutputData/Density')

  CALL WRITE_SLATER_COEFFICIENTS(C_slater, ham_2_size, nstate_2, './OutputData/C_slater.dat')
  CALL WRITE_MULTI_ELECTRON_EXPECTATIONS(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, Nx, dx, './OutputData/Expectations_2.dat')
  CALL WRITE_ENERGIES(Energies_2, nstate_2, './OutputData/Energies2.dat')


  !#################### TIME_DEPENDENT CALCULATION #################
  !Time dependent Schrodinger equation is solved using Crank-Nicholson scheme

  !First calculate all <n|X|m> matrix elements
  !Only upper triangle is needed, since <n|X|m> = CONJG(<m|X|n>)
  OPEN(10, FILE = './OutputData/Nxm.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  !$omp parallel
  !$omp do
  DO n = 1, nstate_2
    DO m = 1, nstate_2
      Nxm_elems(n,m) = many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,&
        & ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, m, Nx, dx, norbs)
      WRITE(10,*) n, m, REAL(Nxm_elems(n,m)), AIMAG(Nxm_elems(n,m))
    END DO
  END DO
  !$omp end do
  !$omp end parallel
  CLOSE(10)

  WRITE(log_string,*) "Running Crank-Nicolson scheme"
  LOG_INFO(log_string)
  !For such energy we have to trigger transition, this is only for a check.
  !omega_ac = Energies_2(2) - Energies_2(1)

  OPEN(10, FILE = './OutputData/C_max_time.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  WRITE(10,*) '#omega_ac [meV] C_max_1 C_max_2 ...'
  !$omp parallel private(omega_ac, C_time, C_max_time, A_crank_nicolson, B_crank_nicolson, ipiv, INFO)
  !$omp do
  DO iomega = 1, N_omega_ac_steps
    omega_ac = iomega*domega_ac
    C_max_time = DCMPLX(0.0d0, 0.0d0)
    C_time = DCMPLX(0.0d0, 0.0d0)
    C_time(1) = DCMPLX(1.0d0, 0.0d0) !Initial condition, assuming full occupation of lowest energy state
    !OPEN(10, FILE = './OutputData/Time_dependent.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
    DO it = 0, N_t_steps
      A_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
      B_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
      ipiv = 0
      DO n = 1, nstate_2
        A_crank_nicolson(n,n) = 1.0d0
        B_crank_nicolson(n,1) = C_time(n)

        DO m = 1, nstate_2
          A_crank_nicolson(n,m) = A_crank_nicolson(n,m) - dt/(2*imag) * energy_phase_offset(Energies_2(n), Energies_2(m), (it + 1)*dt) * v_ac(f_ac, omega_ac, (it + 1)*dt) * Nxm_elems(n,m)
          B_crank_nicolson(n,1) = B_crank_nicolson(n,1) + dt/(2*imag) * C_time(m) * energy_phase_offset(Energies_2(n), Energies_2(m), it*dt) * v_ac(f_ac, omega_ac, it*dt) * Nxm_elems(n,m)
        END DO
      END DO

      CALL ZGESV(nstate_2, 1, A_crank_nicolson, nstate_2, IPIV, B_crank_nicolson, nstate_2, INFO)
      B_crank_nicolson = B_crank_nicolson !/ SUM(ABS(B_crank_nicolson(:,1))**2)
      !IF (INFO /= 0) !PRINT*, "ERROR in ZGESV, INFO = ", INFO
      !WRITE(10,*) it*dt / ns2au, (ABS(B_crank_nicolson(n,1))**2, n = 1, nstate_2)
      C_time(:) = B_crank_nicolson(:,1)


      !Check for maximal value of C_time
      DO n = 1, nstate_2
        IF (ABS(C_time(n))**2 > ABS(C_max_time(n))**2) C_max_time(n) = C_time(n)
      END DO
    END DO

    WRITE(10,*) omega_ac / eV2au * 1e3, (ABS(C_max_time(n))**2, n = 1, nstate_2)
    FLUSH(10)
  END DO
  !$omp end do
  !$omp end parallel
  CLOSE(10)
  !#################################################################

  CALL CLOSE_LOGGER()

  DEALLOCATE (potential)
  DEALLOCATE (Hamiltonian_1_crs)
  DEALLOCATE (row_1_crs)
  DEALLOCATE (column_1_crs)
  DEALLOCATE (Particle_density)
  DEALLOCATE (Psi_1)
  DEALLOCATE (Energies_1)
  DEALLOCATE (C_slater)
  DEALLOCATE (Energies_2)
  DEALLOCATE (Hamiltonian_2_crs)
  DEALLOCATE (column_2_crs)
  DEALLOCATE (row_2_crs)
  DEALLOCATE (N_changed_indeces)
  DEALLOCATE (Combinations)
  DEALLOCATE (Changed_indeces)
  DEALLOCATE (N_ham_2_elems_in_prev_rows)
  DEALLOCATE (Combination_current)
  DEALLOCATE (Nxm_elems)
  DEALLOCATE (C_time)
  DEALLOCATE (C_max_time)
  DEALLOCATE (A_crank_nicolson)
  DEALLOCATE (B_crank_nicolson)
  DEALLOCATE (IPIV)
END PROGRAM MAIN
