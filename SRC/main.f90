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
  USE swap

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
  COMPLEX*16, ALLOCATABLE :: Psi_cross(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
  COMPLEX*16, ALLOCATABLE :: Psi_1_noSO(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction (no SO for bonding and anti-bonding states)
  COMPLEX*16, ALLOCATABLE :: Psi_LR(:, :) !Wavefunction localized in left / right quantum dot
  COMPLEX*16, ALLOCATABLE :: Psi_init(:,:) !Initial wavefunction
  COMPLEX*16, ALLOCATABLE :: Psi_density_init(:) !Initial density
  COMPLEX*16, ALLOCATABLE :: C_slater(:,:) !Coefficients of consecutive Slater determinants (\varphi_i) for many-body hamiltonian such that many-body wavefunction
                                           !is given by \Psi_2 = \sum_{i=1}^{n} [c_i \varphi_i(r_1, r_2, ..., r_n)]
  REAL*8, ALLOCATABLE :: Particle_density(:,:)
  REAL*8, ALLOCATABLE :: Particle_density_up(:,:)
  REAL*8, ALLOCATABLE :: Particle_density_down(:,:)
  REAL*8, ALLOCATABLE :: Energies_1(:) !Single electron energies
  REAL*8, ALLOCATABLE :: Energies_1_noSO(:) !Single electron energies (no SO for bonding and anti-bonding states)
  COMPLEX*16, ALLOCATABLE :: Nxm_single_elems(:,:) ! <n|x|m> elements for single electron wavefunctions
  COMPLEX*16, ALLOCATABLE :: C_single_max_time(:)
  COMPLEX*16, ALLOCATABLE :: C_single_time(:)
  COMPLEX*16, ALLOCATABLE :: Cm(:)
  COMPLEX*16, ALLOCATABLE :: A_single_crank_nicolson(:,:)
  COMPLEX*16, ALLOCATABLE :: B_single_crank_nicolson(:,:)
  INTEGER*4, ALLOCATABLE :: IPIV_single(:)
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
  COMPLEX*16, ALLOCATABLE :: Spin_t(:,:)
  COMPLEX*16, ALLOCATABLE :: Spin_t_single(:,:)
  COMPLEX*16, ALLOCATABLE :: many_body_sigma_z_expected_value_ind(:)
  INTEGER*4 :: INFO
  INTEGER*4 :: it, iomega
  REAL*8 :: omega_ac
  INTEGER :: ix, iy
  INTEGER*4 :: i,n, m
  REAL*8 :: x, y, Rx, Ry, Rb,mu
  INTEGER*4 :: ham_1_size, nonzero_ham_1
  INTEGER*4 :: ham_2_size, nonzero_ham_2
  INTEGER*4 :: t_max_int
  INTEGER*4 :: nmax
  REAL*8 :: slater_norm
  CHARACTER(LEN=100) :: filename
  LOGICAL :: isSO

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
    & 'F = ', f_ac/F2au,&
    & 'V_b = ', Vb/eV2au,&
    & 'V_0 = ', V0/eV2au
  LOG_INFO(log_string)

  max_num_threads = 1 !omp_get_max_threads()
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

  ham_1_size = get_ham_1_size(Nx, Ny, norbs)
  nonzero_ham_1 = get_nonzero_ham_1(Nx, Ny)
  ham_2_size = get_ham_2_size(nstate_1, k_electrons)
  nonzero_ham_2 = get_nonzero_ham_2(ham_2_size, nstate_1, k_electrons)
  WRITE(log_string,*) 'ham_1_size = ', ham_1_size, 'nonzero_ham_1 = ', nonzero_ham_1, 'ham_2_size = ', ham_2_size, 'nonzero_ham_2 = ', nonzero_ham_2
  LOG_INFO(log_string)

  slater_norm = SQRT(GAMMA(k_electrons + 1.0d0))
  nmax = 6
  t_max_int = INT(t_max/dt)
  ! PRINT*, t_max_int
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
  ALLOCATE (Particle_density_up(ham_1_size, nstate_2))
  ALLOCATE (Particle_density_down(ham_1_size, nstate_2))
  ALLOCATE (Psi_1(ham_1_size, nstate_1))
  ALLOCATE (Psi_cross(ham_1_size, 1))
  ALLOCATE (Psi_1_noSO(ham_1_size, nstate_1))
  ALLOCATE (Psi_LR(ham_1_size, 2))
  ALLOCATE (Psi_init(ham_1_size,2))
  ALLOCATE (Psi_density_init(ham_1_size))
  ALLOCATE (Energies_1(nstate_1))
  ALLOCATE (Energies_1_noSO(nstate_1))
  ALLOCATE (Nxm_single_elems(nstate_1, nstate_1))
  ALLOCATE (C_single_time(nstate_1))
  ALLOCATE (C_single_max_time(nstate_1))
  ALLOCATE (A_single_crank_nicolson(nstate_1, nstate_1))
  ALLOCATE (B_single_crank_nicolson(nstate_1, 1))
  ALLOCATE (IPIV_SINGLE(nstate_1))
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
  ALLOCATE (Cm(nstate_1))
  ALLOCATE (C_max_time(nstate_2))
  ALLOCATE (A_crank_nicolson(nstate_2, nstate_2))
  ALLOCATE (B_crank_nicolson(nstate_2, 1))
  ALLOCATE (IPIV(nstate_2))
  ALLOCATE (Spin_t(t_max_int,nmax))
  ALLOCATE (Spin_t_single(t_max_int,nmax))
  ALLOCATE(many_body_sigma_z_expected_value_ind(k_electrons))
  !PRINT*, "Memory check..."


  N_changed_indeces(:,:) = 0
  Changed_indeces(:,:,:,:) = 0

  !PRINT*, "Memory check 2..."
  !#################### ONE-ELECTRON PROBLEM ####################
  Rx = dx * Nx/1.15
  Ry = dx * Ny/1.4
  Rb = Ry / 2.0
  mu = 10.0
  potential = 0.0d0
  DO ix = -Nx, Nx
    DO iy = -Ny, Ny
      x = ix * dx
      y = iy * dx
      potential(ix, iy) = -V0 /((1+(x**2/Rx**2)**mu)*(1+(y**2/Ry**2)**mu)) + Vb /((1+(x**2/Rb**2)**mu)*(1+(y**2/Ry**2)**mu))
      !potential(ix, iy) = 0.5 * m_eff * omega**2 * (x**2 + y**2)
    END DO
  END DO
  CALL WRITE_POTENTIAL(potential, Nx, Ny, './OutputData/potential.dat')

  LOG_INFO(log_string)
  if(.not.(nstate_2 < (4 * (nstate_2 + 1)).and.nstate_2 < ham_2_size.and.(4 * (nstate_2 + 1)) < ham_2_size)) THEN
    WRITE(log_string,*) 'In ARPES NEV < NCV ≤ N, please lower nstate_2. STOP.'
    LOG_INFO(log_string)

    print *, "EXIT: see log file"
    STOP
  endif
  !#################### NO SPIN ORBIT ####################
  isSO = .FALSE.

  CALL CREATE_ONE_ELECTRON_HAMILTONIAN_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Nx, Ny, norbs, potential, isSO)
  CALL DIAGONALIZE_ARPACK_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Psi_1_noSO, Energies_1_noSO, nstate_1)
  CALL GET_LEFT_RIGHT_WAVEFUNCTION_PURE(Psi_1_noSO, nstate_1, ham_1_size, Psi_LR, nx, ny, norbs)

  CALL WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi_1_noSO, ham_1_size, nstate_1, norbs, Nx, Ny, dx, './OutputData/Psi_1_noSO')
  CALL WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi_LR, ham_1_size, 2, norbs, Nx, Ny, dx, './OutputData/Psi_LR')
  CALL WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_LR, ham_1_size, 2, norbs, './OutputData/Expectations_LR.dat')
  CALL WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1_noSO, ham_1_size, nstate_1, norbs, './OutputData/Expectations_noSO_1.dat')
  CALL WRITE_ENERGIES(Energies_1_noSO, nstate_1, './OutputData/Energies1_noSO.dat')
  !#################### WITH SPIN ORBIT ####################
  isSO = .TRUE.
  CALL CREATE_ONE_ELECTRON_HAMILTONIAN_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Nx, Ny, norbs, potential, isSO)

  CALL DIAGONALIZE_ARPACK_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Psi_1, Energies_1, nstate_1)

  ! !###############################################################
  ! !#################### TIME_DEPENDENT CALCULATION FOR SINGLE ELECTRON PROBLEM ####################
  ! !Time dependent Schrodinger equation is solved using Crank-Nicholson scheme
  ! WRITE(log_string,*) "Calculating <n|x|m> elements for single electron"
  ! LOG_INFO(log_string)
  ! !First calculate all <n|X|m> matrix elements
  ! !Only upper triangle is needed, since <n|X|m> = CONJG(<m|X|n>)
  ! OPEN(11, FILE = './OutputData/Nxm_single_electrons.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  ! !$omp parallel
  ! !$omp do
  ! DO n = 1, nstate_1
  !   DO m = 1, nstate_1
  !     Nxm_single_elems(n,m) = single_electron_x_expected_value(Psi_1(:, n), Psi_1(:, m), norbs, Nx, dx, ham_1_size)
  !     WRITE(11,*) n, m, REAL(Nxm_single_elems(n,m)), AIMAG(Nxm_single_elems(n,m))
  !   END DO
  ! END DO
  ! !$omp end do
  ! !$omp end parallel
  ! CLOSE(11)


  ! WRITE(log_string,*) "Running Crank-Nicolson scheme for single electron"
  ! LOG_INFO(log_string)
  ! !For such energy we have to trigger transition, this is only for a check.
  ! !omega_ac = Energies_2(2) - Energies_2(1)

  ! OPEN(11, FILE = './OutputData/C_single_max_time.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  ! WRITE(11,*) '#omega_ac [meV] C_max_1 C_max_2 ...'
  ! !PRINT*, "omega steps ", N_omega_ac_steps
  ! !$omp parallel private(omega_ac, C_single_time, C_single_max_time, A_single_crank_nicolson, B_single_crank_nicolson, IPIV_SINGLE, INFO)
  ! !$omp do
  ! DO iomega = 1, N_omega_ac_steps
  !   !omega_ac = iomega*domega_ac
  !   omega_ac = Energies_1(iomega+1) - Energies_1(1)
  !   C_single_max_time = DCMPLX(0.0d0, 0.0d0)
  !   C_single_time = DCMPLX(0.0d0, 0.0d0)
  !   C_single_time(1) = DCMPLX(1.0d0, 0.0d0) !Initial condition, assuming full occupation of lowest energy state
  !   WRITE(filename, '("./OutputData/Time_dependent", I0, ".dat")') iomega
  !   !OPEN(10, FILE = './OutputData/Time_dependent.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  !   OPEN(10, FILE=filename, ACTION='WRITE', FORM='FORMATTED')
  !   DO it = 0, N_t_steps
  !     A_single_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
  !     B_single_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
  !     IPIV_SINGLE = 0
  !     DO n = 1, nstate_1
  !       A_single_crank_nicolson(n,n) = 1.0d0
  !       B_single_crank_nicolson(n,1) = C_single_time(n)

  !       DO m = 1, nstate_1
  !         A_single_crank_nicolson(n,m) = A_single_crank_nicolson(n,m) - dt/(2*imag) * energy_phase_offset(Energies_1(n), Energies_1(m), (it + 1)*dt) * v_ac(f_ac, omega_ac, (it + 1)*dt) * Nxm_single_elems(n,m)
  !         B_single_crank_nicolson(n,1) = B_single_crank_nicolson(n,1) + dt/(2*imag) * C_single_time(m) * energy_phase_offset(Energies_1(n), Energies_1(m), it*dt) * v_ac(f_ac, omega_ac, it*dt) * Nxm_single_elems(n,m)
  !       END DO
  !     END DO

  !     CALL ZGESV(nstate_1, 1, A_single_crank_nicolson, nstate_1, IPIV_SINGLE, B_single_crank_nicolson, nstate_1, INFO)
  !     !B_single_crank_nicolson = B_single_crank_nicolson !/ SUM(ABS(B_crank_nicolson(:,1))**2)
  !     !IF (INFO /= 0) !PRINT*, "ERROR in ZGESV, INFO = ", INFO
  !     WRITE(10,*) it*dt / ns2au, (ABS(B_single_crank_nicolson(n,1))**2, n = 1, nstate_2)
  !     C_single_time(:) = B_single_crank_nicolson(:,1)


  !     !Check for maximal value of C_time
  !     DO n = 1, nstate_1
  !       IF (ABS(C_single_time(n))**2 > ABS(C_single_max_time(n))**2) C_single_max_time(n) = C_single_time(n)
  !     END DO
  !   END DO
  !   CLOSE(10)
  !   WRITE(11,*) omega_ac / eV2au * 1e3, (ABS(C_single_max_time(n))**2, n = 1, nstate_1)
  !   FLUSH(11)

  ! END DO
  ! !$omp end do
  ! !$omp end parallel
  ! CLOSE(11)
  ! !#################################################################

  !Writing single-electron problem data to a file
  CALL WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi_1, ham_1_size, nstate_1, norbs, Nx, Ny, dx, './OutputData/Psi_1')
  CALL WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1, ham_1_size, nstate_1, norbs, './OutputData/Expectations_1.dat')
  CALL WRITE_ENERGIES(Energies_1, nstate_1, './OutputData/Energies1.dat')


  IF (nstate_2 .eq. 0) THEN
    WRITE(log_string,*) "Program terminated after single electron calculation due to nstate_2 == 0"
    LOG_INFO(log_string)
    STOP
  END IF

  !#################### TWO-ELECTRON PROBLEM #####################
  !Many problem is solved using Galerkin method, where basis
  !is chosen as combinnation of Slater determinants.

  !Initialize combination_row to generate next sets
  CALL INIT_COMBINATION(Combination_current, k_electrons)
  DO i = 1, ham_2_size
    CALL GET_COMBINATION(Combination_current, nstate_1, k_electrons)
    Combinations(i, :) = Combination_current(:)
    !PRINT*, Combinations(i, :)
  END DO
  CALL GET_CHANGED_INDECES(Changed_indeces, Combinations, N_changed_indeces, ham_2_size, k_electrons)
  CALL INIT_PREV_ELEMS(N_ham_2_elems_in_prev_rows, N_changed_indeces, ham_2_size, nonzero_ham_2)

  CALL CREATE_MANY_BODY_HAMILTONIAN_CRS(Hamiltonian_2_crs, column_2_crs, row_2_crs, N_changed_indeces, Changed_indeces, Combinations,&
  & N_ham_2_elems_in_prev_rows, Psi_1, Energies_1, ham_1_size, nstate_1, nonzero_ham_2, ham_2_size, k_electrons, norbs, Nx, Ny, dx, eps_r)

  Energies_2(:) = 0.0d0
  CALL DIAGONALIZE_ARPACK_CRS(Hamiltonian_2_crs, column_2_crs, row_2_crs, nonzero_ham_2, ham_2_size, C_slater, Energies_2, nstate_2)

  CALL WRITE_SLATER_COEFFICIENTS(C_slater, ham_2_size, nstate_2, './OutputData/C_slater.dat')
  CALL WRITE_MULTI_ELECTRON_EXPECTATIONS(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, dx, './OutputData/Expectations_2.dat')
  CALL WRITE_ENERGIES(Energies_2, nstate_2, './OutputData/Energies2.dat')

  CALL GET_INIT_COEFFICIENTS(Psi_LR, Psi_1, C_slater, ham_1_size, ham_2_size,  nstate_1, nstate_1, nstate_2, Combinations,k_electrons, Cm)
  
  CALL TIME_EVOLUTION_SPIN_EXPECTATION(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, nmax, k_electrons, Cm, nstate_1, nstate_2, t_max_int, dt, Energies_2, Spin_t,nx, ny, norbs)
  CALL WRITE_TIME_EVOLUTION(Spin_t, t_max_int, nmax, './OutputData/Spin_time_evolution.dat')
  !CALL GET_INIT_COEFFICIENTS(Psi_LR, Psi_1, C_slater, ham_1_size, ham_2_size,  nstate_1_noSO, nstate_2, Combinations,k_electrons, Cm)
  
  !CALL TIME_EVOLUTION_SPIN_EXPECTATION(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, nmax, k_electrons, Cm, nstate_1, nstate_2, t_max_int, dt, Energies_2, Spin_t,nx, ny, norbs)
  !CALL WRITE_TIME_EVOLUTION(Spin_t, t_max_int, nmax, './OutputData/Spin_time_evolution.dat')

  !CALL TIME_EVOLUTION_SINGLE_SPIN_EXPECTATION(Psi_1, ham_1_size, nmax, Cm, nstate_1, nstate_2, t_max_int, dt, Energies_1, Spin_t_single, Nx, Ny, norbs)
  !CALL WRITE_TIME_EVOLUTION(Spin_t_single, t_max_int, nmax, './OutputData/Spin_time_evolution_single.dat')
  ! !#################### TIME_DEPENDENT CALCULATION FOR MANY BODY PROBLEM ####################
  ! !Time dependent Schrodinger equation is solved using Crank-Nicholson scheme

  ! !First calculate all <n|X|m> matrix elements
  ! !Only upper triangle is needed, since <n|X|m> = CONJG(<m|X|n>)
  ! OPEN(10, FILE = './OutputData/Nxm.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  ! !$omp parallel
  ! !$omp do
  ! DO n = 1, nstate_2
  !   DO m = 1, nstate_2
  !     Nxm_elems(n,m) = many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,&
  !       & ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, m, Nx, dx, norbs)
  !     WRITE(10,*) n, m, REAL(Nxm_elems(n,m)), AIMAG(Nxm_elems(n,m))
  !   END DO
  ! END DO
  ! !$omp end do
  ! !$omp end parallel
  ! CLOSE(10)
  ! !OPEN(10, FILE = './OutputData/Spin_time_evolution.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  ! !Spin_t(:) = Cm(:)
  ! !DO it = 0, N_t_steps
  ! !  DO n = 1, nstate_2
  ! !    DO m = 1, nstate_2
  ! !      Spin_t(n) = Spin_t(n) + Cm(m) * many_body_sigma_z_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, m) * EXP(-imag * Energies_2(m) * it*dt)
  ! !    END DO
  ! !  END DO
  ! !  WRITE(10,*) it*dt / ns2au, (REAL(Spin_t(n)), n = 1, nstate_2)
  ! !END DO
  ! !CLOSE(10)
  ! WRITE(log_string,*) "Running Crank-Nicolson scheme"
  ! LOG_INFO(log_string)
  ! !For such energy we have to trigger transition, this is only for a check.
  ! !omega_ac = Energies_2(2) - Energies_2(1)

  ! OPEN(10, FILE = './OutputData/C_max_time.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  ! WRITE(10,*) '#omega_ac [meV] C_max_1 C_max_2 ...'
  ! !$omp parallel private(omega_ac, C_time, C_max_time, A_crank_nicolson, B_crank_nicolson, ipiv, INFO)
  ! !$omp do
  ! DO iomega = 1, N_omega_ac_steps
  !   omega_ac = iomega*domega_ac
  !   C_max_time = DCMPLX(0.0d0, 0.0d0)
  !   C_time = DCMPLX(0.0d0, 0.0d0)
  !   C_time(1) = DCMPLX(1.0d0, 0.0d0) !Initial condition, assuming full occupation of lowest energy state
  !   !OPEN(10, FILE = './OutputData/Time_dependent.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  !   DO it = 0, N_t_steps
  !     A_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
  !     B_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
  !     ipiv = 0
  !     DO n = 1, nstate_2
  !       A_crank_nicolson(n,n) = 1.0d0
  !       B_crank_nicolson(n,1) = C_time(n)

  !       DO m = 1, nstate_2
  !         A_crank_nicolson(n,m) = A_crank_nicolson(n,m) - dt/(2*imag) * energy_phase_offset(Energies_2(n), Energies_2(m), (it + 1)*dt) * v_ac(f_ac, omega_ac, (it + 1)*dt) * Nxm_elems(n,m)
  !         B_crank_nicolson(n,1) = B_crank_nicolson(n,1) + dt/(2*imag) * C_time(m) * energy_phase_offset(Energies_2(n), Energies_2(m), it*dt) * v_ac(f_ac, omega_ac, it*dt) * Nxm_elems(n,m)
  !       END DO
  !     END DO

  !     CALL ZGESV(nstate_2, 1, A_crank_nicolson, nstate_2, IPIV, B_crank_nicolson, nstate_2, INFO)
  !     B_crank_nicolson = B_crank_nicolson !/ SUM(ABS(B_crank_nicolson(:,1))**2)
  !     !IF (INFO /= 0) !PRINT*, "ERROR in ZGESV, INFO = ", INFO
  !     !WRITE(10,*) it*dt / ns2au, (ABS(B_crank_nicolson(n,1))**2, n = 1, nstate_2)
  !     C_time(:) = B_crank_nicolson(:,1)


  !     !Check for maximal value of C_time
  !     DO n = 1, nstate_2
  !       IF (ABS(C_time(n))**2 > ABS(C_max_time(n))**2) C_max_time(n) = C_time(n)
  !     END DO
  !   END DO

  !   WRITE(10,*) omega_ac / eV2au * 1e3, (ABS(C_max_time(n))**2, n = 1, nstate_2)
  !   FLUSH(10)
  ! END DO
  ! !$omp end do
  ! !$omp end parallel
  ! CLOSE(10)
  !#################################################################

  CALL CLOSE_LOGGER()

  DEALLOCATE (potential)
  DEALLOCATE (Hamiltonian_1_crs)
  DEALLOCATE (row_1_crs)
  DEALLOCATE (column_1_crs)
  DEALLOCATE (Particle_density)
  DEALLOCATE (Psi_1)
  DEALLOCATE (Psi_cross)
  DEALLOCATE (Psi_1_noSO)
  DEALLOCATE (Psi_LR)
  DEALLOCATE (Psi_init)
  DEALLOCATE (Psi_density_init)
  DEALLOCATE (Energies_1)
  DEALLOCATE (Energies_1_noSO)
  DEALLOCATE (Nxm_single_elems)
  DEALLOCATE (C_single_time)
  DEALLOCATE (C_single_max_time)
  DEALLOCATE (A_single_crank_nicolson)
  DEALLOCATE (B_single_crank_nicolson)
  DEALLOCATE (IPIV_single)
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
  DEALLOCATE (Cm)
  DEALLOCATE (C_max_time)
  DEALLOCATE (A_crank_nicolson)
  DEALLOCATE (B_crank_nicolson)
  DEALLOCATE (IPIV)
  DEALLOCATE (Spin_t)
  DEALLOCATE (Spin_t_single)
  DEALLOCATE(many_body_sigma_z_expected_value_ind)
END PROGRAM MAIN
