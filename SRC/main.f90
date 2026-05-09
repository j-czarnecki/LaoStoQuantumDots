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
USE time_dependence
USE potentials
USE omp_lib
USE broydenV2

IMPLICIT NONE

COMPLEX*16, ALLOCATABLE :: Hamiltonian_1_crs(:) !One electron hamiltonian in compressed row storage
INTEGER*4, ALLOCATABLE :: column_1_crs(:)
INTEGER*4, ALLOCATABLE :: row_1_crs(:)

COMPLEX*16, ALLOCATABLE :: Hamiltonian_1_no_potential_crs(:) !One electron hamiltonian in compressed row storage
INTEGER*4, ALLOCATABLE :: column_1_no_potential_crs(:)
INTEGER*4, ALLOCATABLE :: row_1_no_potential_crs(:)

COMPLEX*16, ALLOCATABLE :: Hamiltonian_2_crs(:) !Many-body hamiltonian in compressed row storage
INTEGER*4, ALLOCATABLE :: column_2_crs(:)
INTEGER*4, ALLOCATABLE :: row_2_crs(:)
REAL*8, ALLOCATABLE :: Potential_confinement(:, :)
REAL*8, ALLOCATABLE :: Potential_image(:, :)
REAL*8, ALLOCATABLE :: Potential_image_new(:, :)
REAL*8, ALLOCATABLE :: Potential_image_broyden(:)
REAL*8, ALLOCATABLE :: Potential_image_broyden_new(:)
REAL*8, ALLOCATABLE :: Potential_total(:, :)
REAL*8, ALLOCATABLE :: Potential_zero(:, :)
!COMPLEX*16, ALLOCATABLE :: psi_lapack(:, :, :, :, :)
!REAL*8, ALLOCATABLE :: ev_lapack(:)
COMPLEX*16, ALLOCATABLE :: Psi_1(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
COMPLEX*16, ALLOCATABLE :: Psi_cross(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
COMPLEX*16, ALLOCATABLE :: Psi_1_noSO(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction (no SO for bonding and anti-bonding states)
COMPLEX*16, ALLOCATABLE :: Psi_LR(:, :) !Wavefunction localized in left / right quantum dot
COMPLEX*16, ALLOCATABLE :: Psi_init(:, :) !Initial wavefunction
COMPLEX*16, ALLOCATABLE :: Psi_density_init(:) !Initial density
COMPLEX*16, ALLOCATABLE :: C_slater(:, :) !Coefficients of consecutive Slater determinants (\varphi_i) for many-body hamiltonian such that many-body wavefunction
!is given by \Psi_2 = \sum_{i=1}^{n} [c_i \varphi_i(r_1, r_2, ..., r_n)]
REAL*8, ALLOCATABLE :: Particle_density(:, :)
REAL*8, ALLOCATABLE :: Particle_density_up(:, :)
REAL*8, ALLOCATABLE :: Particle_density_down(:, :)
REAL*8, ALLOCATABLE :: Energies_1(:) !Single electron energies
REAL*8, ALLOCATABLE :: Energies_1_noSO(:) !Single electron energies (no SO for bonding and anti-bonding states)
COMPLEX*16, ALLOCATABLE :: Cm(:)
REAL*8, ALLOCATABLE :: Energies_2(:) !Many-body energies in the same order as multi-body wavefunctions in \Psi_2
INTEGER*4, ALLOCATABLE :: Combination_current(:)
INTEGER*4, ALLOCATABLE :: Combinations(:, :)
INTEGER*1, ALLOCATABLE :: N_changed_indeces(:, :) !Number of changed indeces in given element of multi-body hamiltonian.
!Needed for Slater-Condon rules implementation.
!0 - no indeces chaned,
!1 - one index changed,
!2 - two indeces changed,
!3 - three or more indeces chnged.
INTEGER*4, ALLOCATABLE :: Changed_indeces(:, :, :, :) !Specifies which index has been changed to whic for example 1 -> 2 and 3 -> 4.
INTEGER*4, ALLOCATABLE :: N_ham_2_elems_in_prev_rows(:) !Contains number of non-zero elements in previous rows of multi-body hamiltonian.
!Needed for multithreading, so that at each row we can start filling up new elements independently.
COMPLEX*16, ALLOCATABLE :: Spin_t(:, :)
COMPLEX*16, ALLOCATABLE :: Spin_t_single(:, :)
COMPLEX*16, ALLOCATABLE :: many_body_sigma_z_expected_value_ind(:)
INTEGER*4 :: INFO
INTEGER*4 :: it, iomega
REAL*8 :: omega_ac
INTEGER :: ix, iy
INTEGER*4 :: i, n, m
REAL*8 :: x, y, Rx, Ry, Rb, mu
INTEGER*4 :: ham_1_size, nonzero_ham_1
INTEGER*4 :: ham_2_size, nonzero_ham_2
INTEGER*4 :: t_max_int
INTEGER*4 :: nmax
REAL*8 :: slater_norm
CHARACTER(LEN=100) :: filename
LOGICAL :: isSO
INTEGER*4 :: n_sc_iter
REAL*8 :: max_pot_rel_error

!OMP specific
INTEGER*4 :: max_num_threads

CALL INIT_LOGGER()

max_num_threads = omp_get_max_threads()
WRITE (log_string, *) 'Number of threads: ', max_num_threads
LOG_INFO(log_string)

!read a file
CALL INDATA_GET("./OutputData/quantum_dot.nml")
WRITE (log_string, *) 'Calculation parameters: ',&
  & 'Nx = ', Nx,&
  & 'Ny = ', Ny,&
  & 'dx = ', dx / nm2au,&
  & 'norbs = ', norbs,&
  & 'nstate_1 = ', nstate_1,&
  & 'nstate_2 = ', nstate_2,&
  & 'k_electrons = ', k_electrons
LOG_INFO(log_string)
WRITE (log_string, *) 'Physical parameters: ',&
  & 'th = ', th / eV2au,&
  & 'tl = ', tl / eV2au,&
  & 'td = ', td / eV2au,&
  & 'dso = ', dso / eV2au,&
  & 'drso = ', drso / eV2au,&
  & 'dE = ', dE / eV2au,&
  & 'g = ', g
LOG_INFO(log_string)
WRITE (log_string, *) 'External parameters: ',&
  & 'omega = ', omega / eV2au,&
  & 'Bx = ', Bx / T2au,&
  & 'By = ', By / T2au,&
  & 'Bz = ', Bz / T2au,&
  & 'domega_ac = ', domega_ac / eV2au,&
  & 'omega_ac_max = ', omega_ac_max / eV2au,&
  & 'F = ', f_ac / F2au,&
  & 'V_b = ', Vb / eV2au,&
  & 'V_0 = ', V0 / eV2au
LOG_INFO(log_string)

ham_1_size = get_ham_1_size(Nx, Ny, norbs)
nonzero_ham_1 = get_nonzero_ham_1(Nx, Ny)
ham_2_size = get_ham_2_size(nstate_1, k_electrons)
nonzero_ham_2 = get_nonzero_ham_2(ham_2_size, nstate_1, k_electrons)
WRITE (log_string, *) 'ham_1_size = ', ham_1_size, 'nonzero_ham_1 = ', nonzero_ham_1, 'ham_2_size = ', ham_2_size, 'nonzero_ham_2 = ', nonzero_ham_2
LOG_INFO(log_string)

slater_norm = SQRT(GAMMA(k_electrons + 1.0d0))
nmax = 6
t_max_int = INT(t_max / dt)
! PRINT*, t_max_int
! PRINT*, "Number of combinations: ", ham_2_size
! PRINT*, "Size of one-body hamiltonian: ", ham_1_size
! PRINT*, 'Nonzero elements in hamiltonian 1: ', nonzero_ham_1
! PRINT*, 'Nonzero elements in hamiltonian 2: ', nonzero_ham_2

ALLOCATE (Combination_current(k_electrons))
ALLOCATE (Potential_confinement(-Nx:Nx, -Ny:Ny))
ALLOCATE (Potential_image(-Nx:Nx, -Ny:Ny))
ALLOCATE (Potential_image_new(-Nx:Nx, -Ny:Ny))
ALLOCATE (Potential_total(-Nx:Nx, -Ny:Ny))
ALLOCATE (Potential_zero(-Nx:Nx, -Ny:Ny))
ALLOCATE (Hamiltonian_1_crs(nonzero_ham_1))
ALLOCATE (column_1_crs(nonzero_ham_1))
ALLOCATE (row_1_crs(ham_1_size + 1))
ALLOCATE (Hamiltonian_1_no_potential_crs(nonzero_ham_1))
ALLOCATE (column_1_no_potential_crs(nonzero_ham_1))
ALLOCATE (row_1_no_potential_crs(ham_1_size + 1))
ALLOCATE (Particle_density(ham_1_size, nstate_2))
ALLOCATE (Particle_density_up(ham_1_size, nstate_2))
ALLOCATE (Particle_density_down(ham_1_size, nstate_2))
ALLOCATE (Psi_1(ham_1_size, nstate_1))
ALLOCATE (Psi_cross(ham_1_size, 1))
ALLOCATE (Psi_1_noSO(ham_1_size, nstate_1))
ALLOCATE (Psi_LR(ham_1_size, 2))
ALLOCATE (Psi_init(ham_1_size, 2))
ALLOCATE (Psi_density_init(ham_1_size))
ALLOCATE (Energies_1(nstate_1))
ALLOCATE (Energies_1_noSO(nstate_1))
ALLOCATE (C_slater(ham_2_size, nstate_2))
ALLOCATE (Energies_2(nstate_2))
ALLOCATE (Hamiltonian_2_crs(nonzero_ham_2))
ALLOCATE (column_2_crs(nonzero_ham_2))
ALLOCATE (row_2_crs(ham_2_size + 1))
ALLOCATE (N_changed_indeces(ham_2_size, ham_2_size))
ALLOCATE (Combinations(ham_2_size, k_electrons))
ALLOCATE (Changed_indeces(ham_2_size, ham_2_size, 2, 2))
ALLOCATE (N_ham_2_elems_in_prev_rows(ham_2_size))
ALLOCATE (Cm(nstate_1))
ALLOCATE (Spin_t(t_max_int, nmax))
ALLOCATE (Spin_t_single(t_max_int, nmax))
ALLOCATE (many_body_sigma_z_expected_value_ind(k_electrons))
!PRINT*, "Memory check..."

N_changed_indeces(:, :) = 0
Changed_indeces(:, :, :, :) = 0

!PRINT*, "Memory check 2..."

Potential_confinement = 0.0d0
Potential_image = -1e-8 * eV2au
Potential_image_new = 0.0d0
Potential_zero = 0.0d0
!CALL COMPUTE_GAUSSIAN_POTENTIAL(Potential_confinement, Nx, Ny, V0, Vb)
CALL COMPUTE_HARMONIC_OSCILLATOR_POTENTIAL(Potential_confinement, Nx, Ny, omega, m_eff)

IF (initialize_image) THEN
  CALL COMPUTE_GAUSSIAN_PACKETS_POTENTIAL(Potential_image, d_image, eps_r, shift_packet, sigma_packet, Nx, Ny, dx)
END IF

CALL WRITE_POTENTIAL(Potential_confinement, Nx, Ny, './OutputData/Potential_confinement.dat')
CALL WRITE_POTENTIAL(Potential_image, Nx, Ny, './OutputData/Potential_image_initial.dat')

IF (.not. (nstate_2 < (4 * (nstate_2 + 1)) .and. nstate_2 < ham_2_size .and. (4 * (nstate_2 + 1)) < ham_2_size)) THEN
  WRITE (log_string, *) 'In ARPES NEV < NCV ≤ N, please lower nstate_2. STOP.'
  LOG_ERROR(log_string)
  STOP
END IF

! Computed for checking purposes, so that one can calculate the expectation value of energy of Hamiltonian without potential
CALL CREATE_ONE_ELECTRON_HAMILTONIAN_CRS(Hamiltonian_1_no_potential_crs, column_1_no_potential_crs, row_1_no_potential_crs, nonzero_ham_1, ham_1_size, Nx, Ny, norbs, Potential_zero)

DO n_sc_iter = 1, max_sc_iter
  WRITE (log_string, '(a, I0)') "==== SC_ITER: ", n_sc_iter
  LOG_INFO(log_string)

  Potential_total = Potential_confinement + Potential_image

  !#################### ONE-ELECTRON PROBLEM ####################
  CALL CREATE_ONE_ELECTRON_HAMILTONIAN_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Nx, Ny, norbs, Potential_total)

  CALL DIAGONALIZE_ARPACK_CRS(Hamiltonian_1_crs, column_1_crs, row_1_crs, nonzero_ham_1, ham_1_size, Psi_1, Energies_1, nstate_1)

  CALL SINGLE_ELECTRON_TIME_DEPENDENCE(Psi_1, Energies_1, ham_1_size, nstate_1)

  !Writing single-electron problem data to a file
  IF (n_sc_iter == 1) THEN
    CALL WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi_1, ham_1_size, nstate_1, norbs, Nx, Ny, dx, './OutputData/Psi_1_no_image')
    CALL WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1, ham_1_size, nstate_1, norbs, './OutputData/Expectations_1_no_image.dat')
    CALL WRITE_ENERGIES(Energies_1, nstate_1, './OutputData/Energies1_no_image.dat')
  END IF

  ! Checking whether we have to compute many-body problem
  IF (nstate_2 .ne. 0) THEN
    WRITE (log_string, *) "Solving many-body problem"
    LOG_INFO(log_string)

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

    DO n = 1, nstate_2
      WRITE (log_string, *) "Slater normalization: ", SUM(ABS(C_slater(:, n))**2)
      LOG_INFO(log_string)
    END DO

    CALL CALCULATE_PARTICLE_DENSITY(Particle_density, Psi_1, C_slater, N_changed_indeces, Changed_indeces, Combinations, ham_1_size, ham_2_size, nstate_1, nstate_2, k_electrons, dx)
    CALL WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(DCMPLX(Particle_density), ham_1_size, nstate_2, norbs, Nx, Ny, dx, './OutputData/ParticleDensity')

    ! Write many-electron data to files
    IF (n_sc_iter == 1) THEN
      CALL WRITE_SLATER_COEFFICIENTS(C_slater, ham_2_size, nstate_2, './OutputData/C_slater_no_image.dat')
      CALL WRITE_MULTI_ELECTRON_EXPECTATIONS(Hamiltonian_1_no_potential_crs, column_1_no_potential_crs, row_1_no_potential_crs, nonzero_ham_1, Potential_image, Potential_confinement, Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, dx, './OutputData/Expectations_2_no_image.dat')
      CALL WRITE_ENERGIES(Energies_2, nstate_2, './OutputData/Energies2_no_image.dat')
    END IF

    CALL MANY_BODY_TIME_DEPENDENCE(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, Energies_1, Energies_2, ham_1_size, nstate_1, ham_2_size, nstate_2, k_electrons)

    CALL CALCULATE_IMAGE_POTENTIAL(Potential_image_new, d_image, Particle_density, ham_1_size, nstate_2, Nx, Ny)
  ELSE
    WRITE (log_string, *) "Calculating image potential from 1e wavefunctions"
    LOG_INFO(log_string)
    CALL CALCULATE_IMAGE_POTENTIAL(Potential_image_new, d_image, ABS(Psi_1)**2, ham_1_size, nstate_1, Nx, Ny)

    WRITE (log_string, *) "Skipping two-electron problem due to nstate_2 == 0"
    LOG_INFO(log_string)
  END IF

  ! Convergence condition check
  max_pot_rel_error = MAXVAL(ABS((Potential_image - Potential_image_new) / Potential_image))
  WRITE (log_string, *) "Max relative error is", max_pot_rel_error
  LOG_INFO(log_string)

  IF (max_pot_rel_error < eps_potential) THEN
    WRITE (log_string, *) "Convergence reached!"
    LOG_INFO(log_string)
    EXIT
  END IF

  ! Broyden mixing of image potential
  Potential_image_broyden = RESHAPE(Potential_image, [SIZE(Potential_image)])
  Potential_image_broyden_new = RESHAPE(Potential_image_new, [SIZE(Potential_image_new)])

  CALL MIX_BROYDEN(SIZE(Potential_image_broyden), Potential_image_broyden_new, Potential_image_broyden, sc_alpha, n_sc_iter, 4, .FALSE.)

  Potential_image = RESHAPE(Potential_image_broyden, SHAPE(Potential_image))

  CALL WRITE_POTENTIAL(Potential_image, Nx, Ny, './OutputData/Potential_image_iter.dat')
  CALL WRITE_POTENTIAL(Potential_total, Nx, Ny, './OutputData/Potential_total_iter.dat')

END DO ! End of self-consistent loop

IF (n_sc_iter > 2) THEN
  ! Write 1e data
  CALL WRITE_SINGLE_ELECTRON_WAVEFUNCTIONS(Psi_1, ham_1_size, nstate_1, norbs, Nx, Ny, dx, './OutputData/Psi_1')
  CALL WRITE_SINGLE_ELECTRON_EXPECTATIONS(Psi_1, ham_1_size, nstate_1, norbs, './OutputData/Expectations_1.dat')
  CALL WRITE_ENERGIES(Energies_1, nstate_1, './OutputData/Energies1.dat')

  ! Write 2e data
  CALL WRITE_SLATER_COEFFICIENTS(C_slater, ham_2_size, nstate_2, './OutputData/C_slater.dat')
  CALL WRITE_MULTI_ELECTRON_EXPECTATIONS(Hamiltonian_1_no_potential_crs, column_1_no_potential_crs, row_1_no_potential_crs, nonzero_ham_1, Potential_image, Potential_confinement, Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, norbs, Nx, Ny, dx, './OutputData/Expectations_2.dat')
  CALL WRITE_ENERGIES(Energies_2, nstate_2, './OutputData/Energies2.dat')
END IF

CALL WRITE_POTENTIAL(Potential_image, Nx, Ny, './OutputData/Potential_image.dat')
CALL WRITE_POTENTIAL(Potential_total, Nx, Ny, './OutputData/Potential_total.dat')

CALL CLOSE_LOGGER()

CALL MIX_BROYDEN(SIZE(Potential_image_broyden), Potential_image_broyden_new, Potential_image_broyden, sc_alpha, n_sc_iter, 4, .TRUE.)

DEALLOCATE (Potential_confinement)
DEALLOCATE (Potential_image)
DEALLOCATE (Potential_image_new)
DEALLOCATE (Potential_total)
DEALLOCATE (Potential_zero)
DEALLOCATE (Hamiltonian_1_crs)
DEALLOCATE (row_1_crs)
DEALLOCATE (column_1_crs)
DEALLOCATE (Hamiltonian_1_no_potential_crs)
DEALLOCATE (row_1_no_potential_crs)
DEALLOCATE (column_1_no_potential_crs)
DEALLOCATE (Particle_density)
DEALLOCATE (Particle_density_up)
DEALLOCATE (Particle_density_down)
DEALLOCATE (Psi_1)
DEALLOCATE (Psi_cross)
DEALLOCATE (Psi_1_noSO)
DEALLOCATE (Psi_LR)
DEALLOCATE (Psi_init)
DEALLOCATE (Psi_density_init)
DEALLOCATE (Energies_1)
DEALLOCATE (Energies_1_noSO)
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
DEALLOCATE (Cm)
DEALLOCATE (Spin_t)
DEALLOCATE (Spin_t_single)
DEALLOCATE (many_body_sigma_z_expected_value_ind)
END PROGRAM MAIN
