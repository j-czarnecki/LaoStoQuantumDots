#include "macros_def.f90"
MODULE time_dependence

USE logger
USE utility
USE indata
USE many_body
USE constants
IMPLICIT NONE
CONTAINS

SUBROUTINE SINGLE_ELECTRON_TIME_DEPENDENCE(Psi_1, Energies_1, ham_1_size, nstate_1)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstate_1
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
  REAL*8, INTENT(IN) :: Energies_1(nstate_1) !Single electron energies

  COMPLEX*16, ALLOCATABLE :: Nxm_single_elems(:, :) ! <n|x|m> elements for single electron wavefunctions
  COMPLEX*16, ALLOCATABLE :: C_single_max_time(:)
  COMPLEX*16, ALLOCATABLE :: C_single_time(:)
  COMPLEX*16, ALLOCATABLE :: A_single_crank_nicolson(:, :)
  COMPLEX*16, ALLOCATABLE :: B_single_crank_nicolson(:, :)
  INTEGER*4, ALLOCATABLE :: IPIV_single(:)
  INTEGER*4 :: INFO
  CHARACTER(LEN=100) :: filename
  INTEGER*4 :: n, m, iomega, it
  REAL*8 :: omega_ac

  ALLOCATE (Nxm_single_elems(nstate_1, nstate_1))
  ALLOCATE (C_single_time(nstate_1))
  ALLOCATE (C_single_max_time(nstate_1))
  ALLOCATE (A_single_crank_nicolson(nstate_1, nstate_1))
  ALLOCATE (B_single_crank_nicolson(nstate_1, 1))
  ALLOCATE (IPIV_SINGLE(nstate_1))

  !Time dependent Schrodinger equation is solved using Crank-Nicholson scheme
  WRITE (log_string, *) "Calculating <n|x|m> elements for single electron"
  LOG_INFO(log_string)
  !First calculate all <n|X|m> matrix elements
  !Only upper triangle is needed, since <n|X|m> = CONJG(<m|X|n>)
  OPEN (11, FILE='./OutputData/Nxm_single_electrons.dat', ACTION='WRITE', FORM='FORMATTED')
  !$omp parallel
  !$omp do
  DO n = 1, nstate_1
    DO m = 1, nstate_1
      Nxm_single_elems(n, m) = single_electron_x_expected_value(Psi_1(:, n), Psi_1(:, m), norbs, Nx, dx, ham_1_size)
      WRITE (11, *) n, m, REAL(Nxm_single_elems(n, m)), AIMAG(Nxm_single_elems(n, m))
    END DO
  END DO
  !$omp end do
  !$omp end parallel
  CLOSE (11)

  WRITE (log_string, *) "Running Crank-Nicolson scheme for single electron"
  LOG_INFO(log_string)
  !For such energy we have to trigger transition, this is only for a check.
  !omega_ac = Energies_2(2) - Energies_2(1)

  OPEN (11, FILE='./OutputData/C_single_max_time.dat', ACTION='WRITE', FORM='FORMATTED')
  WRITE (11, *) '#omega_ac [meV] C_max_1 C_max_2 ...'
  !PRINT*, "omega steps ", N_omega_ac_steps
  !$omp parallel private(omega_ac, C_single_time, C_single_max_time, A_single_crank_nicolson, B_single_crank_nicolson, IPIV_SINGLE, INFO)
  !$omp do
  DO iomega = 1, N_omega_ac_steps
    !omega_ac = iomega*domega_ac
    omega_ac = Energies_1(iomega + 1) - Energies_1(1)
    C_single_max_time = DCMPLX(0.0d0, 0.0d0)
    C_single_time = DCMPLX(0.0d0, 0.0d0)
    C_single_time(1) = DCMPLX(1.0d0, 0.0d0) !Initial condition, assuming full occupation of lowest energy state
    WRITE (filename, '("./OutputData/Time_dependent", I0, ".dat")') iomega
    !OPEN(10, FILE = './OutputData/Time_dependent.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
    OPEN (10, FILE=filename, ACTION='WRITE', FORM='FORMATTED')
    DO it = 0, N_t_steps
      A_single_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
      B_single_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
      IPIV_SINGLE = 0
      DO n = 1, nstate_1
        A_single_crank_nicolson(n, n) = 1.0d0
        B_single_crank_nicolson(n, 1) = C_single_time(n)

        DO m = 1, nstate_1
          A_single_crank_nicolson(n, m) = A_single_crank_nicolson(n, m) - dt / (2 * imag) * energy_phase_offset(Energies_1(n), Energies_1(m), (it + 1) * dt) * v_ac(f_ac, omega_ac, (it + 1) * dt) * Nxm_single_elems(n, m)
          B_single_crank_nicolson(n, 1) = B_single_crank_nicolson(n, 1) + dt / (2 * imag) * C_single_time(m) * energy_phase_offset(Energies_1(n), Energies_1(m), it * dt) * v_ac(f_ac, omega_ac, it * dt) * Nxm_single_elems(n, m)
        END DO
      END DO

      CALL ZGESV(nstate_1, 1, A_single_crank_nicolson, nstate_1, IPIV_SINGLE, B_single_crank_nicolson, nstate_1, INFO)
      !B_single_crank_nicolson = B_single_crank_nicolson !/ SUM(ABS(B_crank_nicolson(:,1))**2)
      !IF (INFO /= 0) !PRINT*, "ERROR in ZGESV, INFO = ", INFO
      WRITE (10, *) it * dt / ns2au, (ABS(B_single_crank_nicolson(n, 1))**2, n=1, nstate_1)
      C_single_time(:) = B_single_crank_nicolson(:, 1)

      !Check for maximal value of C_time
      DO n = 1, nstate_1
        IF (ABS(C_single_time(n))**2 > ABS(C_single_max_time(n))**2) C_single_max_time(n) = C_single_time(n)
      END DO
    END DO
    CLOSE (10)
    WRITE (11, *) omega_ac / eV2au * 1e3, (ABS(C_single_max_time(n))**2, n=1, nstate_1)
    FLUSH (11)

  END DO
  !$omp end do
  !$omp end parallel
  CLOSE (11)

  DEALLOCATE (Nxm_single_elems)
  DEALLOCATE (C_single_time)
  DEALLOCATE (C_single_max_time)
  DEALLOCATE (A_single_crank_nicolson)
  DEALLOCATE (B_single_crank_nicolson)
  DEALLOCATE (IPIV_SINGLE)

END SUBROUTINE SINGLE_ELECTRON_TIME_DEPENDENCE

SUBROUTINE MANY_BODY_TIME_DEPENDENCE(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, Energies_1, Energies_2, ham_1_size, nstate_1, ham_2_size, nstate_2, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstate_1, ham_2_size, nstate_2, k_electrons
  COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstate_1) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
  COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstate_2) !Coefficients of consecutive Slater determinants (\varphi_i) for many-body hamiltonian such that many-body wavefunction
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size) !Number of changed indeces in given element of multi-body hamiltonian.
  INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2) !Specifies which index has been changed to whic for example 1 -> 2 and 3 -> 4.

  REAL*8, INTENT(IN) :: Energies_1(nstate_1) !Single electron energies
  REAL*8, INTENT(IN) :: Energies_2(nstate_2) !Single electron energies

  COMPLEX*16, ALLOCATABLE :: Nxm_elems(:, :) !matrix elements of x-position operator <n|X|m> wheren n and m are multi-body eigenstates
  COMPLEX*16, ALLOCATABLE :: C_max_time(:)
  COMPLEX*16, ALLOCATABLE :: C_time(:) !Those are the time-dependent coefficients of the many-body wavefunction, such that
  !\Psi(t) = \sum_{m} [c_m(t) * exp(-i*E_n*t / \hbar) |m> ].
  ! Where m is a multi-body eigenstate.
  ! first index denotes time, second is the number of state. Time is stored such that t(1) > t(2) > ... > t(n)
  COMPLEX*16, ALLOCATABLE :: A_crank_nicolson(:, :)
  COMPLEX*16, ALLOCATABLE :: B_crank_nicolson(:, :)
  INTEGER*4, ALLOCATABLE :: IPIV(:)
  INTEGER*4 :: INFO
  CHARACTER(LEN=100) :: filename
  INTEGER*4 :: n, m, iomega, it
  REAL*8 :: omega_ac

  ALLOCATE (Nxm_elems(nstate_2, nstate_2))
  ALLOCATE (C_time(nstate_2))
  ALLOCATE (C_max_time(nstate_2))
  ALLOCATE (A_crank_nicolson(nstate_2, nstate_2))
  ALLOCATE (B_crank_nicolson(nstate_2, 1))
  ALLOCATE (IPIV(nstate_2))

  !Time dependent Schrodinger equation is solved using Crank-Nicholson scheme

  !First calculate all <n|X|m> matrix elements
  !Only upper triangle is needed, since <n|X|m> = CONJG(<m|X|n>)
  OPEN (10, FILE='./OutputData/Nxm.dat', ACTION='WRITE', FORM='FORMATTED')
  !$omp parallel
  !$omp do
  DO n = 1, nstate_2
    DO m = 1, nstate_2
      Nxm_elems(n, m) = many_body_x_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,&
        & ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, m, Nx, dx, norbs)
      WRITE (10, *) n, m, REAL(Nxm_elems(n, m)), AIMAG(Nxm_elems(n, m))
    END DO
  END DO
  !$omp end do
  !$omp end parallel
  CLOSE (10)
  !OPEN(10, FILE = './OutputData/Spin_time_evolution.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
  !Spin_t(:) = Cm(:)
  !DO it = 0, N_t_steps
  !  DO n = 1, nstate_2
  !    DO m = 1, nstate_2
  !      Spin_t(n) = Spin_t(n) + Cm(m) * many_body_sigma_z_expected_value(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstate_1, nstate_2, n, m) * EXP(-imag * Energies_2(m) * it*dt)
  !    END DO
  !  END DO
  !  WRITE(10,*) it*dt / ns2au, (REAL(Spin_t(n)), n = 1, nstate_2)
  !END DO
  !CLOSE(10)
  WRITE (log_string, *) "Running Crank-Nicolson scheme"
  LOG_INFO(log_string)
  !For such energy we have to trigger transition, this is only for a check.
  !omega_ac = Energies_2(2) - Energies_2(1)

  OPEN (10, FILE='./OutputData/C_max_time.dat', ACTION='WRITE', FORM='FORMATTED')
  WRITE (10, *) '#omega_ac [meV] C_max_1 C_max_2 ...'
  !$omp parallel private(omega_ac, C_time, C_max_time, A_crank_nicolson, B_crank_nicolson, ipiv, INFO)
  !$omp do
  DO iomega = 1, N_omega_ac_steps
    omega_ac = iomega * domega_ac
    C_max_time = DCMPLX(0.0d0, 0.0d0)
    C_time = DCMPLX(0.0d0, 0.0d0)
    C_time(1) = DCMPLX(1.0d0, 0.0d0) !Initial condition, assuming full occupation of lowest energy state
    !OPEN(10, FILE = './OutputData/Time_dependent.dat', ACTION = 'WRITE', FORM = 'FORMATTED')
    DO it = 0, N_t_steps
      A_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
      B_crank_nicolson = DCMPLX(0.0d0, 0.0d0)
      ipiv = 0
      DO n = 1, nstate_2
        A_crank_nicolson(n, n) = 1.0d0
        B_crank_nicolson(n, 1) = C_time(n)

        DO m = 1, nstate_2
          A_crank_nicolson(n, m) = A_crank_nicolson(n, m) - dt / (2 * imag) * energy_phase_offset(Energies_2(n), Energies_2(m), (it + 1) * dt) * v_ac(f_ac, omega_ac, (it + 1) * dt) * Nxm_elems(n, m)
          B_crank_nicolson(n, 1) = B_crank_nicolson(n, 1) + dt / (2 * imag) * C_time(m) * energy_phase_offset(Energies_2(n), Energies_2(m), it * dt) * v_ac(f_ac, omega_ac, it * dt) * Nxm_elems(n, m)
        END DO
      END DO

      CALL ZGESV(nstate_2, 1, A_crank_nicolson, nstate_2, IPIV, B_crank_nicolson, nstate_2, INFO)
      B_crank_nicolson = B_crank_nicolson !/ SUM(ABS(B_crank_nicolson(:,1))**2)
      !IF (INFO /= 0) !PRINT*, "ERROR in ZGESV, INFO = ", INFO
      !WRITE(10,*) it*dt / ns2au, (ABS(B_crank_nicolson(n,1))**2, n = 1, nstate_2)
      C_time(:) = B_crank_nicolson(:, 1)

      !Check for maximal value of C_time
      DO n = 1, nstate_2
        IF (ABS(C_time(n))**2 > ABS(C_max_time(n))**2) C_max_time(n) = C_time(n)
      END DO
    END DO

    WRITE (10, *) omega_ac / eV2au * 1e3, (ABS(C_max_time(n))**2, n=1, nstate_2)
    FLUSH (10)
  END DO
  !$omp end do
  !$omp end parallel
  CLOSE (10)

  DEALLOCATE (Nxm_elems)
  DEALLOCATE (C_time)
  DEALLOCATE (C_max_time)
  DEALLOCATE (A_crank_nicolson)
  DEALLOCATE (B_crank_nicolson)
  DEALLOCATE (IPIV)

END SUBROUTINE MANY_BODY_TIME_DEPENDENCE

END MODULE time_dependence
