#include "macros_def.f90"

MODULE postprocessing
USE indata
USE utility
USE writers
USE utility
USE many_body
USE logger

IMPLICIT NONE
CONTAINS

SUBROUTINE CALCULATE_INTEGRAL_ELEMENTS()
  IMPLICIT NONE

  COMPLEX*16, ALLOCATABLE :: Psi_1(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
  COMPLEX*16, ALLOCATABLE :: V_tilde_upper(:, :)
  COMPLEX*16, ALLOCATABLE :: V_tilde_slice(:)
  INTEGER*4 :: ham_1_size, ham_2_size, v_tilde_elems
  COMPLEX*16 :: integral_value
  INTEGER*4 :: i

  CALL INDATA_GET("/net/ascratch/people/plgjczarnecki/LAO-STO-QD-testing-2e/RUN_Bz_1.0/OutputData/quantum_dot.nml")

  ham_1_size = get_ham_1_size(Nx, Ny, norbs)
  ham_2_size = get_ham_2_size(nstate_1, k_electrons)
  v_tilde_elems = (nstate_1 * (nstate_1 + 1)) / 2 !Number of elements in upper triangle of hermitian matrix V_tilde

  ALLOCATE (Psi_1(ham_1_size, nstate_1))
  ALLOCATE (V_tilde_upper(v_tilde_elems, ham_1_size))
  ALLOCATE (V_tilde_slice(ham_1_size))

  CALL READ_SINGLE_ELECTRON_WAVEFUNCTIONS("/net/ascratch/people/plgjczarnecki/LAO-STO-QD-testing-2e/RUN_Bz_1.0", Psi_1, ham_1_size, nstate_1, norbs)

  CALL CALCULATE_V_TILDE(Psi_1, ham_1_size, nstate_1, V_tilde_upper, v_tilde_elems, norbs, Nx, Ny, dx)

  OPEN (unit=1, FILE='./OutputData/V_integrals.dat', FORM="FORMATTED", ACTION="WRITE")
  WRITE (1, *) '#No. state    V_{iiii}'
  DO i = 1, nstate_1
    WRITE (log_string, *) "V_{iiii} integral for i = ", i
    LOG_INFO(log_string)

    CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper, ham_1_size, nstate_1, v_tilde_elems, i, i)
    CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, i), Psi_1(:, i),&
                                          & Psi_1(:, i), Psi_1(:, i),&
                                          & V_tilde_slice(:), ham_1_size, integral_value, norbs, eps_r)
    WRITE (1, *) i, integral_value
  END DO
  CLOSE (1)

  DEALLOCATE (Psi_1)
  DEALLOCATE (V_tilde_upper)
  DEALLOCATE (V_tilde_slice)

END SUBROUTINE

SUBROUTINE CALCULATE_IMAGE_EXPECTATION_VALUE(path)
  CHARACTER(len=*), INTENT(IN) :: path
  REAL*8, ALLOCATABLE :: Potential(:, :)
  COMPLEX*16, ALLOCATABLE :: Psi_1(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
  INTEGER*4 :: ham_1_size, ham_2_size, v_tilde_elems
  CHARACTER(len=400) :: filename

  WRITE (filename, '(A)') TRIM(path)//"/OutputData/quantum_dot.nml"
  CALL INDATA_GET(filename)

  ham_1_size = get_ham_1_size(Nx, Ny, norbs)
  ham_2_size = get_ham_2_size(nstate_1, k_electrons)
  v_tilde_elems = (nstate_1 * (nstate_1 + 1)) / 2 !Number of elements in upper triangle of hermitian matrix V_tilde

  ALLOCATE (Potential(-Nx:Nx, -Ny:Ny))
  ALLOCATE (Psi_1(ham_1_size, nstate_1))

  WRITE (filename, '(A)') TRIM(path)//"/OutputData/Potential_image.dat"
  CALL READ_POTENTIAL(filename, Potential, Nx, Ny)
  CALL READ_SINGLE_ELECTRON_WAVEFUNCTIONS(path, Psi_1, ham_1_size, 1, norbs)

  WRITE (filename, '(A)') TRIM(path)//"/OutputData/V_image_1e.dat"
  OPEN (11, FILE=filename, ACTION='WRITE', FORM='FORMATTED')
  WRITE (11, '(E20.8)') REAL(potential_expected_value(Psi_1(:, 1), Psi_1(:, 1), Potential, ham_1_size, Nx, Ny, norbs) / eV2au)
  CLOSE (11)

  DEALLOCATE (Potential)
  DEALLOCATE (Psi_1)
END SUBROUTINE CALCULATE_IMAGE_EXPECTATION_VALUE

END MODULE postprocessing
