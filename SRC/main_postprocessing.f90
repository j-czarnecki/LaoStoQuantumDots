#include "macros_def.f90"
PROGRAM MAIN_POSTPROCESSING
  USE indata
  USE utility
  USE many_body
  USE logger

  USE omp_lib

  IMPLICIT NONE

  COMPLEX*16, ALLOCATABLE :: Psi_1(:, :) !Eigenfunction of single-electron hamiltonian, describing real-space wavefunction
  COMPLEX*16, ALLOCATABLE :: V_tilde_upper(:,:)
  COMPLEX*16, ALLOCATABLE :: V_tilde_slice(:)
  INTEGER*4 :: ham_1_size, ham_2_size, v_tilde_elems
  COMPLEX*16 :: integral_value
  INTEGER*4 :: i

  CALL INIT_LOGGER()

  CALL INDATA_GET("/net/ascratch/people/plgjczarnecki/LAO-STO-QD-testing-2e/RUN_Bz_1.0/OutputData/quantum_dot.nml")


  ham_1_size = get_ham_1_size(Nx, Ny, norbs)
  ham_2_size = get_ham_2_size(nstate_1, k_electrons)
  v_tilde_elems = (nstate_1*(nstate_1 + 1))/2 !Number of elements in upper triangle of hermitian matrix V_tilde


  ALLOCATE (Psi_1(ham_1_size, nstate_1))
  ALLOCATE (V_tilde_upper(v_tilde_elems, ham_1_size))
  ALLOCATE (V_tilde_slice(ham_1_size))

  CALL READ_SINGLE_ELECTRON_WAVEFUNCTIONS("/net/ascratch/people/plgjczarnecki/LAO-STO-QD-testing-2e/RUN_Bz_1.0", Psi_1, ham_1_size, nstate_1, norbs)

  CALL CALCULATE_V_TILDE(Psi_1, ham_1_size, nstate_1, V_tilde_upper, v_tilde_elems, norbs, Nx, Ny, dx)

  OPEN(unit = 1, FILE= './OutputData/V_integrals.dat', FORM = "FORMATTED", ACTION = "WRITE")
  WRITE(1,*) '#No. state    V_{iiii}'
  DO i = 1, nstate_1
    WRITE(log_string,*) "V_{iiii} integral for i = ", i
    LOG_INFO(log_string)

    CALL GET_SLICE_FROM_HERMITIAN_MATRIX(V_tilde_slice, V_tilde_upper,  ham_1_size, nstate_1, v_tilde_elems, i, i)
    CALL CALCULATE_INTERACTION_ELEMENTS(Psi_1(:, i), Psi_1(:, i),&
                                          & Psi_1(:, i), Psi_1(:, i),&
                                          & V_tilde_slice(:), ham_1_size, integral_value, norbs, eps_r)
    WRITE(1,*) i, integral_value
  END DO
  CLOSE(1)

  CALL CLOSE_LOGGER()

  DEALLOCATE (Psi_1)
  DEALLOCATE (V_tilde_upper)
  DEALLOCATE (V_tilde_slice)

END PROGRAM