#include "macros_def.f90"
MODULE swap
    USE indata
    USE constants
    USE logger
    USE utility
    USE combinatory
    USE many_body
    
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE GET_LEFT_RIGHT_WAVEFUNCTION(Psi_1_noSO, nstate_1, ham_1_size, Psi_LR, nx, ny, norbs)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nstate_1, ham_1_size, norbs, ny, nx
    COMPLEX*16, INTENT(IN) :: Psi_1_noSO(ham_1_size, nstate_1)

    COMPLEX*16, INTENT(OUT) :: Psi_LR(ham_1_size, 2)
    REAL*8 :: mask_L(ham_1_size), mask_R(ham_1_size)
    INTEGER*4 :: i1, i, j, n

    mask_L(:) = 0.0
    mask_R(:) = 0.0
    n = 0

    DO i1 = 1, ham_1_size
      Psi_LR(i1,1) = 1.0/SQRT(2.0) * (Psi_1_noSO(i1,3) - Psi_1_noSO(i1,1))
      Psi_LR(i1,2) = 1.0/SQRT(2.0) * (Psi_1_noSO(i1,4) + Psi_1_noSO(i1,2))
    END DO
    
    END SUBROUTINE

    SUBROUTINE GET_LEFT_RIGHT_WAVEFUNCTION_PURE(Psi_1_noSO, nstate_1, ham_1_size, Psi_LR, nx, ny, norbs)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nstate_1, ham_1_size, norbs, ny, nx
    COMPLEX*16, INTENT(IN) :: Psi_1_noSO(ham_1_size, nstate_1)

    COMPLEX*16, INTENT(OUT) :: Psi_LR(ham_1_size, 2)
    REAL*8 :: mask_L(ham_1_size), mask_R(ham_1_size)
    REAl*8 :: normL, normR
    INTEGER*4 :: i1, i, j, n

    mask_L(:) = 0.0
    mask_R(:) = 0.0
    n = 0
    DO j = -ny, ny
      DO i = -nx, nx
        IF (i < 0) THEN
          mask_L(n*norbs+1:(n+1)*norbs) = 1.0
          mask_R(n*norbs+1:(n+1)*norbs) = 0.0
        ELSE
          mask_L(n*norbs+1:(n+1)*norbs) = 0.0
          mask_R(n*norbs+1:(n+1)*norbs) = 1.0
        END IF
        n = n + 1
      END DO
    END DO

    DO i1 = 1, ham_1_size,norbs
      Psi_LR(i1,1) = 1.0/SQRT(2.0) * (Psi_1_noSO(i1,3) + Psi_1_noSO(i1,1)) * mask_L(i1)
      Psi_LR(i1+1,2) = 1.0/SQRT(2.0) * (Psi_1_noSO(i1,3) - Psi_1_noSO(i1,1)) * mask_R(i1)
    END DO
    normL = REAL(SUM(ABS(Psi_LR(:,1))**2))
    normR = REAL(SUM(ABS(Psi_LR(:,2))**2))

    Psi_LR(:,1) = Psi_LR(:,1)/SQRT(normL)
    Psi_LR(:,2) = Psi_LR(:,2)/SQRT(normR)
    
  END SUBROUTINE

  SUBROUTINE GET_INIT_COEFFICIENTS(Psi_LR, Psi_1, C_slater, ham_1_size, ham_2_size,  nstate_1_noSO, nstate_1, nstate_2, Combinations,k_electrons, Cm)
    IMPLICIT NONE 

    COMPLEX*16, INTENT(IN) :: Psi_LR(ham_1_size, nstate_1_noSO), C_slater(ham_2_size, nstate_2), Psi_1(ham_1_size, nstate_1)
    INTEGER*4, INTENT(IN) :: ham_1_size, nstate_1_noSO, ham_2_size, nstate_2, k_electrons, nstate_1
    INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
    COMPLEX*16, INTENT(OUT) :: Cm(nstate_2)
    COMPLEX*16 :: li, rj, lj, ri, coeff
    INTEGER*4 :: i, j, im, k
    Cm(:) = DCMPLX(0.0d0, 0.0d0)
    DO im = 1, nstate_2
      coeff = DCMPLX(0.0d0, 0.0d0)
      DO k = 1, ham_2_size 
        i = Combinations(k,1)
        j = Combinations(k,2)

        li = DOT_PRODUCT(Psi_1(:,i), Psi_LR(:,1))
        lj = DOT_PRODUCT(Psi_1(:,j), Psi_LR(:,1))
        ri = DOT_PRODUCT(Psi_1(:,i), Psi_LR(:,2))
        rj = DOT_PRODUCT(Psi_1(:,j), Psi_LR(:,2))
        coeff = coeff + CONJG(C_slater(k,im)) * (li * rj - lj * ri)
      END DO
      Cm(im) = coeff
    END DO

  END SUBROUTINE

  SUBROUTINE TIME_EVOLUTION_SPIN_EXPECTATION(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces, ham_1_size, ham_2_size, nmax, k_electrons, Cm, nstates_1, nstates_2, t_max_int, dt, Energies_2, Spin_t,nx, ny, norbs)
    IMPLICIT NONE

    COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
    COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2), Cm(nstates_2)
    INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
    INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
    INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
    INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size, k_electrons
    INTEGER*4, INTENT(IN) :: nstates_1, nstates_2, nmax, t_max_int, nx, ny, norbs
    REAL*8, INTENT(IN) :: Energies_2(nstates_2)
    REAL*8, INTENT(IN) :: dt 
    INTEGER*4 :: ti, n, m 
    REAL*8 :: t
    COMPLEX*16 :: Sz_L_tab(nstates_2, nstates_2), Sz_R_tab(nstates_2, nstates_2), Sx_L_tab(nstates_2, nstates_2), Sx_R_tab(nstates_2, nstates_2), Sy_L_tab(nstates_2, nstates_2), Sy_R_tab(nstates_2, nstates_2)

    COMPLEX*16, INTENT(OUT) :: Spin_t(t_max_int, nmax)
    Spin_t(:,:) = 0.0
    DO n = 1, nstates_2
      DO m = 1, nstates_2
        Sx_L_tab(n,m) = CONJG(Cm(n)) * Cm(m) * (many_body_sigma_x_expected_value_L(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs))
        Sx_R_tab(n,m) = CONJG(Cm(n)) * Cm(m) * (many_body_sigma_x_expected_value_R(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs))
        Sy_L_tab(n,m) = CONJG(Cm(n)) * Cm(m) * (many_body_sigma_y_expected_value_L(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs))
        Sy_R_tab(n,m) = CONJG(Cm(n)) * Cm(m) * (many_body_sigma_y_expected_value_R(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs))
        Sz_L_tab(n,m) = CONJG(Cm(n)) * Cm(m) * (many_body_sigma_z_expected_value_L(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs))
        Sz_R_tab(n,m) = CONJG(Cm(n)) * Cm(m) * (many_body_sigma_z_expected_value_R(Psi_1, C_slater, Combinations, N_changed_indeces, Changed_indeces,ham_1_size, ham_2_size, k_electrons, nstates_1, nstates_2, n, m, nx, ny, norbs))
      END DO
    END DO

    DO ti = 1, t_max_int
      t = (ti-1) * dt
      DO n = 1, nstates_2
        DO m = 1, nstates_2
          Spin_t(ti, 1) = Spin_t(ti, 1) + Sx_L_tab(n,m) * EXP(-imag * (Energies_2(m) - Energies_2(n)) * t)
          Spin_t(ti, 2) = Spin_t(ti, 2) + Sx_R_tab(n,m) * EXP(-imag * (Energies_2(m) - Energies_2(n)) * t)
          Spin_t(ti, 3) = Spin_t(ti, 3) + Sy_L_tab(n,m) * EXP(-imag * (Energies_2(m) - Energies_2(n)) * t)
          Spin_t(ti, 4) = Spin_t(ti, 4) + Sy_R_tab(n,m) * EXP(-imag * (Energies_2(m) - Energies_2(n)) * t)
          Spin_t(ti, 5) = Spin_t(ti, 5) + Sz_L_tab(n,m) * EXP(-imag * (Energies_2(m) - Energies_2(n)) * t)
          Spin_t(ti, 6) = Spin_t(ti, 6) + Sz_R_tab(n,m) * EXP(-imag * (Energies_2(m) - Energies_2(n)) * t)
        END DO
      END DO
    END DO
  END SUBROUTINE 

    SUBROUTINE GET_SPIN_DENSITY(Psi_1, C_slater, Combinations, &
                        N_changed_indeces, Changed_indeces, &
                        ham_1_size, ham_2_size, k_electrons, &
                        nstates_1, nstates_2, norbs, Nx, Ny, &
                        spin_density)

        IMPLICIT NONE

        ! --- INPUT ---
        INTEGER*4, INTENT(IN) :: ham_1_size, ham_2_size
        INTEGER*4, INTENT(IN) :: k_electrons
        INTEGER*4, INTENT(IN) :: nstates_1, nstates_2
        INTEGER*4, INTENT(IN) :: norbs
        INTEGER*4, INTENT(IN) :: Nx, Ny

        COMPLEX*16, INTENT(IN) :: Psi_1(ham_1_size, nstates_1)
        COMPLEX*16, INTENT(IN) :: C_slater(ham_2_size, nstates_2)

        INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
        INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
        INTEGER*4, INTENT(IN) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)

        ! --- OUTPUT ---
        COMPLEX*16, INTENT(OUT) :: spin_density(-Nx:Nx, -Ny:Ny)

        ! --- LOCAL ---
        INTEGER*4 :: ix, iy
        INTEGER*4 :: n, m,i, j
        INTEGER*4 :: site_index
        INTEGER*4 :: orb_start, orb_end

        INTEGER*4, ALLOCATABLE :: Ordering(:,:)

        ! -----------------------------------------
        spin_density(:,:) = (0.0d0, 0.0d0)

        
        ALLOCATE(Ordering(-nx:nx, -ny:ny))

        !Complete natural ordering matrix
        !.......
        !5 6 7 8
        !1 2 3 4
        Ordering(:,:) = 0
        n = 0
        DO j = -ny, ny
        DO i = -nx, nx
            Ordering(i, j) = n
            n = n + 1
        END DO
        END DO


        
        DO n = 1, nstates_2
            DO m = 1, nstates_2

                DO iy = -Ny, Ny
                    DO ix = -Nx, Nx
                        site_index = (iy + Ny)*(2*Nx + 1) + (ix + Nx)
                        

                        orb_start = Ordering(ix,iy)*norbs + 1
                        orb_end   = orb_start + norbs - 1

                        spin_density(ix,iy) = spin_density(ix,iy) + &
                            many_body_sigma_z_expected_value( &
                                Psi_1(orb_start:orb_end, :), &
                                C_slater, Combinations, &
                                N_changed_indeces, Changed_indeces, &
                                norbs, ham_2_size, k_electrons, &
                                nstates_1, nstates_2, n, m )

                    END DO
                END DO

            END DO
        END DO

    END SUBROUTINE GET_SPIN_DENSITY

END MODULE swap

