#include "macros_def.f90"

MODULE potentials
USE utility
USE constants
USe indata
USE logger
IMPLICIT NONE
CONTAINS

SUBROUTINE COMPUTE_HARMONIC_OSCILLATOR_POTENTIAL(Potential, Nx, Ny, omega, m_eff)
  INTEGER*4, INTENT(IN) :: Nx, Ny
  REAL*8, INTENT(IN) :: omega, m_eff
  REAL*8, INTENT(OUT) :: Potential(-Nx:Nx, -Ny:Ny)
  INTEGER*4 :: ix, iy
  REAL*8 :: x, y

  Potential = 0.0d0
  DO ix = -Nx, Nx
    x = ix * dx
    DO iy = -Ny, Ny
      y = iy * dx
      Potential(ix, iy) = 0.5d0 * m_eff * omega**2 * (x**2 + y**2)
    END DO
  END DO

END SUBROUTINE COMPUTE_HARMONIC_OSCILLATOR_POTENTIAL

SUBROUTINE COMPUTE_GAUSSIAN_POTENTIAL(Potential, Nx, Ny, V0, Vb)
  INTEGER*4, INTENT(IN) :: Nx, Ny
  REAL*8, INTENT(IN) :: V0, Vb
  REAL*8, INTENT(OUT) :: Potential(-Nx:Nx, -Ny:Ny)
  INTEGER*4 :: ix, iy

  REAL*8 :: Rx, Ry, Rb, mu
  REAL*8 :: x, y
  Rx = dx * Nx / 1.15
  Ry = dx * Ny / 1.4
  Rb = Ry / 2.0
  mu = 10.0

  Potential = 0.0d0
  DO ix = -Nx, Nx
    x = ix * dx
    DO iy = -Ny, Ny
      y = iy * dx
      Potential(ix, iy) = -V0 / ((1 + (x**2 / Rx**2)**mu) * (1 + (y**2 / Ry**2)**mu)) + Vb / ((1 + (x**2 / Rb**2)**mu) * (1 + (y**2 / Ry**2)**mu))
    END DO
  END DO

END SUBROUTINE COMPUTE_GAUSSIAN_POTENTIAL

SUBROUTINE CALCULATE_IMAGE_POTENTIAL(V_image, d_image, Particle_density, ham_1_size, nstate, Nx, Ny)
  INTEGER*4, INTENT(IN) :: ham_1_size, nstate, Nx, Ny
  REAL*8, INTENT(IN) :: d_image
  REAL*8, INTENT(IN) :: Particle_density(ham_1_size, nstate)
  REAL*8, INTENT(OUT) :: V_image(-Nx:Nx, -Ny:Ny)
  INTEGER*4 :: i, ix, iy
  REAL*8 :: x, y, x_prime, y_prime, r_relative

  V_image = 0.0d0
  DO ix = -Nx, Nx
    DO iy = -Ny, Ny
      x = ix * dx
      y = iy * dx
      DO i = 1, ham_1_size
        x_prime = get_x_from_psi_index(i, norbs, Nx, dx)
        y_prime = get_y_from_psi_index(i, norbs, Nx, Ny, dx)
        r_relative = SQRT((x - x_prime)**2 + (y - y_prime)**2 + d_image**2)
        ! Minus sign due to opposite charge of the image.
        ! Assuming I only care about the ground state charge density
        V_image(ix, iy) = V_image(ix, iy) - Particle_density(i, 1) / (eps_r * r_relative)
      END DO
    END DO
  END DO

END SUBROUTINE CALCULATE_IMAGE_POTENTIAL

END MODULE potentials
