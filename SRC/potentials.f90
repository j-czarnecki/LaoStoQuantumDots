#include "macros_def.f90"

MODULE potentials
USE utility
USE constants
USe indata
USE logger
IMPLICIT NONE
CONTAINS

SUBROUTINE COMPUTE_GAUSSIAN_PACKETS_POTENTIAL(Potential, d_image, eps_r, x_shift, sigma_x, Nx, Ny, dx)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: Nx, Ny
  REAL*8, INTENT(IN) :: x_shift, sigma_x, dx, d_image, eps_r
  REAL*8, INTENT(OUT) :: Potential(-Nx:Nx, -Ny:Ny)

  REAL*8 :: Charge_1(-Nx:Nx, -Ny:Ny), Charge_2(-Nx:Nx, -Ny:Ny)
  INTEGER*4 :: ix1, iy1, ix2, iy2
  REAL*8 :: x1, y1, x2, y2
  REAL*8 :: r

  Potential = 0.0d0

  DO ix1 = -Nx, Nx
    x1 = ix1 * dx
    DO iy1 = -Ny, Ny
      y1 = iy1 * dx
      Charge_1(ix1, iy1) = calculate_gaussian_particle_density(x1, y1, sigma_x, sigma_x, x_shift, 0.0d0)
      Charge_2(ix1, iy1) = calculate_gaussian_particle_density(x1, y1, sigma_x, sigma_x, -x_shift, 0.0d0)
    END DO
  END DO

  Charge_1 = Charge_1 / SUM(Charge_1)
  Charge_2 = Charge_2 / SUM(Charge_2)

  WRITE (log_string, *) "Charge_1 sum = ", SUM(Charge_1)
  LOG_INFO(log_string)

  WRITE (log_string, *) "Charge_2 sum = ", SUM(Charge_2)
  LOG_INFO(log_string)

  DO ix1 = -Nx, Nx
    x1 = ix1 * dx
    DO iy1 = -Ny, Ny
      y1 = iy1 * dx
      DO ix2 = -Nx, Nx
        x2 = ix2 * dx
        DO iy2 = -Ny, Ny
          y2 = iy2 * dx
          r = SQRT((x1 - x2)**2 + (y1 - y2)**2 + d_image**2)
          Potential(ix1, iy1) = Potential(ix1, iy1) - (Charge_1(ix2, iy2) + Charge_2(ix2, iy2)) / (eps_r * r)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_GAUSSIAN_PACKETS_POTENTIAL

PURE RECURSIVE REAL * 8 FUNCTION calculate_gaussian_particle_density(x, y, sigma_x, sigma_y, x_shift, y_shift)
  !! Compute unnormalized charge density of a gaussian packet.
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x, y
  REAL*8, INTENT(IN) :: sigma_x, sigma_y
  REAL*8, INTENT(IN) :: x_shift, y_shift
  calculate_gaussian_particle_density = EXP(-0.5 * (((x - x_shift) / sigma_x)**2 + ((y - y_shift) / sigma_y)**2))
END FUNCTION

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
