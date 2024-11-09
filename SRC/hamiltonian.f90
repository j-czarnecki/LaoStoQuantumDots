#include "macros_def.f90"
MODULE hamiltonian
  USE indata
  USE constants
  USE logger
  IMPLICIT NONE

CONTAINS

  SUBROUTINE CREATE_ONE_ELECTRON_HAMILTONIAN_CRS(Ham_1, column_crs, row_crs, nonzero, ham_1_size, nx, ny, norbitals, Potential)
    IMPLICIT NONE
    COMPLEX*16, INTENT(OUT) :: Ham_1(nonzero)
    INTEGER*4, INTENT(OUT) :: column_crs(nonzero), row_crs(ham_1_size + 1)
    INTEGER*4, INTENT(IN) :: ham_1_size, nx, ny, norbitals, nonzero
    REAL*8, INTENT(IN) :: Potential(-nx:nx, -ny:ny)

    COMPLEX*16 :: H_0(norbitals, norbitals)
    COMPLEX*16 :: H_so(norbitals, norbitals)
    COMPLEX*16 :: H_B(norbitals, norbitals)
    COMPLEX*16 :: H_x(norbitals, norbitals)
    COMPLEX*16 :: H_x_const(norbitals, norbitals)
    COMPLEX*16 :: H_y(norbitals, norbitals)
    COMPLEX*16 :: H_mix(norbitals, norbitals)
    COMPLEX*16 :: Unity(norbitals, norbitals)
    COMPLEX*16 :: H_onsite(norbitals, norbitals)

    INTEGER*4, ALLOCATABLE :: Ordering(:,:)

    COMPLEX*16 :: orb
    REAL*8 :: y
    INTEGER*4 :: i,j,n, n_row, n_col

    WRITE(log_string, *) 'Creating one-electron Hamiltonian'
    LOG_INFO(log_string)

    H_0 = dcmplx(0.0d0, 0.0d0)
    H_so = dcmplx(0.0d0, 0.0d0)
    H_B = dcmplx(0.0d0, 0.0d0)
    H_x = dcmplx(0.0d0, 0.0d0)
    H_x_const = dcmplx(0.0d0, 0.0d0)
    H_y = dcmplx(0.0d0, 0.0d0)
    H_mix = dcmplx(0.0d0, 0.0d0)
    Unity = dcmplx(0.0d0, 0.0d0)

    ALLOCATE(Ordering(-nx:nx, -ny:ny))

    !Complete natural ordering matrix
    !.......
    !5 6 7 8
    !1 2 3 4
    Ordering(:,:) = 0
    n = 1
    DO j = -ny, ny
      DO i = -nx, nx
        Ordering(i, j) = n
        n = n + 1
      END DO
    END DO


    DO n = 1, norbitals
      Unity(n,n) = 1.0d0
    END DO

    !Compute hamiltonian parts constant over iterations
    !Onsite terms without potential
    !Onsite elements - H_0, H_so, H_B do not need lower triangle,
    !because we will only store upper triangle ofthe sparse hamiltonian.
    !Therefore we do not compute FILL_HERMITIAN on them
    H_0(1,1) = 4*tl - dE
    H_0(2,2) = 4*tl - dE
    H_0(3,3) = 2*tl + 2*th
    H_0(4,4) = 2*tl + 2*th
    H_0(5,5) = 2*tl + 2*th
    H_0(6,6) = 2*tl + 2*th

    !Atomic spin-robit coupling
    H_so(1,4) = imag*dso/3.0d0
    H_so(2,3) = imag*dso/3.0d0
    H_so(1,6) = -dso/3.0d0
    H_so(2,5) = dso/3.0d0
    H_so(3,5) = imag*dso/3.0d0
    H_so(4,6) = -imag*dso/3.0d0

    !Interaction with magnetic field
    !Angular momentum coupling
    H_B(1,3) = imag*mub*Bx
    H_B(2,4) = imag*mub*Bx
    H_B(1,5) = -imag*mub*By
    H_B(2,6) = -imag*mub*By
    H_B(3,5) = imag*mub*Bz
    H_B(4,6) = imag*mub*Bz
    !Zeeman term
    DO n = 1, norbitals
      H_B(n,n) = (-1)**(n+1) * 0.5*g*mub*Bz
      IF (MOD(n,2) .EQ. 1) THEN
        H_B(n, n + 1) = 0.5*g*mub*(Bx - imag*By)
      END IF
    END DO

    !Hopping x
    !Diagonal terms
    H_x_const(1,1) = -tl
    H_x_const(2,2) = -tl
    H_x_const(3,3) = -tl
    H_x_const(4,4) = -tl
    H_x_const(5,5) = -th
    H_x_const(6,6) = -th
    !Off-diagonal terms
    !Upper triangle
    H_x_const(1,5) = -drso/2.0d0
    H_x_const(2,6) = -drso/2.0d0
    !Lower triangle
    H_x_const(5,1) = drso/2.0d0
    H_x_const(6,2) = drso/2.0d0

    !Hopping y
    !Diagonal terms
    H_y(1,1) = -tl
    H_y(2,2) = -tl
    H_y(3,3) = -th
    H_y(4,4) = -th
    H_y(5,5) = -tl
    H_y(6,6) = -tl
    !Off-diagonal terms
    !Upper triangle
    H_y(1,3) = -drso/2.0d0
    H_y(2,4) = -drso/2.0d0
    !Lower triangle
    H_y(3,1) = drso/2.0d0
    H_y(4,2) = drso/2.0d0

    !Hopping diagonal
    !Upper triangle
    H_mix(3,5) = td/2.0d0
    H_mix(4,6) = td/2.0d0
    !Lower triangle
    H_mix(5,3) = td/2.0d0
    H_mix(6,4) = td/2.0d0

    n = 1
    DO j = -ny, ny
      DO i = -nx, nx
        y = j*dx
        orb = EXP(-0.5 * imag * Bz * (dx) * (2 * y))
        H_onsite = H_0 + H_so + H_B + Unity*potential(i,j)
        H_x = H_x_const * orb

        !left facet and left bottom corner
        IF (i .EQ. -nx .AND. j .NE. ny) THEN
          !Here order of hamiltonians is crucial!
          !Compute H_onsite first
          DO n_row = 1, norbitals
            row_crs((Ordering(i,j) - 1)*norbitals + n_row) = n
            DO n_col = n_row, norbitals
              !Top triangle of H_onsite
              IF (ABS(H_onsite(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_onsite(n_row, n_col)
                column_crs(n) = (Ordering(i,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full hopping x matrix
            DO n_col = 1, norbitals
              IF (ABS(H_x(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_x(n_row, n_col)
                column_crs(n) = (Ordering(i + 1,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full hopping y matrix
            DO n_col = 1, norbitals
              IF (ABS(H_y(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_y(n_row, n_col)
                column_crs(n) = (Ordering(i,j + 1) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full x + y hopping matrix
            DO n_col = 1, norbitals
              IF (ABS(H_mix(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = -H_mix(n_row, n_col) !Note the minus
                column_crs(n) = (Ordering(i + 1, j + 1) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO
          END DO

        !right facet and right bottom corner
        ELSE IF (i .EQ. nx .AND. j .NE. ny) THEN
          DO n_row = 1, norbitals
            row_crs((Ordering(i,j) - 1)*norbitals + n_row) = n
            DO n_col = n_row, norbitals
              !Top triangle of H_onsite
              IF (ABS(H_onsite(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_onsite(n_row, n_col)
                column_crs(n) = (Ordering(i,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full -x + y hopping matrix
            DO n_col = 1, norbitals
              IF (ABS(H_mix(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_mix(n_row, n_col) !Note the minus
                column_crs(n) = (Ordering(i - 1, j + 1) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full hopping y matrix
            DO n_col = 1, norbitals
              IF (ABS(H_y(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_y(n_row, n_col)
                column_crs(n) = (Ordering(i,j + 1) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO
          END DO


        !top facet and left top corner
        ELSE IF (j .EQ. ny .AND. i .NE. nx) THEN
          DO n_row = 1, norbitals
            row_crs((Ordering(i,j) - 1)*norbitals + n_row) = n

            !Top triangle of H_onsite
            DO n_col = n_row, norbitals
              IF (ABS(H_onsite(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_onsite(n_row, n_col)
                column_crs(n) = (Ordering(i,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full hopping x matrix
            DO n_col = 1, norbitals
              IF (ABS(H_x(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_x(n_row, n_col)
                column_crs(n) = (Ordering(i + 1,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO
          END DO

        !right-top corner
        ELSE IF (i .EQ. nx .AND. j .EQ. ny) THEN
          DO n_row = 1, norbitals
            row_crs((Ordering(i,j) - 1)*norbitals + n_row) = n
            !Top triangle of H_onsite
            DO n_col = n_row, norbitals
              IF (ABS(H_onsite(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_onsite(n_row, n_col)
                column_crs(n) = (Ordering(i,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO
          END DO

        !Interior and bottom facet
        ELSE
          !Here order of hamiltonians is crucial!
          !Compute H_onsite first
          DO n_row = 1, norbitals
            row_crs((Ordering(i,j) - 1)*norbitals + n_row) = n
            DO n_col = n_row, norbitals
              !Top triangle of H_onsite
              IF (ABS(H_onsite(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_onsite(n_row, n_col)
                column_crs(n) = (Ordering(i,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full hopping x matrix
            DO n_col = 1, norbitals
              IF (ABS(H_x(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_x(n_row, n_col)
                column_crs(n) = (Ordering(i + 1,j) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full -x + y hopping matrix
            DO n_col = 1, norbitals
              IF (ABS(H_mix(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_mix(n_row, n_col) !Note the minus
                column_crs(n) = (Ordering(i - 1, j + 1) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full hopping y matrix
            DO n_col = 1, norbitals
              IF (ABS(H_y(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = H_y(n_row, n_col)
                column_crs(n) = (Ordering(i,j + 1) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO

            !Full x + y hopping matrix
            DO n_col = 1, norbitals
              IF (ABS(H_mix(n_row, n_col)) .NE. 0.0) THEN
                Ham_1(n) = -H_mix(n_row, n_col) !Note the minus
                column_crs(n) = (Ordering(i + 1, j + 1) - 1)*norbitals + n_col
                n = n + 1
              END IF
            END DO
          END DO
        END IF

      END DO !End y-position loop
    END DO !End x-position loop
    row_crs(ham_1_size + 1) = n
    !PRINT*, 'Included non-zero elements in CRS format: ', n

  END SUBROUTINE CREATE_ONE_ELECTRON_HAMILTONIAN_CRS





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!create full matrix containing Hamiltoanian !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE HAMILTONIAN_CREATE(ham, size, nnx, nny, norbitals, potential)
    IMPLICIT NONE
    SAVE
    INTEGER, INTENT(IN) :: norbitals, size
    INTEGER, INTENT(IN) :: nnx, nny
    COMPLEX*16, INTENT(OUT) :: ham(size, size)
    REAL*8, INTENT(OUT) :: potential(-nnx:nnx, -nny:nny)

    COMPLEX*16, ALLOCATABLE :: h_onsite(:, :)
    COMPLEX*16, ALLOCATABLE :: h_onsite_B(:, :)
    COMPLEX*16, ALLOCATABLE :: h_hopping_x(:, :)
    COMPLEX*16, ALLOCATABLE :: h_hopping_y(:, :)
    COMPLEX*16, ALLOCATABLE :: h_hopping_diagonal(:, :)
    COMPLEX*16, ALLOCATABLE :: h_hopping_diagonal_m(:, :)

    INTEGER :: i, j, nn, n_elems
    INTEGER, ALLOCATABLE :: ordering(:, :)
    COMPLEX*16 :: orb, one
    REAL*8 :: pot, yi
    one = dcmplx(0.0, 1.0)

    ALLOCATE (h_onsite(norbitals, norbitals))
    ALLOCATE (h_onsite_B(norbitals, norbitals))
    ALLOCATE (h_hopping_x(norbitals, norbitals))
    ALLOCATE (h_hopping_y(norbitals, norbitals))
    ALLOCATE (h_hopping_diagonal(norbitals, norbitals))
    ALLOCATE (h_hopping_diagonal_m(norbitals, norbitals))

    ALLOCATE (ordering(-nnx:nnx, -nny:nny))

    ordering = 0
    nn = 1
   DO j = -nny, nny
      DO i = -nnx, nnx
        ordering(i, j) = nn
        nn = nn + 1
      END DO
    END DO

    !Julian: Was transposed
    h_hopping_y = RESHAPE((/dcmplx(-tl, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(-tl, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) &
                                      /), SHAPE(h_hopping_y))

    !Julian: Was transposed
    h_onsite_B = RESHAPE((/dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) &
                                     /), SHAPE(h_onsite_B))

    !Hopping one to the right and one down
    !Julian: Was transposed
    h_hopping_diagonal = RESHAPE((/dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) &
                                             /), SHAPE(h_hopping_diagonal))

    !Hopping one to the right and one up
    !Julian: Was transposed
    h_hopping_diagonal_m = RESHAPE((/dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) &
                                               /), SHAPE(h_hopping_diagonal))

    ham = dcmplx(0.0, 0.0)
    n_elems = 0
    DO i = -nnx, nnx
      DO j = -nny, nny

        pot = potential(i, j)
        !pot = 0.5 * 0.286 * (37.378e-3*eV2au)**2 * ((i*dx)**2 + (j*dx)**2)
        !definition of submatrices to create a total Hamiltonian
        !Julian: Was transposed
        h_onsite = RESHAPE((/dcmplx(4.0 * tl - dE + pot + 0.5 * g * mub * Bz, 0.0), dcmplx(0.5 * g * mub * Bx, -0.5 * g * mub * By), dcmplx(0.0, mub * Bx), dcmplx(0.0, dso / 3.0), dcmplx(0.0, -mub * By), dcmplx(-dso / 3.0, 0.0), &
                                       dcmplx(0.5 * g * mub * Bx, 0.5 * g * mub * By), dcmplx(4.0 * tl - dE + pot - 0.5 * g * mub * Bz, 0.0), dcmplx(0.0, dso / 3.0), dcmplx(0.0, mub * Bx), dcmplx(dso / 3.0, 0.0), dcmplx(0.0, -mub * By), &
                                       dcmplx(0.0, -mub * Bx), dcmplx(0.0, -dso / 3.0), dcmplx(2.0 * tl + 2.0 * th + pot + 0.5 * g * mub * Bz, 0.0), dcmplx(0.5 * g * mub * Bx, -0.5 * g * mub * By), dcmplx(0.0, dso / 3.0 + mub * Bz), dcmplx(0.0, 0.0), &
                                       dcmplx(0.0, -dso / 3.0), dcmplx(0.0, -mub * Bx), dcmplx(0.5 * g * mub * Bx, 0.5 * g * mub * By), dcmplx(2.0 * tl + 2.0 * th + pot - 0.5 * g * mub * Bz, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, -dso / 3.0 + mub * Bz), &
                                       dcmplx(0.0, mub * By), dcmplx(dso / 3.0, 0.0), dcmplx(0.0, -dso / 3.0 - mub * Bz), dcmplx(0.0, 0.0), dcmplx(2.0 * tl + 2.0 * th + pot + 0.5 * g * mub * Bz, 0.0), dcmplx(0.5 * g * mub * Bx, -0.5 * g * mub * By), &
                                       dcmplx(-dso / 3.0, 0.0), dcmplx(0.0, mub * By), dcmplx(0.0, 0.0), dcmplx(0.0, dso / 3.0 - mub * Bz), dcmplx(0.5 * g * mub * Bx, 0.5 * g * mub * By), dcmplx(2.0 * tl + 2.0 * th + pot - 0.5 * g * mub * Bz, 0.0) &
                                       /), SHAPE(h_onsite))

        yi = j * dx
        !orb = EXP(-0.5 * one * Bz * (dx) * (2 * yi))
        orb = 1.0
        !Julian: Was transposed
        h_hopping_x = RESHAPE((/dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0) * orb, dcmplx(0.0, 0.0), &
                                          dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0) * orb, &
                                          dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                          dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                          dcmplx(drso / 2.0, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0) * orb, dcmplx(0.0, 0.0), &
                                          dcmplx(0.0, 0.0), dcmplx(drso / 2.0, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0) * orb &
                                          /), SHAPE(h_hopping_x))

        IF (i .EQ. -nnx .AND. j .NE. -nny .AND. j .NE. nny) THEN
          !left facet
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i, j + 1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j + 1), ham, h_hopping_diagonal_m, size, norbitals)

        ELSE IF (i .EQ. nnx .AND. j .NE. -nny .AND. j .NE. nny) THEN
          !right facet
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i, j + 1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i - 1, j + 1), ham, h_hopping_diagonal, size, norbitals) !Julian

        ELSE IF (j .EQ. nny .AND. i .NE. -nnx .AND. i .NE. nnx) THEN
          !top facet
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)
          !CALL put_matrix(ordering(i, j), ordering(i + 1, j - 1), ham, h_hopping_diagonal, size, norbitals) !Julian

        ELSE IF (i .EQ. -nnx .AND. j .EQ. -nny) THEN
          !left-bottom corner
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i, j + 1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j + 1), ham, h_hopping_diagonal_m, size, norbitals)

        ELSE IF (i .EQ. nnx .AND. j .EQ. -nny) THEN
          !right-bottom corner
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i, j + 1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i - 1, j + 1), ham, h_hopping_diagonal, size, norbitals) !Julian !!!

        ELSE IF (i .EQ. -nnx .AND. j .EQ. nny) THEN
          !left-top corner
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)
          !CALL put_matrix(ordering(i, j), ordering(i + 1, j - 1), ham, h_hopping_diagonal, size, norbitals) !!!

        ELSE IF (i .EQ. nnx .AND. j .EQ. nny) THEN
          !right-top corner
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)

        ELSE
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i, j + 1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j + 1), ham, h_hopping_diagonal_m, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i - 1, j + 1), ham, h_hopping_diagonal, size, norbitals)
          n_elems = n_elems + 1
        END IF

      END DO
    END DO

    CALL fill_hermitian(ham, size)

    DEALLOCATE (h_onsite)
    DEALLOCATE (h_onsite_B)
    DEALLOCATE (h_hopping_x)
    DEALLOCATE (h_hopping_y)
    DEALLOCATE (h_hopping_diagonal)
    DEALLOCATE (h_hopping_diagonal_m)
    DEALLOCATE (ordering)
  END SUBROUTINE HAMILTONIAN_CREATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! put submatirx into large matrix !!!!!!!!! !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE put_matrix(ni, nj, matrix_A, matrix_B, sizeA, sizeB)

    INTEGER, INTENT(IN) :: ni, nj
    INTEGER, INTENT(IN) :: sizeA
    INTEGER, INTENT(IN) :: sizeB
    COMPLEX*16, INTENT(INOUT) :: matrix_A(sizeA, sizeA)
    COMPLEX*16, INTENT(IN) :: matrix_B(sizeB, sizeB)

    INTEGER :: i, j

    IF (MOD(sizeA, sizeB) .NE. 0) STOP "put_matrix: sizeA should be an integral multiple of sizeB"

    DO i = 1, sizeB
      DO j = 1, sizeB
        matrix_A((ni - 1) * sizeB + i, (nj - 1) * sizeB + j) = matrix_A((ni - 1) * sizeB + i, (nj - 1) * sizeB + j) + matrix_B(i, j)
      END DO
    END DO
  END SUBROUTINE put_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! make full matrix from the upper part when it is hermitian !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fill_hermitian(matrix_A, sizeA)

    INTEGER, INTENT(IN) :: sizeA
    COMPLEX*16, INTENT(INOUT) :: matrix_A(sizeA, sizeA)
    INTEGER :: i, j

    DO i = 1, sizeA
      DO j = i + 1, sizeA
        matrix_A(j, i) = CONJG(matrix_A(i, j))
      END DO
    END DO

  END SUBROUTINE fill_hermitian

END MODULE hamiltonian
