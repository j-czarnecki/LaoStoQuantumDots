MODULE hamiltonian
  USE indata
  IMPLICIT NONE

CONTAINS

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

    INTEGER :: i, j, nn
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
    DO j = -nnx, nnx
      DO i = -nny, nny
        ordering(i, j) = nn
        nn = nn + 1
      END DO
    END DO

    h_hopping_y = TRANSPOSE(RESHAPE((/dcmplx(-tl, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(-tl, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(drso / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0), dcmplx(0.0, 0.0), &
                                      dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) &
                                      /), SHAPE(h_hopping_y)))

    h_onsite_B = TRANSPOSE(RESHAPE((/dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                     dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) &
                                     /), SHAPE(h_onsite_B)))

    h_hopping_diagonal = TRANSPOSE(RESHAPE((/dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                             dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) &
                                             /), SHAPE(h_hopping_diagonal)))

    h_hopping_diagonal_m = TRANSPOSE(RESHAPE((/dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                               dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-td / 2.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0) &
                                               /), SHAPE(h_hopping_diagonal)))

    ham = dcmplx(0.0, 0.0)
    DO i = -nnx, nnx
      DO j = -nny, nny

        pot = potential(i, j)
        !definition of submatrices to create a total Hamiltonian
        h_onsite = TRANSPOSE(RESHAPE((/dcmplx(4.0 * tl - dE + pot + 0.5 * g * mub * Bz, 0.0), dcmplx(0.5 * g * mub * Bx, -0.5 * g * mub * By), dcmplx(0.0, mub * Bx), dcmplx(0.0, dso / 3.0), dcmplx(0.0, -mub * By), dcmplx(-dso / 3.0, 0.0), &
                                       dcmplx(0.5 * g * mub * Bx, 0.5 * g * mub * By), dcmplx(4.0 * tl - dE + pot - 0.5 * g * mub * Bz, 0.0), dcmplx(0.0, dso / 3.0), dcmplx(0.0, mub * Bx), dcmplx(dso / 3.0, 0.0), dcmplx(0.0, -mub * By), &
                                       dcmplx(0.0, -mub * Bx), dcmplx(0.0, -dso / 3.0), dcmplx(2.0 * tl + 2.0 * th + pot + 0.5 * g * mub * Bz, 0.0), dcmplx(0.5 * g * mub * Bx, -0.5 * g * mub * By), dcmplx(0.0, dso / 3.0 + mub * Bz), dcmplx(0.0, 0.0), &
                                       dcmplx(0.0, -dso / 3.0), dcmplx(0.0, -mub * Bx), dcmplx(0.5 * g * mub * Bx, 0.5 * g * mub * By), dcmplx(2.0 * tl + 2.0 * th + pot - 0.5 * g * mub * Bz, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, -dso / 3.0 + mub * Bz), &
                                       dcmplx(0.0, mub * By), dcmplx(dso / 3.0, 0.0), dcmplx(0.0, -dso / 3.0 - mub * Bz), dcmplx(0.0, 0.0), dcmplx(2.0 * tl + 2.0 * th + pot + 0.5 * g * mub * Bz, 0.0), dcmplx(0.5 * g * mub * Bx, -0.5 * g * mub * By), &
                                       dcmplx(-dso / 3.0, 0.0), dcmplx(0.0, mub * By), dcmplx(0.0, 0.0), dcmplx(0.0, dso / 3.0 - mub * Bz), dcmplx(0.5 * g * mub * Bx, 0.5 * g * mub * By), dcmplx(2.0 * tl + 2.0 * th + pot - 0.5 * g * mub * Bz, 0.0) &
                                       /), SHAPE(h_onsite)))

        yi = j * dx
        orb = EXP(-0.5 * one * Bz * (dx) * (2 * yi))
        h_hopping_x = TRANSPOSE(RESHAPE((/dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0) * orb, dcmplx(0.0, 0.0), &
                                          dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-drso / 2.0, 0.0) * orb, &
                                          dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                          dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-tl, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), &
                                          dcmplx(drso / 2.0, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0) * orb, dcmplx(0.0, 0.0), &
                                          dcmplx(0.0, 0.0), dcmplx(drso / 2.0, 0.0) * orb, dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(0.0, 0.0), dcmplx(-th, 0.0) * orb &
                                          /), SHAPE(h_hopping_x)))

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
          CALL put_matrix(ordering(i, j), ordering(i - 1, j + 1), ham, h_hopping_diagonal, size, norbitals)

        ELSE IF (j .EQ. nny .AND. i .NE. -nnx .AND. i .NE. nnx) THEN
          !top facet
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)

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
          CALL put_matrix(ordering(i, j), ordering(i - 1, j + 1), ham, h_hopping_diagonal, size, norbitals)

        ELSE IF (i .EQ. -nnx .AND. j .EQ. nny) THEN
          !left-top corner
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)

        ELSE IF (i .EQ. nnx .AND. j .EQ. nny) THEN
          !right-top corner
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)

        ELSE
          CALL put_matrix(ordering(i, j), ordering(i, j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j), ham, h_hopping_x, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i, j + 1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i + 1, j + 1), ham, h_hopping_diagonal_m, size, norbitals)
          CALL put_matrix(ordering(i, j), ordering(i - 1, j + 1), ham, h_hopping_diagonal, size, norbitals)

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
