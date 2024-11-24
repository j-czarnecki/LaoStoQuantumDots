test_suite utility


setup
end setup

teardown

end teardown

test test_get_upper_hermitian_index
  INTEGER*4 :: size, i, j, nn

  size = 3
  nn = 1

  DO i = 1, size
    DO j = i , size
      assert_equal(nn, get_upper_hermitian_index(i, j, size))
      nn = nn + 1
    END DO
  END DO
end test

test test_get_upper_hermitian_index_failure
  INTEGER*4 :: size
  size = 3
  assert_equal(-1, get_upper_hermitian_index(2, 1, size))

  !index bigger that the actual size of matrix
  assert_equal(-1, get_upper_hermitian_index(1, 5, size))
  assert_equal(-1, get_upper_hermitian_index(5, 6, size))
end test

test test_get_opposite_r_index
  INTEGER*4 :: Nx, Ny, norbs, current_index
  INTEGER*4 :: i, j, n, count
  INTEGER*4, ALLOCATABLE :: Natural_ordering(:, :, :)

  !Setup
  Nx = 3
  Ny = 3
  norbs = 2

  ALLOCATE(Natural_ordering(-Nx:Nx, -Ny:Ny, norbs))

  count = 1
  !This corresponds to how my 2D grid with multiple spin-orbitals is being mapped to 1D vector
  DO j = -Ny, Ny
    DO i = -Nx, Nx
      DO n = 1, norbs
        Natural_ordering(i,j,n) = count
        count = count + 1
      END DO
    END DO
  END DO

  !So for each position (x,y) I should get position (-x, -y), but the state should not change.
  DO i = -Nx, Nx
    DO j = -Ny, Ny
      DO n = 1, norbs
        assert_equal(Natural_ordering(-i, -j, n), get_opposite_r_index(Nx, Ny, norbs, Natural_ordering(i, j, n)))
      END DO
    END DO
  END DO

  DEALLOCATE(Natural_ordering)
end test

test test_get_x_from_psi_index
  INTEGER*4 :: Nx, Ny, norbs, current_index
  INTEGER*4 :: i, j, n, count
  REAL*8 :: dx, current_x
  INTEGER*4, ALLOCATABLE :: Natural_ordering(:, :, :)

  !Setup
  Nx = 8
  Ny = 3
  norbs = 2
  dx = 1.0d0

  ALLOCATE(Natural_ordering(-Nx:Nx, -Ny:Ny, norbs))

  count = 1
  !This corresponds to how my 2D grid with multiple spin-orbitals is being mapped to 1D vector
  !Mind the order - column index incremented in outer loop
  DO j = -Ny, Ny
    DO i = -Nx, Nx
      DO n = 1, norbs
        Natural_ordering(i,j,n) = count
        count = count + 1
      END DO
    END DO
  END DO

  !So for each position (x,y) I check the position x.
  DO i = -Nx, Nx
    current_x = i * dx
    DO j = -Ny, Ny
      DO n = 1, norbs
        assert_real_equal(current_x, get_x_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, dx))
      END DO
    END DO
  END DO

  DEALLOCATE(Natural_ordering)
end test


test test_get_y_from_psi_index
  INTEGER*4 :: Nx, Ny, norbs, current_index
  INTEGER*4 :: i, j, n, count
  REAL*8 :: dx, current_y
  INTEGER*4, ALLOCATABLE :: Natural_ordering(:, :, :)

  !Setup
  Nx = 3
  Ny = 8
  norbs = 2
  dx = 1.0d0

  ALLOCATE(Natural_ordering(-Nx:Nx, -Ny:Ny, norbs))

  count = 1
  !This corresponds to how my 2D grid with multiple spin-orbitals is being mapped to 1D vector
  !Mind the order - column index incremented in outer loop
  DO j = -Ny, Ny
    DO i = -Nx, Nx
      DO n = 1, norbs
        Natural_ordering(i,j,n) = count
        count = count + 1
      END DO
    END DO
  END DO

  !So for each position (x,y) I check the position y.
  DO i = -Nx, Nx
    DO j = -Ny, Ny
      current_y = j * dx
      DO n = 1, norbs
        assert_real_equal(current_y, get_y_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, Ny, dx))
      END DO
    END DO
  END DO

  DEALLOCATE(Natural_ordering)
end test

test test_get_slice_from_hermitian_matrix
  COMPLEX*16, ALLOCATABLE :: Hermitian_storage(:,:)
  COMPLEX*16, ALLOCATABLE :: Slice(:)
  INTEGER*4 :: storage_size, slice_size, original_size
  INTEGER*4 :: i,j, n, count

  original_size = 6
  storage_size = (original_size * (original_size + 1))/2
  slice_size = 2

  ALLOCATE(Hermitian_storage(storage_size, slice_size))
  ALLOCATE(Slice(slice_size))

  count = 1
  DO i = 1, original_size
    DO j = i, original_size
      IF (i /= j) THEN
        Hermitian_storage(count,:) = DCMPLX(count*1.0d0, count*1.0d0)
      ELSE
        Hermitian_storage(count,:) = DCMPLX(count*1.0d0, 0) !Diagonal elements should be real
      END IF
      count = count + 1
    END DO
  END DO
  assert_equal(storage_size, count - 1)

  !Verification
  count = 1
  DO i = 1, original_size
    DO j = i, original_size
      IF (i /= j) THEN
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(Slice, Hermitian_storage, slice_size, original_size, storage_size, i, j)
        DO n = 1, slice_size
          assert_real_equal(REAL(Slice(n)), count*1.0d0)
          assert_real_equal(AIMAG(Slice(n)), count*1.0d0)
        END DO

        !Test conjugate element
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(Slice, Hermitian_storage, slice_size, original_size, storage_size, j, i)
        DO n = 1, slice_size
          assert_real_equal(REAL(Slice(n)), count*1.0d0)
          assert_real_equal(AIMAG(Slice(n)), -count*1.0d0)
        END DO
      ELSE
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(Slice, Hermitian_storage, slice_size, original_size, storage_size, i, j)
        DO n = 1, slice_size
          assert_real_equal(REAL(Slice(n)), count*1.0d0)
          assert_real_equal(AIMAG(Slice(n)), 0.0d0)
        END DO

      END IF
      count = count + 1

    END DO
  END DO
  assert_equal(storage_size, count - 1)

  DEALLOCATE(Hermitian_storage)
  DEALLOCATE(Slice)
end test

end test_suite