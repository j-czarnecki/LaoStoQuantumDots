! utility_fun.f90 - a unit test suite for utility.f90
!
! funit generated this file from utility.fun

module utility_fun

 use utility

 implicit none

 logical :: noAssertFailed

 public :: test_utility

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0




 contains



 subroutine test_get_upper_hermitian_index

  INTEGER*4 :: size, i, j, nn

  size = 3
  nn = 1

  DO i = 1, size
    DO j = i , size
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(nn== get_upper_hermitian_index(i, j, size))) then
      print *, " *Assert_Equal failed* in test test_get_upper_hermitian_index &
              &[utility.fun:19]"
      print *, "  ", "nn (",nn,") is not",  get_upper_hermitian_index(i, j, size)
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
      nn = nn + 1
    END DO
  END DO

  numTests = numTests + 1

 end subroutine test_get_upper_hermitian_index


 subroutine test_get_upper_hermitian_index_failure

  INTEGER*4 :: size
  size = 3
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(-1== get_upper_hermitian_index(2, 1, size))) then
      print *, " *Assert_Equal failed* in test test_get_upper_hermitian_index_failure &
              &[utility.fun:28]"
      print *, "  ", "-1 (",-1,") is not",  get_upper_hermitian_index(2, 1, size)
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  !index bigger that the actual size of matrix
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(-1== get_upper_hermitian_index(1, 5, size))) then
      print *, " *Assert_Equal failed* in test test_get_upper_hermitian_index_failure &
              &[utility.fun:31]"
      print *, "  ", "-1 (",-1,") is not",  get_upper_hermitian_index(1, 5, size)
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(-1== get_upper_hermitian_index(5, 6, size))) then
      print *, " *Assert_Equal failed* in test test_get_upper_hermitian_index_failure &
              &[utility.fun:32]"
      print *, "  ", "-1 (",-1,") is not",  get_upper_hermitian_index(5, 6, size)
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine test_get_upper_hermitian_index_failure


 subroutine test_get_opposite_r_index

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
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(Natural_ordering(-i, -j, n)== get_opposite_r_index(Nx, Ny, norbs, Natural_ordering(i, j, n)))) then
      print *, " *Assert_Equal failed* in test test_get_opposite_r_index &
              &[utility.fun:62]"
      print *, "  ", "Natural_ordering(-i, -j, n) (",Natural_ordering(-i, -j, n),") is not",  get_opposite_r_index(Nx, Ny, norbs, Natural_ordering(i, j, n))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
      END DO
    END DO
  END DO

  DEALLOCATE(Natural_ordering)

  numTests = numTests + 1

 end subroutine test_get_opposite_r_index


 subroutine test_get_x_from_psi_index

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
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (current_x &
        +2*spacing(real(current_x)) ) &
        .ge. &
        (get_x_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, dx)) &
            .and. &
     (current_x &
      -2*spacing(real(current_x)) ) &
      .le. &
       (get_x_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, dx)) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_x_from_psi_index &
              &[utility.fun:101]"
      print *, "  ", "get_x_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, dx) (", &
 get_x_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, dx), &
  ") is not", &
 current_x,&
 "within", &
  2*spacing(real(current_x))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
      END DO
    END DO
  END DO

  DEALLOCATE(Natural_ordering)

  numTests = numTests + 1

 end subroutine test_get_x_from_psi_index



 subroutine test_get_y_from_psi_index

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
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (current_y &
        +2*spacing(real(current_y)) ) &
        .ge. &
        (get_y_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, Ny, dx)) &
            .and. &
     (current_y &
      -2*spacing(real(current_y)) ) &
      .le. &
       (get_y_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, Ny, dx)) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_y_from_psi_index &
              &[utility.fun:141]"
      print *, "  ", "get_y_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, Ny, dx) (", &
 get_y_from_psi_index(Natural_ordering(i,j,n), norbs, Nx, Ny, dx), &
  ") is not", &
 current_y,&
 "within", &
  2*spacing(real(current_y))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
      END DO
    END DO
  END DO

  DEALLOCATE(Natural_ordering)

  numTests = numTests + 1

 end subroutine test_get_y_from_psi_index


 subroutine test_get_slice_from_hermitian_matrix

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
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(storage_size== count - 1)) then
      print *, " *Assert_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:173]"
      print *, "  ", "storage_size (",storage_size,") is not",  count - 1
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  !Verification
  count = 1
  DO i = 1, original_size
    DO j = i, original_size
      IF (i /= j) THEN
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(Slice, Hermitian_storage, slice_size, original_size, storage_size, i, j)
        DO n = 1, slice_size
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (REAL(Slice(n)) &
        +2*spacing(real(REAL(Slice(n)))) ) &
        .ge. &
        (count*1.0d0) &
            .and. &
     (REAL(Slice(n)) &
      -2*spacing(real(REAL(Slice(n)))) ) &
      .le. &
       (count*1.0d0) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:182]"
      print *, "  ", "count*1.0d0 (", &
 count*1.0d0, &
  ") is not", &
 REAL(Slice(n)),&
 "within", &
  2*spacing(real(REAL(Slice(n))))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (AIMAG(Slice(n)) &
        +2*spacing(real(AIMAG(Slice(n)))) ) &
        .ge. &
        (count*1.0d0) &
            .and. &
     (AIMAG(Slice(n)) &
      -2*spacing(real(AIMAG(Slice(n)))) ) &
      .le. &
       (count*1.0d0) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:183]"
      print *, "  ", "count*1.0d0 (", &
 count*1.0d0, &
  ") is not", &
 AIMAG(Slice(n)),&
 "within", &
  2*spacing(real(AIMAG(Slice(n))))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
        END DO

        !Test conjugate element
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(Slice, Hermitian_storage, slice_size, original_size, storage_size, j, i)
        DO n = 1, slice_size
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (REAL(Slice(n)) &
        +2*spacing(real(REAL(Slice(n)))) ) &
        .ge. &
        (count*1.0d0) &
            .and. &
     (REAL(Slice(n)) &
      -2*spacing(real(REAL(Slice(n)))) ) &
      .le. &
       (count*1.0d0) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:189]"
      print *, "  ", "count*1.0d0 (", &
 count*1.0d0, &
  ") is not", &
 REAL(Slice(n)),&
 "within", &
  2*spacing(real(REAL(Slice(n))))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (AIMAG(Slice(n)) &
        +2*spacing(real(AIMAG(Slice(n)))) ) &
        .ge. &
        (-count*1.0d0) &
            .and. &
     (AIMAG(Slice(n)) &
      -2*spacing(real(AIMAG(Slice(n)))) ) &
      .le. &
       (-count*1.0d0) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:190]"
      print *, "  ", "-count*1.0d0 (", &
 -count*1.0d0, &
  ") is not", &
 AIMAG(Slice(n)),&
 "within", &
  2*spacing(real(AIMAG(Slice(n))))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
        END DO
      ELSE
        CALL GET_SLICE_FROM_HERMITIAN_MATRIX(Slice, Hermitian_storage, slice_size, original_size, storage_size, i, j)
        DO n = 1, slice_size
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (REAL(Slice(n)) &
        +2*spacing(real(REAL(Slice(n)))) ) &
        .ge. &
        (count*1.0d0) &
            .and. &
     (REAL(Slice(n)) &
      -2*spacing(real(REAL(Slice(n)))) ) &
      .le. &
       (count*1.0d0) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:195]"
      print *, "  ", "count*1.0d0 (", &
 count*1.0d0, &
  ") is not", &
 REAL(Slice(n)),&
 "within", &
  2*spacing(real(REAL(Slice(n))))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( (AIMAG(Slice(n)) &
        +2*spacing(real(AIMAG(Slice(n)))) ) &
        .ge. &
        (0.0d0) &
            .and. &
     (AIMAG(Slice(n)) &
      -2*spacing(real(AIMAG(Slice(n)))) ) &
      .le. &
       (0.0d0) )) then
      print *, " *Assert_Real_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:196]"
      print *, "  ", "0.0d0 (", &
 0.0d0, &
  ") is not", &
 AIMAG(Slice(n)),&
 "within", &
  2*spacing(real(AIMAG(Slice(n))))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
        END DO

      END IF
      count = count + 1

    END DO
  END DO
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(storage_size== count - 1)) then
      print *, " *Assert_Equal failed* in test test_get_slice_from_hermitian_matrix &
              &[utility.fun:204]"
      print *, "  ", "storage_size (",storage_size,") is not",  count - 1
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  DEALLOCATE(Hermitian_storage)
  DEALLOCATE(Slice)

  numTests = numTests + 1

 end subroutine test_get_slice_from_hermitian_matrix


 subroutine funit_setup
  noAssertFailed = .true.
 end subroutine funit_setup


 subroutine funit_teardown

 end subroutine funit_teardown


 subroutine test_utility( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call funit_setup
  call test_get_upper_hermitian_index
  call funit_teardown

  call funit_setup
  call test_get_upper_hermitian_index_failure
  call funit_teardown

  call funit_setup
  call test_get_opposite_r_index
  call funit_teardown

  call funit_setup
  call test_get_x_from_psi_index
  call funit_teardown

  call funit_setup
  call test_get_y_from_psi_index
  call funit_teardown

  call funit_setup
  call test_get_slice_from_hermitian_matrix
  call funit_teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_utility

end module utility_fun
