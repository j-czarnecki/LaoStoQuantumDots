! combinatory_fun.f90 - a unit test suite for combinatory.f90
!
! funit generated this file from combinatory.fun

module combinatory_fun

 use combinatory

 implicit none

 logical :: noAssertFailed

 public :: test_combinatory

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0




 contains





 subroutine test_get_combination

  USE utility
  INTEGER*4, ALLOCATABLE :: Combination(:)
  INTEGER*4, ALLOCATABLE :: Combinations_store(:,:)
  INTEGER*4 :: n_states, k_elems, n_combinations
  INTEGER*4 :: i,j, k, l
  INTEGER*4 :: same_indeces

  n_states = 10
  k_elems = 5
  n_combinations = get_ham_2_size(n_states, k_elems)
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(n_combinations== 252)) then
      print *, " *Assert_Equal failed* in test test_get_combination &
              &[combinatory.fun:25]"
      print *, "  ", "n_combinations (",n_combinations,") is not",  252
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  ALLOCATE(Combination(k_elems))
  ALLOCATE(Combinations_store(n_combinations, k_elems))

  CALL INIT_COMBINATION(Combination, k_elems)
  DO i = 1, n_combinations
    CALL GET_COMBINATION(Combination, n_states, k_elems)
    Combinations_store(i,:) =  Combination(:)
    DO k = 1, k_elems
      DO l = 1, k_elems
        !Combinations should consist of unique single-electron wavefunctions
        IF (k /= l) THEN
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(Combination(k) /= Combination(l))) then
      print *, " *Assert_True failed* in test test_get_combination &
              &[combinatory.fun:38]"
      print *, "  ", "Combination(k) /= Combination(l) is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
        END IF
      END DO
    END DO
  END DO

  !Check whether combinations are unique
  DO i = 1, n_combinations
    DO j = 1, n_combinations
      IF (i /= j) THEN
        same_indeces = 0
        DO k = 1, k_elems
          !Checking whether combinations differ by at least one index of single electron wavefunction
          IF (Combinations_store(i,k) == Combinations_store(j,k)) same_indeces = same_indeces + 1
          !Sanity checks - indeces are bigger than 0 and should be at most equal to number of one electron states
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(Combinations_store(i,k) <= n_states)) then
      print *, " *Assert_True failed* in test test_get_combination &
              &[combinatory.fun:53]"
      print *, "  ", "Combinations_store(i,k) <= n_states is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(Combinations_store(i,k) > 0)) then
      print *, " *Assert_True failed* in test test_get_combination &
              &[combinatory.fun:54]"
      print *, "  ", "Combinations_store(i,k) > 0 is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
        END DO
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(same_indeces < k_elems)) then
      print *, " *Assert_True failed* in test test_get_combination &
              &[combinatory.fun:56]"
      print *, "  ", "same_indeces < k_elems is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
      END IF
    END DO
  END DO


  DEALLOCATE(Combination)
  DEALLOCATE(Combinations_store)

  numTests = numTests + 1

 end subroutine test_get_combination


 subroutine funit_setup

  noAssertFailed = .true.
 end subroutine funit_setup


 subroutine funit_teardown

 end subroutine funit_teardown


 subroutine test_combinatory( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call funit_setup
  call test_get_combination
  call funit_teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_combinatory

end module combinatory_fun
