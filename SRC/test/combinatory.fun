test_suite combinatory


setup

end setup

teardown

end teardown



test test_get_combination
  USE utility
  INTEGER*4, ALLOCATABLE :: Combination(:)
  INTEGER*4, ALLOCATABLE :: Combinations_store(:,:)
  INTEGER*4 :: n_states, k_elems, n_combinations
  INTEGER*4 :: i,j, k, l
  INTEGER*4 :: same_indeces

  n_states = 10
  k_elems = 5
  n_combinations = get_ham_2_size(n_states, k_elems)
  assert_equal(n_combinations, 252)

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
          assert_true(Combination(k) /= Combination(l))
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
          assert_true(Combinations_store(i,k) <= n_states)
          assert_true(Combinations_store(i,k) > 0)
        END DO
        assert_true(same_indeces < k_elems)
      END IF
    END DO
  END DO


  DEALLOCATE(Combination)
  DEALLOCATE(Combinations_store)
end test


end test_suite