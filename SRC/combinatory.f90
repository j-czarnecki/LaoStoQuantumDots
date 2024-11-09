#include "macros_def.f90"
MODULE combinatory
  USE logger
  IMPLICIT NONE
  CONTAINS

SUBROUTINE GET_COMBINATION(Combination, N, k)
  !! Returns set of indices that form a k-element combination from set {1, 2, ..., N}.
  !! Should be initialized with an array conataining first combination with last element diminished by 1.
  !! Example: k=3,then in first iteration Combination = [1, 2, 2] should be passed and will return [1,2,3]
  !! It Could be viewed as traversing an upper triangle of N-dimensional tensor and writing all sets of indices.
  IMPLICIT NONE
  INTEGER*4, INTENT(OUT) :: Combination(k) !! Array containg set of indices forming a k-element combination
  INTEGER*4, INTENT(IN) :: N !! Number of elements in the set
  INTEGER*4, INTENT(IN) :: k !! Number of elements in combination

  INTEGER*4 :: i, j, m

  Combination(k) = Combination(k) + 1

  j = k
  i = 0
  DO  WHILE (j .GT. 1 .AND. Combination(j) .GT. N - i)
    Combination(j - 1) = Combination(j - 1) + 1
    !TODO: This probably could be done only once to reduce some spare operations.
    DO m = j, k
      Combination(m) = Combination(m - 1) + 1
    END DO
    j = j - 1
    i = i + 1
  END DO

END SUBROUTINE

SUBROUTINE INIT_COMBINATION(Combination, k)
  !! Initializes COmbination array so that it could be passed to GET_COBINATION subroutine
  IMPLICIT NONE
  INTEGER*4, INTENT(OUT) :: Combination(k) !! Array containg set of indices forming a k-element combination
                                           !! with last eleent shifted by -1
  INTEGER*4, INTENT(IN) :: k !! Number of elements in combination
  INTEGER*4 :: i

  DO i = 1, k - 1
    Combination(i) = i
  END DO
  Combination(k) = Combination(k - 1)

END SUBROUTINE

SUBROUTINE GET_CHANGED_INDECES(Changed_indeces, Combinations, N_changed_indeces, ham_2_size, k_electrons)
  !! Calculates number of changed indeces between two combinations of k_electrons.
  !! For each pair of combinations it also checks which indeces have been changed.
  IMPLICIT NONE
  INTEGER*4, INTENT(OUT) :: Changed_indeces(ham_2_size, ham_2_size, 2, 2)
  INTEGER*1, INTENT(OUT) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4, INTENT(IN) :: Combinations(ham_2_size, k_electrons)
  INTEGER*4, INTENT(IN) :: ham_2_size, k_electrons

  INTEGER*4 :: i,j,n_changed,k,l
  LOGICAL :: same_index

  DO i = 1, ham_2_size
    DO j = 1, ham_2_size
      n_changed = 0

      !Check how many indeces are changed between combinations specifying row and column
      !For each element of row-combination check heter it exists in column-combination.
      !If not, increment N_changed_indexes by 1 and write index without match to Changed_indeces(i,j,n_chnanged,1)
      DO k = 1, k_electrons
        same_index = .FALSE.
        DO l = 1, k_electrons
          IF(Combinations(i,k) == Combinations(j,l)) THEN
            same_index = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. same_index) THEN
          n_changed = n_changed + 1
          IF (n_changed < 3) THEN
            Changed_indeces(i,j,n_changed,1) = Combinations(i,k)
          END IF
        END IF
      END DO
      N_changed_indeces(i,j) = INT(MIN(n_changed, 3), kind = 1)

      n_changed = 0
      !For each element of column-combination check wheter it exists in row-combination.
      !TODO: Consider whther it could be done in a more efficient way
      DO k = 1, k_electrons
        same_index = .FALSE.
        DO l = 1, k_electrons
          IF(Combinations(j,k) == Combinations(i,l)) THEN
            same_index = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. same_index) THEN
          n_changed = n_changed + 1
          IF (n_changed < 3) THEN
            Changed_indeces(i,j,n_changed,2) = Combinations(j,k)
          END IF
        END IF
      END DO

      !!PRINT*, N_changed_indeces(i,j)
      !PRINT*
      !PRINT*
      !WRITE(*,'(I2)', ADVANCE = 'NO') N_changed_indeces(i,j)
      ! !!WRITE(*,*) 'combination(i) = ', (Combinations(i, n), n = 1, k_electrons)
      ! !!WRITE(*,*) 'combination(j) = ', (Combinations(j, n), n = 1, k_electrons)
      ! WRITE(*, *) 'i = ', i, ' j = ', j, ' N_changed  = ', N_changed_indeces(i,j), ' Unpaired row: ', (Changed_indeces(i,j,k,1), k = 1, 2), ' Unparied column: ', (Changed_indeces(i,j,k,2), k = 1, 2)
      ! !!WRITE(*,*)
    END DO
    !!!WRITE(*,*)
  END DO


END SUBROUTINE GET_CHANGED_INDECES

SUBROUTINE INIT_PREV_ELEMS(N_ham_2_elems_in_prev_rows, N_changed_indeces, ham_2_size, nonzero_ham_2)
  IMPLICIT NONE
  INTEGER*4, INTENT(OUT) :: N_ham_2_elems_in_prev_rows(ham_2_size)
  INTEGER*4, INTENT(IN) :: ham_2_size, nonzero_ham_2
  INTEGER*1, INTENT(IN) :: N_changed_indeces(ham_2_size, ham_2_size)
  INTEGER*4 :: i, j, n

  n = 1
  DO i = 1, ham_2_size
    N_ham_2_elems_in_prev_rows(i) = n
    DO j = i, ham_2_size
      IF (N_changed_indeces(i,j) < 3) n = n + 1
    END DO
  END DO

  !Sanity check whether we calculated number of nonzero elems corectly
  IF (n - 1 < nonzero_ham_2) THEN
    WRITE(log_string,*) 'n - 1 = ', n - 1, ' < nonzero_ham_2 = ', nonzero_ham_2
    LOG_ABNORMAL(log_string)
  ELSE IF (n - 1 > nonzero_ham_2) THEN
    STOP 'ERROR IN INIT_PREV_ELEMS: n > nonzero_ham_2. Bad hamiltonian initialization'
  END IF
END SUBROUTINE

END MODULE combinatory