MODULE broydenV2
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32
CONTAINS
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
! Adapted from QE
!
!-----------------------------------------------------------------------
SUBROUTINE mix_broyden(ndim, deltaout, deltain, alphamix, iter, n_iter, conv)
  !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
  !
  !USE constants_epw, ONLY : DP
  !USE mod_indata, ONLY : nsiter
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: conv      !! If true convergence reache
  INTEGER, INTENT(in) :: ndim      !! Dimension of arrays deltaout, deltain
  INTEGER, INTENT(in) :: iter      !! Current iteration number
  INTEGER, INTENT(in) :: n_iter    !! Number of iterations used in the mixing

  REAL(REAL64), INTENT(in) :: alphamix    !! Mixing factor (0 < alphamix <= 1)
  REAL(REAL64), INTENT(inout) :: deltaout(ndim) !! output Delta at current iteration
  REAL(REAL64), INTENT(inout) :: deltain(ndim)  !! Delta at previous iteration

  !   Here the local variables
  ! max number of iterations used in mixing: n_iter must be .le. maxter
  INTEGER(INT32), PARAMETER :: maxter = 8
  !
  INTEGER(INT32) ::  n, i, j, iwork(maxter), info, iter_used, ipos, inext
  ! work space containing info from previous iterations:
  ! must be kept in memory and saved between calls
  REAL(REAL64), ALLOCATABLE, SAVE :: df(:, :), dv(:, :)
  REAL(REAL64), ALLOCATABLE :: deltainsave(:)
  REAL(REAL64) :: beta(maxter, maxter), gammamix, work(maxter), norm
  REAL(REAL64), EXTERNAL :: DDOT, DNRM2

  ! adjustable PARAMETERs as suggested in the original paper
  REAL(REAL64) wg(maxter), wg0
  DATA wg0/0.01d0/, wg/maxter*1.d0/

  !IF ( iter .lt. 1 ) print *, 'mix_broyden','n_iter is smaller than 1'
  !IF ( n_iter .gt. maxter ) print *, 'mix_broyden','n_iter is too big'
  !IF ( ndim .le. 0 ) print *, 'mix_broyden','ndim .le. 0'

  IF (iter .eq. 1) THEN
    IF (.not. ALLOCATED(df)) ALLOCATE (df(ndim, n_iter))
    IF (.not. ALLOCATED(dv)) ALLOCATE (dv(ndim, n_iter))
  END IF
  IF (conv) THEN
    IF (ALLOCATED(df)) DEALLOCATE (df)
    IF (ALLOCATED(dv)) DEALLOCATE (dv)
    RETURN
  END IF
  IF (.not. ALLOCATED(deltainsave)) ALLOCATE (deltainsave(ndim))
  deltainsave(:) = deltain(:)
  !
  ! iter_used = iter-1  IF iter <= n_iter
  ! iter_used = n_iter  IF iter >  n_iter
  !
  iter_used = min(iter - 1, n_iter)
  !
  ! ipos is the position in which results from the present iteraction
  ! are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
  !
  ipos = iter - 1 - ((iter - 2) / n_iter) * n_iter
  !
  DO n = 1, ndim
    deltaout(n) = deltaout(n) - deltain(n)
  END DO
  !
  IF (iter .gt. 1) THEN
    DO n = 1, ndim
      df(n, ipos) = deltaout(n) - df(n, ipos)
      dv(n, ipos) = deltain(n) - dv(n, ipos)
    END DO
    norm = (DNRM2(ndim, df(1, ipos), 1))**2.d0
    norm = sqrt(norm)
    CALL DSCAL(ndim, 1.d0 / norm, df(1, ipos), 1)
    CALL DSCAL(ndim, 1.d0 / norm, dv(1, ipos), 1)
  END IF
  !
  DO i = 1, iter_used
    DO j = i + 1, iter_used
      beta(i, j) = wg(i) * wg(j) * DDOT(ndim, df(1, j), 1, df(1, i), 1)
    END DO
    beta(i, i) = wg0**2.d0 + wg(i)**2.d0
  END DO
  !
  ! DSYTRF computes the factorization of a real symmetric matrix
  !
  CALL DSYTRF('U', iter_used, beta, maxter, iwork, work, maxter, info)
  !print *, "Broyden factorization", info
  !
  ! DSYTRI computes the inverse of a real symmetric indefinite matrix
  !
  CALL DSYTRI('U', iter_used, beta, maxter, iwork, work, info)
  !print *, "broyden DSYTRI", info
  !
  DO i = 1, iter_used
    DO j = i + 1, iter_used
      beta(j, i) = beta(i, j)
    END DO
  END DO
  !
  DO i = 1, iter_used
    work(i) = DDOT(ndim, df(1, i), 1, deltaout, 1)
  END DO
  !
  DO n = 1, ndim
    deltain(n) = deltain(n) + alphamix * deltaout(n)
  END DO
  !
  DO i = 1, iter_used
    gammamix = 0.d0
    DO j = 1, iter_used
      gammamix = gammamix + beta(j, i) * wg(j) * work(j)
    END DO
    !
    DO n = 1, ndim
      deltain(n) = deltain(n) - wg(i) * gammamix * (alphamix * df(n, i) + dv(n, i))
    END DO
  END DO
  !
  inext = iter - ((iter - 1) / n_iter) * n_iter
  df(:, inext) = deltaout(:)
  dv(:, inext) = deltainsave(:)
  !
  IF (ALLOCATED(deltainsave)) DEALLOCATE (deltainsave)
  !
  RETURN
  !
END SUBROUTINE mix_broyden

END MODULE broydenV2

