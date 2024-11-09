#include "macros_def.f90"
MODULE diagonalize
  USE indata
  USE hamiltonian
  USE logger
  IMPLICIT NONE
  SAVE
  CONTAINS

  SUBROUTINE DIAGONALIZE_ARPACK_CRS(ham_csr, ja_csr, ia_csr, nonzero, size, psi, ev, nnstate)
    IMPLICIT NONE

    !COMPLEX*16, INTENT(OUT):: psi(nstate, -nnx:nnx, -nny:nny, norbitals / 2, 2)
    ! CSR args
    INTEGER*4, INTENT(IN) :: nonzero
    COMPLEX*16, INTENT(IN) :: ham_csr(nonzero)
    INTEGER*4, INTENT(IN) :: ja_csr(nonzero)
    INTEGER*4, INTENT(IN) :: ia_csr(size + 1)

    COMPLEX*16, INTENT(OUT) :: psi(size, nnstate)
    REAL*8, INTENT(OUT):: ev(nnstate)
    INTEGER, INTENT(IN) :: nnstate
    !REAL*8, INTENT(INOUT):: potential(-nnx:nnx, -nny:nny)
    !COMPLEX*16, INTENT(IN) :: ham(:, :)
    INTEGER*4, INTENT(IN) :: size
    INTEGER*4, ALLOCATABLE :: state_index(:)
    INTEGER :: i, j
    !COMPLEX*16, ALLOCATABLE :: ham(:, :)

    !For ARPACK ZNAUPD
    INTEGER*4 :: ido  !it must be zero on the first call to znaupd
    INTEGER*4 :: dim !dimension of hamiltonian
    INTEGER*4 :: nev  ! number of eigenvalues obtained (?) from znaupd 0 < nev < N-1
    REAL*8 :: accuracy !this is the tolerance of error in Ritz (?) values
    COMPLEX*16, ALLOCATABLE :: resid(:)  !vector of lenght N
    INTEGER*4 :: ncv !number of columns in V matrix (?) wybierane arbitralnie(?)
    COMPLEX*16, ALLOCATABLE :: Arnoldi_vector(:, :)
    INTEGER*4 :: iparam(11) !parameters fo znaupd
    INTEGER*4 :: ipntr(14)
    COMPLEX*16, ALLOCATABLE :: workd(:)
    COMPLEX*16, ALLOCATABLE :: workl(:)
    INTEGER*4 :: workl_length
    REAL*8, ALLOCATABLE :: rwork(:)
    INTEGER*4 :: info
    INTEGER*4 :: ldv
    INTEGER*4 :: ierr

    !For APRACK ZNEUPD
    LOGICAL*4 :: rvec
    INTEGER*4, ALLOCATABLE :: SELECT(:)
    COMPLEX*16, ALLOCATABLE :: E_eigenvalues(:)
    COMPLEX*16 :: sigma
    COMPLEX*16, ALLOCATABLE :: WORKev(:)
    INTEGER :: n_converged

    INTEGER*4, ALLOCATABLE:: eigen_index(:)
    COMPLEX*16, ALLOCATABLE :: Eigenvectors_normalised(:, :)
    REAL*8 :: probability_norm

    WRITE(log_string, *) 'Diagonalizing Hamiltonian'
    LOG_INFO(log_string)

!Declaration of ARPACK values
    rvec = .TRUE.
    ido = 0
    nev = nnstate
    accuracy = -1 !max machine precision
    ncv = 4 * (nev + 1)
    ldv = size
    iparam(1) = 1 !exact shift strategy
    iparam(2) = 0
    iparam(3) = 400   !number of ARPACK loop iterations allowed
    iparam(4) = 1      !number of blocks used in reccurency (only 1 works, according to documentation)
    ! 5 exit parameter- number of converged eigenvalues, 6 not referenced
    iparam(7) = 1   ! 1 for bulit-in shift-invert strategy (?), 3 for self declared shift-invert (?)
    ! 8 if shift-invert mode is declared separately
    !9 exit parameter
    ipntr = 0
    workl_length = 3 * ncv**2 + 5 * ncv   ! at least 3*ncv**2 + 5*ncv
    info = 0
    dim = size   !Number of points on the grid *4 due to spin-orbit hamiltonian and electron-hole

    !PRINT*, 'NCV = ', ncv, " NEV = ", nev, "size = ", size

    ALLOCATE (resid(dim))
    ALLOCATE (Arnoldi_vector(dim, ncv))
    ALLOCATE (workd(3 * dim))
    ALLOCATE (workl(workl_length))
    ALLOCATE (rwork(ncv))
    ALLOCATE (SELECT(ncv))
    ALLOCATE (E_eigenvalues(nev + 1))
    ALLOCATE (WORKev(2 * ncv))
    ALLOCATE (eigen_index(nev))
    ALLOCATE (Eigenvectors_normalised(dim, nev))
    ALLOCATE (state_index(nev))

    !#################################### ARPACK LOOP #############################################
    DO

      CALL ZNAUPD(ido, 'I', dim, 'SR', nev, accuracy, resid, ncv, Arnoldi_vector, ldv, iparam, ipntr, &
      &  workd, workl, workl_length, rwork, info)
      IF (info .LT. 0) PRINT *, "ERROR in znaupd, info = ", info

      IF (ido .EQ. 99) EXIT
      IF ((ido .NE. -1) .AND. (ido .NE. 1)) THEN
        !prom        !PRINT*, "ido = ", ido
        STOP "ZNAUPD error"
      END IF
      !Hermitian matrix and vector multiplication
      CALL avmult(dim, nonzero, ia_csr, ja_csr, ham_csr, workd(ipntr(1)), workd(ipntr(2)))

    END DO
    ! !##################################################################################################
    !!PRINT*, "ZNAUPD info: ", info
    !!PRINT*,  iparam(5), "eigenvalues found, calling zneupd"
    IF (info .NE. 0) STOP

    !ARPACK eigenvalues extraction
    CALL ZNEUPD(rvec, 'A', SELECT, E_eigenvalues, Arnoldi_vector, ldv, sigma, WORKev, 'I', dim, &
      & 'SR', nev, accuracy, resid, ncv, Arnoldi_vector, ldv, iparam, ipntr, workd, workl, workl_length, rwork, ierr)

    ! !PRINT*, "ZNEUPD ierr: ", ierr
    n_converged = iparam(5)
    IF (n_converged .NE. nev) THEN
      PRINT *, "WARNING: number of eigenvalues converged: ", n_converged, "is different than nstate:", nnstate
      STOP
    END IF

    CALL indexx(n_converged, REAL(E_eigenvalues), state_index)

    !Sorted energies from ARPACK
    DO i = 1, nnstate
      ev(i) = REAL(E_eigenvalues(state_index(i)))
    END DO

    !Julian: consider whether we need to normalize the wavefuctions
    !This may be spare - consult with docummentation
    !Normalisation of wavefuctions
    DO i = 1, nnstate
      DO j = 1, size
        Eigenvectors_normalised(j, i) = Arnoldi_vector(j, state_index(i))
      END DO
      probability_norm = SUM(ABS(Eigenvectors_normalised(:, i))**2)
      !Eigenvectors_normalised(:, i) = Eigenvectors_normalised(:, i) / SQRT(probability_norm)
      psi(:, i) =  Eigenvectors_normalised(:, i) / SQRT(probability_norm)
    END DO


    DEALLOCATE (state_index)
    DEALLOCATE (resid)
    DEALLOCATE (Arnoldi_vector)
    DEALLOCATE (workd)
    DEALLOCATE (workl)
    DEALLOCATE (rwork)
    DEALLOCATE (SELECT)
    DEALLOCATE (E_eigenvalues)
    DEALLOCATE (WORKev)
    DEALLOCATE (eigen_index)
    DEALLOCATE (Eigenvectors_normalised)
  END SUBROUTINE DIAGONALIZE_ARPACK_CRS



  SUBROUTINE DIAGONALIZE_ARPACK(ham, size, psi, ev, nnstate)
    USE indata
    USE hamiltonian
    IMPLICIT NONE

    !COMPLEX*16, INTENT(OUT):: psi(nstate, -nnx:nnx, -nny:nny, norbitals / 2, 2)
    COMPLEX*16, INTENT(OUT) :: psi(size, nnstate)
    REAL*8, INTENT(OUT):: ev(nnstate)
    INTEGER, INTENT(IN) :: nnstate
    !REAL*8, INTENT(INOUT):: potential(-nnx:nnx, -nny:nny)
    COMPLEX*16, INTENT(IN) :: ham(:, :)
    INTEGER*4, INTENT(IN) :: size
    INTEGER*4, ALLOCATABLE :: state_index(:)
    INTEGER :: i, j, flag
    INTEGER :: nn !, size, iorb, is, ix, iy
    !COMPLEX*16, ALLOCATABLE :: ham(:, :)

    ! CSR args
    INTEGER :: nonzero
    COMPLEX*16, ALLOCATABLE :: ham_csr(:)
    INTEGER, ALLOCATABLE :: ja_csr(:)
    INTEGER, ALLOCATABLE :: ia_csr(:)

    !For ARPACK ZNAUPD
    INTEGER*4 :: ido  !it must be zero on the first call to znaupd
    INTEGER*4 :: dim !dimension of hamiltonian
    INTEGER*4 :: nev  ! number of eigenvalues obtained (?) from znaupd 0 < nev < N-1
    REAL*8 :: accuracy !this is the tolerance of error in Ritz (?) values
    COMPLEX*16, ALLOCATABLE :: resid(:)  !vector of lenght N
    INTEGER*4 :: ncv !number of columns in V matrix (?) wybierane arbitralnie(?)
    COMPLEX*16, ALLOCATABLE :: Arnoldi_vector(:, :)
    INTEGER*4 :: iparam(11) !parameters fo znaupd
    INTEGER*4 :: ipntr(14)
    COMPLEX*16, ALLOCATABLE :: workd(:)
    COMPLEX*16, ALLOCATABLE :: workl(:)
    INTEGER*4 :: workl_length
    REAL*8, ALLOCATABLE :: rwork(:)
    INTEGER*4 :: info
    INTEGER*4 :: ldv
    INTEGER*4 :: ierr

    !For APRACK ZNEUPD
    LOGICAL*4 :: rvec
    INTEGER*4, ALLOCATABLE :: SELECT(:)
    COMPLEX*16, ALLOCATABLE :: E_eigenvalues(:)
    COMPLEX*16 :: sigma
    COMPLEX*16, ALLOCATABLE :: WORKev(:)
    INTEGER :: n_converged

    INTEGER*4, ALLOCATABLE:: eigen_index(:)
    COMPLEX*16, ALLOCATABLE :: Eigenvectors_normalised(:, :)
    REAL*8 :: probability_norm

    !size = (2 * nnx + 1) * (2 * nny + 1) * norbitals
    !ALLOCATE (ham(size, size))

    !print*
    !!PRINT*, "Calculations in progress ..."
    !CALL HAMILTONIAN_CREATE(ham, size, nnx, nny, norbitals, potential)

    !delete elements - the matrix become UPPER
    !Julian: this is not needed, we can iterate over upper traingle of matrix
    ! and effect will be the same
    ! DO i = 1, size
    !   DO j = 1, i - 1
    !     ham(i, j) = dcmplx(0.0, 0.0)
    !   END DO
    ! END DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !MKL ham to CSR format
    nonzero = 0
    DO i = 1, size
      DO j = i, size
        IF (ABS(ham(i, j)) .NE. 0.0d0) THEN
          nonzero = nonzero + 1
        END IF
      END DO
    END DO

    !PRINT*, "Trying to allocte compressed storage row hamiltonian..."
    ALLOCATE (ham_csr(nonzero))
    ALLOCATE (ja_csr(nonzero))
    ALLOCATE (ia_csr(size + 1))

    !PRINT*, "Elements:", size*size
    !PRINT*, "Nonzero elements:", nonzero

    nn = 1
    DO i = 1, size
      flag = 0
      DO j = i, size
        IF (ABS(ham(i, j)) .NE. 0.0) THEN
          ham_csr(nn) = ham(i, j)
          ja_csr(nn) = j
          IF (flag .EQ. 0) THEN
            ia_csr(i) = nn
            flag = 1
          END IF
          nn = nn + 1
        END IF
      END DO
    END DO
    ia_csr(size + 1) = nn

!Declaration of ARPACK values
    rvec = .TRUE.
    ido = 0
    nev = nnstate
    accuracy = -1 !max machine precision
    ncv = 4 * (nev + 1)
    ldv = size
    iparam(1) = 1 !exact shift strategy
    iparam(2) = 0
    iparam(3) = 400   !number of ARPACK loop iterations allowed
    iparam(4) = 1      !number of blocks used in reccurency (only 1 works, according to documentation)
    ! 5 exit parameter- number of converged eigenvalues, 6 not referenced
    iparam(7) = 1   ! 1 for bulit-in shift-invert strategy (?), 3 for self declared shift-invert (?)
    ! 8 if shift-invert mode is declared separately
    !9 exit parameter
    ipntr = 0
    workl_length = 3 * ncv**2 + 5 * ncv   ! at least 3*ncv**2 + 5*ncv
    info = 0
    dim = size   !Number of points on the grid *4 due to spin-orbit hamiltonian and electron-hole

    !PRINT*, 'NCV = ', ncv, " NEV = ", nev, "size = ", size

    ALLOCATE (resid(dim))
    ALLOCATE (Arnoldi_vector(dim, ncv))
    ALLOCATE (workd(3 * dim))
    ALLOCATE (workl(workl_length))
    ALLOCATE (rwork(ncv))
    ALLOCATE (SELECT(ncv))
    ALLOCATE (E_eigenvalues(nev + 1))
    ALLOCATE (WORKev(2 * ncv))
    ALLOCATE (eigen_index(nev))
    ALLOCATE (Eigenvectors_normalised(dim, nev))
    ALLOCATE (state_index(nev))

    !#################################### ARPACK LOOP #############################################
    DO

      CALL ZNAUPD(ido, 'I', dim, 'SR', nev, accuracy, resid, ncv, Arnoldi_vector, ldv, iparam, ipntr, &
      &  workd, workl, workl_length, rwork, info)
      IF (info .LT. 0) PRINT *, "ERROR in znaupd, info = ", info

      IF (ido .EQ. 99) EXIT
      IF ((ido .NE. -1) .AND. (ido .NE. 1)) THEN
        !prom        !PRINT*, "ido = ", ido
        STOP "ZNAUPD error"
      END IF
      !Hermitian matrix and vector multiplication
      CALL avmult(dim, nonzero, ia_csr, ja_csr, ham_csr, workd(ipntr(1)), workd(ipntr(2)))

    END DO
    ! !##################################################################################################
    !!PRINT*, "ZNAUPD info: ", info
    !!PRINT*,  iparam(5), "eigenvalues found, calling zneupd"
    IF (info .NE. 0) STOP

    !ARPACK eigenvalues extraction
    CALL ZNEUPD(rvec, 'A', SELECT, E_eigenvalues, Arnoldi_vector, ldv, sigma, WORKev, 'I', dim, &
      & 'SR', nev, accuracy, resid, ncv, Arnoldi_vector, ldv, iparam, ipntr, workd, workl, workl_length, rwork, ierr)

    ! !PRINT*, "ZNEUPD ierr: ", ierr
    n_converged = iparam(5)
    IF (n_converged .NE. nev) THEN
      PRINT *, "WARNING: number of eigenvalues converged: ", n_converged, "is different than nstate:", nnstate
      STOP
    END IF

    CALL indexx(n_converged, REAL(E_eigenvalues), state_index)

    !Sorted energies from ARPACK
    DO i = 1, nnstate
      ev(i) = REAL(E_eigenvalues(state_index(i)))
    END DO

    !Julian: consider whether we need to normalize the wavefuctions
    !This may be spare - consult with docummentation
    !Normalisation of wavefuctions
    DO i = 1, nnstate
      DO j = 1, size
        Eigenvectors_normalised(j, i) = Arnoldi_vector(j, state_index(i))
      END DO
      probability_norm = SUM(ABS(Eigenvectors_normalised(:, i))**2)
      !Eigenvectors_normalised(:, i) = Eigenvectors_normalised(:, i) / SQRT(probability_norm)
      psi(:, i) =  Eigenvectors_normalised(:, i) / SQRT(probability_norm)
    END DO

    ! DO is = 1, nstate
    !   nn = 1
    !   DO iy = -nny, nny
    !   DO ix = -nnx, nnx
    !     DO iorb = 1, norbitals / 2
    !       DO i = 1, 2
    !         psi(is, ix, iy, iorb, i) = Eigenvectors_normalised(nn, is)
    !         nn = nn + 1
    !       END DO
    !     END DO
    !   END DO
    !   END DO
    ! END DO

    !DEALLOCATE (ham)
    DEALLOCATE (ham_csr)
    DEALLOCATE (ja_csr)
    DEALLOCATE (ia_csr)
    DEALLOCATE (state_index)

  END SUBROUTINE DIAGONALIZE_ARPACK

  SUBROUTINE avmult(dimat, maxnonzero, imat, jmat, vmat, vap, wap)
    IMPLICIT NONE
    ! computes  w = M v
    ! for Hermitian matrices M in Compressed Row Storage format
    INTEGER, INTENT(IN) :: dimat, maxnonzero
    INTEGER, INTENT(IN) :: imat(dimat + 1), jmat(maxnonzero)
    COMPLEX*16, INTENT(IN) :: vmat(maxnonzero)
    COMPLEX*16, INTENT(IN) :: vap(dimat)
    COMPLEX*16, INTENT(OUT) :: wap(dimat)

    INTEGER :: ni, nj

    wap(:) = 0.
    !
    DO ni = 1, dimat
      wap(ni) = wap(ni) + vmat(imat(ni)) * vap(jmat(imat(ni)))
      DO nj = imat(ni) + 1, imat(ni + 1) - 1
        wap(ni) = wap(ni) + vmat(nj) * vap(jmat(nj))
        wap(jmat(nj)) = wap(jmat(nj)) + CONJG(vmat(nj)) * vap(ni)
      END DO
    END DO

  END SUBROUTINE avmult

  SUBROUTINE indexx(nelem, arr, indx)
! Indexes an array arr(1:n), i.e., outputs the array indx(1:n)
! such that arr(indx(j)) is in ascending order for j = 1, 2, . . . , N .
! The input quantities n and arr are not changed.
! from NR cap 8
! code taken from
!   www.phys.uu.nl/DU/num_recipes/fortran.208/f77/recipes/indexx.for
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nelem      ! elements in arr and indx
    REAL*8, INTENT(IN) :: arr(nelem)
    INTEGER, INTENT(OUT) :: indx(nelem)

    INTEGER, PARAMETER :: M = 7, NSTACK = 50

    INTEGER :: i, indxt, ir, itemp, j, jstack, k, l
    INTEGER :: istack(NSTACK)
    REAL*8 :: a

    DO j = 1, nelem
      indx(j) = j
    END DO
    jstack = 0
    l = 1
    ir = nelem

1   IF (ir - l .LT. M) THEN
    DO j = l + 1, ir
      indxt = indx(j)
      a = arr(indxt)
      DO i = j - 1, l, -1
        IF (arr(indx(i)) .LE. a) GOTO 2
        indx(i + 1) = indx(i)
      END DO
      i = l - 1
2     indx(i + 1) = indxt
    END DO
    IF (jstack .EQ. 0) RETURN
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
    ELSE
    k = (l + ir) / 2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    IF (arr(indx(l)) .GT. arr(indx(ir))) THEN
      itemp = indx(l)
      indx(l) = indx(ir)
      indx(ir) = itemp
    END IF
    IF (arr(indx(l + 1)) .GT. arr(indx(ir))) THEN
      itemp = indx(l + 1)
      indx(l + 1) = indx(ir)
      indx(ir) = itemp
    END IF
    IF (arr(indx(l)) .GT. arr(indx(l + 1))) THEN
      itemp = indx(l)
      indx(l) = indx(l + 1)
      indx(l + 1) = itemp
    END IF
    i = l + 1
    j = ir
    indxt = indx(l + 1)
    a = arr(indxt)
3   CONTINUE
    i = i + 1
    IF (arr(indx(i)) .LT. a) GOTO 3
4   CONTINUE
    j = j - 1
    IF (arr(indx(j)) .GT. a) GOTO 4
    IF (j .LT. i) GOTO 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    GOTO 3
5   indx(l + 1) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    IF (jstack .GT. NSTACK) STOP 'NSTACK too small in indexx'
    IF (ir - i + 1 .GE. j - l) THEN
      istack(jstack) = ir
      istack(jstack - 1) = i
      ir = j - 1
    ELSE
      istack(jstack) = j - 1
      istack(jstack - 1) = l
      l = i
    END IF
    END IF
    GOTO 1

  END SUBROUTINE indexx

END MODULE diagonalize
