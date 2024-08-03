MODULE diagonalize
  USE indata
  USE hamiltonian
  IMPLICIT  NONE
  SAVE

CONTAINS


SUBROUTINE DIAGONALIZE_LAPACK ( psi, ev , nnx, nny, norbitals, nnstate, potential)
  USE indata 
  USE hamiltonian
  IMPLICIT NONE
  !SAVE

  COMPLEX*16, INTENT(OUT):: psi(nstate, -nnx:nnx, -nny:nny, norbitals/2, 2)
  REAL*8, INTENT(OUT):: ev(nstate)
  INTEGER, INTENT(IN) :: nnx, nny, norbitals, nnstate
  REAL*8, INTENT(INOUT):: potential(-nnx:nnx,-nny:nny)
  
  
  INTEGER :: size, ix, iy, is, iorb, nn
  COMPLEX*16, allocatable :: ham(:,:)
  REAL*8 :: probability_norm
  
  !do rozwiazania problemu wlasnego w LAPACKU
  CHARACTER*1 JOBVL, JOBVR
  INTEGER INFO, LWORK, i, j, ii
  COMPLEX*16, ALLOCATABLE :: W(:)
  COMPLEX*16, ALLOCATABLE :: VL(:,:)
  COMPLEX*16, ALLOCATABLE :: VR(:,:)
  COMPLEX*16, ALLOCATABLE :: WORK(:)
  DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
  COMPLEX*16 :: tmp

  JOBVL='N'
  JOBVR='V'

  size=(2*nnx+1)*(2*nny+1)*norbitals
  allocate(ham(size,size))

  ALLOCATE(VL(size, size))
  ALLOCATE(VR(size, size))
  ALLOCATE(W(size))
  LWORK=4*size+1
  ALLOCATE(WORK(LWORK))
  ALLOCATE(RWORK(2*size))



  print*
  print*, "Calculations in progress ..."
  CALL HAMILTONIAN_CREATE(ham,size,nnx,nny,norbitals,potential)

  
  CALL ZGEEV( JOBVL, JOBVR, size, ham, size, W, VL, size, & 
  VR, size, WORK, LWORK, RWORK, INFO )
  

  !sortowanie Energii
  do i=1,size
    do j=1,size-1
    if(REAL(W(j)).gt.REAL(W(j+1))) then
    !energie
    tmp=W(j)
    W(j)=W(j+1)
    W(j+1)=tmp

    !wartosci wlasne
    do ii=1,size
        tmp=VR(ii,j)
        VR(ii,j)=VR(ii,j+1)
        VR(ii,j+1)=tmp
    enddo

    endif
    enddo
  enddo

  do i=1,nstate
    ev(i)=W(i)
  enddo

  !Normalisation of wavefuctions
  DO i = 1, nstate
    probability_norm = SUM(ABS(VR(:,i))**2)
    VR(:,i) = VR(:,i) / sqrt(probability_norm)
  END DO  

  do is=1,nstate
      nn=1
      do iy=-nny,nny
      do ix=-nnx,nnx
        do iorb=1,norbitals/2
           do i=1,2
              psi(is,ix,iy,iorb,i)=VR(nn,is)
              nn=nn+1
           enddo
        enddo
      enddo
      enddo
  enddo

deallocate(ham)
DEALLOCATE(VL)
DEALLOCATE(VR)
DEALLOCATE(W)
DEALLOCATE(WORK)
DEALLOCATE(RWORK)

END SUBROUTINE DIAGONALIZE_LAPACK




SUBROUTINE DIAGONALIZE_ARPACK ( psi, ev , nnx, nny, norbitals, nnstate, potential)
  USE indata 
  USE hamiltonian
  IMPLICIT NONE

  COMPLEX*16, INTENT(OUT):: psi(nstate, -nnx:nnx, -nny:nny, norbitals/2, 2)
  REAL*8, INTENT(OUT):: ev(nstate)
  INTEGER, INTENT(IN) :: nnx, nny, norbitals, nnstate
  REAL*8, INTENT(INOUT):: potential(-nnx:nnx,-nny:nny)
  
  INTEGER :: i,j, flag
  INTEGER :: size, ix, iy, is, iorb, nn
  COMPLEX*16, allocatable :: ham(:,:)
  
  ! CSR args
  INTEGER :: nonzero
  COMPLEX*16, allocatable :: ham_csr(:)
  INTEGER, allocatable :: ja_csr(:)
  INTEGER, allocatable :: ia_csr(:)
  
  !For ARPACK ZNAUPD 
  INTEGER*4 :: ido  !it must be zero on the first call to znaupd
  INTEGER*4 :: dim !dimension of hamiltonian
  INTEGER*4 :: nev  ! number of eigenvalues obtained (?) from znaupd 0 < nev < N-1
  REAL*8 :: accuracy !this is the tolerance of error in Ritz (?) values
  COMPLEX*16, ALLOCATABLE :: resid(:)  !vector of lenght N
  INTEGER*4 :: ncv !number of columns in V matrix (?) wybierane arbitralnie(?)
  COMPLEX*16, ALLOCATABLE :: Arnoldi_vector(:,:)
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
  INTEGER*4, ALLOCATABLE :: select(:) 
  COMPLEX*16, ALLOCATABLE :: E_eigenvalues(:)
  COMPLEX*16 :: sigma 
  COMPLEX*16, ALLOCATABLE :: WORKev(:)
  INTEGER :: n_converged

  INTEGER*4, ALLOCATABLE:: eigen_index(:)
  COMPLEX*16, ALLOCATABLE :: Eigenvectors_normalised(:,:)
  REAL*8 :: probability_norm

  size=(2*nnx+1)*(2*nny+1)*norbitals
  allocate(ham(size,size))
  
  !print*
  !print*, "Calculations in progress ..."
  CALL HAMILTONIAN_CREATE(ham,size,nnx,nny,norbitals,potential)

  !delete elements - the matrix become UPPER
  do i=1,size
    do j=1,i-1
      ham(i,j)=dcmplx(0.0,0.0)
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !MKL ham to CSR format
  nonzero=0
  do i=1,size
    do j=i,size
      if(abs(ham(i,j)).ne.0.0) then 
        nonzero=nonzero+1
      endif
    enddo
  enddo

  allocate(ham_csr(nonzero))
  allocate(ja_csr(nonzero))
  allocate(ia_csr(size+1))
  
  !print*, "Elements:", size*size
  !print*, "Nonzero elements:", nonzero


  nn=1
  do i=1,size
  flag=0
  do j=i,size
    if(abs(ham(i,j)).ne.0.0) then 
      ham_csr(nn)=ham(i,j)
      ja_csr(nn)=j
      if(flag.eq.0) then
        ia_csr(i)=nn
        flag=1
      endif
    nn=nn+1
    endif
  enddo
  enddo
  ia_csr(size+1)=nn

!Declaration of ARPACK values
  rvec = .TRUE.
  ido = 0
  nev = nnstate
  accuracy = -1 !max machine precision
  ncv = 4*(nev+1)
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
  workl_length = 3*ncv**2 + 5*ncv   ! at least 3*ncv**2 + 5*ncv
  info = 0
  dim = size   !Number of points on the grid *4 due to spin-orbit hamiltonian and electron-hole

  ALLOCATE(resid(dim))
  ALLOCATE(Arnoldi_vector(dim,ncv))
  ALLOCATE(workd(3*dim))
  ALLOCATE(workl(workl_length))
  ALLOCATE(rwork(ncv))
  ALLOCATE(select(ncv))
  ALLOCATE(E_eigenvalues(nev+1))
  ALLOCATE(WORKev(2*ncv))
  ALLOCATE(eigen_index(nev))
  ALLOCATE(Eigenvectors_normalised(dim, nev))


  !#################################### ARPACK LOOP #############################################
  DO 
    
    CALL ZNAUPD(ido, 'I', dim, 'SR', nev, accuracy, resid, ncv, Arnoldi_vector, ldv, iparam, ipntr, &
    &  workd, workl, workl_length, rwork, info)
    IF( info .lt. 0 ) PRINT*, "ERROR in znaupd, info = ", info       
    
    IF (ido .eq. 99) EXIT
    IF((ido .ne. -1) .AND. (ido .ne. 1)) THEN
    !prom        PRINT*, "ido = ", ido 
      STOP "ZNAUPD error"
    END IF 
    !Hermitian matrix and vector multiplication
    CALL avmult(dim, nonzero, ia_csr, ja_csr, ham_csr, workd(ipntr(1)), workd(ipntr(2)))

  END DO 
  ! !##################################################################################################
  !PRINT*, "ZNAUPD info: ", info
  !PRINT*,  iparam(5), "eigenvalues found, calling zneupd"
  IF (info .ne. 0) STOP  
  
  !ARPACK eigenvalues extraction
  CALL ZNEUPD(rvec,'A', select, E_eigenvalues, Arnoldi_vector, ldv ,  sigma, WORKev, 'I', dim, &
    & 'SR', nev, accuracy, resid, ncv, Arnoldi_vector, ldv, iparam, ipntr, workd, workl, workl_length, rwork, ierr) 

  ! PRINT*, "ZNEUPD ierr: ", ierr      
  n_converged = iparam(5)
  IF (n_converged .ne. nev) THEN 
    PRINT*, "WARNING: number of eigenvalues converged: ", n_converged, "is different than nstate:", nnstate
    STOP 
  END IF

  !Sorted energies from ARPACK 
  DO i = 1, nnstate
      ev(i) = REAL(E_eigenvalues(i))
  END DO 
    
  !Normalisation of wavefuctions
  DO i = 1, nnstate
    DO j = 1, size
      Eigenvectors_normalised(j,i) = Arnoldi_vector(j,i)
    END DO 
    probability_norm = SUM(ABS(Eigenvectors_normalised(:,i))**2)
    Eigenvectors_normalised(:,i) = Eigenvectors_normalised(:,i) / sqrt(probability_norm)
  END DO  

  do is=1,nstate
      nn=1
      do iy=-nny,nny
      do ix=-nnx,nnx
        do iorb=1,norbitals/2
           do i=1,2
              psi(is,ix,iy,iorb,i)=Eigenvectors_normalised(nn,is)
              nn=nn+1
           enddo
        enddo
      enddo
      enddo
  enddo

  deallocate(ham)
  deallocate(ham_csr)
  deallocate(ja_csr)
  deallocate(ia_csr)

END SUBROUTINE DIAGONALIZE_ARPACK


SUBROUTINE avmult(dimat, maxnonzero, imat, jmat, vmat, vap, wap)
  IMPLICIT NONE
  ! computes  w = M v
  ! for Hermitian matrices M in Compressed Row Storage format
  INTEGER,    INTENT(IN) :: dimat, maxnonzero
  INTEGER,    INTENT(IN) :: imat(dimat+1), jmat(maxnonzero)
  COMPLEX*16, INTENT(IN) :: vmat(maxnonzero)
  COMPLEX*16, INTENT(IN) :: vap(dimat)
  COMPLEX*16, INTENT(OUT) :: wap(dimat)
  
  INTEGER :: ni, nj

  wap(:)= 0.
  !
  DO ni= 1, dimat
    wap(ni) = wap(ni) + vmat(imat(ni)) * vap(jmat(imat(ni)))
    DO nj = imat(ni)+1, imat(ni+1) - 1
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

INTEGER, PARAMETER :: M= 7 , NSTACK= 50

INTEGER :: i, indxt, ir, itemp, j, jstack, k, l
INTEGER :: istack(NSTACK)
REAL*8 :: a

DO j= 1, nelem
  indx(j)= j
END DO
jstack= 0
l= 1
ir= nelem

1 IF (ir-l .LT. M) THEN
  DO j= l+1, ir
    indxt= indx(j)
    a= arr(indxt)
    DO i= j-1, l, -1
      IF (arr(indx(i)) .LE. a) GOTO 2
      indx(i+1)= indx(i)
    END DO
    i= l-1
2   indx(i+1)= indxt
  END DO
  IF (jstack .EQ. 0) RETURN
  ir= istack(jstack)
  l= istack(jstack-1)
  jstack= jstack-2
ELSE
  k= (l+ir)/2
  itemp= indx(k)
  indx(k)= indx(l+1)
  indx(l+1)= itemp
  IF (arr(indx(l)) .GT. arr(indx(ir))) THEN
    itemp= indx(l)
    indx(l)= indx(ir)
    indx(ir)= itemp
  ENDIF
  IF (arr(indx(l+1)) .GT. arr(indx(ir))) THEN
    itemp= indx(l+1)
    indx(l+1)= indx(ir)
    indx(ir)= itemp
  ENDIF
  IF (arr(indx(l)) .GT. arr(indx(l+1))) THEN
    itemp= indx(l)
    indx(l)= indx(l+1)
    indx(l+1)= itemp
  ENDIF
  i= l+1
  j= ir
  indxt= indx(l+1)
  a= arr(indxt)
3 CONTINUE
  i= i+1
  IF (arr(indx(i)) .LT. a) GOTO 3
4 CONTINUE
  j=j-1
  IF (arr(indx(j)) .GT. a) GOTO 4
  IF (j .LT. i) GOTO 5
  itemp= indx(i)
  indx(i)= indx(j)
  indx(j)= itemp
  GOTO 3
5 indx(l+1)= indx(j)
  indx(j)= indxt
  jstack= jstack + 2
  IF (jstack .GT. NSTACK) STOP 'NSTACK too small in indexx'
  IF (ir-i+1 .GE. j-l) THEN
    istack(jstack)= ir
    istack(jstack-1)= i
    ir= j-1
  ELSE
    istack(jstack)= j-1
    istack(jstack-1)= l
    l= i
  ENDIF
ENDIF
GOTO 1

END SUBROUTINE indexx


END MODULE diagonalize
