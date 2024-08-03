MODULE hamiltonian
  USE indata
  IMPLICIT  NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!create full matrix containing Hamiltoanian !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HAMILTONIAN_CREATE(ham, size, nnx, nny, norbitals, potential)
IMPLICIT NONE
SAVE
    INTEGER, INTENT(IN) :: norbitals, size
    INTEGER, INTENT(IN) :: nnx, nny
    COMPLEX*16, INTENT(OUT) :: ham(size,size)
    REAL*8, INTENT(OUT) :: potential(-nnx:nnx,-nny:nny)

    COMPLEX*16, ALLOCATABLE :: h_onsite(:,:)
    COMPLEX*16, ALLOCATABLE :: h_onsite_B(:,:)
    COMPLEX*16, ALLOCATABLE :: h_hopping_x(:,:)
    COMPLEX*16, ALLOCATABLE :: h_hopping_y(:,:)
    COMPLEX*16, ALLOCATABLE :: h_hopping_diagonal(:,:)
    COMPLEX*16, ALLOCATABLE :: h_hopping_diagonal_m(:,:)

    integer :: i, j, nn
    INTEGER, ALLOCATABLE :: ordering(:,:)
    COMPLEX*16 :: orb, one
    REAL*8 :: pot, yi
    one=dcmplx(0.0,1.0)

    allocate(h_onsite(norbitals,norbitals))
    allocate(h_onsite_B(norbitals,norbitals))
    allocate(h_hopping_x(norbitals,norbitals))
    allocate(h_hopping_y(norbitals,norbitals))
    allocate(h_hopping_diagonal(norbitals,norbitals))
    allocate(h_hopping_diagonal_m(norbitals,norbitals))

    allocate(ordering(-nnx:nnx,-nny:nny))

    ordering = 0
    nn = 1
    DO j = -nnx, nnx
      DO i = -nny, nny
          ordering(i,j) = nn
          nn=nn+1
      END DO 
    END DO 

    h_hopping_y = transpose(reshape((/ dcmplx(-tl,0.0), dcmplx(0.0,0.0), dcmplx(-drso/2.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
    dcmplx(0.0,0.0), dcmplx(-tl,0.0), dcmplx(0.0,0.0), dcmplx(-drso/2.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
    dcmplx(drso/2.0,0.0), dcmplx(0.0,0.0), dcmplx(-th,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
    dcmplx(0.0,0.0), dcmplx(drso/2.0,0.0), dcmplx(0.0,0.0), dcmplx(-th,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
    dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-tl,0.0), dcmplx(0.0,0.0), &
    dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-tl,0.0) &
 /), shape(h_hopping_y))) 

    h_onsite_B = transpose(reshape((/ dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
      dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
      dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
      dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
      dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
      dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0) &
        /), shape(h_onsite_B))) 

    h_hopping_diagonal = transpose(reshape((/ dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(td/2.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(td/2.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(td/2.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(td/2.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0) &
          /), shape(h_hopping_diagonal))) 
    
    h_hopping_diagonal_m = transpose(reshape((/ dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-td/2.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-td/2.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-td/2.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
        dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-td/2.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0) &
        /), shape(h_hopping_diagonal))) 


    ham=dcmplx(0.0,0.0)
    DO i = -nnx, nnx
      DO j = -nny, nny
      
        
        pot=potential(i,j)
        !definition of submatrices to create a total Hamiltonian
        h_onsite = transpose(reshape((/ dcmplx(4.0*tl-dE+pot+0.5*g*mub*Bz,0.0), dcmplx(0.5*g*mub*Bx,-0.5*g*mub*By), dcmplx(0.0,mub*Bx), dcmplx(0.0,dso/3.0), dcmplx(0.0,-mub*By), dcmplx(-dso/3.0,0.0), &
        dcmplx(0.5*g*mub*Bx,0.5*g*mub*By), dcmplx(4.0*tl-dE+pot-0.5*g*mub*Bz,0.0), dcmplx(0.0,dso/3.0), dcmplx(0.0,mub*Bx), dcmplx(dso/3.0,0.0), dcmplx(0.0,-mub*By), &
        dcmplx(0.0,-mub*Bx), dcmplx(0.0,-dso/3.0), dcmplx(2.0*tl+2.0*th+pot+0.5*g*mub*Bz,0.0), dcmplx(0.5*g*mub*Bx,-0.5*g*mub*By), dcmplx(0.0,dso/3.0+mub*Bz), dcmplx(0.0,0.0), &
        dcmplx(0.0,-dso/3.0), dcmplx(0.0,-mub*Bx), dcmplx(0.5*g*mub*Bx,0.5*g*mub*By), dcmplx(2.0*tl+2.0*th+pot-0.5*g*mub*Bz,0.0), dcmplx(0.0,0.0), dcmplx(0.0,-dso/3.0+mub*Bz), &
        dcmplx(0.0,mub*By), dcmplx(dso/3.0,0.0), dcmplx(0.0,-dso/3.0-mub*Bz), dcmplx(0.0,0.0), dcmplx(2.0*tl+2.0*th+pot+0.5*g*mub*Bz,0.0), dcmplx(0.5*g*mub*Bx,-0.5*g*mub*By), &
        dcmplx(-dso/3.0,0.0), dcmplx(0.0,mub*By), dcmplx(0.0,0.0), dcmplx(0.0,dso/3.0-mub*Bz), dcmplx(0.5*g*mub*Bx,0.5*g*mub*By), dcmplx(2.0*tl+2.0*th+pot-0.5*g*mub*Bz,0.0) &
      /), shape(h_onsite)))

        yi=j*dx
        orb=exp(-0.5*one*Bz*(dx)*(2*yi))
        h_hopping_x = transpose(reshape((/ dcmplx(-tl,0.0)*orb, dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-drso/2.0,0.0)*orb, dcmplx(0.0,0.0), &
           dcmplx(0.0,0.0), dcmplx(-tl,0.0)*orb, dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-drso/2.0,0.0)*orb, &
           dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-tl,0.0)*orb, dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
           dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-tl,0.0)*orb, dcmplx(0.0,0.0), dcmplx(0.0,0.0), &
           dcmplx(drso/2.0,0.0)*orb, dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-th,0.0)*orb, dcmplx(0.0,0.0), &
           dcmplx(0.0,0.0), dcmplx(drso/2.0,0.0)*orb, dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(0.0,0.0), dcmplx(-th,0.0)*orb &
        /), shape(h_hopping_x)))  

        IF(i.eq.-nnx.and.j.ne.-nny.and.j.ne.nny) THEN
          !left facet
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j), ham, h_hopping_x, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i,j+1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j+1), ham, h_hopping_diagonal_m, size, norbitals)

        ELSE IF(i.eq.nnx.and.j.ne.-nny.and.j.ne.nny) THEN
          !right facet
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i,j+1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i-1,j+1), ham, h_hopping_diagonal, size, norbitals)

        ELSE IF(j.eq.nny.and.i.ne.-nnx.and.i.ne.nnx) then
          !top facet
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j), ham, h_hopping_x, size, norbitals)
          
        ELSE IF(i.eq.-nnx.and.j.eq.-nny) THEN
          !left-bottom corner
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j), ham, h_hopping_x, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i,j+1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j+1), ham, h_hopping_diagonal_m, size, norbitals)
        
        ELSE IF(i.eq.nnx.and.j.eq.-nny) THEN
          !right-bottom corner
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i,j+1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i-1,j+1), ham, h_hopping_diagonal, size, norbitals)
        
        ELSE IF(i.eq.-nnx.and.j.eq.nny) THEN
          !left-top corner
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j), ham, h_hopping_x, size, norbitals)

        ELSE IF(i.eq.nnx.and.j.eq.nny) THEN
          !right-top corner
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)

        ELSE
          CALL put_matrix(ordering(i,j), ordering(i,j), ham, h_onsite, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j), ham, h_hopping_x, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i,j+1), ham, h_hopping_y, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i+1,j+1), ham, h_hopping_diagonal_m, size, norbitals)
          CALL put_matrix(ordering(i,j), ordering(i-1,j+1), ham, h_hopping_diagonal, size, norbitals)
          
        ENDIF

      ENDDO
    ENDDO

    call fill_hermitian(ham, size)

    deallocate(h_onsite)
    deallocate(h_onsite_B)
    deallocate(h_hopping_x)
    deallocate(h_hopping_y)
    deallocate(h_hopping_diagonal)
    deallocate(h_hopping_diagonal_m)
    deallocate(ordering)
END SUBROUTINE HAMILTONIAN_CREATE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! put submatirx into large matrix !!!!!!!!! !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE put_matrix(ni, nj, matrix_A, matrix_B, sizeA, sizeB)

    INTEGER, INTENT(IN) :: ni,nj
    INTEGER, INTENT(IN) :: sizeA
    INTEGER, INTENT(IN) :: sizeB
    COMPLEX*16, INTENT(INOUT) :: matrix_A(sizeA,sizeA) 
    COMPLEX*16, INTENT(IN) :: matrix_B(sizeB,sizeB) 

    INTEGER :: i,j

    if(mod(sizeA,sizeB).ne.0) STOP "put_matrix: sizeA should be an integral multiple of sizeB"

    DO i = 1, sizeB
      DO j = 1, sizeB
        matrix_A((ni-1)*sizeB+i,(nj-1)*sizeB+j)=matrix_A((ni-1)*sizeB+i,(nj-1)*sizeB+j)+matrix_B(i,j)
      ENDDO
    ENDDO

END SUBROUTINE put_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! make full matrix from the upper part when it is hermitian !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE fill_hermitian(matrix_A, sizeA)

  INTEGER, INTENT(IN) :: sizeA
  COMPLEX*16, INTENT(INOUT) :: matrix_A(sizeA,sizeA) 
  INTEGER :: i,j

  DO i = 1, sizeA
    DO j = i+1, sizeA
      matrix_A(j,i)=conjg(matrix_A(i,j))
    ENDDO
  ENDDO

END SUBROUTINE fill_hermitian

END MODULE hamiltonian
