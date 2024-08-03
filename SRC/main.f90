
PROGRAM MAIN
  USE indata
  USE hamiltonian
  USE diagonalize
    
  IMPLICIT NONE
  
  REAL*8, allocatable :: potential(:,:)
  INTEGER :: i, state_to_write
  COMPLEX*16, allocatable :: psi_lapack(:,:,:,:,:)
  REAL*8, allocatable :: ev_lapack(:)
  COMPLEX*16, allocatable :: psi_arpack(:,:,:,:,:)
  REAL*8, allocatable :: ev_arpack(:)
  
  INTEGER :: ix, iy, iorb
  REAL*8 :: x, y, omega

  !read a file 
  CALL INDATA_GET("../RUNS/quantum_dot.nml")
  write(*,*) "Physical parameters:"
  write(*,*) "th=", th/eV2au
  write(*,*) "tl=", tl/eV2au
  write(*,*) "td=", td/eV2au
  write(*,*) "dso=", dso/eV2au
  write(*,*) "drso=", drso/eV2au
  write(*,*) "dE=", dE/eV2au
  write(*,*) "g=", g
  write(*,*) "Nx=", Nx
  write(*,*) "Ny=", Ny
  write(*,*) "dx=", dx/nm2au
  write(*,*) "norbs=", norbs
  write(*,*) "nstate=", nstate

  write(*,*)
  write(*,*) "External parameters:"
  write(*,*) "Bx=", Bx/T2au
  write(*,*) "By=", By/T2au
  write(*,*) "Bz=", Bz/T2au

 
  allocate(potential(-Nx:Nx,-Ny:Ny))
  allocate(psi_lapack(nstate,-Nx:Nx,-Ny:Ny,norbs/2,2))
  allocate(ev_lapack(nstate))
  allocate(psi_arpack(nstate,-Nx:Nx,-Ny:Ny,norbs/2,2))
  allocate(ev_arpack(nstate))
  
  omega=3.0e-3*eV2au
  potential=0.0
  do ix=-Nx,Nx
      do iy=-Ny,Ny
        x=ix*dx
        y=iy*dx
        potential(ix,iy)=0.5*omega**2*(x**2+y**2)
      enddo
  enddo
  
  CALL DIAGONALIZE_LAPACK( psi_lapack, ev_lapack, Nx, Ny, norbs, nstate, potential)
  CALL DIAGONALIZE_ARPACK( psi_arpack, ev_arpack, Nx, Ny, norbs, nstate, potential)
  
  write(*,*) "Porownanie energii LAPACK i ARPACK"
  do i=1,nstate
   write(*,'(200e20.12)') ev_lapack(i)/eV2au, ev_arpack(i)/eV2au
  enddo

  
  state_to_write=1
  OPEN(1,FILE="psi_lapack.dat")
  do ix=-Nx,Nx
      do iy=-Ny,Ny
      x=ix*dx
      y=iy*dx
      write(1,*) x, y, (abs(psi_lapack(state_to_write,ix,iy,iorb,1))**2+abs(psi_lapack(state_to_write,ix,iy,iorb,2))**2, iorb=1,3)
      enddo
  enddo
  CLOSE(1)
  
  OPEN(1,FILE="psi_arpack.dat")
  do ix=-Nx,Nx
      do iy=-Ny,Ny
      x=ix*dx
      y=iy*dx
      write(1,*) x, y, (abs(psi_arpack(state_to_write,ix,iy,iorb,1))**2+abs(psi_arpack(state_to_write,ix,iy,iorb,2))**2, iorb=1,3)
      enddo
  enddo
  CLOSE(1)


  deallocate(potential)
  deallocate(psi_lapack)
  deallocate(ev_lapack)
  deallocate(psi_arpack)
  deallocate(ev_arpack)
END PROGRAM MAIN

