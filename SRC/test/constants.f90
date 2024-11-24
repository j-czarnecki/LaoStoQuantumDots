  !-----------------------------------------------------------------------
  MODULE constants
  !-----------------------------------------------------------------------
  
  use, intrinsic :: iso_fortran_env
  IMPLICIT NONE
  SAVE

  ! numerical accurancy
  integer, parameter :: SP = REAL32
  integer, parameter :: DP = REAL64
  integer, parameter :: QP = REAL128
  !
  ! Mathematical constants
  ! 
  REAL(DP), PARAMETER :: pi     = 3.141592653589793238462643383279502884197169399375105820974944d0
  COMPLEX(DP), PARAMETER :: ci    = (0.d0, 1.d0)
  COMPLEX(DP), PARAMETER :: cone  = (1.d0, 0.d0)
  COMPLEX(DP), PARAMETER :: czero = (0.d0, 0.d0)
  COMPLEX*16, PARAMETER :: imag = dcmplx(0.0d0, 1.0d0)
  !
  ! Unit conversion factors
  !
  REAL(DP), PARAMETER :: eV2au   = 0.03674932587122423
  REAL(DP), PARAMETER :: nm2au   = 18.89726133921252
  REAL(DP), PARAMETER :: T2au    = 4.254382E-6
  REAL(DP), PARAMETER :: ns2au   = 4.1341373336493 * 1e7
  REAL(DP), PARAMETER :: F2au    = 1.0d0/(5.14220652*1e11) ![V/m] ---> [a.u.]
  REAL(DP), PARAMETER :: mub     = 0.5
  REAL*8, PARAMETER :: m_eff     = 0.286

  !Hamiltonian geometry parameters
  INTEGER*4, PARAMETER :: N_INTERIOR_ELEMENTS = 47 !Number of non-zero elements corresponding to nodes on the interior of grid
  INTEGER*4, PARAMETER :: N_RIGHT_FACET_ELEMENTS = 33 !Number of non-zero elements corresponding to nodes on the right edge of grid
  INTEGER*4, PARAMETER :: N_LEFT_FACET_ELEMENTS = 43 !Number of non-zero elements corresponding to nodes on the left edge of grid
  INTEGER*4, PARAMETER :: N_TOP_FACET_ELEMENTS = 29 !Number of non-zero elements corresponding to nodes on the top edge of grid
  INTEGER*4, PARAMETER :: N_RIGHT_TOP_CORNER_ELEMENTS = 19 !Number of non-zero elements corresponding to node in the right-top corner of grid
  END MODULE constants

