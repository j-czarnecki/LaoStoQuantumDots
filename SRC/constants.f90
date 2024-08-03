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
  !
  ! Unit conversion factors
  !
  REAL(DP), PARAMETER :: eV2au   = 0.03674932587122423
  REAL(DP), PARAMETER :: nm2au   = 18.89726133921252
  REAL(DP), PARAMETER :: T2au    = 4.254382E-6
  REAL(DP), PARAMETER :: mub     = 0.5


  END MODULE constants

