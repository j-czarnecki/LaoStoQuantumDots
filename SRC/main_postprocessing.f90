PROGRAM MAIN_POSTPROCESSING
USE postprocessing
USE logger
USE omp_lib

IMPLICIT NONE

INTEGER*4 :: nargs, i
CHARACTER(len=256) :: arg
CHARACTER(len=256) :: path

nargs = COMMAND_ARGUMENT_COUNT()
WRITE (*, *) "Count of CLI arguments: ", nargs

IF (nargs > 1) THEN
  WRITE (*, *) "ERROR: Number of arguments cannot be > 1!"
END IF

CALL GET_COMMAND_ARGUMENT(1, arg)
path = TRIM(arg)
WRITE (*, *) path

CALL INIT_LOGGER()

!CALL CALCULATE_INTEGRAL_ELEMENTS()
CALL CALCULATE_IMAGE_EXPECTATION_VALUE(TRIM(path))

CALL CLOSE_LOGGER()

END PROGRAM
