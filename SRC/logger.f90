MODULE logger
USE omp_lib
IMPLICIT NONE
SAVE
INTEGER*4, PRIVATE, PARAMETER :: LOGGER_UNIT = 66
INTEGER*4, PRIVATE, PARAMETER :: MAX_LOG_LEN = 2000
REAL*8, PRIVATE :: T_START, T_END
CHARACTER(LEN=MAX_LOG_LEN) :: log_string
CONTAINS

RECURSIVE SUBROUTINE INIT_LOGGER()
    OPEN(unit = LOGGER_UNIT, FILE = "./log.log", FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(LOGGER_UNIT, *) "==== START ===="
    FLUSH(LOGGER_UNIT)
    CALL CPU_TIME(T_START)
END SUBROUTINE

RECURSIVE SUBROUTINE CLOSE_LOGGER()
    CALL CPU_TIME(T_END)
    WRITE(LOGGER_UNIT, *) "Time of simulation (seconds): ", T_END - T_START
    WRITE(LOGGER_UNIT, *) "==== END ===="
    CLOSE(LOGGER_UNIT)

END SUBROUTINE

RECURSIVE SUBROUTINE LOG_STRING_DEBUG(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)
    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
    &'TID = ', omp_get_thread_num(), ' ',&
    &'DEBUG: ' // TRIM(logMsg)
    FLUSH(LOGGER_UNIT)
END SUBROUTINE


RECURSIVE SUBROUTINE LOG_STRING_INFO(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)

    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
    &'TID = ', omp_get_thread_num(), ' ',&
    &"INFO: " // TRIM(logMsg)
    FLUSH(LOGGER_UNIT)
END SUBROUTINE

RECURSIVE SUBROUTINE LOG_STRING_ABNORMAL(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)

    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
    &'TID = ', omp_get_thread_num(), ' ',&
    &"ABNORMAL: " // TRIM(logMsg)
    FLUSH(LOGGER_UNIT)

END SUBROUTINE

RECURSIVE SUBROUTINE LOG_STRING_ERROR(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)

    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
    &'TID = ', omp_get_thread_num(), ' ',&
    &"ERROR: " // TRIM(logMsg)
    FLUSH(LOGGER_UNIT)

END SUBROUTINE

END MODULE