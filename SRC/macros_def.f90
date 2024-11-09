!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file contains macros definitions used in the code
! It should be included in every file that uses those macros
! via #include "macros_def.f90"
! It makes use of logger functions defined in mod_logger.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Debug traces only if the DEBUG flag is active
! One has to specify -fpp flag during compilation to include preprocessor directives
! -DDEBUG to include debug traces
#ifdef DEBUG
#define LOG_DEBUG(logMsg) CALL LOG_STRING_DEBUG(logMsg)
#else
#define LOG_DEBUG(logMsg)
#endif

! Always define info traces
#define LOG_INFO(logMsg) CALL LOG_STRING_INFO(logMsg)
! Always define abnormal traces
#define LOG_ABNORMAL(logMsg) CALL LOG_STRING_ABNORMAL(logMsg)
! Always define error traces
#define LOG_ERROR(logMsg) CALL LOG_STRING_ERROR(logMsg)