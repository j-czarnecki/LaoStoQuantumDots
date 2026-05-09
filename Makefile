TARGET = ./bin/lao_sto_qd.x
POSTPROCESSING_TARGET = ./bin/post_qd.x
SRC_DIR = SRC
OBJ_DIR = OBJ
MOD_DIR = MOD
## the command shell
SHELL = /bin/sh

CXX = g++
CC = gcc
F90 = gfortran #ifort
F77 = gfortran #ifort

CFLAGS =

## name of the program to link the program units
LNK = $(F90)

## flags to use at link time
LFLAGS = $(F90FLAGS)

## LLIBS are libraries to use at link time
LIBS = -L ${SRC_DIR}/lib -larpack
LIB_OMP = -fopenmp

LIBS_MKL = -I${MKLROOT}/include \
	   			 -I/opt/intel/mkl/include \
	   			 -Wl,--start-group \
           ${MKLROOT}/lib/intel64/libmkl_gf_lp64.so \
           ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so \
           ${MKLROOT}/lib/intel64/libmkl_core.so \
           -Wl,--end-group \
           -lgomp -lpthread -lm -ldl


F90FLAGS = -O3 -cpp -m64 -ffree-line-length-none $(LIB_OMP) -J$(MOD_DIR)

# Common debug base flags
F90_DEBUG_BASE = -O0 -cpp -m64 -ffree-line-length-none -g $(LIB_OMP) -J$(MOD_DIR) -DDEBUG
F90_DEBUG_BASE += -Wall -Wextra -Wpedantic -Wconversion
F90_DEBUG_BASE += -fbacktrace

# Standard debug build with runtime checks
F90_DEBUG_FLAGS = $(F90_DEBUG_BASE)
F90_DEBUG_FLAGS += -fcheck=all
F90_DEBUG_FLAGS += -ffpe-trap=invalid,zero,overflow
F90_DEBUG_FLAGS += -finit-real=snan -finit-integer=-2147483647 -finit-logical=true -finit-character=42

# AddressSanitizer build
F90_ASAN_FLAGS = $(F90_DEBUG_BASE)
F90_ASAN_FLAGS += -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer
LIBS_ASAN = -fsanitize=address -fsanitize=leak

F90_TSAN_FLAGS = -O0 -cpp -m64 -ffree-line-length-none -g $(LIB_OMP) -J$(MOD_DIR) -DDEBUG
F90_TSAN_FLAGS += -fsanitize=thread
LIBS_TSAN = -fsanitize=thread

#################################################################
#      SOURCE AND OBJECT FILES
#################################################################

# --- Automatically find all source files recursively ---
SRC_FILES_ALL := $(shell find $(SRC_DIR) -name '*.f90')

# --- Exclude the two main programs from the common source set ---
SRC_COMMON := $(filter-out $(SRC_DIR)/main.f90 $(SRC_DIR)/main_postprocessing.f90, $(SRC_FILES_ALL))

# --- Define two build sets ---
SRC_FILES_MAIN := $(SRC_COMMON) $(SRC_DIR)/main.f90
SRC_FILES_POST := $(SRC_COMMON) $(SRC_DIR)/main_postprocessing.f90

# --- Define corresponding object files ---
OBJS_MAIN := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES_MAIN))
OBJS_POST := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES_POST))

ASMS_MAIN := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.s,$(SRC_FILES_MAIN))
ASMS_POST := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.s,$(SRC_FILES_POST))

# --- Automatically find directories with unit tests ---
UNITTEST_DIRS := $(dir $(shell find $(SRC_DIR) -type d -name test))

#################################################################
#      RULES
#################################################################

.PHONY: all debug profile tsan asan assembly test clean post run_slurm

$(TARGET) : $(OBJS_MAIN)
	$(LNK) -o $(TARGET) $(LFLAGS) $^ $(LIBS) $(LIBS_MKL)

$(POSTPROCESSING_TARGET) : $(OBJS_POST)
	$(LNK) -o $(POSTPROCESSING_TARGET) $(LFLAGS) $^ $(LIBS) $(LIBS_MKL)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(F90FLAGS) -c $< -o $@

$(OBJ_DIR)/%.s : $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(F90FLAGS) -S $< -o $@

all: $(TARGET)

debug: F90FLAGS = $(F90_DEBUG_FLAGS)
debug: $(TARGET)

profile: F90FLAGS = -O3 -cpp -m64 -ffree-line-length-none -Wall -pg -fopt-info $(LIB_OMP)
profile: $(TARGET)

#To avoid Thread Sanitizer error about bad memory mapping
#echo 0 | sudo tee /proc/sys/kernel/randomize_va_space
#To include suppressions run as
#TSAN_OPTIONS="suppressions=thread_suppressions.txt:history_size=7" bin/lao_sto_qd.x
tsan: F90FLAGS = $(F90_TSAN_FLAGS)
tsan: LIBS += $(LIBS_TSAN)
tsan: $(TARGET)

asan: F90FLAGS = $(F90_ASAN_FLAGS)
asan: LIBS += $(LIBS_ASAN)
asan: $(TARGET)

assembly: $(ASMS)

test:
	mkdir -p $(SRC_DIR)/test/$(MOD_DIR)
	cp SRC/*.f90 SRC/test/
	@export FC="$(F90)" && export CC="$(CC)" && export CXX="$(CXX)" && export FSFLAG=-I && export FCFLAGS="$(F90FLAGS)" && cd $(SRC_DIR)/test &&	funit
	cd ../../

post: $(POSTPROCESSING_TARGET)

run_slurm:
	cd Runner && python3 runnerMain.py && cd ..

clean:
	rm -rf $(MOD_DIR)
	rm -rf $(OBJ_DIR)
	rm -f $(TARGET)
	rm -f $(SRC_DIR)/test/*.f90
	rm -rf $(SRC_DIR)/test/$(MOD_DIR)
	rm -rf $(SRC_DIR)/test/*.o
	rm -f *.mod
	rm -rf $(SRC_DIR)/*.i90
# Exclude as funit should be changed to pFUnit
#@export CC="$(CC)" && export CXX="$(CXX)" && cd $(SRC_DIR)/test &&	funit --clean && cd ../../



#################################################################
#      DEPENDENCIES
#################################################################
$(OBJ_DIR)/main.o: $(OBJ_DIR)/indata.o \
									 $(OBJ_DIR)/hamiltonian.o \
									 $(OBJ_DIR)/diagonalize.o \
									 $(OBJ_DIR)/combinatory.o \
									 $(OBJ_DIR)/many_body.o \
									 $(OBJ_DIR)/writers.o \
									 $(OBJ_DIR)/utility.o \
									 $(OBJ_DIR)/constants.o \
									 $(OBJ_DIR)/logger.o \
									 $(OBJ_DIR)/swap.o \
									 $(OBJ_DIR)/time_dependence.o \
									 $(OBJ_DIR)/potentials.o \
									 $(OBJ_DIR)/broydenV2.o

$(OBJ_DIR)/main.o:
$(OBJ_DIR)/combinatory.o:

$(OBJ_DIR)/broydenV2.o:

$(OBJ_DIR)/time_dependence.o: $(OBJ_DIR)/logger.o \
															$(OBJ_DIR)/utility.o \
															$(OBJ_DIR)/many_body.o \
															$(OBJ_DIR)/constants.o \
															$(OBJ_DIR)/indata.o

$(OBJ_DIR)/potentials.o: $(OBJ_DIR)/logger.o \
												 $(OBJ_DIR)/utility.o \
												 $(OBJ_DIR)/constants.o \
												 $(OBJ_DIR)/indata.o

$(OBJ_DIR)/indata.o: $(OBJ_DIR)/constants.o \
										 $(OBJ_DIR)/logger.o

$(OBJ_DIR)/hamiltonian.o: $(OBJ_DIR)/indata.o \
													$(OBJ_DIR)/logger.o

$(OBJ_DIR)/diagonalize.o: $(OBJ_DIR)/indata.o \
													$(OBJ_DIR)/hamiltonian.o \
													$(OBJ_DIR)/logger.o

$(OBJ_DIR)/many_body.o: $(OBJ_DIR)/combinatory.o \
												$(OBJ_DIR)/utility.o \
												$(OBJ_DIR)/logger.o

$(OBJ_DIR)/writers.o: $(OBJ_DIR)/indata.o \
											$(OBJ_DIR)/constants.o \
											$(OBJ_DIR)/utility.o \
											$(OBJ_DIR)/many_body.o \
											$(OBJ_DIR)/logger.o

$(OBJ_DIR)/utility.o: $(OBJ_DIR)/constants.o \
											$(OBJ_DIR)/diagonalize.o

$(OBJ_DIR)/logger.o:

$(OBJ_DIR)/constants.o:

$(OBJ_DIR)/swap.o: $(OBJ_DIR)/indata.o \
    							 $(OBJ_DIR)/constants.o \
    							 $(OBJ_DIR)/logger.o \
    							 $(OBJ_DIR)/utility.o \
    							 $(OBJ_DIR)/combinatory.o \
    							 $(OBJ_DIR)/many_body.o


$(OBJ_DIR)/main_postprocessing.o: $(OBJ_DIR)/indata.o \
																	$(OBJ_DIR)/hamiltonian.o \
																	$(OBJ_DIR)/diagonalize.o \
																	$(OBJ_DIR)/combinatory.o \
																	$(OBJ_DIR)/many_body.o \
																	$(OBJ_DIR)/writers.o \
																	$(OBJ_DIR)/utility.o \
																	$(OBJ_DIR)/constants.o \
																	$(OBJ_DIR)/logger.o

