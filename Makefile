TARGET = ./bin/lao_sto_qd.x
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
LIBS = -L ./lib -larpack
LIB_OMP = -fopenmp

LIBS_MKL = -I${MKLROOT}/include \
	   -I/opt/intel/mkl/include \
	   -Wl,--start-group \
           ${MKLROOT}/lib/intel64/libmkl_gf_lp64.so \
           ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so \
           ${MKLROOT}/lib/intel64/libmkl_core.so \
          -Wl,--end-group \
          -lgomp -lpthread -lm -ldl


F90FLAGS = -Ofast -cpp -m64 -ffree-line-length-none $(LIB_OMP) -J$(MOD_DIR)

#################################################################
#      objects files
#################################################################

OBJS =  $(OBJ_DIR)/main.o \
			  $(OBJ_DIR)/constants.o \
        $(OBJ_DIR)/indata.o \
        $(OBJ_DIR)/hamiltonian.o \
				$(OBJ_DIR)/diagonalize.o \
				$(OBJ_DIR)/combinatory.o \
				$(OBJ_DIR)/many_body.o \
				$(OBJ_DIR)/writers.o \
				$(OBJ_DIR)/utility.o \
				$(OBJ_DIR)/logger.o

ASMS =  $(OBJ_DIR)/main.s \
        $(OBJ_DIR)/constants.s \
        $(OBJ_DIR)/indata.s \
        $(OBJ_DIR)/hamiltonian.s \
        $(OBJ_DIR)/diagonalize.s \
        $(OBJ_DIR)/combinatory.s \
        $(OBJ_DIR)/many_body.s \
        $(OBJ_DIR)/writers.s \
        $(OBJ_DIR)/utility.s \
        $(OBJ_DIR)/logger.s
#################################################################
#      rules
#################################################################

$(TARGET) : $(OBJS)
	$(LNK) -o $(TARGET) $(LFLAGS) $^ $(LIBS) $(LIBS_MKL)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(F90FLAGS) -c $< -o $@

$(OBJ_DIR)/%.s : $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(F90FLAGS) -S $< -o $@

all: $(TARGET)

#To avoid Thread Sanitizer error about bad memory mapping
#echo 0 | sudo tee /proc/sys/kernel/randomize_va_space
#To include suppressions run as
#TSAN_OPTIONS="suppressions=thread_suppressions.txt:history_size=7" bin/lao_sto_qd.x
debug: F90FLAGS = -O0 -cpp -m64 -ffree-line-length-none -Wall -g $(LIB_OMP) -DDEBUG
debug: $(TARGET)

profile: F90FLAGS = -Ofast -cpp -m64 -ffree-line-length-none -Wall -pg -fopt-info $(LIB_OMP)
profile: $(TARGET)

tsan: F90FLAGS = -O0 -cpp -m64 -ffree-line-length-none -Wall -g $(LIB_OMP) -DDEBUG -fsanitize=thread
tsan: $(TARGET)

asan: F90FLAGS = -O0 -cpp -m64 -ffree-line-length-none -Wall -g $(LIB_OMP) -DDEBUG -fsanitize=address
asan: $(TARGET)

assembly: $(ASMS)

test:
	mkdir -p $(SRC_DIR)/test/$(MOD_DIR)
	cp SRC/*.f90 SRC/test/
	@export FC="$(F90)" && export CC="$(CC)" && export CXX="$(CXX)" && export FSFLAG=-I && export FCFLAGS="$(F90FLAGS)" && cd $(SRC_DIR)/test &&	funit
	cd ../../

clean:
	rm -rf $(MOD_DIR)
	rm -rf $(OBJ_DIR)
	rm -f $(TARGET)
	rm -f $(SRC_DIR)/test/*.f90
	rm -rf $(SRC_DIR)/test/$(MOD_DIR)
	rm -rf $(SRC_DIR)/test/*.o
	rm -f *.mod
	@export CC="$(CC)" && export CXX="$(CXX)" && cd $(SRC_DIR)/test &&	funit --clean && cd ../../
#Dependencies
$(OBJ_DIR)/main.o: $(OBJ_DIR)/indata.o \
									 $(OBJ_DIR)/hamiltonian.o \
									 $(OBJ_DIR)/diagonalize.o \
									 $(OBJ_DIR)/combinatory.o \
									 $(OBJ_DIR)/many_body.o \
									 $(OBJ_DIR)/writers.o \
									 $(OBJ_DIR)/utility.o \
									 $(OBJ_DIR)/constants.o \
									 $(OBJ_DIR)/logger.o

$(OBJ_DIR)/main.o:
$(OBJ_DIR)/combinatory.o:

$(OBJ_DIR)/indata.o: $(OBJ_DIR)/constants.o

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
											$(OBJ_DIR)/many_body.o

$(OBJ_DIR)/utility.o:

$(OBJ_DIR)/logger.o:

$(OBJ_DIR)/constants.o:

