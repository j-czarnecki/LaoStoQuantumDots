TARGET = ../bin/lao_sto_qd.x

## the command shell
SHELL = /bin/sh


F90 = gfortran #ifort
F77 = gfortran #ifort

CC =
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


F90FLAGS = -Ofast  -m64 -ffree-line-length-none $(LIB_OMP)
F77FLAGS = -Ofast  -m64 -ffree-line-length-none $(LIB_OMP)

#################################################################
#      objects files
#################################################################

OBJS =  main.o \
			  constants.o \
        indata.o \
        hamiltonian.o \
				diagonalize.o \
				combinatory.o \
				many_body.o \
				writers.o \
				utility.o

ASMS =  main.s \
			  constants.s \
        indata.s \
        hamiltonian.s \
				diagonalize.s \
				combinatory.s \
				many_body.s \
				writers.s \
				utility.s
#################################################################
#      rules
#################################################################

$(TARGET) : $(OBJS)
	$(LNK) -o $(TARGET) $(LFLAGS) $^ $(LIBS) $(LIBS_MKL)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

%.o : %.f
	$(F77) $(F77FLAGS) -c $< -o $@

%.s : %.f90
	$(F90) $(F90FLAGS) -S $< -o $@

all: $(TARGET)

#To avoid Thread Sanitizer error about bad memory mapping
#echo 0 | sudo tee /proc/sys/kernel/randomize_va_space
#To include suppressions run as
#TSAN_OPTIONS="suppressions=thread_suppressions.txt:history_size=7" bin/lao_sto_qd.x
debug: F90FLAGS = -O0 -fsanitize=thread -m64 -ffree-line-length-none -Wall -g $(LIB_OMP)
debug: F77FLAGS = -O0 -fsanitize=thread -m64 -ffree-line-length-none -Wall -g $(LIB_OMP)
debug: $(TARGET)

profile: F90FLAGS = -Ofast  -m64 -ffree-line-length-none -Wall -pg -fopt-info $(LIB_OMP)
profile: F77FLAGS = -Ofast  -m64 -ffree-line-length-none -Wall -pg -fopt-info $(LIB_OMP)
profile: $(TARGET)

assembly: $(ASMS)

clean:
	rm -f $(OBJS) $(OBJS_PLOT) $(OBJS_ANALYZE)
	rm -f *.mod
	rm -f $(TARGET)
	rm -f *.s

#Dependencies
main.o: indata.o \
				hamiltonian.o \
				diagonalize.o \
				combinatory.o \
				many_body.o \
				writers.o \
				utility.o \
				constants.o

combinatory.o:

indata.o: constants.o

hamiltonian.o: indata.o

diagonalize.o: indata.o \
							 hamiltonian.o

many_body.o: combinatory.o \
						 utility.o

writers.o: indata.o \
					 constants.o \
					 utility.o \
					 many_body.o

utility.o:

constants.mod: #Has to be .mod, or else module is not found






	