TARGET = ../bin/lao_sto_qd.x


## the command shell
SHELL = /bin/sh


F90 = gfortran #ifort
F77 = gfortran #ifort
CC =


F90FLAGS = -Ofast  -m64 -ffree-line-length-none
F77FLAGS = -Ofast  -m64 -ffree-line-length-none

#F90FLAGS = -w -O3
#F77FLAGS = -w -O3 


CFLAGS =

## name of the program to link the program units
LNK = $(F90)

## flags to use at link time
LFLAGS = $(F90FLAGS)

## LLIBS are libraries to use at link time
LIBS = -L ./lib -larpack
LIB_OMP = #-fopenmp

LIBS_MKL = -I${MKLROOT}/include \
	   -I/opt/intel/mkl/include \
	   -Wl,--start-group \
           ${MKLROOT}/lib/intel64/libmkl_gf_lp64.so \
           ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so \
           ${MKLROOT}/lib/intel64/libmkl_core.so \
          -Wl,--end-group \
          -lgomp -lpthread -lm -ldl


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
	

#################################################################
#      rules 
#################################################################

$(TARGET) : $(OBJS)
	$(LNK) $(LIB_OMP) -o $(TARGET) $(LFLAGS) $^ $(LIBS) $(LIBS_MKL)

%.o : %.f90
	$(F90) $(LIB_OMP) $(F90FLAGS) -c $< -o $@

%.o : %.f
	$(F77) $(LIB_OMP) $(F77FLAGS) -c $< -o $@

all: $(TARGET)

debug: F90FLAGS = -O0  -fsanitize=thread -m64 -ffree-line-length-none -Wall -g
debug: F77FLAGS = -O0  -fsanitize=thread -m64 -ffree-line-length-none -Wall -g
debug: $(TARGET)

profile: F90FLAGS = -Ofast  -m64 -ffree-line-length-none -Wall -pg -fopt-info
profile: F77FLAGS = -Ofast  -m64 -ffree-line-length-none -Wall -pg -fopt-info
profile: $(TARGET)

clean:
	rm -f $(OBJS) $(OBJS_PLOT) $(OBJS_ANALYZE)
	rm -f *.mod
	rm -f $(TARGET)

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






	