TARGET = fftwpsi3D.exe

#################################################################
#  sp3    Computes the Coulomb integrals from the single-particle
#        wave-functions generated by sp2
# 18 nov 2003        PGI compiler               by Andrea Bertoni
#################################################################
#   Makefile PGIlinux 
#################################################################

## the command shell
SHELL = /bin/bash

## name of the fortran 90, fortran 77, C compilers
F90 = ifort #-autodouble
##pgf90 -Mfreeform -r8 
##-byteswapio
F77 = ifort -autodouble
##pgf77 -r8
CC = cc

## flags to use for fortran 90, fortran 77, C compilation
F90FLAGS = -C -w
##-Mbounds
F77FLAGS = -C -w
CFLAGS = -O

## name of the program to link the program units
LNK = $(F90)

## flags to use at link time
LFLAGS = $(F90FLAGS)

## include dirs
F90_INCLUDE = -I ${MKL_INC} -I ${MKL_INC}/intel64/fftw

## LLIBS are libraries to use at link time
LLIBS = -L ${MKL_LIB} -lmkl_intel_lp64 -lmkl_intel_thread    \
        -lmkl_core -liomp5 -lpthread -lm



#################################################################
#      objects files
#################################################################

OBJS =  main.o              \
        mod_staticdata.o    \
        mod_indata.o        \
        mod_logga.o         \
        mod_fftw3params.o   \

#################################################################
#      rules 
#################################################################

$(TARGET) : $(OBJS)
	$(LNK) -o $(TARGET) $(LFLAGS) $^ $(LLIBS)


%.o : %.f90
	$(F90) $(F90FLAGS) $(F90_INCLUDE) -c $< -o $@

%.o : %.f
	$(F77) $(F77FLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f *.mod
	rm -f $(TARGET)


#################################################################
#      dependencies
#################################################################

main.o :                mod_staticdata.o    \
			mod_indata.o        \
			mod_logga.o         \
			mod_fftw3params.o   \
			initsymci.o         \
			coulombftran.o      \
			outputci.o          \
			quicksortci.o       \
			adapt.o

mod_staticdata.o: 

mod_indata.o:		mod_staticdata.o

mod_logga.o:		mod_indata.o

mod_fftw3params.o:  

initsymci.o:

coulombftran.o:		mod_indata.o     \
			adapt.o

outputci.o:		mod_indata.o

adapt.o: 

quicksortci.o:		mod_indata.o
