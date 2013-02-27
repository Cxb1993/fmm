.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  

LD = gfortran
CF90 = gfortran

AR = ar rv
OPT_FFLAGS = -std=f2008 -march=native -O3 -fcheck=all

FFLAGS =  $(OPT_FFLAGS)
LINKFLAGS = $(OPT_FFLAGS)

LIBS = -DYA_BLAS -DYA_LAPACK -DYA_BLASMULT -lblas -llapack
INCLUDE	= -L/opt/local/lib/

.f90.o:
	$(CF90) $(FFLAGS) $(INCLUDE) -o $*.o -c $*.f90 $(LIBS)

include source.files
include target.mk

OBJECTS = $(f90FILES:.f90=.o)
OMOD = $(MODULES:.f90=.o)

all:	$(TARGET)

$(TARGET): $(OMOD) $(OBJECTS)
	$(LD) $(LINKFLAGS) $(INCLUDE) -o $(TARGET) $(OMOD) $(OBJECTS) $(LIBS)
clean:
	rm -f $(TARGET)
	rm -rf *.o *.mod *.so *.a
