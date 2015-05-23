##################################################
##### MACRO DEFINITIONS ##########################
##################################################

TARGET  = a.out 
# COMMON_MOD1 = file.f
# COMMON_MOD2 = file.for 
# COMMON_MOD3 = file.f90
COMMON_MOD3 = basis.f90 boundary.f90 main.f90 
COMMON_MOD  = $(COMMON_MOD1)       $(COMMON_MOD2)         $(COMMON_MOD3)
OBJECTS     = $(COMMON_MOD1:.f=.o) $(COMMON_MOD2:.for=.o) $(COMMON_MOD3:.f90=.o)
# f90: (ansi/iso or iso standard) modern fortran file
# for: (ansi/iso standard) punchcard fortran file 
# f  : (ibm language) punchcard fortran file

# # gnu compiler: 
# # FORTRAN = gfortran
# FORTRAN = /usr/local/bin/gfortran
# FFLAGS += -fimplicit-none -fbounds-check
# FFLAGS1 = -std=legacy -x f77
# FFLAGS2 = -std=legacy
# FFLAGS3 = 
# intel compiler: 
# FORTRAN = ifort
FORTRAN = /opt/intel/bin/ifort
# FFLAGS += 
FFLAGS1 = -nostand -f66 
FFLAGS2 = -nostand 
# FFLAGS3 = 


# # gnu library
# FFLAGS  += -I/usr/include
# LDFLAGS += -L/usr/lib
# LDFLAGS += -llapack -lblas
# # openblas library
# FFLAGS  += -I/usr/local/opt/openblas/include
# LDFLAGS += -L/usr/local/opt/openblas/lib
# LDFLAGS += -lopenblas
# mkl library
MKLROOT  = /opt/intel/mkl
# MKLROOT  = /opt/intel/composer_xe_2015.0.077/mkl
# FFLAGS  += -i8 -I$(MKLROOT)/include -I$(MKLROOT)/include/ilp64 
FFLAGS  += -i8 -I$(MKLROOT)/include -I$(MKLROOT)/include/intel64/ilp64 
LDFLAGS += $(MKLROOT)/lib/libmkl_blas95_ilp64.a $(MKLROOT)/lib/libmkl_lapack95_ilp64.a
# LDFLAGS += -L$(MKLROOT)/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
LDFLAGS += $(MKLROOT)/lib/libmkl_intel_ilp64.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_sequential.a -lpthread -lm


##################################################
##### INFERENCE RULRS ############################
##################################################

.PHONY:    all check clean
.SUFFIXES: .o .mod .f .for .f90

all: ${TARGET} 
check:  ${OBJECTS}
clean:; rm -f *.o *.mod fort.*

${TARGET}: ${OBJECTS}; ${FORTRAN} -o $@ ${OBJECTS} ${LDFLAGS}
%.o: %.mod
.f90.o:; ${FORTRAN} -c $< ${FFLAGS} ${FFLAGS3}
.for.o:; ${FORTRAN} -c $< ${FFLAGS} ${FFLAGS2}
.f.o:  ; ${FORTRAN} -c $< ${FFLAGS} ${FFLAGS1}


##################################################
##### FILE DEPENDENCIES ##########################
##################################################

${OBJECTS}: ${COMMON_MOD}
# main.o: input.o
# input.o: diag.o
