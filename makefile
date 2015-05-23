##################################################
##### MACRO DEFINITIONS ##########################
##################################################

TARGET  = a.out 
# COMMON_MOD1 = file.f
# COMMON_MOD2 = file.for 
# COMMON_MOD3 = file.f90
COMMON_MOD3 += nrtype.f90 nrutil.f90 nr.f90 plgndr.f90 bessj0.f90 bessj1.f90 bessj.f90 bessy0.f90 bessy1.f90 bessy.f90 
COMMON_MOD3 += const.f90 
COMMON_MOD3 += global.f90 hamiltonian.f90 basis.f90 boundary.f90 inner.f90 outer.f90 main.f90 
COMMON_MOD  = $(COMMON_MOD1)       $(COMMON_MOD2)         $(COMMON_MOD3)
OBJECTS     = $(COMMON_MOD1:.f=.o) $(COMMON_MOD2:.for=.o) $(COMMON_MOD3:.f90=.o)
# f90: (ansi/iso or iso standard) modern fortran file
# for: (ansi/iso standard) punchcard fortran file 
# f  : (ibm language) punchcard fortran file

# # gnu compiler: 
# FORTRAN = gfortran
# FFLAGS += -fimplicit-none -fbounds-check
# FFLAGS1 = -std=legacy -x f77
# FFLAGS2 = -std=legacy
# # FFLAGS3 = 

# intel compiler: 
FORTRAN = ifort
# FFLAGS += 
FFLAGS1 = -nostand -f66 
FFLAGS2 = -nostand 
# FFLAGS3 = 

# # gnu library: 
# FFLAGS  += -I/usr/include
# LDFLAGS += -L/usr/lib
# LDFLAGS += -llapack -lblas

# # openblas library (mac homebrew): 
# FFLAGS  += -I/usr/local/opt/openblas/include
# LDFLAGS += -L/usr/local/opt/openblas/lib
# LDFLAGS += -lopenblas

# mkl library (mac): 
MKLROOT  = /opt/intel/mkl
FFLAGS  += -i8 -I$(MKLROOT)/include -I$(MKLROOT)/include/intel64/ilp64 
LDFLAGS += $(MKLROOT)/lib/libmkl_blas95_ilp64.a $(MKLROOT)/lib/libmkl_lapack95_ilp64.a
LDFLAGS += $(MKLROOT)/lib/libmkl_intel_ilp64.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_sequential.a -lpthread -lm

# # mkl library (linux): 
# MKLROOT  = /opt/intel/mkl
# FFLAGS  += -i8 -I$(MKLROOT)/include -I$(MKLROOT)/include/intel64/ilp64 
# LDFLAGS += $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_ilp64.a 
# LDFLAGS += -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm


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

# ${OBJECTS}: ${COMMON_MOD}
main.o: const.o global.o hamiltonian.o basis.o boundary.o inner.o outer.o
outer.o: const.o global.o hamiltonian.o 
inner.o: const.o global.o hamiltonian.o 
boundary.o: const.o global.o 
basis.o: const.o global.o hamiltonian.o 
hamiltonian.o: const.o global.o
global.o: const.o 
