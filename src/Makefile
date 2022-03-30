include ../makefile.inc

libraries = -L/usr/lib -llapack -lblas -lsymspg

serial_dependencies := spglib_f08 constants error tool type matrix schrodinger band_contribution \
                       spectral_function number_density hqt        
parallel_dependencies := cpu hqt_parallel

dependencies := $(serial_dependencies)
ifeq ($(mpi),true)
	dependencies += $(parallel_dependencies)
endif

target := hqt.o
ifeq ($(mpi),true)
	target := hqt_parallel.o
endif

.PHONEY : debug libhqt $(dependencies)

default : hqt

$(dependencies) :
	@$(compiler) ../src/$@.f90 $(compiler_flags) $(libraries) -c -o $@.o 

libhqt: $(dependencies)
	@ar crv ../libhqt.a $(target)
