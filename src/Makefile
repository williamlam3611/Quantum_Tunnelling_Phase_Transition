include ../make.inc

libraries := -L/usr/lib -llapack -lblas -lsymspg

modules := spglib_f08 constants error tool \
           type import matrix schrodinger  \
           band_contribution spectral_function     \
           number_density potential path

ifeq ($(mpi),true)
	modules += cpu 
	modules += hqt.parallel
else
	modules += hqt.serial
endif

.PHONEY : debug libhqt $(modules)

default : libhqt

$(modules) :
	@($(compiler) ../src/$@.f90 $(compiler_flags) $(libraries) -c)

libhqt: $(modules)
	@(ar cr ../libhqt.a *.o)
	@($(compiler) *.o $(compiler_flags) $(libraries) -shared -o ../libhqt.so)
	@(for module in *.mod; do                    \
		cp -p $$module ../$$module; \
	done)
