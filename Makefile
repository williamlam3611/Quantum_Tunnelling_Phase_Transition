compiler = mpifort
compiler_flags = -O3 -fbackslash -fopenmp
compiler_flags_debug = -Og -fbacktrace -fbackslash -fopenmp

libraries = -L/usr/lib -llapack -lblas -lsymspg
modules = spglib_f08.f90 cpu.f90 timer.f90 import.f90 route.f90 potential.f90 transform.f90 bulk.f90 matrix.f90 sort.f90 spectral.f90 density.f90 export.f90 graph.f90 hqtlib.f90
target = main.f90
output = a.out

.PHONEY : debug clean run

default :
	@$(compiler) $(compiler_flags) $(modules) $(target) -o $(output) $(libraries)
	
run : default
	@mpirun -np 4 ./a.out
	@make clean --no-print-directory
	
debug : 
	@$(compiler) $(compiler_flags_debug) $(modules) $(target) -o $(output) $(libraries)
	@mpirun -np 4 ./a.out
	@make clean --no-print-directory
	
clean :
	@rm --recursive --force *.o *.mod *.swp || true
