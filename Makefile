compiler = mpifort
compiler_flags = -O3 -fbackslash -fopenmp
compiler_flags_debug = -Og -fbacktrace -fbackslash -fopenmp

libraries = -L/usr/lib -llapack -lblas -lsymspg
modules = spglib_f08.f90 cpu.f90 timer.f90 import.f90 route.f90 potential.f90 transform.f90 bulk.f90 matrix.f90 sort.f90 spectral.f90 density.f90 export.f90 graph.f90 hqtlib.f90
output = a.out

.PHONEY : main clean run

main :
	$(compiler) $(compiler_flags) $(modules) main.f90 -o $(output) $(libraries)
	
poisson : 
	$(compiler) $(compiler_flags) $(modules) main_poisson.f90 -o $(output) $(libraries)

run : main
	mpirun -np 4 ./a.out
	make clean

run_poisson : poisson
	mpirun -np 4 ./a.out
	make clean

debug : 
	$(compiler) $(compiler_flags_debug) $(modules) $(target) -o $(output) $(libraries)
	mpirun -np 4 ./a.out
	make clean
	
clean :
	rm *.o *.mod
