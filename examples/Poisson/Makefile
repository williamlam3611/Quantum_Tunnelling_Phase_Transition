include ../../make.inc

libraries := -lhqt

default : all

all : poisson

poisson : 
	@(cd bin && $(compiler) $(compiler_flags) ../poisson.f90 -o ../poisson.out $(libraries))
	
clean :
	@(cd bin && rm -f *.o *.mod *.swp)
	@(rm -f *.a *.so *.mod *.out)
