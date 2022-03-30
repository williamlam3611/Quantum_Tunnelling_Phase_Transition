PREFIX ?= /usr/local

.PHONEY : debug install all libhqt clean

default : all

install: default
	@(install -d "$(DESTDIR)$(PREFIX)/lib/")
	@(install "libhqt.a" "$(DESTDIR)$(PREFIX)/lib/libhqt.a")
	
uninstall: default
	@(rm -f "$(DESTDIR)$(PREFIX)/lib/libhqt.a")

all: libhqt
	
libhqt:
	@(cd bin && $(MAKE) -f ../src/Makefile libhqt)
	
clean :
	@(cd bin && rm -f *.o *.mod *.swp)
	@(rm -f *.a *.so)
