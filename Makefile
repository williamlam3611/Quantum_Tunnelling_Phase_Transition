DESTDIR ?= /usr/local

.PHONEY : debug install all libhqt clean

default : all

install: default
	@(install -d "$(DESTDIR)/lib/")
	@(install "libhqt.a" "$(DESTDIR)/lib/libhqt.a")
	@(install "libhqt.so" "$(DESTDIR)/lib/libhqt.so")
	
uninstall: default
	@(rm -f "$(DESTDIR)$(PREFIX)/lib/libhqt.a")
	@(rm -f "$(DESTDIR)$(PREFIX)/lib/libhqt.so")

all: libhqt
	
libhqt:
	@(cd bin && $(MAKE) --no-print-directory -f ../src/Makefile libhqt)
	
clean :
	@(cd bin && rm -f *.o *.mod *.swp)
	@(rm -f *.a *.so)
