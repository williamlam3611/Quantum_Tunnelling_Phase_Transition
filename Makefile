DESTDIR ?= /usr/local

.PHONEY : debug install all libhqt clean

default : all

install:
	@(install -d "$(DESTDIR)/lib/")
	@(install libhqt.a "$(DESTDIR)/lib/libhqt.a")
	@(install libhqt.so "$(DESTDIR)/lib/libhqt.so")
	@(for module in hqt*.mod; do                    \
		install $$module "$(DESTDIR)/lib/$$module"; \
	done)
	@(echo [hqt] installed)
	
uninstall:
	@(rm -f "$(DESTDIR)/lib/libhqt.a")
	@(rm -f "$(DESTDIR)/lib/libhqt.so")
	@(cd "$(DESTDIR)/lib/" && \
		for module in hqt*.mod; do \
			rm -f "$$module";        \
		done)
	@(echo [hqt] uninstalled)

all: libhqt
	
libhqt:
	@(cd bin && $(MAKE) --no-print-directory -f ../src/Makefile libhqt)
	@(echo [hqt] made)
	
clean :
	@(cd bin && rm -f *.o *.mod *.swp)
	@(rm -f *.a *.so *.mod)
	@(echo [hqt] cleaned)
