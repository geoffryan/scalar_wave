MAKEFILE_IN = /Users/geoff/Coding/CAL/scalar_wave/Makefile.in
include $(MAKEFILE_IN)

CFLAGS ?= -Wall
CURL ?= curl
UNTAR ?= tar -xvf
CD ?= cd
RM ?= rm -f
OS ?= MacOSX10.6.8

LIBS = -lhdf5 -L$(HDF_HOME)/lib
LOCLIBS = cow/libcow.a

ALL = scalar_wave1d scalar_wave2d scalar_wave3d

default: $(ALL)

scalar_wave1d: scalar_wave1d.c $(LOCLIBS)
	$(CC) $(CFLAGS) -o $@ $^ -I./cow/ $(LIBS)
	
scalar_wave2d: scalar_wave2d.c $(LOCLIBS)
	$(CC) $(CFLAGS) -o $@ $^ -I./cow/ $(LIBS)
	
scalar_wave3d: scalar_wave3d.c $(LOCLIBS)
	$(CC) $(CFLAGS) -o $@ $^ -I./cow/ $(LIBS)
	
cow/libcow.a: .FORCE
	$(MAKE) -C cow libcow.a MAKEFILE_IN=$(MAKEFILE_IN)
	
clean:
	$(MAKE) -C cow clean MAKEFILE_IN=$(MAKEFILE_IN)
	$(RM) scalar_wave
	
.FORCE:

