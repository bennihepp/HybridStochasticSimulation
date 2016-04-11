SHELL = /bin/sh
FORTRAN = gfortran
FFLAGS  = -O2 -W -Wunused -Wuninitialized -fPIC -fno-strict-aliasing

SOURCES = $(wildcard src/*.f)
OBJECTS = $(SOURCES:.f=.o)
LOBJECTS = $(SOURCES:.f=.lo)

all: $(LOBJECTS) $(OBJECTS)

%.lo: %.f
	$(LIBTOOL) --mode=compile --tag=FC $(FORTRAN) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(LIBTOOL) --mode=compile --tag=FC $(FORTRAN) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(LOBJECTS) $(OBJECTS)

