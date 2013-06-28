SHELL = /bin/sh
FORTRAN = gfortran
FFLAGS  = -W -Wunused -Wuninitialized -fPIC -fno-strict-aliasing

SOURCES = $(wildcard src/*.f)
OBJECTS = $(SOURCES:.f=.o)

all: $(OBJECTS)
	gcc -fPIC -shared src/lsodar.o src/opkda1.o src/opkda2.o -lgfortran -o libtest.so
	gcc -L`pwd` -o test test.c -ltest -lgfortran

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS)

