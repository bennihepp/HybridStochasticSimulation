SHELL = /bin/sh
CC    = gcc
JAVAH   = javah
CFLAGS = -O3 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
JDKDIR = $(JAVA_HOME)
JDKINCLUDEDIR = $(JDKDIR)/include
JDKPLATFORMDIR = $(JDKINCLUDEDIR)/$(PLATFORM)
LDFLAGS = -shared -lgfortran -fno-omit-frame-pointer
INCLUDES = -I../include -I./include -I$(JDKINCLUDEDIR) -I$(JDKPLATFORMDIR)

SOTARGET = lib/liblsodarjni.$(DLLEXT)
TARGET = lib/liblsodarjni.la
SOURCES = $(wildcard src/*.c)
SOURCES_DEPS = $(wildcard src/*.h)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.c=.o)
LOBJECTS = $(SOURCES:.c=.lo)
EXTERNAL_LOBJECTS = ../src/jni_utils.lo $(wildcard ../odepack/src/*.lo)
JAVACLASSPATH=../../target/classes
#JAVACLASS = ch.ethz.bhepp.ode.lsodar.LsodarSolver ch.ethz.bhepp.ode.lsodar.LsodarDirectSolver
#JAVACLASSFILE = $(JAVACLASSPATH)/ch/ethz/bhepp/ode/lsodar/LsodarSolver.class $(JAVACLASSPATH)/ch/ethz/bhepp/ode/lsodar/LsodarDirectSolver.class

all: $(SOTARGET)

$(SOTARGET): $(TARGET)
	ln -f -s .libs/liblsodarjni.$(DLLEXT) lib/liblsodarjni.$(DLLEXT)

$(TARGET): $(LOBJECTS) $(EXTERNAL_LOBJECTS) lsodar_odepack
	$(LIBTOOL) --tag=CC --mode=link $(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(LOBJECTS) $(EXTERNAL_LOBJECTS) \
		-rpath /usr/local/lib

lsodar_odepack: ../odepack/Makefile
	$(MAKE) -C ../odepack

$(LOBJECTS): $(HEADERS) $(SOURCES_DEPS)

jniheaders:

#jniheaders:
#	CLASSPATH=$(JAVACLASSPATH) $(JAVAH) -d include -jni $(JAVACLASS)

%.lo: %.c $(HEADERS) jniheaders
	$(LIBTOOL) --tag=CC --mode=compile $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

#%.o: %.c $(HEADERS)
#	$(LIBTOOL) --tag=CC --mode=compile $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean: cleanjni cleanodepack

cleanodepack:
	$(MAKE) -C ../odepack clean

cleanjni:
	rm -f $(SOTARGET) $(TARGET) $(LOBJECTS) $(OBJECTS)

