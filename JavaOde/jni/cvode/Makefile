SHELL = /bin/sh
CC    = gcc
#LIBTOOL = glibtool
JAVAH   = javah
CFLAGS = -O2 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
CDEFS = 
JDKDIR = $(JAVA_HOME)
JDKINCLUDEDIR = $(JDKDIR)/include
JDKPLATFORMDIR = $(JDKINCLUDEDIR)/$(PLATFORM)
LDFLAGS = -shared -L$(SUNDIALS_HOME)/lib -lsundials_cvode -lsundials_nvecserial -fno-omit-frame-pointer
INCLUDES = -I../include -I$(SUNDIALS_HOME)/include -I$(JDKINCLUDEDIR) -I$(JDKPLATFORMDIR)

SOTARGET = lib/libcvodejni.$(DLLEXT)
TARGET = lib/libcvodejni.la
SOURCES = $(wildcard src/*.c)
SOURCES_DEPS = $(wildcard src/*.h)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.c=.o)
LOBJECTS = $(SOURCES:.c=.lo)
EXTERNAL_LOBJECTS = ../src/jni_utils.lo
JAVACLASSPATH=../../target/classes
#JAVACLASS = ch.ethz.bhepp.ode.cvode.CVodeSolver
#JAVACLASSFILE = $(JAVACLASSPATH)/ch/ethz/bhepp/ode/cvode/CVodeSolver.class

all: $(SOTARGET)

$(SOTARGET): $(TARGET)
	ln -f -s .libs/libcvodejni.$(DLLEXT) lib/libcvodejni.$(DLLEXT)

$(TARGET): $(LOBJECTS) $(EXTERNAL_LOBJECTS)
	$(LIBTOOL) --tag=CC --mode=link $(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(LOBJECTS) $(EXTERNAL_LOBJECTS) \
		-rpath /usr/local/lib

$(LOBJECTS): $(HEADERS) $(SOURCES_DEPS)

%.lo: %.c $(HEADERS) jniheaders
	$(LIBTOOL) --tag=CC --mode=compile $(CC) $(CFLAGS) $(CDEFS) $(INCLUDES) -c $< -o $@

#%.o: %.c $(HEADERS)
#	$(LIBTOOL) --tag=CC --mode=compile $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(SOTARGET) $(TARGET) $(LOBJECTS) $(OBJECTS)

jniheaders:

#jniheaders: $(JAVACLASSFILE)
#	CLASSPATH=$(JAVACLASSPATH) $(JAVAH) -d include -jni $(JAVACLASS)

