SHELL = /bin/sh
CC    = gcc
LIBTOOL = libtool
JAVAH   = javah
CFLAGS = -O2 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
CDEFS = 
JDKDIR = $(JAVA_HOME)
#JDKDIR = /usr/lib/jvm/java-7-oracle
#JDKDIR = /usr/local/bsse/el6/java/jdk7
JDKINCLUDEDIR = $(JDKDIR)/include
JDKPLATFORMDIR = $(JDKINCLUDEDIR)/linux
#JDKDIR = /c/Program\ Files/Java/jdk1.7.0_21/include
#JDKPLATFORMDIR = /c/Program\ Files/Java/jdk1.7.0_21/include/win32
LDFLAGS = -shared -L$(SUNDIALS_HOME)/lib -lsundials_cvode -lsundials_nvecserial -lgfortran -fno-omit-frame-pointer
INCLUDES = -I../include -I./include/ -I$(SUNDIALS_HOME)/include -I$(JDKINCLUDEDIR) -I$(JDKPLATFORMDIR)

SOTARGET = lib/libcvodejni.so
TARGET = lib/libcvodejni.la
SOURCES = $(wildcard src/*.c)
SOURCES_DEPS = $(wildcard src/*.h)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.c=.o)
LOBJECTS = $(SOURCES:.c=.lo)
EXTERNAL_LOBJECTS = ../src/jni_utils.lo
JAVACLASSPATH=../../target/classes
JAVACLASSPATH=../../bin
JAVACLASS = ch.ethz.bhepp.ode.cvode.CVodeSolver
JAVACLASSFILE = $(JAVACLASSPATH)/ch/ethz/bhepp/ode/cvode/CVodeSolver.class

all: $(SOTARGET)

$(SOTARGET): $(TARGET)
	ln -f -s .libs/libcvodejni.so lib/libcvodejni.so

$(TARGET): $(LOBJECTS) $(EXTERNAL_LOBJECTS)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(LOBJECTS) $(EXTERNAL_LOBJECTS) \
		-rpath /usr/local/lib

$(LOBJECTS): $(HEADERS) $(SOURCES_DEPS)

%.lo: %.c $(HEADERS) jniheaders
	$(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(CDEFS) $(INCLUDES) -c $< -o $@

#%.o: %.c $(HEADERS)
#	$(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(SOTARGET) $(TARGET) $(LOBJECTS) $(OBJECTS)

jniheaders: $(JAVACLASSFILE)
	CLASSPATH=$(JAVACLASSPATH) $(JAVAH) -d include -jni $(JAVACLASS)

