SHELL = /bin/sh
CC    = gcc
CXX   = g++
#LIBTOOL = glibtool
JAVAH   = javah
CFLAGS = -O3 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
#CFLAGS = -g -O0 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
CDEFS = 
JDKDIR = $(JAVA_HOME)
PLATFORM=darwin
#JDKDIR = /usr/lib/jvm/java-7-oracle
#JDKDIR = /usr/local/bsse/el6/java/jdk7
JDKINCLUDEDIR = $(JDKDIR)/include
JDKPLATFORMDIR = $(JDKINCLUDEDIR)/$(PLATFORM)
#JDKDIR = /c/Program\ Files/Java/jdk1.7.0_21/include
#JDKPLATFORMDIR = /c/Program\ Files/Java/jdk1.7.0_21/include/win32
LDFLAGS = -shared -fno-omit-frame-pointer
INCLUDES = -I../include -I$(JDKINCLUDEDIR) -I$(JDKPLATFORMDIR)

SOTARGET = lib/libboostjni.dylib
TARGET = lib/libboostjni.la
SOURCES = $(wildcard src/*.cpp)
SOURCES_DEPS = $(wildcard src/*.h)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.cpp=.o)
LOBJECTS = $(SOURCES:.cpp=.lo)
JAVACLASSPATH=../../target/classes
JAVACLASSPATH=../../bin
JAVACLASS = ch.ethz.bhepp.ode.boost.FixedBoostOdeintSolver
JAVACLASSFILE = $(JAVACLASSPATH)/ch/ethz/bhepp/ode/boost/FixedBoostOdeintSolver.class

all: $(SOTARGET)

$(SOTARGET): $(TARGET)
	ln -f -s .libs/libboostjni.dylib lib/libboostjni.dylib

$(TARGET): $(LOBJECTS)
	$(LIBTOOL) --tag=CXX --mode=link $(CXX) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(LOBJECTS) $(EXTERNAL_LOBJECTS) \
		-rpath /usr/local/lib

$(LOBJECTS): $(HEADERS) $(SOURCES_DEPS)

%.lo: %.cpp $(HEADERS) jniheaders
	$(LIBTOOL) --tag=CXX --mode=compile $(CXX) $(CFLAGS) $(CDEFS) $(INCLUDES) -c $< -o $@

%.lo: %.c $(HEADERS) jniheaders
	$(LIBTOOL) --tag=CC --mode=compile $(CC) $(CFLAGS) $(CDEFS) $(INCLUDES) -c $< -o $@

#%.o: %.c $(HEADERS)
#	$(LIBTOOL) --tag=CXX --mode=compile $(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(SOTARGET) $(TARGET) $(LOBJECTS) $(OBJECTS)

jniheaders:

#jniheaders: $(JAVACLASSFILE)
#	CLASSPATH=$(JAVACLASSPATH) $(JAVAH) -d include -jni $(JAVACLASS)

