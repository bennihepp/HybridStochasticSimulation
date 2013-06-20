SHELL = /bin/sh
CC    = gcc
LIBTOOL = libtool
CFLAGS = -O2 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
JDKDIR = /usr/lib/jvm/java-7-oracle/include/
JDKPLATFORMDIR = /usr/lib/jvm/java-7-oracle/include/linux
#JDKDIR = /c/Program\ Files/Java/jdk1.7.0_21/include
#JDKPLATFORMDIR = /c/Program\ Files/Java/jdk1.7.0_21/include/win32
LDFLAGS = -shared -lgfortran -fno-omit-frame-pointer
INCLUDES = -I./include/ -I$(JDKDIR) -I$(JDKPLATFORMDIR)

SOTARGET = lib/liblsodarjni.so
TARGET = lib/liblsodarjni.la
SOURCES = $(wildcard src/*.c)
SOURCES_DEPS = $(wildcard src/*.h)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.c=.o)
LOBJECTS = $(SOURCES:.c=.lo)
LSODAR_LOBJECTS = $(wildcard ../odepack/src/*.lo)
JAVACLASS = ch.ethz.khammash.nativeode.lsodar.LsodarSolver
JAVACLASSFILE = ../bin/ch/ethz/khammash/nativeode/lsodar/LsodarSolver.class

all: $(SOTARGET)

$(SOTARGET): $(TARGET)
	ln -f -s .libs/liblsodarjni.so lib/liblsodarjni.so

$(TARGET): $(LOBJECTS) lsodar_odepack
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(LOBJECTS) $(LSODAR_LOBJECTS) \
		-rpath /usr/local/lib

lsodar_odepack: ../odepack/Makefile
	$(MAKE) -C ../odepack

$(LOBJECTS): $(HEADERS) $(SOURCES_DEPS)

$(HEADERS): jniheaders

jniheaders: $(JAVACLASSFILE)
	CLASSPATH=../bin javah -d include -jni $(JAVACLASS)

%.lo: %.c $(HEADERS)
	$(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.c $(HEADERS)
	$(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean: cleanjni cleanodepack

cleanodepack:
	$(MAKE) -C ../odepack clean

cleanjni:
	rm -f $(SOTARGET) $(TARGET) $(LOBJECTS) $(OBJECTS)

