SHELL = /bin/sh
CC    = gcc
LIBTOOL = libtool
CFLAGS = -O2 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
CDEFS = 
JDKDIR = $(JAVA_HOME)
#JDKDIR = /usr/lib/jvm/java-7-oracle
#JDKDIR = /usr/local/bsse/el6/java/jdk7
JDKINCLUDEDIR = $(JDKDIR)/include
JDKPLATFORMDIR = $(JDKINCLUDEDIR)/linux
#JDKDIR = /c/Program\ Files/Java/jdk1.7.0_21/include
#JDKPLATFORMDIR = /c/Program\ Files/Java/jdk1.7.0_21/include/win32
INCLUDES = -I./include/ -I$(JDKINCLUDEDIR) -I$(JDKPLATFORMDIR)

SOURCES = $(wildcard src/*.c)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.c=.o)
LOBJECTS = $(SOURCES:.c=.lo)

TARGET1 = lib/liblsodarjni.so
TARGET2 = lib/libcvodejni.so

all: $(TARGET1) $(TARGET2)

liblsodar: $(TARGET1)

libcvode: $(TARGET2)

$(TARGET1): lsodar
	ln -f -s ../lsodar/lib/liblsodarjni.so lib/liblsodarjni.so

$(TARGET2): cvode
	ln -f -s ../cvode/lib/libcvodejni.so lib/libcvodejni.so

lsodar: $(LOBJECTS)
	$(MAKE) -C lsodar/

cvode: $(LOBJECTS)
	$(MAKE) -C cvode/

clean: lsodar_clean cvode_clean
	rm -f src/jni_utils.o src/jni_utils.lo

lsodar_clean:
	$(MAKE) -C lsodar/ clean

cvode_clean:
	$(MAKE) -C cvode/ clean

$(LOBJECTS): $(HEADERS) $(SOURCES)

%.lo: %.c $(HEADERS)
	$(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(CDEFS) $(INCLUDES) -c $< -o $@

