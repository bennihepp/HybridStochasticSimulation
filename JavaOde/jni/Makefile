SHELL = /bin/sh
CC    = gcc
CFLAGS = -O2 -fno-strict-aliasing -D_JNI_IMPLEMENTATION_ -Wall -fPIC
CDEFS = 
JDKDIR = $(JAVA_HOME)
PLATFORM=darwin
JDKDIR = $(JAVA_HOME)
JDKINCLUDEDIR = $(JDKDIR)/include
JDKPLATFORMDIR = $(JDKINCLUDEDIR)/$(PLATFORM)
INCLUDES = -I./include/ -I$(JDKINCLUDEDIR) -I$(JDKPLATFORMDIR)

SOURCES = $(wildcard src/*.c)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.c=.o)
LOBJECTS = $(SOURCES:.c=.lo)

TARGET1 = lib/liblsodarjni.$(DLLEXT)
TARGET2 = lib/libcvodejni.$(DLLEXT)
TARGET3 = lib/libboostjni.$(DLLEXT)

.PHONY: liblsodar
liblsodar: $(TARGET1)

.PHONY: libcvode
libcvode: $(TARGET2)

.PHONY: libboostjni
libboostjni: $(TARGET3)

.PHONY: all
all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): lsodar
	ln -f -s ../lsodar/lib/liblsodarjni.$(DLLEXT) lib/liblsodarjni.$(DLLEXT)

$(TARGET2): cvode
	ln -f -s ../cvode/lib/libcvodejni.$(DLLEXT) lib/libcvodejni.$(DLLEXT)

$(TARGET3): boost
	ln -f -s ../boost/lib/libboostjni.$(DLLEXT) lib/libboostjni.$(DLLEXT)

lsodar: $(LOBJECTS)
	$(MAKE) -C lsodar/

cvode: $(LOBJECTS)
	$(MAKE) -C cvode/

boost: $(LOBJECTS)
	$(MAKE) -C boost/

clean: lsodar_clean cvode_clean boost_clean
	rm -f src/jni_utils.o src/jni_utils.lo

lsodar_clean:
	$(MAKE) -C lsodar/ clean
	rm -f lib/liblsodarjni.$(DLLEXT)

cvode_clean:
	$(MAKE) -C cvode/ clean
	rm -f lib/libcvodejni.$(DLLEXT)

boost_clean:
	$(MAKE) -C boost/ clean
	rm -f lib/libboostjni.$(DLLEXT)

$(LOBJECTS): $(HEADERS) $(SOURCES)

%.lo: %.c $(HEADERS)
	$(LIBTOOL) --tag=CC --mode=compile $(CC) $(CFLAGS) $(CDEFS) $(INCLUDES) -c $< -o $@

