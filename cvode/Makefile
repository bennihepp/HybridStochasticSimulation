SHELL = /bin/sh
CC    = gcc
CFLAGS = -fno-strict-aliasing -Wall -fPIC -I/usr/local/include/
LDFLAGS = -fno-omit-frame-pointer
CVODELIBS = /usr/local/lib/libsundials_cvode.la /usr/local/lib/libsundials_nvecserial.la -lm
INCLUDES = 
LIBTOOL = libtool

TARGET = pdmp_example
#SOURCES = $(wildcard src/*.c)
SOURCES = src/prob_oscillator.c src/pdmp_example.c
SOURCES_DEPS = $(wildcard src/*.h)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(SOURCES:.c=.o)
LOBJECTS = $(SOURCES:.c=.lo)

all: $(TARGET)

$(TARGET): $(LOBJECTS)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) $(LDFLAGS) $(LOBJECTS) -o $(TARGET) $(CVODELIBS)

$(LOBJECTS): $(HEADERS) $(SOURCES_DEPS)

%.lo: %.c
	$(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS) $(LOBJECTS)

