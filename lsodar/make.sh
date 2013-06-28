#!/usr/bin/bash

g77-3 -fPIC -c lsodar.f opkda1.f opkda2.f
gcc-3 -fPIC -shared lsodar.o opkda1.o opkda2.o -lg2c -o libtest.so
#gcc -fPIC -c test.c
gcc-3 -L`pwd` -o test test.c -ltest -lg2c

gcc-3 -c -D_REENTRANT -fPIC LsodarOdeSolverJNI.c -I/cygdrive/c/Program\ Files/Java/jdk1.7.0_21/include/ -I/cygdrive/c/Program\ Files/Java/jdk1.7.0_21/include/win32
#gcc -fPIC -shared LsodarOdeSolverJNI.o -o liblsodarjni.dll

#gcc -Wall -D_JNI_IMPLEMENTATION_ -Wl,--kill-at  -mno-cygwin -I/cygdrive/c/Program\ Files/Java/jdk1.7.0_21/include -I/cygdrive/c/Program\ Files/Java/jdk1.7.0_21/include/win32 -shared LsodarOdeSolverJNI.c -o liblsodarjni.dll

gcc-3 -fPIC -shared lsodar.o opkda1.o opkda2.o -lg2c -o libtest.so

