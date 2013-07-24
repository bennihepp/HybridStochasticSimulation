#!/bin/sh

JAVA_PREFIX=$HOME/java/lib
LOCAL_PREFIX=$HOME/local

export CLASSPATH="../bin:$CLASSPATH"
export CLASSPATH="../../JavaOde/bin:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/commons-configuration-1.9.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/commons-lang-2.6.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/commons-logging-1.1.3.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/commons-lang3-3.1.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/commons-math3-3.2.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/guava-14.0.1.jar:$CLASSPATH"
#export CLASSPATH="$JAVA_PREFIX/JaCoP-3.2.jar:$CLASSPATH"
#export CLASSPATH="$JAVA_PREFIX/Jama-1.0.1.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/jamtio.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/ejml-0.23.jar:$CLASSPATH"
#export CLASSPATH="$JAVA_PREFIX/jcommon-1.0.17.jar:$CLASSPATH"
#export CLASSPATH="$JAVA_PREFIX/jfreechart-1.0.14.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/jgrapht-jdk1.6.jar:$CLASSPATH"
#export CLASSPATH="$JAVA_PREFIX/matlabcontrol-4.1.0.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/aopalliance.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/guice-3.0.jar:$CLASSPATH"
export CLASSPATH="$JAVA_PREFIX/javax.inject.jar:$CLASSPATH"

export LD_LIBRARY_PATH=../../JavaOde/jni/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LOCAL_PREFIX/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/4.7.3/:$LD_LIBRARY_PATH
#echo "$CLASSPATH"
#echo "$LD_LIBRARY_PATH"

java ch.ethz.khammash.hybridstochasticsimulation.Main $1

