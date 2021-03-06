#!/bin/sh

export CLASSPATH="./bin:$CLASSPATH"
export CLASSPATH="../JavaOde/bin:$CLASSPATH"
export CLASSPATH="lib/commons-lang3-3.1.jar:$CLASSPATH"
export CLASSPATH="lib/commons-lang-2.6.jar:$CLASSPATH"
export CLASSPATH="lib/commons-configuration-1.9.jar:$CLASSPATH"
export CLASSPATH="lib/commons-logging-1.1.3.jar:$CLASSPATH"
export CLASSPATH="lib/commons-math3-3.2.jar:$CLASSPATH"
export CLASSPATH="lib/guava-14.0.1.jar:$CLASSPATH"
export CLASSPATH="lib/JaCoP-3.2.jar:$CLASSPATH"
export CLASSPATH="lib/Jama-1.0.1.jar:$CLASSPATH"
export CLASSPATH="lib/jamtio.jar:$CLASSPATH"
export CLASSPATH="lib/jcommon-1.0.17.jar:$CLASSPATH"
export CLASSPATH="lib/jfreechart-1.0.14.jar:$CLASSPATH"
export CLASSPATH="lib/jgrapht-jdk1.6.jar:$CLASSPATH"
export CLASSPATH="lib/matlabcontrol-4.1.0.jar:$CLASSPATH"

export LD_LIBRARY_PATH=../JavaOde/jni/lib
#export LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/4.7.3/:$LD_LIBRARY_PATH
#echo "$CLASSPATH"
#echo "$LD_LIBRARY_PATH"

java ch.ethz.khammash.hybridstochasticsimulation.Main

