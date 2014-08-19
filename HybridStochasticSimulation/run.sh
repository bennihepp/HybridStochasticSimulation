#!/bin/sh

export CLASSPATH="./target/classes:$CLASSPATH"
export CLASSPATH="../JavaOde/target/classes:$CLASSPATH"
export CLASSPATH="../lib/commons-lang3-3.1.jar:$CLASSPATH"
export CLASSPATH="../lib/commons-math3-3.3.jar:$CLASSPATH"
export CLASSPATH="../lib/guava-14.0.1.jar:$CLASSPATH"
export CLASSPATH="../lib/javax.inject.jar:$CLASSPATH"
export CLASSPATH="../lib/slf4j-api-1.7.5.jar:$CLASSPATH"
export CLASSPATH="../lib/slf4j-jdk14-1.7.5.jar:$CLASSPATH"
export CLASSPATH="../lib/jcommon-1.0.17.jar:$CLASSPATH"
export CLASSPATH="../lib/jfreechart-1.0.14.jar:$CLASSPATH"
export CLASSPATH="../lib/jgrapht-jdk1.6.jar:$CLASSPATH"
export CLASSPATH="../lib/matlabcontrol-4.1.0.jar:$CLASSPATH"
export CLASSPATH="../lib/jmatio.jar:$CLASSPATH"

export DYLD_LIBRARY_PATH=../JavaOde/jni/lib
export LD_LIBRARY_PATH=../JavaOde/jni/lib

java ch.ethz.bhepp.hybridstochasticsimulation.GUI

