#!/bin/bash

if [ "$3" = "" ]; then
	HOURS=94
else
	HOURS=$3
fi

CONFIGFILE=$1
JOBNAME=$2

CONFIGFILENAME=${CONFIGFILE##*/}
CONFIGBASENAME=${CONFIGFILENAME%.xml}
STDOUTFILENAME=${CONFIGBASENAME}.out
STDERRFILENAME=${CONFIGBASENAME}.err
STDOUTPATH=logs/new/$STDOUTFILENAME
STDERRPATH=logs/new/$STDERRFILENAME

TOTALNCPUS=16

bsub -oo $STDOUTPATH -eo $STDERRPATH -J $JOBNAME -R "span[ptile=${TOTALNCPUS}]" -n ${TOTALNCPUS} -W $HOURS:00 bash runjob.sh $CONFIGFILE

