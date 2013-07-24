#!/bin/bash

MACHINESFILE=$1
if [ -f "$MACHINESFILE" ]; then
	rm $MACHINESFILE
fi
touch $MACHINESFILE

#for host in "$LSB_HOSTS"; do
#	echo "$host" >> $MACHINESFILE
#done

total_ncpus=0
host=""
for entry in $LSB_MCPU_HOSTS; do
	#echo "entry: $entry"
	if [ "$host" == "" ]; then
		host=$entry
		#echo "host: $host"
	else
		ncpus=$entry
		#echo "ncpus: $ncpus"
		for i in `seq 1 $ncpus`; do
			#echo "i: $i"
			echo "$host" >> $MACHINESFILE
		done
		total_ncpus=$(( total_ncpus + ncpus ))
		host=""
	fi
done
echo $total_ncpus

