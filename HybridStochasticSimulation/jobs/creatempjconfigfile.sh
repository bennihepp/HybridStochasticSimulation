#!/bin/bash

CONFIGFILE=$1
if [ -f "$CONFIGFILE" ]; then
	rm $CONFIGFILE
fi
touch $CONFIGFILE

hostlist=""
ncpuslist=""
total_ncpus=0
listlength=0

host=""
for entry in $LSB_MCPU_HOSTS; do
	#echo "entry: $entry"
	if [ "$host" == "" ]; then
		host=$entry
		#echo "host: $host"
	else
		ncpus=$entry
		#echo "ncpus: $ncpus"
		#echo "listlength: $listlength"
		hostlist[$listlength]=$host
		ncpuslist[$listlength]=$ncpus
		total_ncpus=$(( total_ncpus + ncpus ))
		listlength=$(( listlength + 1 ))
		host=""
	fi
done

#echo "hostlist: ${hostlist[1]}"
echo "ncpuslist: ${ncpuslist[1]}"

echo "$total_ncpus" >> $CONFIGFILE
echo "131072" >> $CONFIGFILE

rank=0
for i in `seq 0 $(( listlength - 1 ))`; do
	host=${hostlist[$i]}
	ncpus=${ncpuslist[$i]}
	port=20000
	for j in `seq 1 $ncpus`; do
		echo "$host@$port@$rank" >> $CONFIGFILE
		rank=$(( rank + 1 ))
		port=$(( port + 2 ))
	done
done
echo $total_ncpus

