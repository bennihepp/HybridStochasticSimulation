#!/bin/bash

BASEDIR=$1
LS=/bin/ls

DIRS=`/bin/ls "$BASEDIR"`

L=0
for d in $DIRS; do
	PATH="$BASEDIR/$d"
    #echo "checking $PATH"
	if [[ "." == "$d" ]]; then
		#echo "found ."
        L=$L
	elif [[ ".." == "$d" ]]; then
		#echo "found .."
        L=$L
	elif [[ -d "$PATH" ]]; then
		#echo entering "$PATH"
		q=`/bin/bash countjavalines.sh "$PATH"`
        echo "q=$q"
		#/bin/bash countjavalines.sh $PATH
        w=1
		L=$((L + w))
	elif [[ -f "$PATH" ]]; then
        #echo "found a file"
		if [[ "$PATH" =~ ".*\.java$" ]]; then
			WCOUT=`/bin/wc -l "$PATH"`
            echo "wcout=$WCOUT"
			for w in $WCOUT; do
				L=$((L+w))
				break
			done
		fi
    else
        #echo "neither file nor directory"
        L=$L
	fi
done
echo $L

