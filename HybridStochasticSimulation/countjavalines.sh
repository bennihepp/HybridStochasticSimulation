#!/bin/bash

BASEDIR=$1
LS=/bin/ls

DIRS=`/bin/ls "$BASEDIR"`

L=0
for d in "$DIRS"; do
	PATH="$BASEDIR/$d"
	if [[ "." == "$d" ]]; then
		echo "found ."
	elif [[ ".." == "$d" ]]; then
		echo "found .."
	elif [[ -d "$PATH" ]]; then
		#echo entering "$PATH"
		q=`/bin/bash countjavalines.sh $PATH`
		L=$((L + q))
	elif [[ -f "$PATH" ]]; then
		if [[ "$PATH" =~ ".*\.java$" ]]; then
			WCOUT=`/bin/wc -l "$PATH"`
			for w in $WCOUT; do
				L=$((L+w))
				break
			done
		fi
	fi
done
echo $L

