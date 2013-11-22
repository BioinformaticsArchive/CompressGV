#!/bin/bash

v="data/outcome/$1.verbose";
l="data/outcome/$1";

del="data/humvar/deleterious/";
neut="data/humvar/neutral/";

if [ -f "${l}.log" ]
then
	echo "$1 - DONE";
	exit 0;
fi

if [ -f "${l}.err" ]
then
	echo "$1 - DONE - ERROR";
	exit 0;
fi

if [ -z $2 ]
then
	d=`< "${del}$1" wc -l`;
	if [ $d -lt 3 ]
	then
		echo "Skipping $1 - only $d deleterious";
		exit 0;
	fi
	n=`< "${del}$1" wc -l`;
	if [ $n -lt 3 ]
	then
		echo "Skipping $1 - only $n neutral";
		exit 0;
	fi
fi

../bin/grantham data/alignments/combined/$1.aln "${del}$1" "${neut}$1" /dev/null 2>"${v}.tmp" > "${l}.tmp";

if [ $? -eq 0 ]
then
	echo "$1 - OK";
	mv "${v}.tmp" "${v}.log";
	mv "${l}.tmp" "${l}.log";
else
	echo "$1 - FAIL";
	mv "${v}.tmp" "${v}.err";
	mv "${l}.tmp" "${l}.err";
fi
