#!/bin/bash

if [ -z $1 ]
then
	echo "Specify input file as argument";
	exit 1;
fi;

# 1. sed - convert species designation line to > only
# 2. tr - remove all new line characters (Windows + Unix)
# 3. cut - remove > from first species
# 4. tr - place species on new lines using > from (1)
< $1 sed 's|^>.*$|>|' | tr -d '\n\r' | cut -b 2- | tr '>' '\n'
