#!/bin/bash

#---------------------------------------------------------------------------------------
#
#Copyright 2013 Arran Schlosberg.
#
#This file is part of https://github.com/aschlosberg/SNP (SNP)
#
#    SNP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SNP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SNP. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------------------

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
