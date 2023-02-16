#!/bin/zsh

# This script takes the output of Jfilament tracking as its argument and coverts it into a CSV file with two fileds:
#	--> frame : frame number in time series
#	--> length : filament length in pixels

# output file name
name=`echo $1 | gsed 's/txt$/csv/g'`

# output data as CSV
tail -n +14 $1 | tr '\t' ',' | cut -d',' -f1,2 | gsed '1i\frame,length' > $name

# Ankit Roy
# 8th December, 2022
