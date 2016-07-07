#!/bin/bash

#if [ ! -f 1111.dat ]; then
	#cat chi_1111_w1* > temp_dat
	#LC_ALL=C sort -g -k1 -k2 -k3 temp_dat > 1111.dat
#fi

if [ ! -f 1010.dat ]; then
	cat chi_1010_w1* > temp_dat
	LC_ALL=C sort -g -k1 -k2 -k3 temp_dat > 1010.dat
fi

#if [ ! -f 1001.dat ]; then
	#cat chi_1001_w1* > temp_dat
	#LC_ALL=C sort -g -k1 -k2 -k3 temp_dat > 1001.dat
#fi

mv chi_*.dat dat/ 2> /dev/null
mv *.pom dat/ 2> /dev/null
rm temp_dat 2> /dev/null

exit 0
