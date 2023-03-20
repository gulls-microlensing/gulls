#!/bin/bash

for i in l*.header; do
    rootname=${i%.header}
    echo $i $rootname
    if [ ! -e kgriz/$rootname.KgrizVRIw ]; then
	good=`grep 'TOTAL NUMBER OF STARS' $i | wc -l | awk '{if($1==1) print 1; else print 0}'`
	if [ $good -eq 1 ]; then
	    ./apply_extinction.pl $rootname
	    mv kgriz/$rootname.KgrizVRIw kgriz/$rootname.tmp
	    ./kgrizcoordconvert kgriz/$rootname.tmp > kgriz/$rootname.KgrizVRIw
	    rm kgriz/$rootname.tmp
	fi
    fi
done
