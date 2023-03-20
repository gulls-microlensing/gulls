#!/bin/bash

#This script performs all the housekeeping type setup required to run an gulls simulation
#It makes sure that the directory structure that gulls expects for output exists

#Load the common preamble
source ${SCRIPTDIR}gullsPreamble.sh
#set location of executables, must be kept in this dir? Or give full dir in 
export SCRIPTDIR=/Users/penny/gulls/scripts/
#export SRCDIR='/home/stingray/johnson.7080/gulls/gulls_sj/bin/'

#echo first line
#Remove any files remaining on compute nodes
#ssh node52 "rm -r $outputdir$runname/" 
#ssh node53 "rm -r $outputdir$runname/" 
#ssh node54 "rm -r $outputdir$runname/" 
#ssh node55 "rm -r $outputdir$runname/" 

#Setup output directory structure

#echo $finaldir, $SCRIPTDIR

#Make sure there is an output root directory
if [ ! -d $finaldir ]; then
   
   mkdir $finaldir
   echo Created $finaldir
fi

#echo Hello

if [ ! -d $finaldir$runname ]; then
    mkdir $finaldir$runname
    echo Created $finaldir$runname
fi

#Make sure there are output subdirectories
for rt in $finaldir; do
    for direc in raw logs lc analysis fm; do 
	diraim=$rt$runname/$direc
	if [ ! -d $diraim ]; then
	    mkdir $diraim
	    echo Created $diraim
	fi
    done
done

#Copy across the parameter files that were used

#the observatory list
cp -i -n $obslist $finaldir$runname/logs/

grep -v '^#' $obslist | grep -v '^$' | while read obs; do
    #each observatory
    cp -i -n $obsdir$obs $finaldir$runname/logs/
    #each observatory's weather profile
    cp -i -n $weatherdir`grep WEATHER_PROFILE $obsdir$obs | awk '{print $2}'` $finaldir$runname/logs/
    #each observatory's field list
    cp -i -n $obsdir`grep FIELDCENTRES $obsdir$obs | awk '{print $2}'` $finaldir$runname/logs/
    #each observatory's observation sequence list
    cp -i -n $obsdir`grep OBSERVATION_SEQUENCE $obsdir$obs | awk '{print $2}'` $finaldir$runname/logs/
    #each observatory's detector parameters
    cp -i -n $obsdir`grep DETECTOR $obsdir$obs | awk '{print $2}'` $finaldir$runname/logs/
    #Executable, path in paramfile needs to be relative to GULLS_BASE_DIR
    cp -i -n $SRCDIR$executable $finaldir$runname/logs/
    
done

#the parameter file
#echo 'Here'
#echo $paramdir$paramfile
#echo $finaldir$runname/logs/
#original line
#cp -i --reply=no $paramdir$paramfile $finaldir$runname/logs/
#cp -i --reply=no $base_dir$param_file_dir$paramfile $finaldir$runname/logs/
#had to remove --reply for an -f, could be dangerous.
cp -i -n  $base_dir$param_file_dir$paramfile $finaldir$runname/logs/
