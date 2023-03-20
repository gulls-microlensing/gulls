#!/bin/bash

#export SCRIPTDIR=/home/stingray/johnson.7080/gulls/gulls_sj/scripts/
export SCRIPTDIR=/Users/penny/gulls/scripts/

if [ $# -ne 1 ]; then
  echo "Usage:"
  echo "./gullsRun.sh <parameterFile>"
  exit
fi

#get the env variable for the base_dir                                                                                                                                                                    
base_dir=$GULLS_BASE_DIR
# the name of the paramfile is an arguement                                                                                                                                                               
paramfile=$1
# right now defaulting to have a parameterFiles dir in base_dir, may not be the best but working withit for now                                                                                           
param_file_dir="parameterFiles/"


prm=$1
shortprm=${prm##*/}

if [ `grep RUN_NAME $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'` != "${shortprm%.*}" ]; then
    echo "Warning parameter file name does not match run name. Are you sure you want to continue?"
    echo "$prm"
    grep RUNNAME $prm
    read tmp
fi

#setup
#echo $SCRIPTDIR
#echo Here 
echo Setup run...
${SCRIPTDIR}gullsSetupRun.sh $prm

#run
echo Launch run...
${SCRIPTDIR}gullsLaunch.sh $prm

#tidy up
#${SCRIPTDIR}gullsTidy.sh $prm &

#analyze
#${SCRIPTDIR}gullsAnalyze.pl $prm
