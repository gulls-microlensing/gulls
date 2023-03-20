#!/bin/bash

export SCRIPTDIR=/Users/penny/gulls/scripts/
export SCRIPTDIR=/Users/penny/gulls/scripts/

if [ $# -lt 1 ]; then
  echo "Usage:"
  echo "./gullsRun.sh <parameterFile> {-subruns=<list> {-subruns=<list> ... } }"
  exit
fi

prm=$1
shortprm=${prm##*/}

if [ `grep RUN_NAME $prm | awk -v FS='=' '{print $2}'` != "${shortprm%.*}" ]; then
    echo "Warning parameter file name does not match run name. Are you sure you want to continue?"
    echo "$prm"
    grep RUNNAME $prm
    read tmp
fi

#setup
#echo Setup run...
#${SCRIPTDIR}gullsSetupRun.sh $prm

#run
echo Launch run...
${SCRIPTDIR}gullsLaunch.sh $@

#tidy up
#${SCRIPTDIR}gullsTidy.sh $prm &

#analyze
#${SCRIPTDIR}gullsAnalyze.pl $prm
