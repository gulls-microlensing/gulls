#!/bin/bash

#This script creates scripts to launch each sub-run of the simulation, and sets them running

if [ $# -eq 0 ]; then
	echo "Usage: ./gullsLaunch.sh <parameter file> {-fields=<list>} {-subruns=<list>}"
        exit

fi

#Load the common preamble
source ${SCRIPTDIR}gullsPreamble.sh

#if [ -e ~/.QUEUE ]; then
#	queue=`head -n1 ~/.QUEUE | awk '{print $1}'`
#else
#	queue=mfun
#fi

#fieldlist="for i in {0..168}; do if [ `echo $i | awk '{if(int($1/13)>=11 || ($1%13)<=1 || ($1%13)==12) print 1; else print 0}'` -eq 0 ]; then echo $i fi done"


fieldlist=`for i in {0..168}; do echo -n "$i "; done`
subrunlist=`for ((i=0;i<$nsubruns;i++)); do echo -n "$i "; done`

fcommands=0  #Number of field commands issued
scommands=0 #Number of subrun commands issued

if [ $# -gt 1 ]; then
        args=( "$@" )
        for ((i=1;i<$#;i++)); do
                if [[ ${args[$i]} =~ "-fields=" ]]; then
                        if [ $fcommands -eq 0 ]; then
                                fieldlist=${args[$i]#-fields=}
                                fcommands=1
                        else
                                fieldlist="$fieldlist ${args[$i]#-fields=}"
                        fi
                fi
                if [[ ${args[$i]} =~ "-subruns=" ]]; then
                        if [ $scommands -eq 0 ]; then
				subrunlist=${args[$i]#-subruns=}
                                scommands=1
                        else
                                subrunlist="$subrunlist ${args[$i]#-subruns=}"
                        fi
                fi
        done
fi

if [ $# -eq 0 ] || [ $# -gt 1 ] && [ $scommands -eq 0 ] && [ $fcommands -eq 0 ]; then
	echo "Usage: ./gullsLaunch.sh <parameter file> {-fields=<list>} {-subruns=<list>}"
	exit
fi

#startrun=0
#endrun=$nsubruns

#if [ $# -eq 2 ]; then
#  startrun=$2;
#  endrun=`echo $2 | awk '{print $1+1}'`;
#fi

#if [ $# -eq 3 ]; then
#  startrun=$2;
#  endrun=$3;
#fi

#done


#for ((n=$startrun;n<$endrun;n++)); do
for n in $subrunlist; do
for field in $fieldlist; do

    subrun=$n

    #planetfile= $planetdir$planetroot$field.$subrun
    #echo "$planetfile"
    if [ ! -e $planetdir$planetroot$field.$subrun ] ; then
	echo 'No planet file '$planetdir$planetroot$field.$subrun	
	continue;
    fi

    #Check if we have enough free compute nodes
    ncont=0
    if [ $ncont -eq 1 ]; then
	continue;
    fi

    waitfornodes


    
    launchfile=$finaldir$runname/logs/run_${runname}_${subrun}_$field.sh
    logfilebase=$finaldir$runname/logs/run_${runname}_${subrun}_$field.sh

    #Check to see if this run has been attempted yet
    if [ ! -e $launchfile ]; then
    
        #Write the launch file script for this sub-run
	if [ 1 -eq 1 ]; then
	    echo "#!/bin/bash"
	    echo "#PBS -e $logfilebase.e"
	    echo "#PBS -o $logfilebase.o"
	    echo "#PBS -N ${runname}_${subrun}_${field}"
	    echo "#PBS -l walltime=6:00:00"
	    echo "#PBS -l nodes=1:ppn=1,mem=1GB"
	    echo "#PBS -m a"
	    echo "#PBS -M johnson.7080@osu.edu"
	    
	    #echo "#PBS -q $QUEUE"

   	    #Setup directories on the compute node
	    echo "if [ ! -d $outputdir ]; then"
	    echo "mkdir -p $outputdir"
	    echo "echo Created $outputdir"
	    echo "fi"

	    echo "if [ ! -d $outputdir$runname ]; then"
	    echo "mkdir -p $outputdir$runname"
	    echo "echo Created $outputdir$runname"
	    echo "fi"
	    
	    echo "cd $outputdir"

  	    #echo "${SCRIPTDIR}copyback.sh $paramfile $subrun &"

	    #echo "export LD_LIBRARY_PATH=/share/apps/lib:$LD_LIBRARY_PATH"
	    #echo "$finaldir$runname/logs/$executable -i $paramfile -s $subrun -f $field"
	    #added full path here. maybe it's best to require full path in calling? I need to think about this a bit...
	    echo "$finaldir$runname/logs/$executable -i $base_dir$param_file_dir$paramfile -s $subrun -f $field"
	    echo "mv $outputdir$runname/${runname}_${subrun}_$field.??? $finaldir$runname/raw/"
	    echo "for iii in $outputdir$runname/${runname}_${subrun}_${field}_*.???.fm*; do mv \$iii $finaldir$runname/fm/; done"
	    if [ `echo $OUTPUTLC | awk '{if($1>0) print 1; else print 0}'` -ne 0 ]; then 
		echo "for iii in $outputdir$runname/${runname}_${subrun}_${field}_*.???.lc; do mv \$iii $finaldir$runname/lc/; done"
	    fi

	    if [ $OUTPUTIMAGES -ne 0 ]; then 
		echo "for iii in $outputdir$runname/${runname}_${subrun}_${field}_*.fits; do mv \$iii $finaldir$runname/lc/; done"
	    fi

	fi > $launchfile
	
        #Launch the sub-run
	qsub  $launchfile
	#qsub $launchfile

	sleep 0.1
    fi
done
done
#done

