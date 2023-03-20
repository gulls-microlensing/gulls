#!/bin/bash

# I think this code is all updated to be good for the new paramfiles path stucture

#This piece of code is common to all exigere scripts, and is called at the start of each

# This script is just preparing other analysis scripts? None of this is used in and of the src files, but still needs to be updated to use the paths.txt file
#We need the file paths.txt in this directory to establish some structure. Currently, everything is relative to base_dir.

# should have access to GULLS_BASE_DIR env variable

#get the env variable for the base_dir
base_dir=$GULLS_BASE_DIR
# the name of the paramfile is an arguement
paramfile=$1
# right now defaulting to have a parameterFiles dir in base_dir, may not be the best but working withit for now
param_file_dir="parameterFiles/"

#The parameter file is the only option needed - it contains all the other run-by-run options
if [ $# -lt 1 ]; then
  echo "Usage:"
  echo "   $0 <parameterFile>"
  exit 1
fi

#paramfile=$1
#paramdir=""


#test exit

if [ ! -e $base_dir$param_file_dir$paramfile ]; then
  echo "Parameter file '$base_dir$param_file_dir$paramfile' not found"
  exit
fi

########################################################################
#    Less commonly changed options
########################################################################

#General options:
NNODES=200     #The maximum number of compute nodes to use
QUEUE=mfun

#exigereSetup.sh options:

#exigereLaunch.sh options:

#exigereTidy.sh options:


#GENERAL STRUCTURE: grep VAR_NAME in the base_dir/param_file_dir/paramfile
#TODO do i need to set up the output directory too? this should require attention of the user

base_dir=$GULLS_BASE_DIR

# all of these varaibles are path dependent (or better grouped here), so we give them paths
lensdir=$base_dir`grep LENS_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
lenslist=`grep LENS_LIST $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
#echo  $base_dir
#echo $param_file_dir
#echo $paramfile
#echo $lensdir
srcdir=$base_dir`grep SOURCE_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
srclist=`grep SOURCE_LIST $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
sfdir=$base_dir`grep STARFIELD_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
planetdir=$base_dir`grep PLANET_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
planetroot=`grep PLANET_ROOT $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
obsdir=$base_dir`grep OBSERVATORY_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
obslist=$obsdir`grep OBSERVATORY_LIST $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
weatherdir=$base_dir`grep WEATHER_PROFILE_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`

#these varaibles have no path dependence, so they just get read in from the paramfile
#the output and final directory are independent of the location of this file tree
runname=`grep RUN_NAME $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
outputdir=`grep OUTPUT_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
#finaldir expects full path in parameterfile

finaldir=`grep FINAL_DIR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
band=`grep PRIMARYCOLOUR $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
subrunsize=`grep SUBRUNSIZE $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
nsubruns=`grep NSUBRUNS $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
simlength=`grep SIMULATION_LENGTH $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'` #Length of time over which events could occur in the sim in years
    #(Not the length of the data collection)

#XXX see how this is used, may need a path before depsnding on how it is called later
executable=`grep EXECUTABLE $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
WAITTIME=`grep MAXTIME $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
OUTPUTLC=`grep OUTPUT_LC $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`
OUTPUTIMAGES=`grep OUTPUT_IMAGES $base_dir$param_file_dir$paramfile | awk -v FS='=' '{print $2}'`


logscaleColumns="10 26 28 33 34 35 36 37 42 43 46 47 48 50 79"
linearColumns="2 3 4 9 13 14 15 18 22 23 24 25 30 31 38 39 44 45 54 56 57 58 60 64 66 69 71"
nHistBins=20

########################################################################
#Function definitions
########################################################################

messagecount=0
function waitfornodes ()
{
    if [ -e ~/.NNODES ]; then
	NNODES=`head -n1 ~/.NNODES | awk '{print $1}'`
	QUEUE=`head -n1 ~/.NNODES | awk '{print $2}'`
    fi

    #This function waits for a node to become free
    while [ `qstat | grep $USER | grep $QUEUE | wc -l` -ge $NNODES ]; do
	if [ $messagecount -eq 0 ]; then
	    echo -n Waiting for nodes to become free
	else
	    echo -n .
	fi
	sleep 60
	messagecount=1
	if [ -e ~/.NNODES ]; then
		NNODES=`head -n1 ~/.NNODES | awk '{print $1}'`
		QUEUE=`head -n1 ~/.NNODES | awk '{print $2}'`
	fi
    done

    if [ $messagecount -ne 0 ]; then echo; fi
    
    messagecount=0
}

########################################################################
#    Conventions
########################################################################

#Listed are the conventions the gulls pipeline relies on in order to 
#run successfully. Changes to these could require changes to the
#pipeline scripts, the underlying simulation code and the analysis
#scripts

###  The results and logs of each simulation are stored in a directory
#    named <outputdir>/<runname>/raw

###  The scripts used to run each sub-run, and their stdout and stderr
#    output are put in the directory <outputdir>/<runname>/logs

###  Files associated with the analysis of runs will be placed in
#    the <outputdir>/<runname>/analysis directory

###  The output files of each sub-run are named 
#    <runname>_<subRunNo>.out and <runname>_<subRunNo>.log
#    Both files are necessary for subsequent analysis. The files in 
#    the logs/ directory are not needed for analysis.

###  The lightcurves, and any subsequent analysis files, are stored in 
#    the outputdir/runname/lc directory, and are named
#    lc_<runname>_<subRunNo>_<idNo>.txt

###  The gulls pipeline scripts only take one argument, any other
#    options should be passed via the parameter file. Each run should
#    have a separate parameter file.

###  A copy of the parameter file, observatory list and files and 
#    weather files will be copied to the output files for future 
#    reference.

###  The SYNTHGALROOT option in the parameter file is expected to
#    end with a dot. i.e. gal_H_0.dat.

###  The gulls simulator is launched as
#    ./gulls -i <parameterFile> -s <subRunNo>

###  If you have a short user name, you may have to adapt the scripts 
#    to reliably select only jobs you have submitted from the qstat 
#    running jobs list

###  Directories are listed in the parameter file with their trailing
#    slash. In fact any variable that holds a directory, and whos 
#    name ends in dir or DIR, will include the trailing slash

### Runname must be 10 characters or less, and must not contain an
#   underscore
