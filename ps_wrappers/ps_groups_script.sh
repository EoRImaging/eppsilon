#!/bin/bash
# Script to integrate/PS in a specific way.
# Assumption: subcubes have already been integrated. 
# This script will take those subcubes and integrate them in groups and PS them.
# names of cubes to be included should be in $filename with the following format:
#   Healpix/<name>  (do not include "even_cube" or "odd_cube", but do include Healpix/)

#Parse inputs
while getopts ":d:f:" option
do
    case $option in
	d) FHDdir="$OPTARG";;
	f) filename="$OPTARG";;
	\?) echo "Unknown option: Accepted flags are -d (directory), -f (filename with cube names)"
            exit 1;;
        :) echo "Missing option argument for input flag"
           exit 1;;
   esac
done

#Throw error if dir doesn't exist
if [ ! -d $FHDdir ]; then
    echo FHD directory $dir does not exist.
    exit 1
fi

#Error if filename is not set
if [ -z $filename ]; then
    echo Need to specify filename with option -f
    exit 1
fi
#Error if filename does not exist
if [ ! -e $filename ]; then
    echo Filename $filename does not exist
    exit 1
fi

PSpath=$(idl -e 'print,rootdir("ps")') ### NOTE this only works if idlstartup doesn't have any print statements (e.g. healpix check)
mem=4G ## per core. this should be an option to the script
nslots=10 # cores to use

version=$(basename $filename) # get filename
version="${version%.*}" # strip extension
echo Version is $version

chunk=0
outfile=${FHDdir}/Healpix/${version}_out.log
errfile=${FHDdir}/Healpix/${version}_err.log
message=$(qsub -l h_vmem=$mem,h_stack=512k -V -v file_path_cubes=$FHDdir,obs_list_path=$filename,version=$version,chunk=$chunk,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/integrate_job.sh)
message=($message)
int_id=${message[2]}

outfile=${FHDdir}/ps/${version}_ps_out.log
errfile=${FHDdir}/ps/${version}_ps_err.log

if [ ! -d ${FHDdir}/ps ]; then
    mkdir ${FHDdir}/ps
fi

qsub -hold_jid $int_id -l h_vmem=$mem,h_stack=512k -V -v file_path_cubes=$FHDdir,obs_list_path=$filename,version=$version,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh
# note that obs_list_path is only used to calculate nobs. So this will not get the right answer, but it will be close-ish (off by ~20x)