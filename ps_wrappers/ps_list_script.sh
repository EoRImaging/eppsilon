#!/bin/bash
# Top level script to integrate healpix cubes and run power spectrum code.
# This version will take a file with a list of obsids to include in the integration

#Parse flags for inputs
while getopts ":f:o:e:" option
do
   case $option in
        f) file_path_cubes="$OPTARG";;
        o) obs_list_path="$OPTARG";;
        \?) echo "Unknown option: Accepted flags are -f (file path to cubes), -o (obs list path)"
            exit 1;;
        :) echo "Missing option argument for input flag"
           exit 1;;
   esac
done

#Manual shift to the next flag
shift $(($OPTIND - 1))

#Throw error if no file path to filenames
if [ -z ${file_path_cubes} ]
then
   echo "Need to specify a file path to individual cubes: Example /nfs/complicated_path/fhd_mine/"
   exit 1
fi

#Throw error if file path does not exist
if [ ! -d "$file_path_cubes" ]
then
   echo "Argument after flag -f is not a real directory. Argument should be the file path to the location of cubes to integrate."
   exit 1
fi

#Error if obs_list_path is not set
if [ -z ${obs_list_path} ]; then
    echo Need to specify obs list file path with option -o
    exit 1
fi
#Error if obs list filename does not exist
if [ ! -e "$obs_list_path" ]; then
    echo Obs list file does not exist.
    exit 1
fi

PSpath=$(idl -e 'print,rootdir("ps")') ### NOTE this only works if idlstartup doesn't have any print statements (e.g. healpix check)
mem=4G ## per core. this should be an option to the script
nslots=10 # cores to use

version=$(basename $obs_list_path) # get filename
version="${version%.*}" # strip extension
echo Version is $version

# read in obses and divide into chunks
obs=0
rm ${file_path_cubes}/Healpix/${version}_int_chunk*.txt # remove any old chunk files lying around
while read line
do
    ((chunk=obs/100+1))
    echo $line >> ${file_path_cubes}/Healpix/${version}_int_chunk${chunk}.txt
    ((obs++))
done < $obs_list_path
nchunk=$chunk # number of chunks we ended up with

unset idlist
if [ "$nchunk" -gt "1" ]; then
    # set up files for master integration
    sub_cubes_list=${file_path_cubes}/Healpix/${version}_sub_cubes.txt
    rm $sub_cubes_list # remove any old lists
    # launch separate chunks
    for chunk in $(seq 1 $nchunk); do
	chunk_obs_list=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}.txt
	outfile=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_out.log
	errfile=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_err.log
	message=$(qsub -l h_vmem=$mem,h_stack=512k -V -v file_path_cubes=$file_path_cubes,obs_list_path=$chunk_obs_list,version=$version,chunk=$chunk,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/integrate_job.sh)
	message=($message)
	if [ "$chunk" -eq 1 ]; then idlist=${message[2]}; else idlist=${idlist},${message[2]}; fi
	echo Healpix/Combined_obs_${version}_int_chunk${chunk} >> $sub_cubes_list # trick it into finding our sub cubes
    done
    # master integrator
    chunk=0
    outfile=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_out.log
    errfile=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_err.log
    message=$(qsub -hold_jid $idlist -l h_vmem=$mem,h_stack=512k -V -v file_path_cubes=$file_path_cubes,obs_list_path=$sub_cubes_list,version=$version,chunk=$chunk,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/integrate_job.sh)
    message=($message)
    master_id=${message[2]}
else
    # Just one integrator
    chunk=0
    outfile=${file_path_cubes}/${version}_int_chunk${chunk}_out.log
    errfile=${file_path_cubes}/${version}_int_chunk${chunk}_err.log
    message=$(qsub -l h_vmem=$mem,h_stack=512k -V -v file_path_cubes=$file_path_cubes,obs_list_path=$obs_list_path,version=$version,chunk=$chunk,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/integrate_job.sh)
    message=($message)
    master_id=${message[2]}
fi


outfile=${file_path_cubes}/${version}_ps_out.log
errfile=${file_path_cubes}/${version}_ps_err.log

qsub -hold_jid $master_id -l h_vmem=$mem,h_stack=512k -V -v file_path_cubes=$file_path_cubes,obs_list_path=$obs_list_path,version=$version,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh
