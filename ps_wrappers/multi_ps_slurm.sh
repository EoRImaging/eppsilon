#!/bin/bash

# Batch run eppsilon on a set of ObsIDs with minimal file production.

module load git/2.10.2
module load idl
shopt -s expand_aliases; source $IDL/envi53/bin/envi_setup.bash
PSpath=$(idl -e 'print,rootdir("eppsilon")')

while getopts ":d:f:s:e:w:n:q:m:" option
do
   case $option in
        d) FHDdir="$OPTARG";;			#file path to fhd directory with cubes
        f) obs_file_name="$OPTARG";;		#txt file of obs ids or a single obsid
        s) starting_obs=$OPTARG;;       #starting observation in text file for choosing a range
        e) ending_obs=$OPTARG;;         #ending observation in text file for choosing a range
        w) wallclock_time=$OPTARG;;     	#Time for execution in slurm
        n) ncores=$OPTARG;;             	#Number of cores for slurm
        q) qos=$OPTARG;;        # Quality of service for Oscar
        m) mem=$OPTARG;;                	#Total memory for slurm
        \?) echo "Unknown option."
            exit 1;;
        :) echo "Missing option argument for input flag"
           exit 1;;
   esac
done


# SLURM accounting
if [ -z ${qos} ]
then
    acct='jpober-condo'
    qos='jpober-condo'
fi

if [ ${qos} == 'pri-alanman' ]
then
    acct="default"
fi


if [ ! -e "$obs_file_name" ]
then
    echo "Obs list is either not a file or the file does not exist!"
    echo "Assuming the obs list is a single observation id."

    obs_id_array=("$obs_file_name")
#    version=$obs_file_name  #Currently assuming that the integrate list is a single obsid
else
    first_line=$(head -n 1 $obs_file_name)
    i=0
    while read line
    do
       if [ ! -z "$line" ]; then
          obs_id_array[$i]=$line
          i=$((i + 1))
       fi
    done < "$obs_file_name"

    max=${obs_id_array[$((i-1))]}
    min=${obs_id_array[0]}
#    version=$(basename $obs_file_name) # get filename
#    version="${version%.*}" # strip extension
fi

#Set typical wallclock_time for standard PS with obs ids if not set.
if [ -z ${wallclock_time} ]; then wallclock_time=10:00:00; fi

#Set typical slots needed for standard PS with obs ids if not set.
if [ -z ${ncores} ]; then ncores=10; fi

#Set typical memory needed for standard PS with obs ids if not set.
if [ -z ${mem} ]; then mem=20G; fi



#If minimum not specified, start at minimum of obs_file
if [ -z ${starting_obs} ]
then
   echo "Starting observation not specified: Starting at minimum of "$(basename $obs_file_name)""
   starting_obs=$min
fi

#If maximum not specified, end at maximum of obs_file
if [ -z ${ending_obs} ]
then
   echo "Ending observation not specified: Ending at maximum of "$(basename $obs_file_name)""
   ending_obs=$max
fi

#Create a list of observations using the specified range, or the full observation id file. 

unset good_obs_list
startflag=0
endflag=0

for obs_id in "${obs_id_array[@]}"; do 
     if [ $obs_id == $starting_obs ]; then 
     	startflag=1
     fi
     if [ $startflag -eq 1 ] && [ $endflag -eq 0 ]; then
     	good_obs_list+=($obs_id)
     fi
     if [ $obs_id == $ending_obs ]; then
    	endflag=1
     fi
done	

#echo ${good_obs_list[@]}

nobs=${#obs_id_array[@]}
outfile=${FHDdir}/ps/ps-%A_%a.out

sbatch -A $acct --qos=$qos --mem=$mem -t ${wallclock_time} -n ${ncores} --array=0-$(( $nobs - 1 ))%20 --export=file_path_cubes=$FHDdir,nobs=$nobs,ncores=$ncores -o ${outfile}  ${PSpath}ps_wrappers/PS_list_slurm_job.sh ${good_obs_list[@]}
