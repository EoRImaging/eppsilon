#!/bin/bash
# Top level script to integrate healpix cubes and run power spectrum code.

#Parse flags for inputs
while getopts ":f:s:e:" option
do
   case $option in
        f) file_path_cubes="$OPTARG";;
        s) starting_obs=$OPTARG;;
        e) ending_obs=$OPTARG;;
        \?) echo "Unknown option: Accepted flags are -f (file path to cubes), -s (starting_obs) and -e (ending obs)"
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

#generate the file names in the folder for even and odd

name_array_even_unsorted=($(find $file_path_cubes -name "[[:digit:]]*even_cube.sav")) # ensure the first character is a numeral (to avoid Combined... getting lumped in.)
name_array_odd_unsorted=($(find $file_path_cubes -name "[[:digit:]]*odd_cube.sav"))

#arrange them in order of increasing obs id...which isn't done automatically

name_array_even=($(for a in "${name_array_even_unsorted[@]}";do echo "$a"; done | sort -n))

name_array_odd=($(for a in "${name_array_odd_unsorted[@]}";do echo "$a"; done | sort -n))

i=0
for full_file_name in "${name_array_even[@]}"
do
   file_name=${full_file_name##*/}
   obsid_in_name_even[$i]=${file_name%%_*}
   i=$((i+1))
done

i=0
for full_file_name in "${name_array_odd[@]}"
do
   file_name=${full_file_name##*/}
   obsid_in_name_odd[$i]=${file_name%%_*}
   i=$((i+1))
done

#Start and end integration at specified observation ids if specified
#Tested on the even cubes arbitrarily

if [ ! -z $starting_obs ]
then
   i=0
   for obs_id in "${name_array_even[@]}"
   do
	if [[ "$obs_id" =~ "$starting_obs" ]]
	then
	   starting_index=$i
	fi
	i=$((i+1))	
   done
   
   if [ -z $starting_index ]
   then
	echo "Specified starting observation id not found"
	exit 1
   fi
else
  starting_index=0
fi

if [ ! -z $ending_obs ]
then
   i=0
   for obs_id in "${name_array_even[@]}"
   do
	if [[ "$obs_id" =~ "$ending_obs" ]]
	then
	   ending_index=$i
	fi
	i=$((i+1))	
   done
   
   if [ -z $ending_index ]
   then
	echo "Specified ending observation id not found"
	exit 1
   fi
else
  len=${#name_array_even[*]}
  ending_index=$(expr $len - 1)
fi

min_obs=${obsid_in_name_even[$starting_index]}
max_obs=${obsid_in_name_even[$ending_index]}

PSpath=$(idl -e 'print,rootdir("ps")') ### NOTE this only works if idlstartup doesn't have any print statements (e.g. healpix check)
mem=4G ## per core. this should be an option to the script
nslots=10 # cores to use

outfile=${file_path_cubes}/${obsid_in_name_even[$starting_index]}_${obsid_in_name_even[$ending_index]}_out.log
errfile=${file_path_cubes}/${obsid_in_name_even[$starting_index]}_${obsid_in_name_even[$ending_index]}_err.log

qsub -l h_vmem=$mem,h_stack=512k -V -v file_path_cubes=$file_path_cubes,min_obs=$min_obs,max_obs=$max_obs,starting_index=$starting_index,ending_index=$ending_index,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/PS_job.sh

## The rest is done in PS_job.sh now.

#Use the command to integrate healpix cubes through idl

#input_file="$file_path_cubes""*even_cube.sav"

#idl -e integrate_healpix_cubes -args "$input_file" $starting_index $ending_index

#input_file="$file_path_cubes""*odd_cube.sav"

#idl -e integrate_healpix_cubes -args "$input_file" $starting_index $ending_index

#Make power spectra through a ps wrapper in idl

#input_file="$file_path_cubes""Healpix/"

#idl -e mit_wrapper -args $input_file '${obsid_in_name_even[$starting_index]}-${obsid_in_name_even[$ending_index]}'


