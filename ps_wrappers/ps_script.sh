#!/bin/bash

######################################################################################
# Top level script to integrate healpix cubes and run power spectrum code.
#
# A file path to the fhd directory is needed.
# 
# A file path to a text file listing observation ids OR preintegrated subcubes is
# needed.
# 
# If a text file of observation ids to be used in integration is specified, the obs 
# ids are assumed to be seperated by newlines.
#
# If a text file of preintegrated subcubes is specified, the format should be
# the name of the save file seperated by newlines.  "even_cube.sav" and "odd_cube.sav"
# is not necessary to include, as both will be used anyways.  The subcubes are
# assumed to be in <fhd_directory>/Healpix/. If elsewhere in the FHD directory, the 
# name of the subcubes must specify this in the text file as Other_than_Healpix/<name>.
#
# Set -ps to 1 to skip integration and make cubes only.
# 
# NOTE: print statements must be turned off in idl_startup file (e.g. healpix check)
######################################################################################

#Parse flags for inputs
while getopts ":d:f:p:w:n:m:o:l:h:t:" option
do
   case $option in
        d) FHDdir="$OPTARG";;			#file path to fhd directory with cubes
        f) integrate_list="$OPTARG";;		#txt file of obs ids or subcubes or a single obsid
        p) priority=$OPTARG;;           	#priority level for grid engine qsub
        w) wallclock_time=$OPTARG;;     	#Time for execution in grid engine
        n) nslots=$OPTARG;;             	#Number of slots for grid engine
        m) mem=$OPTARG;;                	#Memory per core for grid engine
	o) ps_only=$OPTARG;;			#Flag for skipping integration to make PS only
        l) legacy=$OPTARG;;                     #Use legacy directory structure. hacky solution to a silly problem.
        h) hold=$OPTARG;;                       #Hold for a job to finish before running. Useful when running immediately after firstpass
	t) tukey_filter=$OPTARG;;               #Apply a Tukey image window filter during eppsilon
        \?) echo "Unknown option: Accepted flags are -d (file path to fhd directory with cubes), -f (obs list or subcube path or single obsid), "
	    echo "-p (priority for grid engine), -w (wallclock time), -n (number of slots), -m (memory allocation), "
	    echo "-o (make ps only without integration), and -l (legacy flag for old file structure),"
	    echo "-h (hold int/ps script on a running job id), and -t (apply a tukey window filter during ps),"
            exit 1;;
        :) echo "Missing option argument for input flag"
           exit 1;;
   esac
done

#Manual shift to the next flag
shift $(($OPTIND - 1))

#Throw error if no file path to FHD directory
if [ -z ${FHDdir} ]
then
   echo "Need to specify a file path to a FHD directory with cubes: Example /nfs/complicated_path/fhd_mine/"
   exit 1
fi

#Throw error if file path does not exist
if [ ! -d "$FHDdir" ]
then
   echo "Argument after flag -d is not a real directory. Argument should be the file path to the location of cubes to integrate."
   exit 1
fi

#Remove extraneous / on FHD directory if present
if [[ $FHDdir == */ ]]; then FHDdir=${FHDdir%?}; fi

#Error if integrate_list is not set
if [ -z ${integrate_list} ]
then
    echo "Need to specify obs list file path or preintegrated subcubes list file path with option -f"
    exit 1
fi

#Warning if integrate list filename does not exist
if [ ! -e "$integrate_list" ]
then
    echo "Integrate list is either not a file or the file does not exist!"
    echo "Assuming the integrate list is a single observation id."

    if [ -z ${ps_only} ]
    then
        echo "ps_only flag must be set if integrate list is a single observation id. Set -o 1 if desired function"
        exit 1
    fi 
    version=$integrate_list  #Currently assuming that the integrate list is a single obsid
else
    version=$(basename $integrate_list) # get filename
    version="${version%.*}" # strip extension
fi

#Default priority if not set.
if [ -z ${priority} ]; then priority=0; fi

#Set typical wallclock_time for standard PS with obs ids if not set.
if [ -z ${wallclock_time} ]; then wallclock_time=40:00:00; fi

#Set typical slots needed for standard PS with obs ids if not set.
if [ -z ${nslots} ]; then nslots=10; fi

#Set typical memory needed for standard PS with obs ids if not set.
if [ -z ${mem} ]; then mem=4.5G; fi

#Set default to do integration
if [ -z ${ps_only} ]; then ps_only=0; fi

#Set default to use current file structure
if [ -z ${legacy} ]; then legacy=0; fi

# create hold string
if [ -z ${hold} ]; then hold_str=""; else hold_str="-hold_jid ${hold}"; fi

# create hold string
if [[ -n ${tukey_filter} ]]; then tukey_filter="Tukey"; fi

### NOTE this only works if idlstartup doesn't have any print statements (e.g. healpix check)
PSpath=$(idl -e 'print,rootdir("ps")')

#Versions made during integrate list logic check above
echo Version is $version

if [ ! -e "$integrate_list" ]
then
    first_line=$integrate_list
else
    first_line=$(head -n 1 $integrate_list)
fi

first_line_len=$(echo ${#first_line})

rm -f ${FHDdir}/Healpix/${version}_int_chunk*.txt # remove any old chunk files lying around

exit_flag=0

#Check that cubes or integrated cubes are present, print and error if they are not
if [ "$ps_only" -ne "1" ]; then 	#only if we're integrating
while read line
do
   if [ "$first_line_len" == 10 ]; then
      if [ "$legacy" -ne "1" ]; then
	  if ! ls $FHDdir/Healpix/$line*cube*.sav &> /dev/null; then
              echo Missing cube for obs $line
	      if [ -z "$hold" ]; then
		  exit_flag=1
	      fi
	  fi
      else
	  if ! ls $FHDdir/$line*cube*.sav &> /dev/null; then
	      echo Missing cube for obs $line
	      exit_flag=1
	  fi
      fi
   else
      if [[ "$first_line" != */* ]]; then
	 check=$FHDdir/Healpix/$line*.sav
      else
	 check=$FHDdir/$line*.sav
      fi
      if ! ls $check &> /dev/null; then
	    echo Missing save file for $line
	    exit_flag=1
      fi
   fi
done < $integrate_list
fi

if [ "$exit_flag" -eq 1 ]; then exit 1; fi

if [ "$first_line_len" == 10 ]; then
    
    # Just PS if flag has been set
    #if [ "$ps_only" -eq "1" ]; then
    #    outfile=${FHDdir}/ps/${version}_ps_out.log
    #    errfile=${FHDdir}/ps/${version}_ps_err.log
    #    if [ ! -d ${FHDdir}/ps ]; then mkdir ${FHDdir}/ps; fi
#	echo "Running only ps code"
        #qsub ${hold_str} -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh
        #exit $?
	
    #fi

    # read in obs ids 100 at a time and divide into chunks to integrate in parallel mode
    obs=0   

    while read line
    do
        ((chunk=obs/100+1))		#integer division results in chunks labeled 0 (first 100), 1 (second 100), etc
        echo $line >> ${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt	#put that obs id into the right txt file
        ((obs++))			#increment obs for the next run through
    done < $integrate_list
    nchunk=$chunk 			#number of chunks we ended up with

else

    if [[ "$first_line" != */* ]]; then
   
        chunk=0 
        while read line
        do
            echo Healpix/$line >> ${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt        #put that obs id into the right txt file
        done < $integrate_list
        nchunk=$chunk                       #number of chunks we ended up with
    
    else

        chunk=0 
        while read line
        do
            echo $line >> ${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt        #put that obs id into the right txt file
        done < $integrate_list
        nchunk=$chunk                       #number of chunks we ended up with

    fi

fi

unset idlist
if [ "$ps_only" -ne "1" ]; then   
    if [ "$nchunk" -gt "1" ]; then

        # set up files for master integration
        sub_cubes_list=${FHDdir}/Healpix/${version}_sub_cubes.txt
        rm $sub_cubes_list # remove any old lists

        # launch separate chunks
        for chunk in $(seq 1 $nchunk); do
	    chunk_obs_list=${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt
	    outfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_out.log
	    errfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_err.log
	    for evenodd in even odd; do
		for pol in XX YY; do 
	    	    message=$(qsub ${hold_str} -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$chunk_obs_list,version=$version,chunk=$chunk,nslots=$nslots,legacy=$legacy,evenodd=$evenodd,pol=$pol -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/integrate_job.sh)
	    	    message=($message)
	    	    if [ "$chunk" -eq 1 ] && [[ "$evenodd" = "even" ]] && [[ "$pol" = "XX" ]]; then idlist=${message[2]}; else idlist=${idlist},${message[2]}; fi
		done
	    done
	    echo Combined_obs_${version}_int_chunk${chunk} >> $sub_cubes_list # trick it into finding our sub cubes
        done

        # master integrator
        chunk=0
        outfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_out.log
        errfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_err.log
	for evenodd in even odd; do
	    for pol in XX YY; do 
        	message=$(qsub -hold_jid $idlist -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$sub_cubes_list,version=$version,chunk=$chunk,nslots=$nslots,legacy=$legacy,evenodd=$evenodd,pol=$pol -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/integrate_job.sh)
        	message=($message)
		if [[ "$evenodd" = "even" ]] && [[ "$pol" = "XX" ]]; then idlist_master=${message[2]}; else idlist_master=${idlist_master},${message[2]}; fi
	    done
	done
	hold_str="-hold_jid ${idlist_master}"

    else

        # Just one integrator
        mv ${FHDdir}/Healpix/${version}_int_chunk1.txt ${FHDdir}/Healpix/${version}_int_chunk0.txt
        chunk=0
        chunk_obs_list=${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt
        outfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_out.log
        errfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_err.log
	for evenodd in even odd; do
	    for pol in XX YY; do
        	message=$(qsub ${hold_str} -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$chunk_obs_list,version=$version,chunk=$chunk,nslots=$nslots,legacy=$legacy,evenodd=$evenodd,pol=$pol -e $errfile -o $outfile -pe chost $nslots ${PSpath}ps_wrappers/integrate_job.sh)
       		message=($message)
		if [[ "$evenodd" = "even" ]] && [[ "$pol" = "XX" ]]; then idlist_int=${message[2]}; else idlist_int=${idlist_int},${message[2]}; fi
	    done
	done
        hold_str="-hold_jid ${idlist_int}"

    fi
else
    echo "Running only ps code" # Just PS if flag has been set
fi

outfile=${FHDdir}/ps/logs/${version}_ps_out
errfile=${FHDdir}/ps/logs/${version}_ps_err


if [ ! -d ${FHDdir}/ps ]; then
    mkdir ${FHDdir}/ps
fi
if [ ! -d ${FHDdir}/ps/logs ]; then
    mkdir ${FHDdir}/ps/logs
fi

###XX, even
pol='xx'
evenodd='even'

#weights/variance
cube_type='weights'
message=$(qsub ${hold_str} -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_even_weights.log -o ${outfile}_xx_even_weights.log -N PS_xe_weights -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id=${message[2]}

#dirty
cube_type='dirty'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_even_dirty.log -o ${outfile}_xx_even_dirty.log -N PS_xe_dirty -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${message[2]}

#model
cube_type='model'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_even_model.log -o ${outfile}_xx_even_model.log -N PS_xe_model -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#res
cube_type='res'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_even_res.log -o ${outfile}_xx_even_res.log -N PS_xe_res -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

###XX, odd
pol='xx'
evenodd='odd'

#weights/variance
cube_type='weights'
message=$(qsub ${hold_str} -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_odd_weights.log -o ${outfile}_xx_odd_weights.log -N PS_xo_weights -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id=${message[2]}

#dirty
cube_type='dirty'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_odd_dirty.log -o ${outfile}_xx_odd_dirty.log -N PS_xo_dirty -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#model
cube_type='model'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_odd_model.log -o ${outfile}_xx_odd_model.log -N PS_xo_model -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#res
cube_type='res'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_xx_odd_res.log -o ${outfile}_xx_odd_res.log -N PS_xo_res -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

###YY, even
pol='yy'
evenodd='even'

#weights/variance
cube_type='weights'
message=$(qsub ${hold_str} -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_even_weights.log -o ${outfile}_yy_even_weights.log -N PS_ye_weights -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id=${message[2]}

#dirty
cube_type='dirty'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_even_dirty.log -o ${outfile}_yy_even_dirty.log -N PS_ye_dirty -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#model
cube_type='model'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_even_model.log -o ${outfile}_yy_even_model.log -N PS_ye_model -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#res
cube_type='res'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_even_res.log -o ${outfile}_yy_even_res.log -N PS_ye_res -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

###YY, odd
pol='yy'
evenodd='odd'

#weights/variance
cube_type='weights'
message=$(qsub ${hold_str} -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_odd_weights.log -o ${outfile}_yy_odd_weights.log -N PS_yo_weights -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id=${message[2]}

#dirty
cube_type='dirty'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_odd_dirty.log -o ${outfile}_yy_odd_dirty.log -N PS_yo_dirty -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#model
cube_type='model'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_odd_model.log -o ${outfile}_yy_odd_model.log -N PS_yo_model -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#res
cube_type='res'
message=$(qsub -hold_jid $id -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$wallclock_time -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_yy_odd_res.log -o ${outfile}_yy_odd_res.log -N PS_yo_res -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh)
message=($message)
id_list=${id_list},${message[2]}

#final plots
if [[ -n ${tukey_filter} ]]; then plot_walltime=10:00:00; else plot_walltime=00:20:00; fi
qsub -hold_jid $id_list -p $priority -l h_vmem=$mem,h_stack=512k,h_rt=$plot_walltime -V -v file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,nslots=$nslots,image_filter_name=$tukey_filter -e ${errfile}_plots.log -o ${outfile}_plots.log -N PS_plots -pe chost $nslots ${PSpath}ps_wrappers/PS_list_job.sh
