#! /bin/bash
#$ -V
#$ -N Cube_integrator
#$ -S /bin/bash

#This script is an extra layer between Grid Engine and IDL commands because
#Grid Engine runs best on bash scripts.

#inputs needed: file_path_cubes, obs_list_path, version, chunk, nslots, evenodd, pol
#chunk is the chunk number when the list was broken up. 0 for "master" or only chunk

echo JOBID ${JOB_ID}

if [ "$legacy" -ne "1" ]; then
    #Create a name for the obs txt file based off of inputs
    evenoddpol_file_paths=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_${evenodd}${pol}_list.txt
    #clear old file paths
    rm $evenoddpol_file_paths

    #***Fill the obs text file with the obsids to integrate
    nobs=0
    while read line
    do
	evenoddpol_file=${file_path_cubes}/Healpix/${line}_${evenodd}_cube${pol}.sav
	echo $evenoddpol_file >> $evenoddpol_file_paths
	((nobs++))
    done < "$obs_list_path"
    #***

    unset int_pids

    #***If the integration has been split up into chunks, name the save file specifically off of that.
    if [ "$chunk" -gt "0" ]; then
	save_file_evenoddpol="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_${evenodd}_cube${pol}.sav
    else
	save_file_evenoddpol="$file_path_cubes"/Healpix/Combined_obs_${version}_${evenodd}_cube${pol}.sav
    fi
    #***

    #***Run the integration IDL script
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e integrate_healpix_cubes -args "$evenoddpol_file_paths" "$save_file_evenoddpol" &
    int_pids+=( $! )
    #***

    wait ${int_pids[@]} # Wait for integration to finish before making PS

else # legacy

    # Get list of obs to integrate, turn them into file names to integrate
    even_file_paths=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_even_list.txt
    odd_file_paths=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_odd_list.txt
    # clear file paths
    rm $even_file_paths
    rm $odd_file_paths
    nobs=0
    while read line
    do
	even_file=${file_path_cubes}/${line}_even_cube.sav
	odd_file=${file_path_cubes}/${line}_odd_cube.sav
	echo $even_file >> $even_file_paths
	echo $odd_file >> $odd_file_paths
	((nobs++))
    done < $obs_list_path

#Integrate cubes
    unset int_pids

    if [ "$chunk" -gt "0" ]; then
	save_file_even="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_even_cube.sav
	save_file_odd="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_odd_cube.sav
    else
	save_file_even="$file_path_cubes"/Healpix/Combined_obs_${version}_even_cube.sav
	save_file_odd="$file_path_cubes"/Healpix/Combined_obs_${version}_odd_cube.sav
    fi
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e integrate_healpix_cubes -args "$even_file_paths" "$save_file_even" &
    int_pids+=( $! )
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e integrate_healpix_cubes -args "$odd_file_paths" "$save_file_odd" &
    int_pids+=( $! )
    wait ${int_pids[@]} # Wait for integration to finish before making PS
fi

# TODO: check that they completed - exit with code 100 if they errored.
