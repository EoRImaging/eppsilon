#! /bin/bash
#integrate cubes

#SBATCH -J cube_integrate
# # SBATCH --mail-type=ALL
# # SBATCH --mail-user=adam_lanman@brown.edu

#inputs needed: file_path_cubes, obs_list_path, version, chunk, ncores
# chunk is the chunk number when the list was broken up. 0 for "master" or only chunk

#echo JOBID ${JOB_ID}

(( nperint = $ncores/4 ))

echo ${file_path_cubes}' , '${version}' , '${chunk}' , '$obs_list_path


if [ "$legacy" -ne "1" ]; then
    # Get list of obs to integrate, turn them into file names to integrate
    evenXX_file_paths=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_evenXX_list.txt
    evenYY_file_paths=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_evenYY_list.txt
    oddXX_file_paths=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_oddXX_list.txt
    oddYY_file_paths=${file_path_cubes}/Healpix/${version}_int_chunk${chunk}_oddYY_list.txt
    # clear file paths
    rm $evenXX_file_paths
    rm $evenYY_file_paths
    rm $oddXX_file_paths
    rm $oddYY_file_paths
    nobs=0
#    echo $obs_list_path
    while read line
    do
	evenXX_file=${file_path_cubes}/Healpix/${line}_even_cubeXX.sav
	evenYY_file=${file_path_cubes}/Healpix/${line}_even_cubeYY.sav
	oddXX_file=${file_path_cubes}/Healpix/${line}_odd_cubeXX.sav
	oddYY_file=${file_path_cubes}/Healpix/${line}_odd_cubeYY.sav
	echo $evenXX_file >> $evenXX_file_paths
	echo $evenYY_file >> $evenYY_file_paths
	echo $oddXX_file >> $oddXX_file_paths
	echo $oddYY_file >> $oddYY_file_paths
	((nobs++))
    done < $obs_list_path

#Integrate cubes
    unset int_pids

    if [ "$chunk" -gt "0" ]; then
	save_file_evenXX="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_even_cubeXX.sav
	save_file_evenYY="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_even_cubeYY.sav
	save_file_oddXX="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_odd_cubeXX.sav
	save_file_oddYY="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_odd_cubeYY.sav
    else
	save_file_evenXX="$file_path_cubes"/Healpix/Combined_obs_${version}_even_cubeXX.sav
	save_file_evenYY="$file_path_cubes"/Healpix/Combined_obs_${version}_even_cubeYY.sav
	save_file_oddXX="$file_path_cubes"/Healpix/Combined_obs_${version}_odd_cubeXX.sav
	save_file_oddYY="$file_path_cubes"/Healpix/Combined_obs_${version}_odd_cubeYY.sav
    fi
# $! = PID of most recent background process
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$evenXX_file_paths" "$save_file_evenXX" &
    int_pids+=( $! )
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$evenYY_file_paths" "$save_file_evenYY" &
    int_pids+=( $! )
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$oddXX_file_paths" "$save_file_oddXX" &
    int_pids+=( $! )
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$oddYY_file_paths" "$save_file_oddYY" &
    int_pids+=( $! )
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
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$even_file_paths" "$save_file_even" &
    int_pids+=( $! )
    /usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$odd_file_paths" "$save_file_odd" &
    int_pids+=( $! )
    wait ${int_pids[@]} # Wait for integration to finish before making PS
fi

# TODO: check that they completed - exit with code 100 if they errored.
