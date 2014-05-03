#! /bin/bash
#integrate cubes
#$ -V
#$ -N Cube_integrator
#$ -S /bin/bash

#inputs needed: file_path_cubes, obs_list_path, version, chunk, nslots
# chunk is the chunk number when the list was broken up. 0 for "master" or only chunk

echo JOBID ${JOB_ID}

(( nperint = $nslots/2 ))

# Get list of obs to integrate, turn them into file names to integrate
even_file_paths=${file_path_cubes}/${version}_int_chunk${chunk}_even_list.txt
odd_file_paths=${file_path_cubes}/${version}_int_chunk${chunk}_odd_list.txt
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
    save_file="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_even_cube.sav
else
    save_file="$file_path_cubes"/Healpix/Combined_obs_${version}_even_cube.sav
fi
/usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$even_file_paths" "$save_file" &
int_pids+=( $! )
if [ "$chunk" -gt "0" ]; then
    save_file="$file_path_cubes"/Healpix/Combined_obs_${version}_int_chunk${chunk}_odd_cube.sav
else
    save_file="$file_path_cubes"/Healpix/Combined_obs_${version}_odd_cube.sav
fi
/usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$odd_file_paths" "$save_file" &
int_pids+=( $! )
wait ${int_pids[@]} # Wait for integration to finish before making PS

# TODO: check that they completed - exit with code 100 if they errored.
