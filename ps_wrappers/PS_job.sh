#! /bin/bash
#integrate cubes and run Bryna's PS code
#$ -V
#$ -S /bin/bash

#inputs needed: file_path_cubes, obsid_in_name_even, starting_index, ending_index, nslots

echo JOBID ${JOB_ID}

(( nperint = $nslots/2 ))

#Integrate cubes
unset int_pids

input_file="$file_path_cubes"/"*even_cube.sav"
/usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$input_file" $starting_index $ending_index &
int_pids+=( $! )
input_file="$file_path_cubes"/"*odd_cube.sav"
/usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nperint -e integrate_healpix_cubes -args "$input_file" $starting_index $ending_index &
int_pids+=( $! )
wait ${int_pids[@]} # Wait for integration to finish before making PS

#Make power spectra through a ps wrapper in idl

input_file="$file_path_cubes""Healpix/"

idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e mit_ps_job -args $input_file '${obsid_in_name_even[$starting_index]}-${obsid_in_name_even[$ending_index]}'

