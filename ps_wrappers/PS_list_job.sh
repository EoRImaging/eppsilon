#! /bin/bash
#run Bryna's PS code on cubes that have already been integrated
#$ -V
#$ -N PS_B
#$ -S /bin/bash

#inputs needed: file_path_cubes, obs_list_path, version, nslots

echo JOBID ${JOB_ID}

nobs=0
while read line
do
  ((nobs++))
done < $obs_list_path

#Make power spectra through a ps wrapper in idl

input_file=${file_path_cubes}/

idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e mit_ps_job -args $input_file $version

