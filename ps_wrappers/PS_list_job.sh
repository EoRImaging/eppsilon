#! /bin/bash
#run Bryna's PS code on cubes that have already been integrated
#$ -V
#$ -S /bin/bash

#inputs needed: file_path_cubes, obs_list_path, version, nslots
#inputs optional: cube_type, pol, evenodd

echo JOBID ${JOB_ID}

nobs=0
while read line
do
  ((nobs++))
done < $obs_list_path

input_file=${file_path_cubes}/

#Default priority if not set.
if [ -z ${cube_type} ] && [ -z ${pol} ] && [ -z ${evenodd} ]; then
    #Make power spectra through a ps wrapper in idl, running all cubes in succession
    idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e mit_ps_job -args $input_file $version
else
    if [ ! -z ${cube_type} ] && [ ! -z ${pol} ] && [ ! -z ${evenodd} ]; then
        #Make power spectra through a ps wrapper in idl, running all cubes in parallel
        idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e mit_ps_job -args $input_file $version $cube_type $pol $evenodd
    else
        echo "Need to specify cube_type, pol, and evenodd altogether"
        exit 1
    fi
fi
