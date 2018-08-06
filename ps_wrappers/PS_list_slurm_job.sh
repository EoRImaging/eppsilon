#! /bin/bash
#run Bryna's PS code on cubes that have already been integrated (in Slurm)
#SBATCH -J PS_B
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=adam_lanman@brown.edu

#inputs needed: file_path_cubes, nobs, version, ncores

#echo "Job ID: " $SLURM_ARRAY_JOB_ID

#nobs=0
#while read line
#do
#  ((nobs++))
#done < $obs_list_path

#Make power spectra through a ps wrapper in idl

module load ghostscript
module load git/2.10.2
module load imagemagick
module load idl
shopt -s expand_aliases; source $IDL/envi53/bin/envi_setup.bash


if [  -z ${SLURM_ARRAY_TASK_ID} ]
then
    ### This is being run from ps_slurm.sh, and should only be a single ObsID
    echo "Integrated mode"
    input_file=${file_path_cubes}/
    idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $ncores -e slurm_ps_job -args $input_file $version # $nobs
else
    ### This is run as an array job from multi_ps_slurm.sh, and each task should choose a single ObsID
    echo "Separated mode"
    obsids=("$@")
    obsid=${obsids[$SLURM_ARRAY_TASK_ID]}
    echo $obsid
    idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $ncores -e slurm_reduced_ps_job -args ${file_path_cubes} $obsid
    ### Delete unneeded files:
    psdir=${file_path_cubes}"/ps"
    rm $psdir"/data/1d_binning/*90deg*"
fi
