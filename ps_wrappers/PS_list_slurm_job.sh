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
module load git/2.2.1
module load imagemagick

input_file=${file_path_cubes}/
idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $ncores -e slurm_ps_job -args $input_file $version # $nobs
