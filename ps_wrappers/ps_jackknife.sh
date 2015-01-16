#! /bin/bash
# A script to make a jackknife cut, integrate, and make the PS
# Read a config file to make sql query and select obsids
# Must pass in a config file - full list of accepted variables is at the end of this file.

config_file=$1
if [ "${config_file}" == "" ]; then
  echo Must supply a config file. Exiting.
  exit 1
fi
filename=$(basename "${config_file}")
obs_name="${filename%.*}"

# config can be passed as full path, relative path, or assumed to be in PS repo.
if [ ! -f "${config_file}" ]; then
  # find the PSpath
  PSpath=$(idl -e 'print,rootdir("ps")') ### NOTE this only works if idlstartup doesn't have any print statements (e.g. healpix check)
  if [ -f "${PSpath}ps_wrappers/configs/${config_file}" ]; then
    config_file=${PSpath}ps_wrappers/configs/${config_file}
  elif [ -f "${PSpath}ps_wrappers/configs/${config_file}.cfg" ]; then
    config_file=${PSpath}ps_wrappers/configs/${config_file}.cfg
  else
    echo Cannot find config file. Exiting.
    exit 1
  fi
fi
echo Using config file ${config_file}

# Load in config file and test things
. ${config_file}
# set defaults
if [ -z ${FHDdir} ]; then FHDdir=$(pwd); fi
if [ -z ${ps_only} ]; then ps_only=0; fi
if [ -z ${priority} ]; then priority=0; fi
if [ -z ${wallclock_time} ]; then wallclock_time=07:00:00; fi
if [ -z ${nslots} ]; then nslots=10; fi
if [ -z ${mem} ]; then mem=4G; fi
if [ -z ${hold} ]; then hold=1; fi

if [ ! -z ${integrate_list} ]; then
  # obsid list was passed in. Just launch the old script.
  ${PSpath}ps_wrappers/ps_script.sh -d ${FHDdir} -f ${integrate_list} -p ${priority} -w ${wallclock_time} -n ${nslots} -m ${mem} -o ${ps_only} -h ${hold}
  exit $?
fi
# obsid list was not specified. Time to make one from the cubes that exist and db querires.
cube_list=$(ls ${FHDdir}/Healpix/1*odd_cubeYY.sav)
# full_obs_list is formatted to use in an "in" sql statement. It has the form (obs1,obs2,...,obsN)
full_obs_list="("
for cube in ${cube_list}; do
  full_obs_list=${full_obs_list}$(basename $cube _odd_cubeYY.sav),
done
full_obs_list="$(echo ${full_obs_list} | rev | cut -c 2- | rev))"

## Now we cut down the list based on database options from config file. For now we have two tables to select from. We could probably do this in one go, but I'm not that good at SQL. So start with qs table.
export PGPASSWORD=BowTie
psql_call="select obsid from qs where obsid in ${full_obs_list}"
if [ ! -z ${window_x_min} ]; then psql_call="${psql_call} and window_x >= ${window_x_min}"; fi
if [ ! -z ${window_x_max} ]; then psql_call="${psql_call} and window_x <= ${window_x_max}"; fi
if [ ! -z ${window_y_min} ]; then psql_call="${psql_call} and window_y >= ${window_y_min}"; fi
if [ ! -z ${window_y_max} ]; then psql_call="${psql_call} and window_y >= ${window_y_max}"; fi
if [ ! -z ${window_tot_min} ]; then psql_call="${psql_call} and window_x+window_y >= ${window_tot_min}"; fi
if [ ! -z ${window_tot_max} ]; then psql_call="${psql_call} and window_x+window_y <= ${window_tot_max}"; fi
if [ ! -z ${wedge_res_x_min} ]; then psql_call="${psql_call} and wedge_res_x >= ${wedge_res_x_min}"; fi
if [ ! -z ${wedge_res_x_max} ]; then psql_call="${psql_call} and wedge_res_x <= ${wedge_res_x_max}"; fi
if [ ! -z ${wedge_res_y_min} ]; then psql_call="${psql_call} and wedge_res_y >= ${wedge_res_y_min}"; fi
if [ ! -z ${wedge_res_y_max} ]; then psql_call="${psql_call} and wedge_res_y >= ${wedge_res_y_max}"; fi
if [ ! -z ${wedge_res_tot_min} ]; then psql_call="${psql_call} and wedge_res_x+wedge_res_y >= ${wedge_res_tot_min}"; fi
if [ ! -z ${wedge_res_tot_max} ]; then psql_call="${psql_call} and wedge_res_x+wedge_res_y <= ${wedge_res_tot_max}"; fi
if [ ! -z ${gal_wedge_x_min} ]; then psql_call="${psql_call} and gal_wedge_x >= ${gal_wedge_x_min}"; fi
if [ ! -z ${gal_wedge_x_max} ]; then psql_call="${psql_call} and gal_wedge_x <= ${gal_wedge_x_max}"; fi
if [ ! -z ${gal_wedge_y_min} ]; then psql_call="${psql_call} and gal_wedge_y >= ${gal_wedge_y_min}"; fi
if [ ! -z ${gal_wedge_y_max} ]; then psql_call="${psql_call} and gal_wedge_y >= ${gal_wedge_y_max}"; fi
if [ ! -z ${gal_wedge_tot_min} ]; then psql_call="${psql_call} and gal_wedge_x+gal_wedge_y >= ${gal_wedge_tot_min}"; fi
if [ ! -z ${gal_wedge_tot_max} ]; then psql_call="${psql_call} and gal_wedge_x+gal_wedge_y <= ${gal_wedge_tot_max}"; fi
if [ ! -z ${ptsrc_wedge_x_min} ]; then psql_call="${psql_call} and ptsrc_wedge_x >= ${ptsrc_wedge_x_min}"; fi
if [ ! -z ${ptsrc_wedge_x_max} ]; then psql_call="${psql_call} and ptsrc_wedge_x <= ${ptsrc_wedge_x_max}"; fi
if [ ! -z ${ptsrc_wedge_y_min} ]; then psql_call="${psql_call} and ptsrc_wedge_y >= ${ptsrc_wedge_y_min}"; fi
if [ ! -z ${ptsrc_wedge_y_max} ]; then psql_call="${psql_call} and ptsrc_wedge_y >= ${ptsrc_wedge_y_max}"; fi
if [ ! -z ${ptsrc_wedge_tot_min} ]; then psql_call="${psql_call} and ptsrc_wedge_x+ptsrc_wedge_y >= ${ptsrc_wedge_tot_min}"; fi
if [ ! -z ${ptsrc_wedge_tot_max} ]; then psql_call="${psql_call} and ptsrc_wedge_x+ptsrc_wedge_y <= ${ptsrc_wedge_tot_max}"; fi
# use psql_call to trim the list
# if you didn't make any wedge cuts, skip this part because obsids aren't in the table until wedgie has done its thing
if [ "${psql_call}" != "select obsid from qs where obsid in ${full_obs_list}" ]; then
  temp=`psql -h eor-00 mwa_qc -U mwa -c "${psql_call};" -t -A`
  if [ -z "${temp}" ]; then
    echo Quality statistics cut removed all obses. Quitting.
    exit 1
  fi
  # repackage again
  full_obs_list="("
  for obs in ${temp}; do
    full_obs_list=${full_obs_list}${obs},
  done
  full_obs_list="$(echo ${full_obs_list} | rev | cut -c 2- | rev))"
else
  echo No quality statistics cut to be made. Moving on.
fi

# Now cut things based on observation options
psql_call="select obsid from observation_info where obsid in ${full_obs_list}"
# Julian date
if [ ! -z ${jd} ]; then
  psql_call="${psql_call} and round(mjd) = ${jd}-2400000";
else
  if [ ! -z ${jd_min} ]; then psql_call="${psql_call} and round(mjd) >= ${jd_min}-2400000"; fi
  if [ ! -z ${jd_max} ]; then psql_call="${psql_call} and round(mjd) >= ${jd_max}-2400000"; fi
fi
# Pointing
if [ ! -z ${pointing} ]; then
  # need to translate between notations
  if [ "${pointing}" -ge "0" ]; then (( pointing=2*pointing )); else (( pointing=-2*pointing-1 )); fi
  psql_call="${psql_call} and gridnum = ${pointing}"
elif [ ! -z ${pointing_min} ] || [ ! -z ${pointing_max} ]; then
  # construct list of allowed pointings to pass into query
  if [ -z ${pointing_max} ]; then pointing_max=10; fi
  if [ -z ${pointing_min} ]; then pointing_min=-10; fi
  points="("
  for i in $(seq ${pointing_min} ${pointing_max}); do
    if [ "${i}" -ge "0" ]; then (( j=2*i )); else (( j=-2*i-1 )); fi
    points=${points}${j},
  done
  points="$(echo ${points} | rev | cut -c 2- | rev))"
  psql_call="${psql_call} and gridnum in ${points}"
fi
# TODO: Field
# use psql_call to trim the list
temp=`psql -h eor-00 mwa_qc -U mwa -c "${psql_call} order by obsid asc;" -t -A`
if [ -z "${temp}" ]; then
  echo Observation info cut removed all obses. Quitting.
  exit 1
fi

# Create a text file with the remaining observations
obsfile=${FHDdir}/Healpix/${obs_name}.txt
if [ -f ${obsfile} ]; then rm ${obsfile}; fi
for obs in $temp; do
    echo $obs >> $obsfile
done
echo Obs list located at $obsfile
echo Running power spectrum code with `wc -l < ${obsfile}` observations 

# Launch ps script
${PSpath}ps_wrappers/ps_script.sh -d ${FHDdir} -f ${obsfile} -p ${priority} -w ${wallclock_time} -n ${nslots} -m ${mem} -o ${ps_only} -h ${hold}



#new_list=`psql -h eor-00 mwa_qc -U mwa -c "select obsid from qs where obsid in ${full_obs_list};" -t -A`

#### Full list of accepted variables in config files
### Script and execution options
# integrate_list		: txt file of obs ids for subcubes. If this is provided, script checks that all obses have cubes, and skips the database options.
# FHDdir (cwd)			: The path to the fhd directory with cubes
# priority (0)			: priority level for grid engine
# wallclock_time (7:00:00)	: Time for execution in grid engine
# nslots (10)			: Number of slots for grid engine
# mem	(4G)			: Memory per core for grid engine
# ps_only (0)			: Flag for skipping integration to make PS only
# (legacy)			: (obsolete) in previous scripts this handled the old FHD directory structure
# hold				: Hold for a job to finish before running. Useful when running immediately after firstpass
### database options
## Danny stats options
# (window/wedge_res/gal_wedge/ptsrc_wedge)_(x/y/tot)_(min/max)    : various combinations of power in different places, polarizations, and min/max.
## Observation options
# pointing_min                  : Minimum pointing to include (counting negative before zenith, zenith=0, positive after zen)
# pointing_max                  : Maximum pointing
# pointing                      : Only pointing to include. Supercedes the previous two options
# jd_min                        : Minimum julian day to include
# jd_max                        : Maximum julian day to include
# jd                            : Only jd to include. Supercedes previous two options
# field                         : EoR field (NOT YET SUPPORTED)


