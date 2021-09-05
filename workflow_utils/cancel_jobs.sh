#!/bin/sh

## use this script to cancel jobs on spartan that are > than job id you specify
## I found it here  https://unix.stackexchange.com/questions/424871/how-to-cancel-jobs-on-slurm-with-job-idjob-number-bigger-than-a-certain-number
## it is not optimal but for the moment it is okay

if [ -z "$1" ] ; then
    echo "Minimum Job Number argument is required.  Run as '$0 jobnum'"
    exit 1
fi

minjobnum="$1"

myself="$(id -u -n)"

for j in $(squeue --user="dvespasiani" --noheader --format='%i') ; do
  if [ "$j" -gt "$minjobnum" ] ; then
    scancel "$j"
  fi
done