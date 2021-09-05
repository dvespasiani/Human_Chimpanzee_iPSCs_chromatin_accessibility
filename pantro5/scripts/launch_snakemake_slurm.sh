#!/bin/bash
# Partition for the job:
#SBATCH --partition=mig

# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="pantro5_atac_pipeline"

# The project ID which this job should run under:
#SBATCH --account="punim0586"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=400000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=10-23:00:00

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=dvespasiani@student.unimelb.edu.au

# output
#SBATCH --output=./slurm_report/slurm.out

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
source /usr/local/module/spartan_new.sh
module load web_proxy
module load r/4.0.0 
module load anaconda3/2020.07

## create tmp directory and export it to store all the tmp files generated during the pipeline
if [ ! -d /data/scratch/projects/punim0586/dvespasiani/tmp ]; then
  mkdir -p /data/scratch/projects/punim0586/dvespasiani/tmp;
fi
export TMPDIR=/data/scratch/projects/punim0586/dvespasiani/tmp

# The job command(s):
## I assume you've already created the conda atac environment

source activate atac
snakemake -j 999 --cluster-config env/cluster.yaml --cluster "sbatch -A {cluster.account} -t {cluster.time} \
 -p {cluster.partition} --nodes {cluster.nodes} --ntasks {cluster.ntasks} \
  --mem {cluster.mem} --mail-user {cluster.mail_user} --mail-type {cluster.mail_type} \
  --output {cluster.output}"