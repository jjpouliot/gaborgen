#!/bin/sh
#SBATCH --job-name=af_model014
#SBATCH --mail-user=andrew.farkas@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output /blue/akeil/andrew.farkas/logs/af_model014_slurm_array_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=3
#SBATCH --ntasks-per-node=4
#SBATCH --qos=akeil # 18 or 180 cpu cores with boost flag
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=2-00:00:00 # max is 31 days or 4 days with boost
#SBATCH --constraint='rome|milan'
#SBATCH --array=1

module purge

#export LIBRARY_PATH=/apps/opencl/intel/system_studio_2019/opencl/SDK/lib64:$LIBRARY_PATH

ml gcc/12.2.0
ml openmpi/4.1.6
ml opencl/18.1.0.013
ml stan

export OPENCL_VENDOR_PATH=/apps/opencl/pyopencl/2018.2.5/etc/OpenCL/vendors

clinfo -l          # prints the platforms / devices

CHAIN=$((SLURM_ARRAY_TASK_ID))
CHAIN_TEXT="$SLURM_ARRAY_TASK_ID"

JOBID_INT=$((SLURM_ARRAY_JOB_ID))
JOBID_TEXT="$SLURM_ARRAY_JOB_ID"

# pocl controls
#export POCL_MAX_PTHREAD_COUNT=3    # exactly 10 workers per rank
#export POCL_AFFINITY=1              # pin each worker

# optional visibility
#export POCL_DEBUG=timing
#export POCL_MSG_LEVEL=2

# quick visibility
export POCL_MAX_PTHREAD_COUNT=3
export POCL_AFFINITY=1
export POCL_DEBUG=all          # ← loudest setting
export POCL_MSG_LEVEL=3
echo "Batch shell POCL vars:"
env | grep ^POCL_

srun --mpi=${HPC_PMIX} \
     --export=ALL,POCL_DEBUG,POCL_MSG_LEVEL,POCL_MAX_PTHREAD_COUNT,POCL_AFFINITY \
     stdbuf -oL -eL \
     /blue/akeil/andrew.farkas/repository/gaborgen/stan_models/fMRI/Model014 \
     sample \
        num_warmup=1000 \
        num_samples=1000 \
        num_chains=1 \
        save_warmup=true \
    data file=/blue/akeil/andrew.farkas/data/fmri_stan_list.json \
opencl platform=0 device=0 \
random seed=$(((CHAIN*100) + JOBID_INT)) \
output file=/blue/akeil/andrew.farkas/chains/model014_chain_${JOBID_TEXT}_${CHAIN_TEXT}.csv \
refresh=10
