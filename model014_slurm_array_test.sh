#!/bin/sh
#SBATCH --job-name=af_model014
#SBATCH --mail-user=andrew.farkas@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output /blue/akeil/andrew.farkas/logs/af_model014_slurm_array_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=4
#SBATCH --qos=akeil # 18 or 180 cpu cores with boost flag
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=24:00:00 # max is 31 days or 4 days with boost
#SBATCH --array=1

ml purge

export LIBRARY_PATH=/apps/opencl/intel/system_studio_2019/opencl/SDK/lib64:$LIBRARY_PATH

ml stan
ml gcc/12.2.0
ml openmpi/4.1.6
ml opencl/18.1.0.013


CHAIN=$((SLURM_ARRAY_TASK_ID))
CHAIN_TEXT="$SLURM_ARRAY_TASK_ID"

JOBID_INT=$((SLURM_ARRAY_JOB_ID))
JOBID_TEXT="$SLURM_ARRAY_JOB_ID"


srun --mpi=${HPC_PMIX} \
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
