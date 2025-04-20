#!/bin/sh
#SBATCH --job-name=af_model016
#SBATCH --mail-user=andrew.farkas@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output /blue/akeil/andrew.farkas/logs/af_model016_slurm_array_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=3
#SBATCH --ntasks-per-node=15
#SBATCH --qos=akeil-b # 18 or 180 cpu cores with boost flag
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=3-20:00:00 # max is 31 days or 4 days with boost
#SBATCH --constraint='rome|milan'
#SBATCH --array=1-4%1

ml stan
ml gcc/12.2.0
ml openmpi/4.1.6

CHAIN=$((SLURM_ARRAY_TASK_ID))
CHAIN_TEXT="$SLURM_ARRAY_TASK_ID"

JOBID_INT=$((SLURM_ARRAY_JOB_ID))
JOBID_TEXT="$SLURM_ARRAY_JOB_ID"


srun --mpi=${HPC_PMIX} /blue/akeil/andrew.farkas/repository/gaborgen/stan_models/fMRI/Model016 \
    sample \
        num_warmup=1000 \
        num_samples=1000 \
        num_chains=1 \
        save_warmup=true \
    data file=/blue/akeil/andrew.farkas/data/fmri_stan_list.json \
    opencl platform=0 device=0 \
    random seed=$(((CHAIN*100) + JOBID_INT)) \
    output file=/blue/akeil/andrew.farkas/chains/model016_chain_${JOBID_TEXT}_${CHAIN_TEXT}.csv \
    refresh=10
