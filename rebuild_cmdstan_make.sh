#!/bin/sh
#SBATCH --job-name=rebuild_cmd_make
#SBATCH --mail-user=andrew.farkas@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output /blue/akeil/andrew.farkas/logs/rebuild_cmd_make%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=00:15:00

ml purge

export LIBRARY_PATH=/apps/opencl/intel/system_studio_2019/opencl/SDK/lib64:$LIBRARY_PATH

ml gcc/12.2.0
ml openmpi/4.1.6
ml opencl/18.1.0.013
ml stan

cd /blue/akeil/andrew.farkas/cmdstan-2.36.0/

make clean-all

make build

