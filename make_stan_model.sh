#!/bin/sh
#SBATCH --job-name=af_make_model
#SBATCH --mail-user=andrew.farkas@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output /blue/akeil/andrew.farkas/logs/af_make_model%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=00:05:00

ml purge

ml gcc/12.2.0
ml openmpi/4.1.6
ml stan

cd /blue/akeil/andrew.farkas/cmdstan-2.36.0/

make /blue/akeil/andrew.farkas/repository/gaborgen/stan_models/fMRI/Model014


#make /blue/akeil/andrew.farkas/repository/gaborgen/stan_models/fMRI/Model013

#make /blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/repository/gaborgen/stan_models/fMRI/Model012


