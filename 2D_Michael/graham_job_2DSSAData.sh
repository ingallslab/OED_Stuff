#!/bin/bash 
#SBATCH --ntasks=1                # number of tasks
#SBATCH --cpus-per-task=1         # number of cores per task (adjust this when using PCT)
#SBATCH --mem-per-cpu=2096M       # memory; default unit is megabytes 
#SBATCH --time=0-00:20            # time (DD-HH:MM) 
#SBATCH --input=matlab_test.m     # input file
#SBATCH --output=matlab_test.log  # output file
module load matlab/2017a 
srun matlab -nodisplay -nosplash -nojvm -singleCompThread