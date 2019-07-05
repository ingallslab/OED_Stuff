#!/bin/bash -l
#SBATCH --job-name=SSA_Fitting_Test_1
#SBATCH --account=def-bingalls # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-01:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=4000             # adjust this according to your the memory requirement per node you need
#SBATCH --mail-user=mrastwoo@edu.uwaterloo.ca # adjust this to match your email address
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2018a
# Remove -singleCompThread below if you are using parallel commands:
srun matlab -nodisplay -singleCompThread -r "TwoD_Fit_V2"