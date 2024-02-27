#!/bin/bash

#SBATCH --job-name= mvgam-trial-job
#SBATCH --output= output.log
#SBATCH --error= error.log
#SBATCH --partition= 
#SBATCH --nodes= 1
#SBATCH --ntasks-per-node= 1
#SBATCH --time= 00:30:00

module load R

Rscript dm-sample.R