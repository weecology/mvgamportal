#!/bin/bash

# Job name and who to send updates to
#SBATCH --job-name=portal_sliding_window
#SBATCH --mail-user=ethanwhite@ulf.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --account=ewhite
#SBATCH --qos=ewhite

# Where to put the outputs: %j expands into the job number (a unique identifier for this job)
#SBATCH --output my_job%j.out
#SBATCH --error my_job%j.err

# Number of nodes to use
#SBATCH --nodes=1

# Number of tasks (usually translate to processor cores) to use: important! this means the number of mpi ranks used, useless if you are not using Rmpi)
#SBATCH --ntasks=1

#number of cores to parallelize with:
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000

# Job run time in [DAYS]
# HOURS:MINUTES:SECONDS
# [DAYS] are optional, use when it is convenient
#SBATCH --time=96:00:00

# Save some useful information to the "output" file
date;hostname;pwd

# Load R and run a script named my_R_script.R
ml R
Rscript mvgam_install.R
Rscript sliding_window_fits.R
