#!/bin/bash

# Job name and who to send updates to
#SBATCH --job-name=4casttime
#SBATCH --mail-user=ethanwhite@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --account=ewhite
#SBATCH --qos=ewhite

# Where to put the outputs: %j expands into the job number (a unique identifier for this job)
#SBATCH --output 4casttime%j.out
#SBATCH --error 4casttime%j.err

# Number of nodes to use
#SBATCH --nodes=1

# Number of tasks (usually translate to processor cores) to use: important! this means the number of mpi ranks used, useless if you are not using Rmpi)
#SBATCH --ntasks=26

# Memory
#SBATCH --mem=4gb

# Job run time in
# [DAYS]:HOURS:MINUTES:SECONDS
# [DAYS] are optional, use when it is convenient
#SBATCH --time=24:00:00

# Save some useful information to the "output" file
date;hostname;pwd

# Load R and run a script named my_R_script.R
ml R
ml pandoc
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) {install.packages("renv")}'
Rscript -e 'renv::restore(prompt=FALSE)'
Rscript -e 'cmdstanr::check_cmdstan_toolchain()'
Rscript -e 'cmdstanr::install_cmdstan()'
Rscript sliding_window_fits.R
