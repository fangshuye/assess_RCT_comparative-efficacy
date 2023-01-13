#!/bin/bash

#SBATCH --time=02:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out 
#SBATCH --mail-user=fangshu.stat@gmail.com   # email address
#SBATCH --mail-type=END
module load r 
module load gcc
R CMD BATCH check_sim_added_effect.R
exit

