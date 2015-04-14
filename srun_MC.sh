#!/bin/bash
#
# srun command for CSCI-6360 group project.
# USAGE: /full/path/to/./q_MC.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --mail-type=END
#SBATCH --mail-user=tany3@rpi.edu
#SBATCH -D /gpfs/u/home/MMPM/MMPMtany/scratch/nCompare/data
#SBATCH --partition debug
#SBATCH -t 60
#SBATCH -N 16
#SBATCH -n 1024
#SBATCH --overcommit
#SBATCH -o /gpfs/u/home/MMPM/MMPMtany/scratch/nCompare/prMC_16_1024.log

srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/MMPM/MMPMtany/scratch/nCompare/q_MC.out --nonstop 2 voronmc.0000.dat 10000 5000 0 2 773 673

