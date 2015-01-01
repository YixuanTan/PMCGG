#!/bin/bash
#
# srun command for CSCI-6360 group project.
# USAGE: /full/path/to/./q_GG.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --mail-type=END
#SBATCH --mail-user=tany3@rpi.edu
#SBATCH -D /gpfs/u/home/MMPM/MMPMtany/scratch/1000GridData/
#SBATCH --partition small
#SBATCH -t 360
#SBATCH -N 64
#SBATCH -n 1024
#SBATCH --overcommit
#SBATCH -o /gpfs/u/home/MMPM/MMPMtany/scratch/thinfilm/proj_64_1024.log

srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/MMPM/MMPMtany/barn/MCgraingrowth/thinfilm/q_GG.out --nonstop 2 ./voronoi.000.dat 100 1 100 4

