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
#SBATCH -D /gpfs/u/home/MMPM/MMPMtany/scratch/thinfilm/stg2/000c
#SBATCH --partition small
#SBATCH -t 360
#SBATCH -N 32
#SBATCH -n 128
#SBATCH --overcommit
#SBATCH -o /gpfs/u/home/MMPM/MMPMtany/scratch/thinfilm/stg2/000c/prMC_32_128.log

srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/MMPM/MMPMtany/scratch/thinfilm/q_MC.out --nonstop 2 voronmc.00000.dat 250000 10000 2 
