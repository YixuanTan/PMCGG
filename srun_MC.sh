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
#SBATCH -D /gpfs/u/home/ACME/ACMEtany/scratch/HalfHalf/AlCuFilm/
#SBATCH --partition small
#SBATCH -t 120
#SBATCH -N 32
#SBATCH -n 1024
#SBATCH --overcommit
#SBATCH -o /gpfs/u/home/ACME/ACMEtany/scratch/HalfHalf/AlCuFilm/Al_Cu.log

srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/ACME/ACMEtany/barn/thinfilm/q_MC.out --nonstop 2 voronmc.0000.dat 1500 10 0 1 673 723

