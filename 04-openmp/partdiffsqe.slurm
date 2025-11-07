#!/bin/bash

#SBATCH --job-name=parseq
#SBATCH --partition=west
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --output=output_%j.txt

srun ./partdiff-seq 1 2 512 2 2 500
