#!/bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:30:00
#SBATCH --output=output_%j.txt
#SBATCH --partition=west

echo "=== 512 Interlines, 800 Iterationen, 24 Prozesse über 2 Knoten==="
srun --mpi=pmi2 ./partdiff.x 1 2 512 2 2 800
