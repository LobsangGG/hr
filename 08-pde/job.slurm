#!/bin/bash
#SBATCH --job-name=partdiff_mpi
#SBATCH --output=partdiff_%j.out
#SBATCH --error=STDOUT
#SBATCH --nodes=2
#SBATCH --ntasks=24
#SBATCH --partition=west
#SBATCH --time=01:30:00

echo "=== Abbruchbedingung: Iterationen ==="

echo "=== RUN 1: 1 Prozess (sequenziell) ==="
mpirun -np 1 ./partdiff.x 1 2 512 2 2 800

echo "=== RUN 2: 12 Prozesse (parallel) ==="
mpirun -np 12 ./partdiff.x 1 2 512 2 2 800

echo "=== RUN 3: 24 Prozesse (parallel) ==="
mpirun -np 24 ./partdiff.x 1 2 512 2 2 800

echo "=== Abbruchbedingung: Präzision ==="

echo "=== RUN 4: 1 Prozess (sequenziell) ==="
mpirun -np 1 ./partdiff.x 1 2 32 2 1 1e-7

echo "=== RUN 5: 12 Prozesse (parallel) ==="
mpirun -np 12 ./partdiff.x 1 2 32 2 1 1e-7 

echo "=== RUN 6: 24 Prozesse (parallel) ==="
mpirun -np 24 ./partdiff.x 1 2 32 2 1 1e-7