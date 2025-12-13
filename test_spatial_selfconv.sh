#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=test_self_conv
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=2g
#SBATCH --time=01:30:00
#SBATCH --account=ners570f25_class
#SBATCH --partition=standard
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd build &&
cmake .. -DENABLE_GCOV=OFF &&
cmake --build . -j18 &&

# 3) Concentric sinusoidal channels, check self-convergence.
for n in {4..1..-1}; do # Nelem
    for m in {32..4..-4}; do # FourierOrder
        	mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/SpatialSolver_selfconv $n 10 $m
    	done
done

