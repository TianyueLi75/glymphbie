#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=convergence
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3g
#SBATCH --time=01:30:00
#SBATCH --account=ners570f25_class
#SBATCH --partition=standard
#SBATCH --exclusive
#SBATCH --output=out/concentric_poiseuille/convergence.out

# cpus per task for OpenMP originally 18
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))
module load gcc/13.2.0  mkl/2023.2.1  openmpi/5.0.7
cd build &&
cmake --build . -j18 &&

# ===================== MPI tests ===========================================
# mpirun -n $NTASK --map-by slot:pe=${OMP_NUM_THREADS} ./bin/Annular_mpi_tests 

# ===================== Convergence tests ===========================================
# 1) Fix Rin, Rout to be straight cylindrical channels, check convergence to true solution.
for n in {1..3}; do # Nelem
    for m in {4..40..4}; do # FourierOrder
        	mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/FuncWall_SpatialSolver_ex $n 10 $m 
    	done
done
