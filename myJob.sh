#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=mpi_test
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=3g
#SBATCH --time=00:05:00
#SBATCH --account=ners570f25_class
#SBATCH --partition=standard
#SBATCH --exclusive
#SBATCH --output=out/MPI_cylindrical.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd build &&
cmake --build . -j18 &&

# ===================== MPI tests ===========================================
# mpirun -n $NTASK --map-by slot:pe=${OMP_NUM_THREADS} ./bin/Annular_mpi_tests 

# # ===================== Convergence tests ===========================================
# 1) Fix Rin, Rout to be straight cylindrical channels, check convergence to true solution. -- and flux
for n in {1..2}; do # Nelem
    for m in {4..40..4}; do # FourierOrder
        	mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/SpatialSolver_tests $n 10 $m
    	done
done
# 2) Eccentric cylindrical channels, check convergence to true Flux.
for n in {1..3}; do # Nelem
    for m in {4..40..4}; do # FourierOrder
        	mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/SpatialSolver_tests $n 10 $m
    	done
done

# # 3) Concentric sinusoidal channels, check self-convergence.
# for n in {4..1..-1}; do # Nelem
#     for m in {32..4..-4}; do # FourierOrder
#         	mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/SpatialSolver_tests $n 10 $m
#     	done
# done
