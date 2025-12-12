#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=eccentricity_sweep
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=3g
#SBATCH --time=00:40:00
#SBATCH --account=ners570f25_class
#SBATCH --partition=standard
#SBATCH --exclusive
#SBATCH --output=out/eccentricity_sweep/ecc_sweep_CC_full.out

# cpus per task for OpenMP originally 18
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))
module load gcc/13.2.0 mkl/2023.2.1 openmpi/5.0.7
cd build &&
cmake --build . -j18 &&

# ===================== Eccentricity Sweep Test (Constant/Constant Walls) =====================
# The C++ executable arguments are:
# ./bin/FuncWall_SpatialSolver_ex <Nelem> <ElemOrder> <FourierOrder> <eccentricity_ratio> <inner_func> <outer_func>

# --- Fixed parameters for the sweep ---
NELEM=4       # Number of elements (N) - Keep moderate
ELEM_ORDER=10 # Element Order (O) - Keep high for accuracy
FOURIER_ORDER=24 # Fourier Order (F) - Keep high for accuracy
INNER_FUNC=0  # Inner Wall: Constant Radius (0)
OUTER_FUNC=0  # Outer Wall: Constant Radius (0)
# parameters are defined in FuncWall_Eccentric_ex.cpp

# eccentricity_ratio is e / (Rout - Rin). Max safe ratio < 1.0.
# error checking for sinusoidal walls is done in FuncWall_Eccentric_ex.cpp

echo "Starting Eccentricity Sweep for Constant/Constant (CC) Walls..."

# 1) Sweep over eccentricity_ratio from 0.0 to 0.9 in steps of 0.1
eccentricity_ratios=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95)
#eccentricity_ratios=(0.0 0.5)

for ratio in "${eccentricity_ratios[@]}"; do
    echo "Running with eccentricity_ratio = $ratio"
    
    # Run with 1 MPI rank (for simpler I/O and as discussed, it's safer for the Q calculation)
    mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} \
    ./bin/FuncWall_Eccentric_ex $NELEM $ELEM_ORDER $FOURIER_ORDER "$ratio" $INNER_FUNC $OUTER_FUNC
    
    # Note: Using quotes for "$ratio" in case it's a fractional value (e.g., 0.9)
done

echo "Eccentricity Sweep complete."