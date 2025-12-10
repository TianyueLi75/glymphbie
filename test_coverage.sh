#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=test_coverage
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=2g
#SBATCH --time=00:30:00
#SBATCH --account=ners570f25_class
#SBATCH --partition=standard
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd build &&
cmake .. -DENABLE_GCOV=ON &&
cmake --build . -j18 &&
ctest -V &&
cd CMakeFiles/GLYMPHBIE.dir/src/ &&
gcov Annular.cpp.gcno