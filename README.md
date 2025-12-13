# GlymphBIE
A simulation platform for glymphatic flow in the perivascular spaces of penetrating arterioles coupled to neuronal activity.

## SETUP notes
### Getting the submodules
Either clone the whole repo with recurse-submodules:

```bash
git clone --recurse-submodules https://github.com/TianyueLi75/glymphbie.git
```

Or if repo is already cloned, use submodule update:
```bash
(git clone https://github.com/TianyueLi75/glymphbie.git)
cd glymphbie
git submodule update --init --recursive
```

To stay up-to-date with the linked submodules: 
```bash
git pull --recurse-submodules
```

### Compile PVFMM
Requires: C++17+ compiler with OpenMP; mpirun; autotools

If on Great Lakes: load the following modules: 
```bash
module list
    1) gcc/13.2.0   2) mkl/2023.2.1   3) openmpi/5.0.7
```

Build PVFMM:
```bash
cd glymphbie/extern/pvfmm
./autogen.sh
./configure
make -j
cd ../..
```

### Compile GlymphBIE
```bash
mkdir build && cd build
cmake ..
cmake --build .
```

You can find the example executables in build/bin, see bash scripts in the main directory for examples on how to call them using MPI while enabling OpenMP. In particular, you can find a convergence study for a set of concentric cylinders under constant pressure in `examples/FuncWall_SpatialSolver_ex.cpp`


