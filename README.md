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

### Compile GlymphBIE tests
In the main directory, go through CMakeLists.txt to update source files being compiled: 
- under "# Create the library", add in source files to GLYMPHBIE_SOURCES; 
- under "# Source discovery, build tests", add in source files to TEST_SOURCES.
```bash
mkdir build && cd build
cmake ..
cmake --build .
```

Resolve any compilation errors from source code, then find the executables in build/bin, for example:
```bash
./bin/Annular_tests
```



Note: Makefile in main folder is NOT updated!