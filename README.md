# HPS-Cholesky

The source code of paper: "HPS Cholesky: Hierarchical Parallelized Supernodal Cholesky with Adaptive Parameters"

## Code Environment

1. In the NUMA architecture, the **libnuma** library is required
2. **OpenBLAS-0.3.9** is recommended. You can obtain it from https://github.com/xianyi/OpenBLAS. 
3. **LAPACKE** library is also required.

## Setup

- Need thread scheduling framework TPSM we implemented 

  ```
  cd TPSM_v2.6
  make
  cd ..
  ```

  

- HPS Cholesky can be compiled with minimal dependencies in the usual CMake-way, e.g.:


```
mkdir build && cd build
cmake ..
make
```
Then the test routine is in `build/test/`










