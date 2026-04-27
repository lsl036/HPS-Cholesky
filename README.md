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


## Cite this work
If you find our work useful in your research, please consider citing:
```
@article{10.1145/3630051,
author = {Lin, Shengle and Yang, Wangdong and Hu, Yikun and Cai, Qinyun and Dai, Minlu and Wang, Haotian and Li, Kenli},
title = {HPS Cholesky: Hierarchical Parallelized Supernodal Cholesky with Adaptive Parameters},
year = {2024},
issue_date = {March 2024},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {11},
number = {1},
issn = {2329-4949},
url = {https://doi.org/10.1145/3630051},
doi = {10.1145/3630051},
journal = {ACM Trans. Parallel Comput.},
month = mar,
articleno = {3},
numpages = {22},
keywords = {Supernodal Cholesky factorization, graph convolutional network, hierarchical parallelization, task stream processing, multi-NUMA architecture}
}
```










