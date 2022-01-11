#!/bin/bash
nvcc -c orientCUDA.cu -I$CONDA_PREFIX/include -I$CONDA_PREFIX/include/python3.9 -I./eigen -I./cutlass/include -I./cutlass/tools/util/include -L$CONDA_PREFIX/lib -L$CONDA_PREFIX/lib/python3.9/config-3.9-x86_64-linux-gnu/ -l python3.9 -lboost_python39 -lboost_system -lboost_numpy39  --compiler-options=-lstdc++,-fPIC,-O3 -o orient_cuda.o --expt-relaxed-constexpr

nvcc -c orient_ext.cpp -I$CONDA_PREFIX/include -I$CONDA_PREFIX/include/python3.9  -I./eigen -I./cutlass/include -I./cutlass/tools/util/include -L$CONDA_PREFIX/lib -L$CONDA_PREFIX/lib/python3.9/config-3.9-x86_64-linux-gnu/ -lpython3.9 -lboost_python39 -lboost_system  -lboost_numpy39 --compiler-options=-lstdc++,-fPIC,-O3 -o orient.o --expt-relaxed-constexpr

g++ -shared orient.o orient_cuda.o -L$CONDA_PREFIX/lib -L$CONDA_PREFIX/lib/python3.9/config-3.9-x86_64-linux-gnu/ -L/usr/local/cuda/lib64 -lpython3.9 -lboost_python39  -lboost_numpy39 -lcudart -o orient.so

