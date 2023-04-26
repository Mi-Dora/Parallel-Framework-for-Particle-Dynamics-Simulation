# Parallel-Framework-for-Particle-Dynamics-Simulation
CSE 6230 Course Project @ Georgia Tech
# Prepare the environment
## For the MPI + OpenMP environment
'''
module load gcc/6.2.0
module load cuda
module load cmake
'''
## For the CUDA environment
'''
module load gcc/8.3.0 openblas/0.3.13 cmake/3.20.3 openmpi/4.1.2 cuda/11.1
'''


## Checkout the source code
Clone the distfw2d, semiring-gemm and cutlass code as follows
```
git clone --recurse-submodules https://code.ornl.gov/gordon-bell-apsp/distfw2d.git
```

## To Compile
```bash
#Make the build directory
mkdir build;
cd build;
#Run cmake
cmake ../
# Run make to compile 
make 
```

In `cori.nersc.gov`, I should use the following `cmake` command.

```
cmake -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment ../
## Or
cmake  -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DCMAKE_BUILD_TYPE=Debug ..
```

## Download the data files

1. Download the dataset from sparse tamu dataset
2. Download the dataset from semantic medline db
3. Follow the procedure in https://code.ornl.gov/gordon-bell-apsp/data.git for preparing the matrix out of semantic medline DB.

## Processing Input for the run
1. Under src/utils/, build and run permute-sparse-struct.cpp to randomly permute the input real world matrix. Random permutation helps balance the load by distorting the structure of the matrix.
2. 2D partition the permuted matrix with src/utils/partition2d

## To Validate
The name of the main test executable is `testBenchDistFW` in the `build` directory.
We use ctest as the test runner to run all tests in `build/tests/`

## To Run
The name of the executable will be `benchDistFW` in the `build` directory.

Running the benchmark in `Pr X Pc` 2D process-grid with `blockSize` as paramater for block cyclic data distribution, and `problemSize`

```
bash
cd build
mpirun -n PrXPc ./benchDistFW -r Pr -c Pc -b blockSize problemSize
```

For running on Summit, sample LSF files have been shipped under run directory. Please use the scripts and modify for your run. 

## Version of the code (-x contain versions: )

1. -x 0 (CPU)
2. -x 1 (CPU) with with lookahead
3. -x 2 (GPU) with lookahead
4. -x 3 (GPU) with lookahead + unified memory
5. -x 4 (GPU) with lookahead + RDMA
6. -x 5 (GPU) with lookahead + RDMA + Stream
7. -x 6 (GPU) with lookahead + RDMA + Stream + ring bcast

## Debug option (-y have different debug option)

1. -y 0 (no debug)
2. -y 1 (strong debug: compare with serial code)
3. -y 2 (medium debug: compare with distributed CPU code)
4. -y 3 (weak debug: compare with itself on 2 different runs)
5. -y 9 (no debug, but run all GPU benchmark)
******* !! ******** Becareful !! ALL debug flag except 0 and 1 will use double the GPU
memory

## Benchmarking

To actually benchmarking the benifit, please remove cudaDeviceSynchronize()
in your timer

To use the ERF_file, go to gen_rankfile.cpp, compile it with g++

Run ./a.out pr pc v b x #ofTotalNode nr nc proc_per_gpu

nr and nc are node grid. that is if nr=3 nc=2, then each node will grab process in this chunck
pr MUST BE DEVISIBE by nr and pc%nc MUST = 0

proc_per_gpu is currently support 1, and 2
