# Parallel-Framework-for-Particle-Dynamics-Simulation
## Solution I: MPI + OpenMP
### Prepare the environment
```
module load gcc cmake openmpi valgrind
module load gcc cmake openmpi valgrind
module load gcc cmake openmpi valgrind
module load gcc cmake openmpi valgrind
```
### Checkout the source code
```
git clone --branch neo https://github.com/Mi-Dora/Parallel-Framework-for-Particle-Dynamics-Simulation.git
```

### To Compile
The CMakeLists.txt is already provided in the home directory. Under the "Parallel-Framework-for-Particle-Dynamics-Simulation" directory, run:
```
make
```

### Download the data files
The particles are generated with randomly distributed positions, velocity, mass and acceleration. Under the "pysrc" directory, we have a "generate.py" script to do this task. Execute this scripte with the following command:
```
python generate.py --n_particle 1000 --n_dim 3 --pos_range [-100, 100, -100, 100] --mass_range [0, 0.01] --output_filename [output.txt]
```
and replace [output.txt] with your own output file name. Note that the "--output_filename" argument is required to run this script, while all other arguments are optional. If they are not provided, the parameters as listed in the command above will be the default setting. For a quick execution, run
```
python generate.py --output_filename output.txt
```
Note: If this does not work, you may try replace "python" with "python3" or "py" as per the running environment.

### Processing Input for the run
We have provided test scripts to do this part. One can refer to "test_dense_mpi_neo.cpp" (for the dense scheduler) or "test_cutoff_mpi_neo.cpp" (for the cutoff scheduler) or "test_centroid_mpi_neo.cpp" (for the centroid scheduler). Each of them will call the processing functions under src/neo, in the "dense_scheduler.cpp" or "cutoff_scheduler.cpp" or "centroid_scheduler.cpp" scripts respectively. For each type of scheduler, one can refer to the corresponding test scripts on the code required to process inputs.

We always start with calling the loadParticles() function, then the padding() function, then the update() function iteratively, and end with calling the gather() and shrink() functions. They are overloded in different scripts under src/neo. They will scatter particles to a 2D-grid processor topology in a similar manner.

### To Validate
### To Run
### Version of the code (-x contain versions: )
### Debug option (-y have different debug option)
### Benchmarking

## Solution II: CUDA
### Prepare the environment
```
module load gcc/8.3.0 
module load openblas/0.3.13 
module load cmake/3.20.3
module load openmpi/4.1.2 
module load cuda/11.1
```


### Checkout the source code
Clone the distfw2d, semiring-gemm and cutlass code as follows
```
git clone --branch cuda https://github.com/Mi-Dora/Parallel-Framework-for-Particle-Dynamics-Simulation.git
```

### To Compile
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

### Download the data files

1. Download the dataset from sparse tamu dataset
2. Download the dataset from semantic medline db
3. Follow the procedure in https://code.ornl.gov/gordon-bell-apsp/data.git for preparing the matrix out of semantic medline DB.

### Processing Input for the run
1. Under src/utils/, build and run permute-sparse-struct.cpp to randomly permute the input real world matrix. Random permutation helps balance the load by distorting the structure of the matrix.
2. 2D partition the permuted matrix with src/utils/partition2d

## To Validate
The name of the main test executable is `testBenchDistFW` in the `build` directory.
We use ctest as the test runner to run all tests in `build/tests/`

### To Run
The name of the executable will be `benchDistFW` in the `build` directory.

Running the benchmark in `Pr X Pc` 2D process-grid with `blockSize` as paramater for block cyclic data distribution, and `problemSize`

```
bash
cd build
mpirun -n PrXPc ./benchDistFW -r Pr -c Pc -b blockSize problemSize
```

For running on Summit, sample LSF files have been shipped under run directory. Please use the scripts and modify for your run. 

### Version of the code (-x contain versions: )

1. -x 0 (CPU)
2. -x 1 (CPU) with with lookahead
3. -x 2 (GPU) with lookahead
4. -x 3 (GPU) with lookahead + unified memory
5. -x 4 (GPU) with lookahead + RDMA
6. -x 5 (GPU) with lookahead + RDMA + Stream
7. -x 6 (GPU) with lookahead + RDMA + Stream + ring bcast

### Debug option (-y have different debug option)

1. -y 0 (no debug)
2. -y 1 (strong debug: compare with serial code)
3. -y 2 (medium debug: compare with distributed CPU code)
4. -y 3 (weak debug: compare with itself on 2 different runs)
5. -y 9 (no debug, but run all GPU benchmark)
******* !! ******** Becareful !! ALL debug flag except 0 and 1 will use double the GPU
memory

### Benchmarking

To actually benchmarking the benifit, please remove cudaDeviceSynchronize()
in your timer

To use the ERF_file, go to gen_rankfile.cpp, compile it with g++

Run ./a.out pr pc v b x #ofTotalNode nr nc proc_per_gpu

nr and nc are node grid. that is if nr=3 nc=2, then each node will grab process in this chunck
pr MUST BE DEVISIBE by nr and pc%nc MUST = 0

proc_per_gpu is currently support 1, and 2
