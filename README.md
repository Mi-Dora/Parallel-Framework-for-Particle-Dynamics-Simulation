# Parallel-Framework-for-Particle-Dynamics-Simulation
All the code on [This repository](https://github.com/Mi-Dora/Parallel-Framework-for-Particle-Dynamics-Simulation/edit/main/)
## Solution I: MPI + OpenMP
### Prepare the environment
```shell
module load gcc cmake openmpi valgrind
module load gcc cmake openmpi valgrind
module load gcc cmake openmpi valgrind
module load gcc cmake openmpi valgrind
```
(During practice, we had to load multiple times to avoid error)
### Checkout the source code
We need to clone the "neo" branch for this solution:
```shell
git clone --branch neo https://github.com/Mi-Dora/Parallel-Framework-for-Particle-Dynamics-Simulation.git
```

### To Compile
The CMakeLists.txt is already provided in the home directory. Under the "Parallel-Framework-for-Particle-Dynamics-Simulation" directory, run:
```shell
make
```

### Download the data files
The particles are generated with randomly distributed positions, velocity, mass and acceleration. Under the "pysrc" directory, we have a "generate.py" script to do this task. Execute this scripte with the following command:
```shell
python generate.py --n_particle 1000 --n_dim 3 --pos_range [-100, 100, -100, 100] --mass_range [0, 0.01] --output_filename [output.txt]
```
and replace [output.txt] with your own output file name. Note that the "--output_filename" argument is required to run this script, while all other arguments are optional. If they are not provided, the parameters as listed in the command above will be the default setting. For a quick execution, run
```shell
python generate.py --output_filename output.txt
```
Note: If this does not work, you may try replace "python" with "python3" or "py" as per the running environment.

### Processing Input for the run
The "generate.py" script is designed to minimize additional processing actions. So there is not much to process. We just directly run the test scripts, in which the data are loaded and parsed for computations.

### To Validate
One can compare the output files generated by different schedulers above and by the sequential version to verify that the sequential algorithm produces the same result as the parallel versions. For the sequential version, just set --np 1 in the .pbs running script as described in the section below.

### To Run
We have provided test scripts to do this part. One can refer to "test_dense_mpi_neo.cpp" (for the dense scheduler) or "test_cutoff_mpi_neo.cpp" (for the cutoff scheduler) or "test_centroid_mpi_neo.cpp" (for the centroid scheduler). Each of them will call the processing functions under src/neo, in the "dense_scheduler.cpp" or "cutoff_scheduler.cpp" or "centroid_scheduler.cpp" scripts respectively. For each type of scheduler, one can refer to the corresponding test scripts on the code required to process inputs.

We always start with calling the loadParticles() function, then the padding() function, then the update() function iteratively, and end with calling the gather() and shrink() functions. They are overloded in different scripts under src/neo. They will scatter particles to a 2D-grid processor topology in a similar manner.

In the home directory, there is a "run_mpi.pbs" script to submit work to the coc-ice-multi cluster. (It also prepare the environment). In the PACE environment, run
```shell
qsub run_mpi.pbs
```
And we're done!

Note: Double check the commands in the .pbs file. We have prepared some input data files for execution and you may need to reproduce them following the instructions above.

### Version of the code (-x contain versions: )
```shell
-x neo
```
This is the final MPI+OpenMP version.

### Debug option (-y have different debug option)
We don't have debug option for the MPI+OpenMP solution. It will be identical to the RELEASE version.

### Benchmarking
We have different branches for this repository; however, all schedulers are merged in the --neo branch so one branch is sufficient to validate out work in the first solution.

## Solution II: CUDA

### Prepare the environment
```shell
module load gcc/8.3.0 
module load openblas/0.3.13 
module load cmake/3.20.3
module load openmpi/4.1.2 
module load cuda/11.1
```

### Checkout the source code
We need to clone the "cuda" branch for this solution:
```shell
git clone --branch cuda https://github.com/Mi-Dora/Parallel-Framework-for-Particle-Dynamics-Simulation.git
```

### To Compile
Under the "Parallel-Framework-for-Particle-Dynamics-Simulation" home directory, there is the CMakeList.txt file, and a "build.config" file that indicates configuration parameters and reffered to by the "build.sh" script. To compile the code, run:
```shell
./build.sh
```
Note: If we encounter the "Permission Denied" error, try this line before running the command above:
```shell
chmod +x build.sh
```

### Download the data files
This part is identical to the MPI+OpenMP solution as we use the same script to generate input data. One can refer to the previous section.

### Processing Input for the run
Similar to the MPI+OpenMP solution, since we are using the same "generate.py" script, there is not much to process. We just read in and parse each line of data.

### To Validate
We have a sequential version algorithm, and a test script in "/test/test_seq.cpp" that runs this algorithm. To execute it, run:
```shell
./build/test/TEST_seq 4096 1 1
```
(more of this will be discussed in the following section)
We also have a test script that checks whether the sequential version and CUDA versions outputs the same computations, as in "/test/test_correctness.cpp". To execute it, run:
```shell
./build/test/TEST_correctness 4096 0
```
Where the first argument（4096）indicates number of particles, and the second argument (0) indicates number of iterations. Since iteration index starts with 0, an input of 0 result in one single iteration. We believe it is sufficient to test correctness by only looking at the first iteration, but you can run additional iterations to double check.

IMPORTANT: These test scripts will look for input data file themselves. You might need to change the directory to where you actually store the data input files in the "test_correctness.cpp" before compilation. For example, right now the path is named as 
```c++
const std::string inputFilename_cuda = "/storage/home/hcocice1/ydu340/Particle-new/data/vis/cuda_" + n_par + "_" + n_iter + ".txt";
const std::string inputFilename_seq = "/storage/home/hcocice1/ydu340/Particle-new/data/vis/seq_" + n_par + "_" + n_iter + ".txt";
```
Where we assume that the filename is related to number of particles.

### To Run
The "run.pbs" script provides a full set of commands to run different test scripts. As long as we locate input data files in the correct directory as described in the above section, run:
```shell
qsub run.pbs
```
and wait for the output! The job will be submmitted to the coc-ice-gpu cluster.

### Version of the code (-x contain versions: )
We don't have a version argument. Different algorithms are stored with different branches of the repository. This CUDA version is stored in the "cuda" branch.

### Debug option (-y have different debug option)
To change to debug option, alter the first line of the "build.config" file in the home directory as follows:
```shell
-DCMAKE_BUILD_TYPE=DEBUG
```

### Benchmarking
The full CUDA solution is in one single "cuda" branch.
