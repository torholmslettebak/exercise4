# exercise4
Exercise from TMA4280, see p4.pdf for info

To build: create build folder in dir, enter build and enter: cmake .. 
To compile, enter make.

This compiles several exe-files.
To run program with MPI: mpirun -np 4 ./sumvec_MPI , or sumvec_comb, -np 4 runs program with 4 processes.
Number of processes must be a power of 2.
To set number of threads when running OpenMP file:  export OMP_NUM_THREADS=2 
