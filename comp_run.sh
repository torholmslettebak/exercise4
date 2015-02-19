cd _build;
cmake .. && make && mpirun -np 4 ./sumvec_MPI 
