cd _build
cmake .. && make
# ./gensumserial 
mpirun -np 4 ./sumvec_MPI 
mpirun -np 4./sumvec_comb
export OMP_NUM_THREADS=4
./sumvec_openmp