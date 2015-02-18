#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include <mpi.h>
#include <string.h>


// Is run with: mpirun -np 4 ./sumvec_MPI
// if you write 2 instead of 4 -> 2 processors
#define M_PI 3.14159265358979323846

Vector generateVector(int length)
{
	Vector vec = createVector(length);
	for (int i = 0; i < length; i++)
	{
		vec -> data[i] = (1.0 / ((i + 1.0)*(i + 1.0)));
	}
	return vec;
}

double sumVector(Vector vec, int length)
{
	double sum = 0;
	#pragma omp parallel for schedule(static) reduction(+:sum)
	for (int i = 0; i < length; i++)
	{
		sum += vec -> data[i];	
	}
	return sum;
}

double difference(double sum)
{
	double actualSum = M_PI * M_PI / 6;
	return actualSum - sum;
}

// used for testing
void printVector2(int length, Vector vec)
{
	printf("Here comes the vector: \n");
	for (int i = 0; i < length; i++)
	{
		printf(" %lf ", vec -> data[i]);
	}
	printf("\n");
}


int isPowerOfTwo(int x)
{
	if ((x%2 == 0) && x>1)
	{
		return 1;
	}
	else 
		return 0;
}

int main(int argc, char **argv)
{
	int length, myid, nproc;
	MPI_Status status;
	double global_sum;
	// From this point on every process executes a seperate copy of the program
	MPI_Init (&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (isPowerOfTwo(nproc) != 1)
	{
		printf("The number of processes needs to be a power of two.\n");
		#ifdef HAVE_MPI
			MPI_Finalize();
		#endif
		return 1;
	}

	for (int k = 3; k < 15; k++)
	{
		double t1, t2, dt;
		t1 = WallTime();
		length = (int) pow(2, k);
		Vector vec;
		int elements_per_proc = (int) length/nproc;
		Vector result = createVector(nproc);
		Vector recVec = createVector(elements_per_proc);
		if (myid == 0) // If I am root, generate vector
		{
			vec = generateVector(length);
		}
		MPI_Scatter(vec -> data, elements_per_proc, MPI_DOUBLE, recVec -> data, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		double local_sum = sumVector(recVec, elements_per_proc);

		MPI_Gather(&local_sum, 1, MPI_DOUBLE, result -> data, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if(myid == 0)
		{
			global_sum = sumVector(result, nproc);
			printf("The global_sum: %lf\n", global_sum);
		}
		if(myid == 0)
		{
			double error = (double) ((M_PI * M_PI / 6) - global_sum);
			printf("The error for 2^k (k = %d)  	  = %d, elements is error = %lf\n", k,length, error);
			t2 = WallTime();
			dt = t2 - t1;
			freeVector(vec);
			printf("The time elapsed for 2^k (k = %d) = %d  elements: dt = %lf\n", k, length, dt);
		}
		freeVector(result);
		freeVector(recVec);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}