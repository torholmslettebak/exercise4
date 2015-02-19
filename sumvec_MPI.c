#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include <mpi.h>
#include <string.h>


// Is run with: mpirun -np 4 ./sumvec_MPI
// if you write 2 instead of 4 -> 2 processors
#define M_PI 3.14159265358979323846



void fillVectorWithSeries(Vector vec)
{
	for (int i = 0; i < vec->len; i++)
	{
		vec -> data[i] = (1.0 / ((i + 1.0)*(i + 1.0)));
	}
}

double sumVector(Vector vec, int start, int stop)
{
	double sum = 0;
	for (int i = start; i < stop; i++)
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

// used for testing of vectors
void printVector2(int length, Vector vec)
{
	printf("Here comes the vector: \n");
	for (int i = 0; i < length; i++)
	{
		printf(" %lf ", vec -> data[i]);
	}
	printf("\n");
}


int isPowerOfTwo (unsigned int x) 
{
	return ((x != 0) && !(x & (x - 1))); 
}

int main(int argc, char **argv)
{
	int length, myid, nproc, tag1, tag2, tag3;
	Vector vec;
	MPI_Status status;
	
	//test();
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
		int *len, *displ;
		t1 = WallTime();
		length = (int) pow(2, k);

		tag1 = 1;
		tag2 = 2;
		tag3 = 3;
		int elements_per_proc = (int) length/nproc;
		if (myid == 0) // If I am root, generate vector
		{
			/* Initialize vector */
			vec = createVector(length);
			double global_sum = 0;
			fillVectorWithSeries(vec);
			splitVector(length, nproc, &len, &displ);
			/* Distribute vector to the other processes */
			for (int i = 1; i < nproc; i++)
			{
				MPI_Send(&displ[i], 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
				MPI_Send((vec -> data), length, MPI_DOUBLE, i, tag3, MPI_COMM_WORLD);
			}

			/* Receives the local sums from each process*/
			double receivedSum = 0;
			for (int i = 1; i < nproc; i++)
			{
				MPI_Recv(&receivedSum, 1, MPI_DOUBLE, i, tag2, MPI_COMM_WORLD, &status);
				global_sum += receivedSum;
			}

			/*sum own elements*/
			double localSum = sumVector(vec, 0, elements_per_proc);
			global_sum = global_sum + localSum;
			double error = (double) ((M_PI * M_PI / 6) - global_sum);
			t2 = WallTime();
			dt = t2 - t1;
			printf("Number of elements: %d\t\t", length);
			printf("time:   %0.10f\t", dt);
			printf("sum:	%0.10f\t", global_sum);
			printf("error:	%0.10f\n", error);
		}

		else
		{
			vec = createVector(length);
			int localIndex;
			/* Each process receives local index, where it is to begin summing the vector*/
			MPI_Recv(&localIndex, elements_per_proc, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);
			/* Each process receives a copy of the vector which is stored in vec*/			
			MPI_Recv((vec -> data), length, MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD, &status);
			/* Computes the sum for each subvector*/
			double localSum = sumVector(vec, localIndex, localIndex + elements_per_proc);
			/* Returns the computed local sum to master process*/
			MPI_Send(&localSum, 1, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD);
		}
		freeVector(vec);
	}

	MPI_Finalize();
	return 0;
}