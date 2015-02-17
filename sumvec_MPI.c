#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include <mpi.h>
#include <string.h>
	
#define M_PI 3.14159265358979323846

/*int main(int argc, char **argv)
{
	int rank, size, tag, i;
	MPI_Status status;
	char message[20];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("size: %d\n", size);
	printf("rank: %d\n", rank);
	tag = 100;
	if (rank == 0) // I am thread number 0 aka rank 0
	{
		strcpy(message, "Hai, world!"); 
		for (i = 1; i < size; i++) // thread 0 send  the message to all the other threads
		{
			printf("I am %d, and I am sending a message to i: %d\n", rank, i);
			MPI_Send(message, 13, MPI_CHAR, i, tag, MPI_COMM_WORLD);
			// thread i will now perform this program for itself
			// but this thread will not have rank 0, and therefor it will performonly receive message
			// and print itself and message
		}
	}
	else
	{
		MPI_Recv(message, 13, MPI_CHAR, 0, tag, MPI_COMM_WORLD, &status);
	}
	printf("process: %d: %s\n", rank, message);
	MPI_Finalize();
	return 0;
}
*/
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
	for (int i = 0; i < length; i++)
	{
		sum += vec -> data[i];
		// printf("I am process: %d\n", myid);
		// printf("Out of: %d\n", nproc);
		
	}
	// freeVector(vec);
	return sum;
}

double difference(double sum)
{
	double actualSum = M_PI * M_PI / 6;
	return actualSum - sum;
}

void printVector2(int length, Vector vec)
{
	printf("Here comes the vector: \n");
	for (int i = 0; i < length; i++)
	{
		printf(" %lf ", vec -> data[i]);
	}
	printf("\n");
}
void printVector3(int length, double* arr)
{
	printf("Here comes the vector: \n");
	for (int i = 0; i < length; i++)
	{
		printf(" %lf ", arr[i]);
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
	int rank, size, tag, i, ierr;
	int length, myid, nproc;
	MPI_Status status;
	double sumfinal;
	// From this point on every process executes a seperate copy of the program
	MPI_Init (&argc, &argv);
	int root = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	// printf("I am process: %d\n", myid);
	// printf("Out of: %d\n", nproc);
	if (isPowerOfTwo(nproc) != 1)
	{
		printf("The number of processes must be a power of two.\n");
		#ifdef HAVE_MPI
			MPI_Finalize();
		#endif
		return 1;
	}

	for (int k = 3; k < 6; k++)
	{
		double t1, t2, dt;
		t1 = WallTime();
		length = (int) pow(2, k);
		Vector vec;
		int elements_per_proc = length/nproc;
		Vector result = createVector(nproc);
		Vector recVec = createVector(elements_per_proc);
		if (myid == 0) // If I am root, generate vector
		{
			result = createVector(length);
			vec = generateVector(length);
			printVector2(length, vec);
		}
		// MPI_Barrier(MPI_COMM_WORLD);
	//	MPI_Scatter(vec -> data, elements_per_proc, MPI_DOUBLE, sub_rand_nums, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter(vec -> data, elements_per_proc, MPI_DOUBLE, recVec -> data, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		// printf("I am process: %d\n", myid);
		// printf("Out of: %d\n", nproc);
		double local_sum = sumVector(recVec, elements_per_proc);
		// printf("The local sum: %lf\n", local_sum);
		// printf("The recVec\n");
	//	printVector2(elements_per_proc, recVec);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(&local_sum, 1, MPI_DOUBLE, result -> data, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (myid == 0)
		{
			printf("The resultvector: %d\n", nproc);
			printVector2(nproc, result);
			//double global_sum;
			double global_sum = sumVector(result, nproc);
			printf("The global_sum: %lf\n", global_sum);
		}
		// if(myid == 0)
		// {
		// 	double global_sum = sumVector(result, length);
		// 	printf("the sum_ %lf\n", global_sum);
		// }
		// a barrier for master thread to wait for result
		// double error = difference(sum);
		// printf("The error for 2^k (k = %d) = %d,  elements is error = %lf\n", k,length, error);
		// t2 = WallTime();
		// dt = t2 - t1;
		// printf("The time elapsed for 2^k (k = %d) = %d elements: dt = %lf\n", k, length, dt);

	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}