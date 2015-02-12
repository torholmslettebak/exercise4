#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

#define M_PI 3.14159265358979323846

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

Vector generateVector(int length)
{
	Vector sumVector = createVector(length);
	for (int i = 0; i < length; i++)
	{
		sumVector -> data[i] = 1.0/((i+1)*(i+1));
	}
	return sumVector;
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

double sumVector(Vector vec)
{
	double sum = 0;
	#pragma omp parallel for schedule(static) reduction(+:sum)
	for (unsigned int i = 0; i < vec -> len; i++)
	{
		// printf("Hai, i'm thread number %d and I am working on i %d\n", omp_get_thread_num(), i);
		sum += vec -> data[i];
	}
	return sum;
}

int main(void)
{
	double sum, t1, t2, dt;
	double actualSum = M_PI*M_PI/6;
	t1 = WallTime();
	for (int k = 3; k < 14; k++)
	{
		int length = pow(2.0, 1.0*k);
		Vector vec = generateVector(length);
		sum = sumVector(vec);
		printf("The computed sum for sum: %.13lf\n", sum);
		printf("The sum as number of elements goes -> inf: %.13lf\n", actualSum);
		printf("The difference S - Sn for n = 2^k for k = %d,  error = %.14lf\n", k, actualSum - sum);		
	}
	t2 = WallTime();
	//printVector2(length, vec);
	printf("Time elapsed summing vectors: %lf \n", t2-t1);
	return 0;
}
//1.6448340718481
//1.6449340668482

