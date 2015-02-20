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
		sum += vec -> data[i];
	}
	return sum;
}

int main(void)
{
	double sum, t1, t2, dt;
	double actualSum = M_PI*M_PI/6;

	for (int k = 3; k < 15; k++)
	{
		t1 = WallTime();
		int length = pow(2.0, 1.0*k);
		Vector vec = generateVector(length);
		sum = sumVector(vec);
		double error = (double) ((M_PI * M_PI / 6) - sum);
		t2 = WallTime();
		dt = t2 - t1;
		printf("Number of elements: \t%d\t", length);
		printf("sum:	%0.10f\t", sum);
		printf("time:   %0.10f\t", dt);
		printf("error:	%0.10f\n", error);		
	}

	return 0;
}
//1.6448340718481
//1.6449340668482

