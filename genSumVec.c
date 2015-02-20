#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

#define M_PI 3.14159265358979323846

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
	for (int i = 0; i < vec -> len; i++)
	{
		sum += vec -> data[i];
	}
	return sum;
}
int main(void)
{
	double sum, t1,t2;
	for (int i = 3; i < 15; i++)
	{
		t1 = WallTime();
		int length = pow(2, i);
		Vector vec = generateVector(length);
		sum = sumVector(vec);
		t2 = WallTime();
		double dt = t2 - t1;
		double error = (double) ((M_PI * M_PI / 6) - sum);
		printf("Number of elements: \t%d\t", length);
		printf("sum:	%0.10f\t", sum);
		printf("time:   %0.10f\t", dt);
		printf("error:	%0.10f\n", error);
	}

	return 0;
}

