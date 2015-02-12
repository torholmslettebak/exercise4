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
	int length = pow(2, 14);
	Vector vec = generateVector(length);
	t1 = WallTime();
	sum = sumVector(vec);
	t2 = WallTime();
	//printVector2(length, vec);
	printf("The sum: %.13lf calculated in: , with serial computation. %lf\n", sum, t2-t1);
	printf("The sum as number of elements -> inf: %.13lf\n", M_PI*M_PI/6);
	return 0;
}
//1.6448340718481
//1.6449340668482

