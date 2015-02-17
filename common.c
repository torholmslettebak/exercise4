#include "common.h"
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>

#ifdef HAVE_MPI 
    MPI_Comm WorldComm;
    MPI_Comm SelfComm;
#endif

Vector createVector(int len)
{
    Vector result = (Vector)calloc(1, sizeof(vector_t));
    result->data = calloc(len, sizeof(double));
    result->len = result->glob_len = len;
    result->stride = 1;
    #ifdef HAVE_MPI
        result->comm = &SelfComm;
    #endif
    result->comm_size = 1;
    result->comm_rank = 0;
    result->displ = NULL;
    result->sizes = NULL;
    return result;
}

// Vector createVectorMPI(int glob_len, MPI_Comm* comm, int allocdata)
// {
    
// }
void freeVector(Vector vec)
{
    free(vec->data);
    free(vec->sizes);
    free(vec->displ);
    free(vec);
}

int getMaxThreads()
{
    #ifdef HAVE_OPENMP
        return omp_get_max_threads();
    #else
        return 1;
    #endif
}

int getCurrentThread()
{
    #ifdef HAVE_OPENMP
        return omp_get_thread_num();
    #else
        return 0;
    #endif
}

double WallTime ()
{
    #ifdef HAVE_MPI
        return MPI_Wtime();
    #elif defined(HAVE_OPENMP)
        return omp_get_wtime();
    #else
        struct timeval tmpTime;
        gettimeofday(&tmpTime,NULL);
        return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
    #endif
}

/*void copyVector(Vector y, const Vector x)
{
  dcopy(&y->len, x->data, &x->stride, y->data, &y->stride);
}
*/

void fillVector(Vector x, double alpha)
{
    int i;
    for (i=0;i<x->len;++i)
    x->data[i*x->stride] = alpha;
}

void splitVector(int globLen, int size, int** len, int** displ)
{
  int i;
  *len = calloc(size,sizeof(int));
  *displ = calloc(size,sizeof(int));
  for (i=0;i<size;++i) {
    (*len)[i] = globLen/size;
    if (globLen % size && i >= (size - globLen % size))
      (*len)[i]++;
    if (i < size-1)
      (*displ)[i+1] = (*displ)[i]+(*len)[i];
  }
}