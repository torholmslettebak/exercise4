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
  return result;
}
void freeVector(Vector vec)
{
  free(vec->data);
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