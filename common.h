#ifndef COMMON_H
#define COMMON_H

#ifdef HAVE_MPI
#include <mpi.h>
extern MPI_Comm WorldComm;
extern MPI_Comm SelfComm;
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

// Value of pi
#define M_PI            3.14159265358979323846

// A structure describing a vector
typedef struc 
{
	double* data;	// Array with data
	int len;        // The local length of the vector
	int glob_len;   // The global length of the vector
	int stride;     // The distance in memory between vector elements
}
vector_t;

//brief Convenience typedef
typedef vector_t* Vector;

// brief Create a vector
// param len The length of the vector
// return The new vector
Vector createVector(int len);

// brief Free up memory allocated to a vector
// param vec The vector to free
void freeVector(Vector vec);

//! \brief Copy a vector: \f$y = x\f$
//! \param y The y vector
//! \param[in] x The x vector
void copyVector(Vector y, const Vector x);

//! \brief Fill a vector with a constant: \f$u(i) = \alpha\,\forall\,i\f$
//! \param u The u vector
//! \param[in] alpha The fill constant
void fillVector(Vector u, double alpha);

//! \brief Print a vector to the terminal for inspection
//! \param u The vector to print
void printVector(const Vector u);

//! \brief Get the maximum number of available threads
//! \return Number of available threads
int getMaxThreads();

//! \brief Get current thread ID
//! \return Current thread ID
int getCurrentThread();

//! \brief Get current wall-clock time
//! \return The current wall time in seconds
double WallTime();
#endif