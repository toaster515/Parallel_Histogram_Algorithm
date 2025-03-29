#ifndef BASIC_DEFS_H
#define BASIC_DEFS_H

#include <cstdint>

#include <mpi.h>

typedef uint64_t dist_sort_t;
typedef uint64_t dist_sort_size_t;

#define DIST_SORT_MAX UINT64_MAX
#define DIST_SORT_SIZE_MAX UINT64_MAX

#define MPI_TYPE_DIST_SORT_T MPI_UINT64_T
#define MPI_TYPE_DIST_SORT_SIZE_T MPI_UINT64_T

#endif /*BASIC_DEFS_H*/
