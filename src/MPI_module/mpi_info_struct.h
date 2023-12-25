#ifndef MPI_INFO_STRUCT_H__
#define MPI_INFO_STRUCT_H__

#include <stdbool.h>

struct MpiInfo
{
    int rank, size;
    bool has_neighbor_below, has_neighbor_above;
};

#endif // MPI_INFO_STRUCT_H__