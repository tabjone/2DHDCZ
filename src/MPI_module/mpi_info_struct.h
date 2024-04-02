#ifndef MPI_INFO_STRUCT_H__
#define MPI_INFO_STRUCT_H__

#include <stdbool.h>

struct MpiInfo
{
    int rank, size;
    bool has_neighbor_below, has_neighbor_above;
    int soft_wall_end_process, soft_wall_end_index;
};

#endif // MPI_INFO_STRUCT_H__