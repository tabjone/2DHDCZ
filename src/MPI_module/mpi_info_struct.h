#ifndef MPI_INFO_STRUCT_H__
#define MPI_INFO_STRUCT_H__

#include <stdbool.h>

struct MpiInfo
{
    int rank, size;
    bool has_neighbor_below, has_neighbor_above;
};

void initialize_mpi_info_struct(struct MpiInfo **mpi_info_ptr);

#endif // MPI_INFO_STRUCT_H__