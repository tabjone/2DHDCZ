#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi_info_struct.h"

void initialize_mpi_info_struct(struct MpiInfo **mpi_info_ptr)
{
    *mpi_info_ptr = malloc(sizeof(struct MpiInfo));
    if (*mpi_info_ptr == NULL) {
        // Handle memory allocation failure
        fprintf(stderr, "Memory allocation failed for mpi_info\n");
        return;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &(*mpi_info_ptr)->size);
    MPI_Comm_rank(MPI_COMM_WORLD, &(*mpi_info_ptr)->rank);

    (*mpi_info_ptr)->has_neighbor_below = true;
    (*mpi_info_ptr)->has_neighbor_above = true;

    if ((*mpi_info_ptr)->rank == 0)
    {
        (*mpi_info_ptr)->has_neighbor_below = false;
    }
    if ((*mpi_info_ptr)->rank == (*mpi_info_ptr)->size - 1)
    {
        (*mpi_info_ptr)->has_neighbor_above = false;
    }
}
