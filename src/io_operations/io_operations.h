#ifndef IO_OPERATIONS_H__
#define IO_OPERATIONS_H__

#include "./background_io/background_io.h"
#include "./general_io/general_io.h"
#include "MPI_module/mpi_info_struct.h"

void save_mpi_info(struct MpiInfo *mpi_info);
void save_simulation_info(struct MpiInfo *mpi_info);

#endif // IO_OPERATIONS_H__