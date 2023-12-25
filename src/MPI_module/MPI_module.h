#ifndef MPI_MODULE_H__
#define MPI_MODULE_H__

#include "mpi_info_struct.h"
#include "./communication/communication.h"

void initialize_mpi_info_struct(struct MpiInfo **mpi_info_ptr);

#endif // MPI_MODULE_H__