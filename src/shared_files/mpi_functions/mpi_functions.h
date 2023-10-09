#ifndef MPI_FUNCTIONS_H__
#define MPI_FUNCTIONS_H__

#include <mpi.h>
#include "shared_files.h"

void communicate_above_below(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info, struct MpiInfo *mpi_info);

#endif // MPI_FUNCTIONS_H__