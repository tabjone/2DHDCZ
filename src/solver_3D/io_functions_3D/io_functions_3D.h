#ifndef IO_FUNCTIONS_H__
#define IO_FUNCTIONS_H__

#include "global_parameters.h"
#include "shared_files.h"

void save_foreground_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, int snap_number, FLOAT_P time);

#endif // IO_FUNCTIONS_H__