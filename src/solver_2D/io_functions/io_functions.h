#ifndef IO_FUNCTIONS_H__
#define IO_FUNCTIONS_H__

#include "hdf5.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "shared_files.h"
#include "global_parameters.h"
#include "global_constants.h"
#include "../precalculated_variables/precalculated_variables.h"

void save_foreground(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, int snap_number, FLOAT_P time);

void load_grid_info(struct GridInfo2D **grid_info, const char *file_path);
void load_background(struct BackgroundVariables *bg, struct GridInfo2D *grid_info, const char *file_path);
FLOAT_P load_foreground(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, const char *file_path);

void read_dataset_line(hid_t dataset, const char* name, hid_t datatype, FLOAT_P *data);

void save_rhs(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc, int save_nr);
void save_info(struct MpiInfo *mpi_info);

#endif // IO_FUNCTIONS_H__