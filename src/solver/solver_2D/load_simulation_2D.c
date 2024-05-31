#include "global_parameters.h"
#include "io_operations/foreground_io/foreground_io_2D/foreground_io_2D.h"
#include "data_structures/background_data/background_data.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_data_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_2D.h"
#include "initialization/initialize_background/solar_s_initialization/solar_s_initialization.h"
#include "MPI_module/mpi_info_struct.h"

void load_simulation_2D(struct BackgroundVariables **bg, struct ForegroundVariables2D **fg, struct ForegroundVariables2D **fg_previous, struct GridInfo2D **grid_info, struct PrecalculatedVariables2D **precalc, struct MpiInfo *mpi_info, int *save_nr, FLOAT_P *t, FLOAT_P *first_t)
{
    // Allocating grid info struct and initializing it
    allocate_grid_info_struct_2D(grid_info);
    initialize_grid_info_2D(*grid_info, mpi_info);

    // Allocating background, foreground and precalculated data structs
    allocate_background_struct(bg, (*grid_info)->nz_full);
    allocate_foreground_struct_2D(fg, *grid_info);
    allocate_foreground_struct_2D(fg_previous, *grid_info);
    allocate_precalculate_data_struct_2D(precalc, (*grid_info)->nz_full, (*grid_info)->ny);

    // Initializing background
    solar_s_initialization(*bg, mpi_info, (*grid_info)->nz_full, (*grid_info)->nz_ghost, (*grid_info)->dz, (*grid_info)->z0);
    initialize_precalculated_data_2D(*precalc, *bg, *grid_info, mpi_info);

    // Loading foreground
    char file_path_foreground[150];
    snprintf(file_path_foreground, sizeof(file_path_foreground), "%s%s/snap%d_%d.h5", SAVE_DIR, RUN_NAME, LOAD_SNAP_NUMBER, mpi_info->rank);
    *t = load_foreground_2D(*fg_previous, *grid_info, file_path_foreground);
    
    
    *first_t = *t;

    *save_nr = LOAD_SNAP_NUMBER + 1;
}