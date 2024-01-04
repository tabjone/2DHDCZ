#include "data_structures/background_data/background_data.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_data_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_1D.h"
#include "initialization/initialize_background/solar_s_initialization/solar_s_initialization.h"
#include "io_operations/io_operations.h"
#include "initialization/initialize_foreground/initialize_foreground_1D/initialize_foreground_1D.h"

#include <stdio.h>

void initialize_simulation_1D(struct BackgroundVariables **bg, struct ForegroundVariables1D **fg, struct ForegroundVariables1D **fg_previous, struct GridInfo1D **grid_info, struct PrecalculatedVariables1D **precalc, struct MpiInfo *mpi_info, int *save_nr, FLOAT_P *t, FLOAT_P *first_t)
{
    *save_nr = 0;
    // Allocating grid info struct and initializing it
    allocate_grid_info_struct_1D(grid_info);
    initialize_grid_info_1D(*grid_info, mpi_info);

    // Allocating background, foreground and precalculated data structs
    allocate_background_struct(bg, (*grid_info)->nz_full);
    allocate_foreground_struct_1D(fg, *grid_info);
    allocate_foreground_struct_1D(fg_previous, *grid_info);
    allocate_precalculate_data_struct_1D(precalc, (*grid_info)->nz_full);

    // Initializing background, foreground and precalculated data structs
    solar_s_initialization(*bg, mpi_info, (*grid_info)->nz_full, (*grid_info)->nz_ghost, (*grid_info)->dz, (*grid_info)->z0);
    initialize_foreground_1D(*fg_previous, *bg, *grid_info, mpi_info);
    initialize_precalculated_data_1D(*precalc, *bg, *grid_info);
    // Saving parameters to file
    save_background(*bg, mpi_info, (*grid_info)->nz_full, (*grid_info)->nz, (*grid_info)->nz_ghost, (*grid_info)->dz, (*grid_info)->z0, (*grid_info)->z1);
    save_mpi_info(mpi_info);
    save_simulation_info(mpi_info);
    save_foreground_1D(*fg_previous, *grid_info, mpi_info, 0, 0.0);
    
    // Initializing time variables
    *save_nr += 1;
    *t = 0.0;
    *first_t = 0.0;

    if (*bg == NULL) {
        printf("Error: Memory allocation for bg failed\n");
    }

    if (*fg == NULL) {
        printf("Error: Memory allocation for fg failed\n");
    }

    if (*fg_previous == NULL) {
        printf("Error: Memory allocation for fg_previous failed\n");
    }

    if (*grid_info == NULL) {
        printf("Error: Memory allocation for grid_info failed\n");
    }

    if (*precalc == NULL) {
        printf("Error: Memory allocation for precalc failed\n");
    }
}
