#include "data_structures/background_data/background_data.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_data_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_2D.h"
#include "initialization/initialize_background/solar_s_initialization/solar_s_initialization.h"
#include "io_operations/io_operations.h"
#include "initialization/initialize_foreground/initialize_foreground_2D/initialize_foreground_2D.h"
#include "array_utilities/array_memory_management/array_memory_management.h"

#include "global_parameters.h"

#include <math.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include <stdio.h>

void initialize_simulation_2D(struct BackgroundVariables **bg, struct ForegroundVariables2D **fg, struct ForegroundVariables2D **fg_previous, struct GridInfo2D **grid_info, struct PrecalculatedVariables2D **precalc, struct MpiInfo *mpi_info, int *save_nr, FLOAT_P *t, FLOAT_P *first_t)
{
    *save_nr = 0;
    // Allocating grid info struct and initializing it
    allocate_grid_info_struct_2D(grid_info);
    initialize_grid_info_2D(*grid_info, mpi_info);

    // Allocating background, foreground and precalculated data structs
    allocate_background_struct(bg, (*grid_info)->nz_full);
    allocate_foreground_struct_2D(fg, *grid_info);
    allocate_foreground_struct_2D(fg_previous, *grid_info);
    allocate_precalculate_data_struct_2D(precalc, (*grid_info)->nz_full);

    // Initializing background, foreground and precalculated data structs
    solar_s_initialization(*bg, mpi_info, (*grid_info)->nz_full, (*grid_info)->nz_ghost, (*grid_info)->dz, (*grid_info)->z0);
    initialize_precalculated_data_2D(*precalc, *bg, *grid_info, mpi_info);
    initialize_foreground_2D(*fg_previous, *bg, *grid_info, mpi_info);
    // Saving parameters to file
    save_background(*bg, mpi_info, (*grid_info)->nz_full, (*grid_info)->nz, (*grid_info)->nz_ghost, (*grid_info)->dz, (*grid_info)->z0, (*grid_info)->z1);
    save_mpi_info(mpi_info);
    save_simulation_info(mpi_info);
    save_foreground_2D(*fg_previous, *grid_info, mpi_info, 0, 0.0);
    
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


    // Printing the amount of resolution pr. pressure scale height

    // Getting the grid info
    int nz_full = (*grid_info)->nz_full;
    FLOAT_P *p0 = (*bg)->p0;
    FLOAT_P dz = (*grid_info)->dz;
    FLOAT_P dy = (*grid_info)->dy;

    // Allocating the Hp array
    FLOAT_P *Hp;
    FLOAT_P dp0_dr;
    allocate_1D_array(&Hp, nz_full);

    dp0_dr = (p0[nz_full-1] - p0[nz_full-2]) / (dz);
    Hp[nz_full-1] = - p0[nz_full-1] / dp0_dr;

    dp0_dr = (p0[1] - p0[0]) / (dz);
    Hp[0] = - p0[0] / dp0_dr;

    for (int i = 1; i < nz_full-1; i++)
    {
        dp0_dr = (p0[i+1] - p0[i-1]) / (2.0*dz);
        Hp[i] = - p0[i] / dp0_dr;
    }

    FLOAT_P my_min_Hp = Hp[0];
    for (int i = 1; i < nz_full; i++)
    {
        if (Hp[i] < my_min_Hp)
        {
            my_min_Hp = Hp[i];
        }
    }
    FLOAT_P min_Hp;
    MPI_Allreduce(&my_min_Hp, &min_Hp, 1, MPI_FLOAT_P, MPI_MIN, MPI_COMM_WORLD);

    // Finding the biggest resolution if polar coordinates
    #if COORDINATES == 0
        FLOAT_P max_resolution_dr = dz;
        FLOAT_P max_resolution_dy = dy;

        if (mpi_info->rank == 0)
        {
            printf("Resolution vertical: %.3f [# of pressure scale height]\n", min_Hp/dz);
            printf("Resolution horizontal: %.3f [# of pressure scale height]\n", min_Hp/dy);
        }
    #elif COORDINATES == 1
        FLOAT_P max_resolution_dr = dz;
        FLOAT_P dtheta = dy;

        int nz_ghost = (*grid_info)->nz_ghost;
        // Finding biggest arc length
        FLOAT_P *r = (*bg)->r;
        FLOAT_P my_max_arc = 2.0 * M_PI * r[nz_full-nz_ghost-1] * dtheta;

        FLOAT_P max_arc;
        MPI_Allreduce(&my_max_arc, &max_arc, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);

        if (mpi_info->rank == 0)
        {
            printf("Resolution vertical: %.3f [# of pressure scale height]\n", min_Hp/dz);
            printf("Resolution horizontal: %.3f [# of pressure scale height]\n", min_Hp/max_arc);
        }
    #endif // COORDINATES

    deallocate_1D_array(Hp);
}
