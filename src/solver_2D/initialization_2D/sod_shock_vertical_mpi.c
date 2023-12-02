#include "initialization_2D.h"
#include <math.h>

void sod_shock_vertical_mpi(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct with a Sod Shock Tube test.

    Parameters
    ----------
    fg : ForegroundVariables2D
        A pointer to the ForegroundVariables2D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo2D
        A pointer to the GridInfo2D struct.
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    int ny = grid_info->ny;

    FLOAT_P dz = grid_info->dz;
    FLOAT_P z0 = grid_info->z0;

    initialize_foreground_zeros(fg, grid_info);

    FLOAT_P z, tanh_value;
    FLOAT_P width = (R_END - R_START) * R_SUN * 0.05;
    FLOAT_P steepness = 1.0;
    FLOAT_P midpoint = (R_START + (R_END - R_START) / 2.0) * R_SUN;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        z = (i-nz_ghost) * dz + z0;
        
        tanh_value = tanh(steepness * (z - midpoint) / width);
        for (int j = 0; j < ny; j++)
        {
        
            fg->p1[i][j] = 0.5 * ((UPPER_PRESSURE_BOUNDARY - LOWER_PRESSURE_BOUNDARY) * tanh_value + UPPER_PRESSURE_BOUNDARY + LOWER_PRESSURE_BOUNDARY);
            fg->s1[i][j] = 0.5 * ((UPPER_ENTROPY_BOUNDARY - LOWER_ENTROPY_BOUNDARY) * tanh_value + UPPER_ENTROPY_BOUNDARY + LOWER_ENTROPY_BOUNDARY);

        }
    }
    update_vertical_boundary_ghostcells_2D(fg->p1, grid_info, mpi_info);
    update_vertical_boundary_ghostcells_2D(fg->s1, grid_info, mpi_info);

    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);
}