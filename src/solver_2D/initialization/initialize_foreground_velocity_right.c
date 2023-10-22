#include "initialization.h"

void initialize_foreground_velocity_right(struct ForegroundVariables *fg, struct GridInfo *grid_info)
{
    /*
    Sets the velocity to 1000.0 in the y-direction and 0.0 in the x- and z-directions.

    Parameters
    ----------
    fg : ForegroundVariables
        A pointer to the ForegroundVariables struct.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    mpi_info : MpiInfo
        A pointer to the MpiInfo struct.
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    initialize_foreground_zeros(fg, grid_info);
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->vy[i][j] = 1000.0; // cm/s
        }
    }

    // Setting velocity to zero at the boundaries
    // Top boundary
    for (int j = 0; j < ny; j++)
    {
        fg->vy[nz_full-nz_ghost-1][j] = 0.0;
    }
    extrapolate_2D_array_constant_up(fg->vy, grid_info);
    // Bottom boundary
    for (int j = 0; j < ny; j++)
    {
        fg->vy[nz_ghost][j] = 0.0;
    }
    extrapolate_2D_array_constant_down(fg->vy, grid_info);
}