#include "boundary_3D.h"

void periodic_boundary_3D(FLOAT_P ***array, struct GridInfo3D *grid_info)
{
    /*
    Applies periodic boundary conditions to a 2D array.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    grid_info : struct
        A pointer to the GridInfo2D struct.
    */

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;

    // Periodic boundary conditions
    for (int i = 0; i < nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                array[i][j][k] = array[nz + i][j][k];
                array[nz_ghost + nz + i][j][k] = array[i+nz_ghost][j][k];
            }
        }
    }
}