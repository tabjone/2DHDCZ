#include "initialization_2D.h"

void initialize_foreground_zeros(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info)
{
    /*
    Sets all foreground variables to zero.

    Parameters
    ----------
    fg : ForegroundVariables2D
        A pointer to the ForegroundVariables2D struct.
    grid_info : GridInfo2D
        A pointer to the GridInfo2D struct.
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = 0.0;
            fg->rho1[i][j] = 0.0;
            fg->T1[i][j] = 0.0;
            fg->s1[i][j] = 0.0;
            fg->vy[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
        }
    }

    // No need to extrapolate since all values are zero
}