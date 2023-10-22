#include "initialization.h"

void initialize_foreground_zeros(struct ForegroundVariables *fg, struct GridInfo *grid_info)
{
    /*
    Sets all foreground variables to zero.

    Parameters
    ----------
    fg : ForegroundVariables
        A pointer to the ForegroundVariables struct.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
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