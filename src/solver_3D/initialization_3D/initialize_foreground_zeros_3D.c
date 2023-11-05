#include "initialization_3D.h"

void initialize_foreground_zeros_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info) 
{
    /*
    Initializes the foreground variables to zero.

    Parameters
    ----------
    fg : ForegroundVariables3D
        A pointer to the ForegroundVariables3D struct.
    grid_info : GridInfo3D
        A pointer to the GridInfo3D struct.
    */
    
    // Getting grid info
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;
    int nx = grid_info->nx;

    // Setting everything to zero
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->s1[i][j][k] = 0.0;
                fg->p1[i][j][k] = 0.0;
                fg->rho1[i][j][k] = 0.0;
                fg->T1[i][j][k] = 0.0;
                fg->vx[i][j][k] = 0.0;
                fg->vy[i][j][k] = 0.0;
                fg->vz[i][j][k] = 0.0;
            }
        }
    }
}
