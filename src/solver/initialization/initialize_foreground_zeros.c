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

    #if DIMENSIONS == 1
        // Getting grid info
        int nz_full = grid_info->nz_full;

        for (int i = 0; i < nz_full; i++)
        {
            fg->p1[i] = 0.0;
            fg->rho1[i] = 0.0;
            fg->T1[i] = 0.0;
            fg->s1[i] = 0.0;
            fg->vz[i] = 0.0;
            #if BFIELD_ON == 1
                fg->Bx[i] = 0.0;
                fg->By[i] = 0.0;
                fg->Bz[i] = 0.0;
            #endif // BFIELD_ON
        }
    #elif DIMENSIONS == 2
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
                #if BFIELD_ON == 1
                    fg->Bx[i][j] = 0.0;
                    fg->By[i][j] = 0.0;
                    fg->Bz[i][j] = 0.0;
                #endif // BFIELD_ON
            }
        }
    #elif DIMENSIONS == 3
        // Getting grid info
        int nx = grid_info->nx;
        int ny = grid_info->ny;
        int nz_full = grid_info->nz_full;

        for (int i = 0; i < nz_full; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    fg->p1[i][j][k] = 0.0;
                    fg->rho1[i][j][k] = 0.0;
                    fg->T1[i][j][k] = 0.0;
                    fg->s1[i][j][k] = 0.0;
                    fg->vx[i][j][k] = 0.0;
                    fg->vy[i][j][k] = 0.0;
                    fg->vz[i][j][k] = 0.0;
                    #if BFIELD_ON == 1
                        fg->Bx[i][j][k] = 0.0;
                        fg->By[i][j][k] = 0.0;
                        fg->Bz[i][j][k] = 0.0;
                    #endif // BFIELD_ON
                }
            }
        }
    #endif // DIMENSIONS

    // No need to extrapolate since all values are zero
}