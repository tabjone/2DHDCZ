#include "equations.h"

void equation_of_state(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    /*
    Calculates the foreground density from the equation of state.

    Parameters
    ----------
    fg : ForegroundVariables
        A pointer to the ForegroundVariables struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    #if DIMENSIONS == 1
        // Getting grid info
        int nz_full = grid_info->nz_full;

        for (int i = 0; i < nz_full; i++)
        {
            fg->rho1[i] = (fg->p1[i]/bg->p0[i] - fg->T1[i]/bg->T0[i]) * bg->rho0[i];
        }

    #elif DIMENSIONS == 2
        // Getting grid info
        int ny = grid_info->ny;
        int nz_full = grid_info->nz_full;

        for (int i = 0; i < nz_full; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                fg->rho1[i][j] = (fg->p1[i][j]/bg->p0[i] - fg->T1[i][j]/bg->T0[i]) * bg->rho0[i];
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
                    fg->rho1[i][j][k] = (fg->p1[i][j][k]/bg->p0[i] - fg->T1[i][j][k]/bg->T0[i]) * bg->rho0[i];
                }
            }
        }
    #endif // DIMENSIONS
}