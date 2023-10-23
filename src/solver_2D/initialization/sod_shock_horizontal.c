#include "initialization.h"

void sod_shock_horizontal(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info)
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
    int ny = grid_info->ny;

    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->vy[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
            if (j < ny/2)
            {
                fg->p1[i][j] = bg->p0[nz_full/2]/10.0;
                fg->rho1[i][j] = bg->rho0[nz_full/2]/10.0;
            }
            else
            {
                fg->p1[i][j] = bg->p0[nz_full/2]/100.0;
                fg->rho1[i][j] = fg->rho1[i][j] = bg->rho0[nz_full/2]/100.0;;
                fg->T1[i][j] = 0.0;
                fg->s1[i][j] = 0.0;
            }
            // Getting T1 from equation of state
            fg->T1[i][j] = bg->T0[i] * (fg->p1[i][j]/bg->p0[i] - fg->rho1[i][j]/bg->rho0[i]);
            // Getting entropy from first law of thermodynamics
            fg->s1[i][j] = c_p * (fg->T1[i][j]/bg->T0[i] - fg->p1[i][j]/bg->p0[i]);
        }
    }
    periodic_boundary_2D(fg->p1, grid_info);
    periodic_boundary_2D(fg->rho1, grid_info);
    periodic_boundary_2D(fg->T1, grid_info);
    periodic_boundary_2D(fg->s1, grid_info);
    periodic_boundary_2D(fg->vy, grid_info);
    periodic_boundary_2D(fg->vz, grid_info);
}