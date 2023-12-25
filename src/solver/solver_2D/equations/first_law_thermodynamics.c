#include "equations.h"

void first_law_thermodynamics(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info)
{
    /*
    Calculates the foreground temperature from the first law of thermodynamics.

    Parameters
    ----------
    fg : ForegroundVariables2D
        A pointer to the ForegroundVariables2D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo2D
        A pointer to the GridInfo2D struct.
    */

    // Calculating spesific gas constant
    FLOAT_P r_star = K_B / (MU * M_U);
    // Spesific heat under constant pressure
    FLOAT_P c_p = r_star /(1.0-1.0/GAMMA);

    // Getting grid info
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->T1[i][j] = bg->T0[i]*(fg->s1[i][j]/c_p + (GAMMA-1)/GAMMA * fg->p1[i][j]/bg->p0[i]);
        }
    }
}