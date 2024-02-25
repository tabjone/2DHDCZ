#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "global_parameters.h"
#include "global_constants.h"

void first_law_of_thermodynamics_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info)
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

    // Spesific heat under constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);

    // Getting grid info
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->T1[i][j] = bg->T0[i]*(fg->s1[i][j]/c_p + (GAMMA-1.0)/GAMMA * fg->p1[i][j]/bg->p0[i]);
        }
    }
}