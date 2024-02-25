#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "global_parameters.h"
#include "global_constants.h"

void first_law_of_thermodynamics_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info)
{
    /*
    Calculates the foreground temperature from the first law of thermodynamics.

    Parameters
    ----------
    fg : ForegroundVariables1D
        A pointer to the ForegroundVariables1D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo1D
        A pointer to the GridInfo1D struct.
    */

    // Spesific heat under constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);

    // Getting grid info
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        fg->T1[i] = bg->T0[i] * (fg->s1[i]/c_p + (GAMMA-1.0)/GAMMA * fg->p1[i]/bg->p0[i]);
    }
}