#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"

void initialize_foreground_zeros_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info)
{
    /*
    Sets all foreground variables to zero.

    Parameters
    ----------
    fg : ForegroundVariables1D
        A pointer to the ForegroundVariables1D struct.
    grid_info : GridInfo1D
        A pointer to the GridInfo1D struct.
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        fg->p1[i] = 0.0;
        fg->rho1[i] = 0.0;
        fg->T1[i] = 0.0;
        fg->s1[i] = 0.0;
        fg->vz[i] = 0.0;
    }

    // No need to extrapolate since all values are zero
}