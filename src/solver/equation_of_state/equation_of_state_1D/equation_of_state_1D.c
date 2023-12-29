#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"

void equation_of_state_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info)
{
    // Getting grid info
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        fg->rho1[i] = (fg->p1[i]/bg->p0[i] - fg->T1[i]/bg->T0[i]) * bg->rho0[i];
    }
}