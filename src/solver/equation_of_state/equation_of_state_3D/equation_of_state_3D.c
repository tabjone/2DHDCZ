#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"

void equation_of_state_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info)
{
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
}