#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "global_parameters.h"
#include "global_constants.h"

void first_law_of_thermodynamics_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info)
{
    // Calculating spesific gas constant
    FLOAT_P r_star = K_B / (MU * M_U);
    // Spesific heat under constant pressure
    FLOAT_P c_p = r_star /(1.0-1.0/GAMMA);

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
                fg->T1[i][j][k] = bg->T0[i]*(fg->s1[i][j][k]/c_p + (GAMMA-1)/GAMMA * fg->p1[i][j][k]/bg->p0[i]);
            }
        }
    }
}