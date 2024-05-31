#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_2D.h"
#include "global_float_precision.h"
#include "global_parameters.h"
#include "global_constants.h"

void calculate_vy_vertical_average_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *pv)
{
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    FLOAT_P dz = grid_info->dz;
    FLOAT_P Lz = R_SUN * (R_END - R_START);

    for (int j = 0; j < ny; j++)
    {
        pv->vy_vertical_average[j] = 0.0;
        for (int i = nz_ghost; i < nz_full-nz_ghost-1; i++)
        {
            pv->vy_vertical_average[j] += fg->vy[i][j];
        }
        pv->vy_vertical_average[j] *= dz / Lz;
    }
}