#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_2D.h"
#include "global_float_precision.h"

void calculate_vz_horizontal_average_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *pv)
{
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P Ly = grid_info->y1 - grid_info->y0;

    for (int i = 0; i < nz_full; i++)
    {
        pv->vz_horizontal_average[i] = 0.0;
        for (int j = 0; j < ny; j++)
        {
            pv->vz_horizontal_average[i] += fg->vz[i][j];
        }
        pv->vz_horizontal_average[i] *= dy / Ly;
    }
}