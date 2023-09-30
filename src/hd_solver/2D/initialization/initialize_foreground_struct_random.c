#include <stdlib.h>
#include "initialization.h"
#include "global_parameters.h"

void initialize_foreground_struct_random(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg , struct GridInfo *grid_info)
{
    int nx = grid_info->nx;
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;


    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->p1[i][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e9;  // random between -1e9 to 1e9
            fg->s1[i][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
            
            fg->vx[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
        }
    }
    extrapolate_2D_array(fg->p1, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, nx);

    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);


    extrapolate_2D(fg, grid_info);

}
