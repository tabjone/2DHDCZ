#include "structs.h"

void deep_copy_foreground_2D(struct ForegroundVariables2D *destination, struct ForegroundVariables2D *source, struct GridInfo *grid_info) 
{
    int nz_full = grid_info->nz_full;
    int nx = grid_info->nx;

    allocate_foreground_struct_2D(&destination, nz_full, nx);

    for (int i = 0; i < nz_full; i++) 
    {
        for (int j = 0; j < nx; j++) 
        {
            destination->p1[i][j] = source->p1[i][j];
            destination->rho1[i][j] = source->rho1[i][j];
            destination->T1[i][j] = source->T1[i][j];
            destination->s1[i][j] = source->s1[i][j];
            destination->vx[i][j] = source->vx[i][j];
            destination->vz[i][j] = source->vz[i][j];
        }
    }

}