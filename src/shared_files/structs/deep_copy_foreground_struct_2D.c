#include "structs.h"

void deep_copy_foreground_2D(struct ForegroundVariables2D *destination, struct ForegroundVariables2D *source) 
{
    int nz_full = source->nz_full;
    int nx = source->nx;
    double dz = source->dz;
    double dx = source->dx;

    allocate_foreground_struct_2D(nz_full, nx, dz, dx, &destination);

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
    
    destination->nz = source->nz;
    destination->nz_ghost = source->nz_ghost;
    destination->nz_full = source->nz_full;
    destination->nx = source->nx;
    destination->dx = source->dx;
    destination->dz = source->dz;
}