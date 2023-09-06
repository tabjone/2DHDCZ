#include "../../../shared_files/structs/structs.h"
#include "./initialization.h"
#include "../../../shared_files/extrapolation/extrapolation.h"

void initialize_foreground_struct_zeros(struct ForegroundVariables2D *foreground_variables)
{
    int nz = foreground_variables->nz;
    int nz_ghost = foreground_variables->nz_ghost;
    int nx = foreground_variables->nx;

    int i, j;

    for (i = nz_ghost; i < nz+nz_ghost; i++)
    {
        for (j = 0; j < nx; j++)
        {
            foreground_variables->p1[i][j] = 0.0;
            foreground_variables->rho1[i][j] = 0.0;
            foreground_variables->T1[i][j] = 0.0;
            foreground_variables->s1[i][j] = 0.0;
            foreground_variables->vx[i][j] = 0.0;
            foreground_variables->vz[i][j] = 0.0;
        }
    }

    extrapolate_2D(foreground_variables);
}