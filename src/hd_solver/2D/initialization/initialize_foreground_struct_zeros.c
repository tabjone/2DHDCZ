#include "../../../shared_files/structs/structs.h"
#include "./initialization.h"

void initialize_foreground_struct_zeros(struct ForegroundVariables2D *foreground_variables)
{
    int nz = foreground_variables->nz;
    int nx = foreground_variables->nx;

    int i, j;

    for (i = 0; i < nz; i++)
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
}