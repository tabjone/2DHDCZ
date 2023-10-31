#include "boundary.h"

void apply_vertical_boundary_damping(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt)
{
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    // Pre-Calculate this in main
    // Calculating damping factor
    FLOAT_P damping_factor[nz_full];
    #if MPI_ON == 0
        calculate_damping(damping_factor, bg, grid_info);
    #elif MPI_ON == 1
        calculate_damping_mpi(damping_factor, bg, grid_info, mpi_info);
    #endif // MPI_ON

    // Finding mean of s1
    FLOAT_P s1_mean = 0.0;
    FLOAT_P tau = 60.0*60.0*4;
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            s1_mean += fg_prev->s1[i][j];
        }
    }
    s1_mean /= (nz_full - 2*nz_ghost)*ny;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[i][j] = damping_factor[i]*fg->s1[i][j] + s1_mean*(1-damping_factor[i])*(1-dt/tau);
            fg->vz[i][j] = damping_factor[i]*fg->vz[i][j];
        }
    }


}