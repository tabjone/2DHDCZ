#include "boundary_3D.h"

void apply_vertical_boundary_damping_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt)
{
    
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    int nx = grid_info->nx;

    // Pre-Calculate this in main
    // Calculating damping factor
    //FLOAT_P damping_factor[nz_full];

    FLOAT_P *damping_factor;
    allocate_1D_array(&damping_factor, nz_full);

    #if MPI_ON == 0
        calculate_damping_3D(damping_factor, bg, grid_info);
    #elif MPI_ON == 1
        calculate_damping_mpi_3D(damping_factor, bg, grid_info, mpi_info);
    #endif // MPI_ON
    
    // Finding mean of s1
    FLOAT_P s1_mean = 0.0;
    FLOAT_P tau = 60.0*60.0*4;
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                s1_mean += fg->s1[i][j][k];
            }
        }
    }
    s1_mean /= (nz_full - 2*nz_ghost)*ny*nx;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->s1[i][j][k] = damping_factor[i]*fg->s1[i][j][k];
                fg->vz[i][j][k] = damping_factor[i]*fg->vz[i][j][k];
            }
        }
    }
    deallocate_1D_array(damping_factor);
}