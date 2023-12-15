#include "boundary.h"

void apply_vertical_boundary_damping(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt)
{
    
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    // Pre-Calculate this in precalc
    // Calculating damping factor
    //FLOAT_P damping_factor[nz_full];

    FLOAT_P *damping_factor;
    allocate_1D_array(&damping_factor, nz_full);

    #if MPI_ON == 0
        calculate_damping(damping_factor, bg, grid_info);
    #elif MPI_ON == 1
        calculate_damping_mpi(damping_factor, bg, grid_info, mpi_info);
    #endif // MPI_ON

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[i][j] = damping_factor[i]*fg->s1[i][j];
            fg->vz[i][j] = damping_factor[i]*fg->vz[i][j];
        }
    }

    // On boundary
    if (!mpi_info->has_neighbor_above)
    {
        for (int j = 0; j < ny; j++)
        {
            *fg->s1[nz_full-nz_ghost-1][j] = UPPER_ENTROPY_BOUNDARY;
            *fg->p1[nz_full-nz_ghost-1][j] = UPPER_PRESSURE_BOUNDARY;
            *fg->vz[nz_full-nz_ghost-1][j] = 0.0;
        }
    }
    if (!mpi_info->has_neighbor_below)
    {
        for (int j = 0; j < ny; j++)
        {
            *fg->s1[nz_ghost][j] = LOWER_ENTROPY_BOUNDARY;
            *fg->p1[nz_ghost][j] = LOWER_PRESSURE_BOUNDARY;
            *fg->vz[nz_ghost][j] = 0.0;
        }
    }

    #if OSCILLATING_S1_BC == 1
        FLOAT_P s1_mean = 0.0;
        for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                s1_mean += fg->s1[i][j];
            }
        }
        s1_mean /= (nz_full - 2*nz_ghost)*ny;
        for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                fg->s1[i][j] += + s1_mean*(1-damping_factor[i])*(1-dt/TAU_BC);
            }
        }
    #endif // OSCILLATING_S1_BC
    deallocate_1D_array(damping_factor);
    
}