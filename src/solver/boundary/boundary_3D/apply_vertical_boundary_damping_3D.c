#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "MPI_module/mpi_info_struct.h"
#include "global_float_precision.h"
#include "global_boundary.h"
#include "boundary_3D.h"
#include "array_utilities/array_memory_management/array_memory_management.h"

void apply_vertical_boundary_damping_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt)
{
    
    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    int nx = grid_info->nx;

    // This should be pre-calculated and put in precalc struct
    FLOAT_P *damping_factor;
    allocate_1D_array(&damping_factor, nz_full);
    calculate_damping_3D(damping_factor, bg, grid_info, mpi_info);

    // Damping factor for soft wall, hard wall
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

    // Boundary restrictions on vz, s1, p1
    if (!mpi_info->has_neighbor_above)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->s1[nz_full-nz_ghost-1][j][k] = UPPER_ENTROPY_BOUNDARY;
                fg->p1[nz_full-nz_ghost-1][j][k] = UPPER_PRESSURE_BOUNDARY;
                fg->vz[nz_full-nz_ghost-1][j][k] = 0.0;

            }
        }
    }
    if (!mpi_info->has_neighbor_below)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->s1[nz_ghost][j][k] = LOWER_ENTROPY_BOUNDARY;
                fg->p1[nz_ghost][j][k] = LOWER_PRESSURE_BOUNDARY;
                fg->vz[nz_ghost][j][k] = 0.0;
            }
        }
    }
    /*
    // Boundary restrictions on vy
    if (!mpi_info->has_neighbor_below)
    {
        for (int i = 0; i < 2*nz_ghost; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                fg->vy[i][j] = fg->vy[2*nz_ghost][j];

            }
        }
    }
    if (!mpi_info->has_neighbor_above)
    {
        for (int i = nz_full - 2*nz_ghost; i < nz_full; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                fg->vy[i][j] = fg->vy[nz_full - 2*nz_ghost - 1][j];
            }
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
    */
    deallocate_1D_array(damping_factor);
}