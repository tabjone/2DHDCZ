#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"
#include "MPI_module/mpi_info_struct.h"
#include "global_float_precision.h"
#include "global_boundary.h"
#include "boundary_2D.h"
#include "array_utilities/array_memory_management/array_memory_management.h"

#include "global_constants.h"
#include "global_parameters.h"
#include <math.h>

// Structure to handle both the sum and the compensation
typedef struct {
    FLOAT_P sum;
    FLOAT_P c;  // Compensation variable
} KahanAccumulator;

void kahan_add(KahanAccumulator *acc, FLOAT_P value) {
    FLOAT_P y = value - acc->c;  // Subtract the compensation
    FLOAT_P t = acc->sum + y;    // Add the value to the sum
    acc->c = (t - acc->sum) - y; // Update the compensation
    acc->sum = t;                // Update the sum
}

void apply_vertical_boundary_damping_2D(struct ForegroundVariables2D *fg, struct ForegroundVariables2D *fg_prev, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *pv, FLOAT_P dt)
{
    
    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    
    /*
    Finding the mass density of the grid and setting the top boundary to be the average of the grid.
    */
    KahanAccumulator my_density_total = {0.0, 0.0};
    FLOAT_P density_grid;

    int end_index = nz_full - nz_ghost;
    if (!mpi_info->has_neighbor_above) {
        end_index -= 1;
    }

    // Sum densities avoiding ghost zones
    for (int i = nz_ghost; i < end_index; i++) {
        for (int j = 0; j < ny; j++) {
            kahan_add(&my_density_total, fg_prev->rho1[i][j]);
        }
    }

    // Use MPI to sum densities across all processes
    MPI_Allreduce(&my_density_total.sum, &density_grid, 1, MPI_FLOAT_P, MPI_SUM, MPI_COMM_WORLD);

    FLOAT_P rho1_add = -density_grid/(1.0*ny);

    // If rho1_add is positive it will create negative entropy, then we should put it at the top of the grid. If it is negative it will create positive entropy, then we should put it at the bottom of the grid.
    
    FLOAT_P rho1top, rho1bot;
    if (rho1_add >= 0.0)
    {
        rho1top = rho1_add;
        rho1bot = 0.0;
    }
    else
    {
        rho1top = 0.0;
        rho1bot = rho1_add;
    }

    FLOAT_P rho0top, rho0bot;
    if (!mpi_info->has_neighbor_above)
    {
        rho0top = bg->rho0[nz_full-nz_ghost-1];
        for (int j = 0; j < ny; j++)
        {   
            fg_prev->rho1[nz_full-nz_ghost-1][j] = rho1top;
        }      
    }
    if (!mpi_info->has_neighbor_below)
    {
        rho0bot = bg->rho0[nz_ghost];
        for (int j = 0; j < ny; j++)
        {
            fg_prev->rho1[nz_ghost][j] = rho1bot;
        }
    }
    MPI_Bcast(&rho0top, 1, MPI_FLOAT_P, mpi_info->size - 1, MPI_COMM_WORLD);
    MPI_Bcast(&rho0bot, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);

    // Finding the entropy at the boundaries from the equation of state
    FLOAT_P c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);
    FLOAT_P s1_top = - rho1top / rho0top * c_p;
    FLOAT_P s1_bot = - rho1bot / rho0bot * c_p;

    // Top boundary
    if (!mpi_info->has_neighbor_above)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[nz_full-nz_ghost-1][j] += s1_top;
            fg->vz[nz_full-nz_ghost-1][j] = 0.0; 
        }
    }

    // Bottom boundary
    if (!mpi_info->has_neighbor_below)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[nz_ghost][j] += s1_bot;
            fg->vz[nz_ghost][j] = 0.0;
        }
    }



    /*
    FLOAT_P c_p, s1_bot;
    FLOAT_P p1_top, p1_bot, p0_top, p0_bot;
    FLOAT_P rho1_bot, rho0_top, rho0_bot;
    
    // Same for all processes
    c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);
    s1_bot = LOWER_ENTROPY_BOUNDARY * c_p;
    p1_top = UPPER_PRESSURE_BOUNDARY;
    p1_bot = LOWER_PRESSURE_BOUNDARY;

    // Top process finds top values and broadcasts them
    if (!mpi_info->has_neighbor_above)
    {
        p0_top = bg->p0[nz_full-nz_ghost-1];
        rho0_top = bg->rho0[nz_full-nz_ghost-1];
    }
    MPI_Bcast(&p0_top, 1, MPI_FLOAT_P, mpi_info->size-1, MPI_COMM_WORLD);
    MPI_Bcast(&rho0_top, 1, MPI_FLOAT_P, mpi_info->size-1, MPI_COMM_WORLD);

    // Bottom process finds bottom values and broadcasts them
    if (!mpi_info->has_neighbor_below)
    {
        p0_bot = bg->p0[nz_ghost];
        rho0_bot = bg->rho0[nz_ghost];
    }
    MPI_Bcast(&p0_bot, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rho0_bot, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);

    rho1_bot = rho0_bot * (-s1_bot/c_p + 1.0/GAMMA*p1_bot/p0_bot);

    s1_top = (rho1_bot/rho0_top + 1.0/GAMMA*p1_top/p0_top)*c_p;
    
    
    //TRY TO MAKE THE BOUNDARY CONVECTIVE: s1 = (s1_boundary - s1)
    
    

    // Top boundary
    if (!mpi_info->has_neighbor_above)
    {
        c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);
        for (int j = 0; j < ny; j++)
        {
            //fg->p1[nz_full-nz_ghost-1][j] = UPPER_PRESSURE_BOUNDARY;
            fg->vz[nz_full-nz_ghost-1][j] = 0.0;
            fg->s1[nz_full-nz_ghost-1][j] = s1_top;
            //fg->s1[nz_full-nz_ghost-2][j] = -0.02 * c_p;
        }
    }
    if (mpi_info->rank == mpi_info->soft_wall_end_process)
    {
        c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);
        for (int j = 0; j < ny; j++)
        {
            //fg->s1[mpi_info->soft_wall_end_index][j] = -0.005 * c_p;
        }
    }

    // Bottom boundary
    if (!mpi_info->has_neighbor_below)
    {
        c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);
        for (int j = 0; j < ny; j++)
        {
            fg->s1[nz_ghost][j] = s1_bot;
            //fg->s1[nz_ghost+1][j] = 0.02 * c_p;
            //fg->p1[nz_ghost][j] = LOWER_PRESSURE_BOUNDARY;
            fg->vz[nz_ghost][j] = 0.0;
            //fg->vy[nz_ghost][j] = 0.0;
        }
    }*/
}