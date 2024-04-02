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

void apply_vertical_boundary_damping_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *pv, FLOAT_P dt)
{
    
    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    // This should be pre-calculated and put in precalc struct
    FLOAT_P *damping_factor = pv->damping_factor;

    FLOAT_P vz_damping;
    FLOAT_P full_damping;
    FLOAT_P width_vz = 1.0e1;
    
    FLOAT_P r, r_range;
    r_range = (R_END - R_START)*R_SUN;
    // Damping factor for soft wall, hard wall
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        r = bg->r[i] - R_START*R_SUN;
        for (int j = 0; j < ny; j++)
        {   
            if (r > r_range/2.0)
            {
                vz_damping = 0.5 * (tanh(fg->vz[i][j]/width_vz)+1.0);
            }
            else 
            {
                vz_damping = 0.5 * (tanh(-fg->vz[i][j]/width_vz)+1.0);
            }
            full_damping = 1.0 + vz_damping * (damping_factor[i]-1.0);

            //fg->vz[i][j] = full_damping*fg->vz[i][j];
            //fg->s1[i][j] = damping_factor[i]*fg->s1[i][j];
            
            //fg->vy[i][j] = damping_factor[i]*fg->vy[i][j];
            //fg->p1[i][j] = damping_factor[i]*fg->p1[i][j];
        }
    }

    FLOAT_P c_p;

    // Top boundary
    if (!mpi_info->has_neighbor_above)
    {
        c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);
        for (int j = 0; j < ny; j++)
        {
            //fg->p1[nz_full-nz_ghost-1][j] = UPPER_PRESSURE_BOUNDARY;
            fg->vz[nz_full-nz_ghost-1][j] = 0.0;
            fg->s1[nz_full-nz_ghost-1][j] = UPPER_ENTROPY_BOUNDARY * c_p;
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
            fg->s1[nz_ghost][j] = LOWER_ENTROPY_BOUNDARY * c_p;
            //fg->s1[nz_ghost+1][j] = 0.02 * c_p;
            //fg->p1[nz_ghost][j] = LOWER_PRESSURE_BOUNDARY;
            fg->vz[nz_ghost][j] = 0.0;
            //fg->vy[nz_ghost][j] = 0.0;
        }
    }
}