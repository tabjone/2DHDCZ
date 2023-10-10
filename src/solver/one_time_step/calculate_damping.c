#include "one_time_step.h"

void calculate_damping(FLOAT_P *damping_factor, struct GridInfo *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Calculates the damping factor for the given boundary conditions.

    Parameters
    ----------
    damping_factor : FLOAT_P*
        A pointer to the array where the damping factor will be stored.
    grid_info : struct
        A pointer to the GridInfo struct.
    */

    // Getting grid info
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    // Setting damping factor
    #if VERTICAL_BOUNDARY_TYPE == 0
        // No damping
        for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
        {
            damping_factor[i] = 1.0;
        }
        // But zero boundary
        if (mpi_info->has_neighbor_below)
        {
            damping_factor[nz_ghost] = 1.0;
        }
        else
        {
            damping_factor[nz_ghost] = 0.0;
        }
        if (mpi_info->has_neighbor_above)
        {
            damping_factor[nz_full-nz_ghost-1] = 1.0;
        }
        else
        {
            damping_factor[nz_full-nz_ghost-1] = 0.0;
        }

    #elif VERTICAL_BOUNDARY_TYPE == 1
        // Damping at the top and bottom
        FLOAT_P damping_coeffs[5] = {0.0, 0.5, 0.75, 0.875, 0.9375};

        // No damping in the middle
        for (int i = nz_ghost+5; i < nz_full-nz_ghost-5+1; i++)
        {
            damping_factor[i] = 1.0;
        }

        // Top damping if no neighbor above
        if (mpi_info->has_neighbor_above)
        {
            for (int i = nz_ghost; i < nz_ghost+5; i++)
            {
                damping_factor[nz_full-i-1] = 1.0;
            }
        }
        else
        {
            for (int i = nz_ghost; i < nz_ghost+5; i++)
            {
                damping_factor[nz_full-i-1] = damping_coeffs[i-nz_ghost];
            }
        }

        // Bottom damping if no neighbor below
        if (mpi_info->has_neighbor_below)
        {
            for (int i = nz_full-nz_ghost-5; i < nz_full-nz_ghost; i++)
            {
                damping_factor[i] = 1.0;
            }
        }
        else
        {
            for (int i = nz_full-nz_ghost-5; i < nz_full-nz_ghost; i++)
            {
                damping_factor[i] = damping_coeffs[i-nz_ghost];
            }
        }
    #endif // VERTICAL_BOUNDARY_TYPE
}