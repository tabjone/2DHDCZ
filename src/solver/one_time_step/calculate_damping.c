#include "one_time_step.h"
 
void calculate_damping(FLOAT_P *damping_factor, struct GridInfo *grid_info)
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
    int nz_full = grid_info->nz_full;

    // Setting damping factor
    #if VERTICAL_BOUNDARY_TYPE == 0
        int nz_ghost = grid_info->nz_ghost;

        printf("Hard wall.\n");
        // No damping
        for (int i = 0; i < nz_full; i++)
        {
            damping_factor[i] = 1.0;
        }
        // Boundaries
        for (int i = 0; i < nz_ghost+1; i++)
        {
            damping_factor[i] = 0.0;
        }
        for (int i = nz_full-nz_ghost-1; i < nz_full; i++)
        {
            damping_factor[i] = 0.0;
        }

    #elif VERTICAL_BOUNDARY_TYPE == 1
        int nz_ghost = grid_info->nz_ghost;

        printf("Soft wall.\n");
        // Damping at the top and bottom
        FLOAT_P damping_coeffs[5] = {0.0, 0.5, 0.75, 0.875, 0.9375};

        for (int i = 0; i < nz_full; i++)
        {
            damping_factor[i] = 1.0;
        }

        for (int i = 0; i < nz_ghost; i++)
        {
            damping_factor[i] = 0.0;
        }
        for (int i = nz_ghost; i < nz_ghost+5; i++)
        {
            damping_factor[i] = damping_coeffs[i-nz_ghost];
        }
        
        for (int i = nz_full-nz_ghost; i < nz_full; i++)
        {
            damping_factor[i] = 0.0;
        }
        for (int i = nz_full - nz_ghost - 5; i < nz_full - nz_ghost; i++) 
        {
            damping_factor[i] = damping_coeffs[nz_full - nz_ghost - 1 - i];
        }
    #elif VERTICAL_BOUNDARY_TYPE == 2
        printf("Periodic.\n");
        // No damping
        for (int i = 0; i < nz_full; i++)
        {
            damping_factor[i] = 1.0;
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    #if DEBUG == 1
        for (int i = 0; i < nz_full; i++)
        {
            //printf("damping_factor[%d] = %f\n", i, damping_factor[i]);
        }
    #endif // DEBUG
}