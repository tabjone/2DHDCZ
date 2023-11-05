#include "boundary_3D.h"
 
void calculate_damping_3D(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo3D *grid_info)
{
    /*
    Calculates the damping factor for the given boundary conditions.

    Parameters
    ----------
    damping_factor : FLOAT_P*
        A pointer to the array where the damping factor will be stored.
    grid_info : struct
        A pointer to the GridInfo2D struct.
    mpi_info : struct
        A pointer to the MpiInfo struct.
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;

    // Setting damping factor
    #if VERTICAL_BOUNDARY_TYPE == 0
        int nz_ghost = grid_info->nz_ghost;
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
        // Soft wall
        // Getting grid info
        int nz_ghost = grid_info->nz_ghost;
        // Calculate what SOFT_WALL_HEIGHT_PERCENTAGE of the domain is
        FLOAT_P damped_domain = (bg->r[nz_full-nz_ghost-1] - bg->r[nz_ghost]) * SOFT_WALL_HEIGHT_PERCENTAGE;
        // Calculate z0 for the bottom and top
        FLOAT_P z0_bot = bg->r[nz_ghost] + damped_domain;
        FLOAT_P z0_top = bg->r[nz_full-nz_ghost-1] - damped_domain;
        // Calculate zb for the bottom and top
        FLOAT_P zb_bot = bg->r[nz_ghost];
        FLOAT_P zb_top = bg->r[nz_full-nz_ghost-1];

        // First setting all to 1 for the inside of the grid
        for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
        {
            damping_factor[i] = 1.0;
        }
        // Then 0 for the ghost cells
        for (int i = 0; i < nz_ghost+1; i++)
        {
            damping_factor[i] = 0.0;
            damping_factor[nz_full-1-i] = 0.0;
        }
        
        // Then the soft wall for the bottom boundary
        int i = nz_ghost;
        while (bg->r[i]<= z0_bot)
        {
            // Calculating the weight
            FLOAT_P w = pow((bg->r[i]-zb_bot)/(z0_bot-zb_bot),ALPHA);
            // Setting the damping factor
            damping_factor[i] = w;
            i++;
        }
        
        // Soft wall for top boundary
        i = nz_full - nz_ghost - 1;
        while (bg->r[i]>= z0_top)
        {
            // Calculating the weight
            FLOAT_P w = pow((zb_top-bg->r[i])/(zb_top-z0_top),ALPHA);
            // Setting the damping factor
            damping_factor[i] = w;
            i--;
        }
         
    #elif VERTICAL_BOUNDARY_TYPE == 2
        // No damping, periodic
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