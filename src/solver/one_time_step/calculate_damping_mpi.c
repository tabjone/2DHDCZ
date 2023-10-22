#include "one_time_step.h"
 
void calculate_damping_mpi(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Calculates the damping factor for the given boundary conditions.

    Parameters
    ----------
    damping_factor : FLOAT_P*
        A pointer to the array where the damping factor will be stored.
    grid_info : struct
        A pointer to the GridInfo struct.
    mpi_info : struct
        A pointer to the MpiInfo struct.
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
        // Soft wall
        // Getting grid info
        int nz_ghost = grid_info->nz_ghost;
        int nz_grid = grid_info->nz;
        // Calculate what SOFT_WALL_HEIGHT_PERCENTAGE of the domain is
        // Bottom r
        FLOAT_P bottom_r;
        if (mpi_info->rank == 0)
        {
            // Broadcast the bottom r to all processes
            MPI_Bcast(&bg->r[nz_ghost], 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);
            bottom_r = bg->r[nz_ghost];
        }
        else
        {
            // Receive the bottom r from process 0
            MPI_Bcast(&bottom_r, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);
        }
        // Top r
        FLOAT_P top_r;
        if (mpi_info->rank == mpi_info->size-1)
        {
            // Broadcast the top r to all processes
            MPI_Bcast(&bg->r[nz_full-nz_ghost-1], 1, MPI_FLOAT_P, mpi_info->size-1, MPI_COMM_WORLD);
            top_r = bg->r[nz_full-nz_ghost-1];
        }
        else
        {
            // Receive the top r from process 0
            MPI_Bcast(&top_r, 1, MPI_FLOAT_P, mpi_info->size-1, MPI_COMM_WORLD);
        }
        



        FLOAT_P damped_domain = (top_r - bottom_r) * SOFT_WALL_HEIGHT_PERCENTAGE;
        // Calculate z0 for the bottom and top
        FLOAT_P z0_bot = bottom_r + damped_domain;
        FLOAT_P z0_top = top_r - damped_domain;
        // Calculate zb for the bottom and top
        FLOAT_P zb_bot = bottom_r;
        FLOAT_P zb_top = top_r;

        // First setting all to 1 for the inside of the grid
        for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
        {
            damping_factor[i] = 1.0;
        }
        // Then 0 for the ghost cells if there is no neighbor above or below
        if (!mpi_info->has_neighbor_below)
        {
            for (int i = 0; i < nz_ghost; i++)
            {
                damping_factor[i] = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < nz_ghost; i++)
            {
                damping_factor[i] = 1.0;
            }
        }
        if (!mpi_info->has_neighbor_above)
        {
            for (int i = nz_full-nz_ghost; i < nz_full; i++)
            {
                damping_factor[i] = 0.0;
            }
        }
        else
        {
            for (int i = nz_full-nz_ghost; i < nz_full; i++)
            {
                damping_factor[i] = 1.0;
            }
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
        printf("Periodic.\n");
        // No damping
        for (int i = 0; i < nz_full; i++)
        {
            damping_factor[i] = 1.0;
        }
    #endif // VERTICAL_BOUNDARY_TYPE
}