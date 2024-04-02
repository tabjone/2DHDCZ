#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "MPI_module/mpi_info_struct.h"
#include "global_float_precision.h"
#include "global_parameters.h"
#include "global_constants.h"
#include "global_boundary.h"
#include <mpi.h>
#include <math.h>
#include <stdbool.h>

#include <stdio.h>

void calculate_damping_2D(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
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

    #if HARD_WALL_VERTICAL == 1
        int nz_ghost = grid_info->nz_ghost;

        // No damping inside grid
        for (int i = 0; i < nz_full; i++)
        {
            damping_factor[i] = 1.0;
        }
        if (mpi_info->rank == 0) // Bottom boundary
        {   
            //damping_factor[nz_ghost+1] = 0.5;
            for (int i = 0; i < nz_ghost+1; i++)
            {
                damping_factor[i] = 0.0;
            }
        }
        if (mpi_info->rank == mpi_info->size-1) // Top boundary
        {
            for (int i = nz_full-nz_ghost-1; i < nz_full; i++)
            {
                damping_factor[i] = 0.0;
            }
            //printf("damping_factor[%d] = %f\n", nz_full-nz_ghost-18, damping_factor[nz_full-nz_ghost-18]);
            //printf("damping_factor[%d] = %f\n", nz_full-nz_ghost-17, damping_factor[nz_full-nz_ghost-17]);
            // Top 10 cells should go linearly from 0 to 1
            for (int i = nz_full-nz_ghost-16; i < nz_full-nz_ghost-1; i++)
            {
                damping_factor[i] = (FLOAT_P)(nz_full-nz_ghost-1-i)/15.0;
                //printf("damping_factor[%d] = %f\n", i, damping_factor[i]);
            }
            //printf("damping_factor[%d] = %f\n", nz_full-nz_ghost-1, damping_factor[nz_full-nz_ghost-1]);
        }
    #endif // HARD_WALL_VERTICAL

    #if SOFT_WALL_VERTICAL == 1
        // Getting grid info
        int nz_ghost = grid_info->nz_ghost;
        // Calculate what SOFT_WALL_HEIGHT_PERCENTAGE of the domain is
        // Bottom r

        FLOAT_P bottom_r, top_r;
        if (mpi_info->rank == 0)
        {
            // Broadcast the bottom r to all processes
            
            bottom_r = bg->r[nz_ghost];
        }
        if (mpi_info->rank == mpi_info->size-1)
        {
            // Broadcast the top r to all processes
            
            top_r = bg->r[nz_full-nz_ghost-1];
        }

        MPI_Bcast(&bottom_r, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);
        MPI_Bcast(&top_r, 1, MPI_FLOAT_P, mpi_info->size-1, MPI_COMM_WORLD);


        FLOAT_P damped_domain = (top_r - bottom_r) * SOFT_WALL_HEIGHT_PERCENTAGE_VERTICAL;
        FLOAT_P damped_domain_bot = damped_domain/2.0;
        // Calculate z0 for the bottom and top
        FLOAT_P z0_bot = bottom_r + damped_domain_bot;
        FLOAT_P z0_top = top_r - damped_domain;

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

        FLOAT_P damping_bottom, damping_top;
        FLOAT_P width = SOFT_WALL_WIDTH*R_SUN;
        // Going trough domain, including boundaries and calculating the damping factor
        for (int i = nz_ghost+1; i < nz_full-nz_ghost; i++)
        {
            damping_bottom = (tanh((bg->r[i] - z0_bot) / width) + 1.0) / 2.0;
            //damping_bottom = 1.0;
            damping_top = (tanh((z0_top - bg->r[i]) / width) + 1.0) / 2.0;

            damping_factor[i] = damping_bottom * damping_top;
        }

        // Setting boundaries to zero if no neighbor above or below
        if (!mpi_info->has_neighbor_below)
        {
            damping_factor[nz_ghost] = 0.0;
        }
        if (!mpi_info->has_neighbor_above)
        {
            damping_factor[nz_full-nz_ghost-1] = 0.0;
        }

    // Finding where tanh(x) > 0.8 so we save this index in the grid and use it in apply_vertical_boundary_damping_2D.c
    bool hit = false;
    int my_hit_index = 0;
    for (int i = nz_full-nz_ghost-1; i >= nz_ghost; i--)
    {
        if ((tanh((z0_top - bg->r[i]) / width) + 1.0) / 2.0 >= 1.0)
        {
            hit = true;
            my_hit_index = i;
            break;
        }
    }
    // We find the top most process that hit
    int my_hit = hit * mpi_info->rank;
    if (hit)
    {
        if (my_hit_index < nz_ghost || my_hit_index > nz_full-nz_ghost)
        {
            my_hit = 0;
        }
    }

    int global_hit;
    MPI_Allreduce(&my_hit, &global_hit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    mpi_info->soft_wall_end_process = global_hit;
    mpi_info->soft_wall_end_index = my_hit_index;
    #endif // SOFT_WALL_VERTICAL

    #if PERIODIC_BOUNDARY_VERTICAL == 1 || NO_BOUNDARY_VERTICAL == 1
        for (int i = 0; i < nz_full; i++)
        {
            damping_factor[i] = 1.0;
        }
    #endif // PERIODIC_BOUNDARY_VERTICAL


    for (int i = 0; i < nz_full; i++)
    {
        //printf("damping_factor[%d] = %f\n", i, damping_factor[i]);
    }

}