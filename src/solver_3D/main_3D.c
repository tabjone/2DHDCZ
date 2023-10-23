#include "functions.h"
#include <mpi.h>

int main_3D(int argc, char *argv[], struct MpiInfo *mpi_info)
{
    
    FLOAT_P t, dt, dt_last, t_since_save;
    FLOAT_P first_t;
    int save_nr;

    // Declare the background variables, foreground variables and grid info
    struct BackgroundVariables *bg;
    struct ForegroundVariables3D *fg, *fg_previous, *tmp_ptr;
    struct GridInfo3D *grid_info;
    
    #if LOAD == 1
    #elif LOAD == 0
        // Initializing the simulation
        save_nr = 0;

        
        #if MPI_ON == 0
            calculate_grid_info(&grid_info, mpi_info);
        #elif MPI_ON == 1
            calculate_grid_info_mpi(mpi_info, &grid_info);
        #endif // MPI_ON

        // Allocating memory for the background and foreground variables
        allocate_background_struct(&bg, grid_info);
        allocate_foreground_struct_3D(&fg, grid_info);
        allocate_foreground_struct_3D(&fg_previous, grid_info);

        // Initialize the background variables and saving it to file
        solar_s_background_initialization(bg, mpi_info, grid_info->nz_full, grid_info->nz_ghost, grid_info->dz, grid_info->z0, grid_info->z1, grid_info->nz);
        communicate_background_ghost_above_below(bg, grid_info, mpi_info);
        save_background(bg, grid_info, mpi_info);
        save_mpi_info(mpi_info);

        
        // Initialize foreground to type set in parameter file
        initialize_foreground(fg_previous, bg, grid_info, mpi_info);
        
        // Saving the foreground variables to file
        save_foreground(fg_previous, grid_info, mpi_info, 0, 0.0);
        save_nr ++;
    
        t = 0.0;
        first_t = 0.0;
    #endif // LOAD
    /*
    t_since_save = 0.0;
    dt_last = 0.0;
    while (t < T)
    { 
        dt = one_time_step(bg, fg_previous, fg, grid_info, mpi_info, dt_last, first_t == t);
        t += dt;

        t_since_save += dt;
        dt_last = dt;
        
        if (t_since_save > SAVE_INTERVAL && SAVE_ALL == 0)
        {
            save_foreground(fg, grid_info, mpi_info, save_nr, t);
            save_nr++;
            t_since_save = 0.0;
        }
        else if (SAVE_ALL == 1)
        {
            save_foreground(fg, grid_info, mpi_info, save_nr, t);
            save_nr++;
        }

        // pointer swap
        tmp_ptr = fg_previous;
        fg_previous = fg;
        fg = tmp_ptr;
        if (mpi_info->rank == 0)
            printf("t = %.2f\n", t);
        break;
    }

    // Save last time step
    save_foreground(fg_previous, grid_info, mpi_info, save_nr, t);
    
    deallocate_grid_info_struct(grid_info);
    deallocate_background_struct(bg);
    deallocate_foreground_struct(fg_previous);
    deallocate_foreground_struct(fg);
    */
    return 1;
}