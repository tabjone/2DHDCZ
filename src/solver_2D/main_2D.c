#include "functions.h"
#include <mpi.h>

int main_2D(int argc, char *argv[], struct MpiInfo *mpi_info)
{
    /*
    This is a GOD function. It should be simplified and split into smaller functions.
    */


    FLOAT_P t, dt, dt_last, t_since_save;
    FLOAT_P first_t;
    int save_nr;

    printf("0\n");
    // Declare the background variables, foreground variables and grid info
    struct BackgroundVariables *bg;
    struct ForegroundVariables2D *fg, *fg_previous, *tmp_ptr;
    struct GridInfo2D *grid_info = NULL;
    
    // Function load_snapshot_2D()
    #if LOAD == 1
        // Loading snapshot

        char file_path_foreground[150];
        char file_path_background[150];
        snprintf(file_path_foreground, sizeof(file_path_foreground), "data/%s/snap%d.h5", RUN_NAME, LOAD_SNAP_NUMBER);
        snprintf(file_path_background, sizeof(file_path_background), "data/%s/background.h5", RUN_NAME);
                
        load_grid_info(&grid_info, file_path_foreground);
        allocate_background_struct(&bg, grid_info->nz_full);
        allocate_foreground_struct_2D(&fg_previous, grid_info);
        allocate_foreground_struct_2D(&fg, grid_info);

        load_background(bg, grid_info, file_path_background);

        t = load_foreground(fg_previous, grid_info, file_path_foreground);
        first_t = t;

        printf("time = %f\n", t);
        save_nr = LOAD_SNAP_NUMBER + 1;

    // Function initialize_simulation_2D()
    #elif LOAD == 0
        // Initializing the simulation
        save_nr = 0;
        printf("1\n");
        calculate_grid_info_mpi(mpi_info, &grid_info);
        printf("2\n");
        // Allocating memory for the background and foreground variables
        allocate_background_struct(&bg, grid_info->nz_full);
        allocate_foreground_struct_2D(&fg, grid_info);
        allocate_foreground_struct_2D(&fg_previous, grid_info);

        printf("3\n");
        // Initialize the background variables and saving it to file
        solar_s_background_initialization(bg, mpi_info, grid_info->nz_full, grid_info->nz_ghost, grid_info->dz, grid_info->z0, grid_info->z1, grid_info->nz);

        save_background(bg, mpi_info, grid_info->nz_full, grid_info->nz, grid_info->nz_ghost, grid_info->dz, grid_info->z0, grid_info->z1);

        printf("neig_above = %d, neig_below = %d\n", mpi_info->has_neighbor_above, mpi_info->has_neighbor_below);

        save_mpi_info(mpi_info);
        save_info(mpi_info);
        
        // Initialize foreground to type set in parameter file
        initialize_foreground(fg_previous, bg, grid_info, mpi_info);
        
        // Saving the foreground variables to file
        save_foreground(fg_previous, grid_info, mpi_info, 0, 0.0);
        save_nr ++;
    
        t = 0.0;
        first_t = 0.0;
    #endif // LOAD
    
    printf("Initialization done\n");

    printf("Precalculating\n");
    struct PrecalculatedVariables *precalc;

    allocate_calculate_precalc_struct(&precalc, bg, grid_info);
    
    #if SAVE_RHS == 1
        save_rhs(fg_previous, bg, grid_info, mpi_info, precalc, 0);
        save_elliptic_vars(fg_previous, bg, grid_info, mpi_info, precalc, 0);
    #endif // SAVE_RHS

    t_since_save = 0.0;
    dt_last = 0.0;
    while (t < T)
    {
        dt = one_time_step(bg, fg_previous, fg, grid_info, mpi_info, precalc, dt_last, first_t == t);
        t += dt;


        if (mpi_info->rank == 0)
            printf("t = %.2f, dt=%.2f\n", t, dt);
        
        t_since_save += dt;
        dt_last = dt;
        
        if (dt_last < 0.1)
        {
            break;
        }
        
        if (t_since_save > SAVE_INTERVAL && SAVE_ALL == 0)
        {
            save_foreground(fg, grid_info, mpi_info, save_nr, t);
            #if SAVE_RHS == 1
                save_rhs(fg, bg, grid_info, mpi_info, precalc, save_nr);
                save_elliptic_vars(fg, bg, grid_info, mpi_info, precalc, save_nr);
            #endif // SAVE_RHS
            save_nr++;
            t_since_save = 0.0;
        }
        else if (SAVE_ALL == 1)
        {
            save_foreground(fg, grid_info, mpi_info, save_nr, t);
            #if SAVE_RHS == 1
                save_rhs(fg, bg, grid_info, mpi_info, precalc, save_nr);
                save_elliptic_vars(fg, bg, grid_info, mpi_info, precalc, save_nr);
            #endif // SAVE_RHS
            save_nr++;
        }

        // pointer swap
        tmp_ptr = fg_previous;
        fg_previous = fg;
        fg = tmp_ptr;
    }

    // Save last time step
    save_foreground(fg_previous, grid_info, mpi_info, save_nr, t);
    #if SAVE_RHS == 1
        save_elliptic_vars(fg_previous, bg, grid_info, mpi_info, precalc, save_nr);
        save_rhs(fg_previous, bg, grid_info, mpi_info, precalc, save_nr);
    #endif // SAVE_RHS
    
    deallocate_grid_info_struct_2D(grid_info);
    deallocate_background_struct(bg);
    deallocate_precalc_struct(precalc);
    deallocate_foreground_struct_2D(fg_previous);
    deallocate_foreground_struct_2D(fg);

    return 0;
}