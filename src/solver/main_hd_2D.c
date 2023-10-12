#include "functions.h"
#include <mpi.h>

int main_hd_2D(int argc, char *argv[])
{
    FLOAT_P t, dt, dt_last, t_since_save;
    FLOAT_P first_t;
    int save_nr;

    // Declare the background variables, foreground variables and grid info
    struct BackgroundVariables *bg;
    struct ForegroundVariables *fg, *fg_previous, *tmp_ptr;
    struct GridInfo *grid_info;

    #if LOAD == 1
        // Loading snapshot

        char file_path_foreground[150];
        char file_path_background[150];
        snprintf(file_path_foreground, sizeof(file_path_foreground), "data/%s/snap%d.h5", RUN_NAME, LOAD_SNAP_NUMBER);
        snprintf(file_path_background, sizeof(file_path_background), "data/%s/background.h5", RUN_NAME);
                
        load_grid_info(&grid_info, file_path_foreground);
        allocate_background_struct(&bg, grid_info);
        allocate_foreground_struct(&fg_previous, grid_info);
        allocate_foreground_struct(&fg, grid_info);

        load_background(bg, grid_info, file_path_background);

        t = load_foreground(fg_previous, grid_info, file_path_foreground);
        first_t = t;

        printf("time = %f\n", t);
        save_nr = LOAD_SNAP_NUMBER + 1;

    #elif LOAD == 0
        // Initializing the simulation
        save_nr = 0;

        // Calculating the size of the grid
        FLOAT_P L_z = (R_END - R_START)*R_SUN;
        FLOAT_P L_y = Y_SIZE*R_SUN;

        // Calculating the size of the grid cells
        FLOAT_P dy = L_y/(NY - 1);
        FLOAT_P dz = L_z/(NZ - 1);

        // Calculating the number of ghost cells
        int nz_ghost;
        if (UPWIND_ORDER >= CENTRAL_ORDER)
        {
            nz_ghost = UPWIND_ORDER;
        }
        else
        {
            nz_ghost = CENTRAL_ORDER;
        }

        // Calculating the number of cells in the full grid
        int nz_full = NZ + 2*nz_ghost;


        FLOAT_P z0 = R_SUN * R_START;
        FLOAT_P z1 = R_SUN * R_END;
        FLOAT_P y0 = 0.0;
        FLOAT_P y1 = R_SUN * Y_SIZE;

        allocate_grid_info_struct(&grid_info, NZ, nz_ghost, nz_full, NY, dz, dy, z0, z1, y0, y1);
        // Allocating memory for the background and foreground variables
        allocate_background_struct(&bg, grid_info);
        allocate_foreground_struct(&fg, grid_info);
        allocate_foreground_struct(&fg_previous, grid_info);

        // Initialize the background variables and saving it to file
        solar_s_background_initialization(bg, grid_info);
        save_background(bg, grid_info);
        
        // Initialize foreground to type set in parameter file
        initialize_foreground(fg_previous, bg, grid_info);
        
        // Saving the foreground variables to file
        save_foreground(fg_previous, grid_info, 0, 0.0);
        save_nr ++;
    
        t = 0.0;
        first_t = 0.0;
    #endif // LOAD

    t_since_save = 0.0;
    dt_last = 0.0;
    while (t < T)
    {
        break;
        dt = one_time_step(bg, fg_previous, fg, grid_info, dt_last, first_t == t);

        t += dt;
        t_since_save += dt;
        dt_last = dt;
        
        if (t_since_save > SAVE_INTERVAL && SAVE_ALL == 0)
        {
            save_foreground(fg, grid_info, save_nr, t);
            save_nr++;
            t_since_save = 0.0;
        }
        else if (SAVE_ALL == 1)
        {
            save_foreground(fg, grid_info, save_nr, t);
            save_nr++;
        }

        // pointer swap
        tmp_ptr = fg_previous;
        fg_previous = fg;
        fg = tmp_ptr;
        
        printf("t = %.2f\n", t);
    }

    // Save last time step
    save_foreground(fg_previous, grid_info, save_nr, t);
    
    deallocate_grid_info_struct(grid_info);
    deallocate_background_struct(bg);
    deallocate_foreground_struct(fg_previous);
    deallocate_foreground_struct(fg);

    return 0;
}