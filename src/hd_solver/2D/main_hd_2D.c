#include "functions.h"

int main_hd_2D(int argc, char *argv[])
{
    
    FLOAT_P t, dt, t_since_save;
    int save_nr;

    // Declare the background variables, foreground variables and grid info
    struct BackgroundVariables *bg;
    struct ForegroundVariables2D *fg, *fg_previous, *tmp_ptr;
    struct GridInfo *grid_info;

    #if LOAD == 1
        // Loading snapshot

        char file_path_foreground[150];
        char file_path_background[150];
        snprintf(file_path_foreground, sizeof(file_path_foreground), "data/%s/snap%d.h5", RUN_NAME, LOAD_SNAP_NUMBER);
        snprintf(file_path_background, sizeof(file_path_background), "data/%s/background.h5", RUN_NAME);

        load_grid_info(&grid_info, file_path_foreground);
        allocate_background_struct(&bg, grid_info->nz_full);

        printf("grid_info->nz_full = %d\n", grid_info->nz_full);
        printf("grid_info->nx = %d\n", grid_info->nx);
        allocate_foreground_struct_2D(&fg_previous, grid_info->nz_full, grid_info->nx);
        allocate_foreground_struct_2D(&fg, grid_info->nz_full, grid_info->nx);
        
        
        load_background(bg, grid_info, file_path_background);

        
        
        t = load_foreground(fg_previous, grid_info, file_path_foreground);

        printf("time = %f\n", t);

        //save_background(bg, grid_info);
        //save_foreground(fg_previous, grid_info, 12, t);
    

    #elif LOAD == 0
        // BLA BLA BLA
        save_nr = 0;

        // Calculating the size of the grid
        FLOAT_P L_z = (R_END - R_START)*R_SUN;
        FLOAT_P L_x = X_SIZE*R_SUN;

        // Calculating the size of the grid cells
        FLOAT_P dx = L_x/(NX - 1);
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
        FLOAT_P x0 = 0.0;
        FLOAT_P x1 = R_SUN * X_SIZE;

        allocate_grid_info_struct(&grid_info, NZ, nz_ghost, nz_full, NX, dz, dx, z0, z1, x0, x1);

        // Allocating memory for the background and foreground variables
        allocate_background_struct(&bg, nz_full);
        allocate_foreground_struct_2D(&fg, nz_full, NX);
        allocate_foreground_struct_2D(&fg_previous, nz_full, NX);

        // Initialize the background variables and saving it to file
        solar_s_background_initialization(bg, grid_info);
        save_background(bg, grid_info);
        
        // Initialize the foreground variables
        if (false)
        {
            initialize_foreground_struct_zeros(fg_previous, grid_info);
        }
        for (int i = grid_info->nz_ghost; i < grid_info->nz_full-grid_info->nz_ghost; i++)
        {
            for (int j = 0; j < grid_info->nx; j++)
            {
                fg_previous->vz[i][j] = 10.0;
                fg_previous->vx[i][j] = 0.0;
            }
        }
        //solve_elliptic_equation(bg, fg_previous, fg_previous, grid_info);


        if (false)
        {
            initialize_foreground_struct_random(fg_previous, bg, grid_info);
        }

        if (false)
        {
            initialize_foreground_struct_ones(fg_previous, grid_info);
        }
        if (false)
        {
            initialize_velocity_right(fg_previous, grid_info);
        }
        if (true)
        {
            initialize_foreground_struct_density_pertubation(fg_previous, bg, grid_info);
        }

        
        // Saving the foreground variables to file
        save_foreground(fg_previous, grid_info, 0, 0.0);
        save_nr ++;
        
        t_since_save = 0.0;
        t = 0.0;
        while (t < T)
        {
            dt = one_time_step(bg, fg_previous, fg, grid_info, t==0);
            t += dt;
            t_since_save += dt;
            
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
        
    #endif // LOAD
        deallocate_grid_info_struct(grid_info);
        deallocate_background_struct(bg);
        deallocate_foreground_struct_2D(fg_previous);
        deallocate_foreground_struct_2D(fg);

    return 0;
}