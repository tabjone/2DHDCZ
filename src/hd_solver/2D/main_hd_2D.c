#include "functions.h"

int main_hd_2D(int argc, char *argv[])
{
    /*
    double t, dt, t_since_save;
    int save_nr = 0;

    // Calculating the size of the grid
    double L_z = (R_END - R_START)*R_SUN;
    double L_x = X_SIZE*R_SUN;

    // Calculating the size of the grid cells
    double dx = L_x/(NX - 1);
    double dz = L_z/(NZ - 1);

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

    double z0 = R_SUN * R_START;
    double z1 = R_SUN * R_END;
    double x0 = 0.0;
    double x1 = R_SUN * X_SIZE;

    // Declare the background and foreground variables
    struct BackgroundVariables *bg;
    struct ForegroundVariables2D *fg, *fg_previous, *tmp_ptr;
    struct GridInfo *grid_info;

    allocate_grid_info_struct(&grid_info, NZ, nz_ghost, nz_full, NX, dz, dx, z0, z1, x0, x1);

    // Allocating memory for the background and foreground variables
    allocate_background_struct(&bg, nz_full);
    allocate_foreground_struct_2D(&fg, nz_full, NX);
    allocate_foreground_struct_2D(&fg_previous, nz_full, NX);

    // Initialize the background variables and saving it to file
    solar_s_background_initialization(bg, grid_info);
    save_background(bg, grid_info);
    
    // Initialize the foreground variables
    if (true)
    {
        initialize_foreground_struct_zeros(fg_previous, grid_info);
    }
    if (false)
    {
        initialize_foreground_struct_ones(fg_previous, grid_info);
    }
    if (false)
    {
        initialize_velocity_right(fg_previous, grid_info);
    }
    if (false)
    {
        initialize_foreground_struct_density_pertubation(fg_previous, grid_info);
    }

    
    // Saving the foreground variables to file
    save_foreground(fg_previous, grid_info, 0, 0.0);
    save_nr ++;
    
    t_since_save = 0.0;
    t = 0.0;
    while (t < T)
    {
        break;
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
        break;
    }

    // Save last time step
    save_foreground(fg_previous, grid_info, save_nr, t);
    
    // Deallocating memory
    deallocate_background_struct(bg);
    deallocate_foreground_struct_2D(fg);
    deallocate_foreground_struct_2D(fg_previous);
    deallocate_grid_info_struct(grid_info);
    return 0;

    */


    
   // Example of loading snap


    struct BackgroundVariables *bg;
    //struct ForegroundVariables2D *fg, *fg_previous, *tmp_ptr;
    struct ForegroundVariables2D *fg_previous;
    struct GridInfo *grid_info;

    load_grid_info(&grid_info, "./data/save_test/snap0.h5");
    allocate_background_struct(&bg, grid_info->nz_full);

    printf("grid_info->nz_full = %d\n", grid_info->nz_full);
    printf("grid_info->nx = %d\n", grid_info->nx);
    allocate_foreground_struct_2D(&fg_previous, grid_info->nz_full, grid_info->nx);
    
    
    load_background(bg, grid_info, "./data/save_test/background.h5");

    save_background(bg, grid_info);

    double time;
    
    time = load_foreground(fg_previous, grid_info, "./data/save_test/snap0.h5");

    printf("time = %f\n", time);


    save_foreground(fg_previous, grid_info, 11, time);


    deallocate_grid_info_struct(grid_info);
    deallocate_background_struct(bg);
    deallocate_foreground_struct_2D(fg_previous);

    return 0;
}