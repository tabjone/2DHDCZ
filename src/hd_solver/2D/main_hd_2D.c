#include "functions.h"

int main_hd_2D(int argc, char *argv[])
{
    printf("Running 2D HD solver\n");
    int nx = x_size/dx - 1; // Number of grid points in x-direction
    int nz = (R_END - R_START)/dz - 1; // Number of grid points in z-direction 

    struct BackgroundVariables *bg;
    struct ForegroundVariables2D *fg, *fg_previous, *tmp_ptr;

    allocate_background_struct(nz, &bg);
    allocate_foreground_struct_2D(nz, nx, &fg);
    allocate_foreground_struct_2D(nz, nx, &fg_previous);

    solar_s_background_initialization(bg);
    save_background(bg);
    

    initialize_foreground_struct_zeros(fg);
    // Copy foreground_variables to fg_previous
    deep_copy_foreground_2D(fg, fg_previous);

    double T = 1.0;
    double dt = 0.1;

    save_foreground(fg_previous, 0);
    double t = 0.0;
    while (t < T)
    {
        one_time_step(bg, fg_previous, fg);
        // Then extrapolate again

        // pointer swap
        tmp_ptr = fg_previous;
        fg_previous = fg;
        fg = tmp_ptr;
        t += dt;
    }
    save_foreground(fg_previous, 1);

    deallocate_background_struct(bg);
    deallocate_foreground_struct_2D(fg);
    deallocate_foreground_struct_2D(fg_previous);
    return 0;
}