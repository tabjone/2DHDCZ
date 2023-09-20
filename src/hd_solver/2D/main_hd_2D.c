#include "functions.h"

int main_hd_2D(int argc, char *argv[])
{
    // Fix this so that we instead chose nx, nz
    double L_z = (R_END - R_START)*R_SUN;
    double L_x = X_SIZE*R_SUN;

    printf("L_z = %.10f\n", L_z/R_SUN);
    printf("L_x = %.10f\n", L_x/R_SUN);

    #if RESOLUTION == 0
        printf("Error: Resolution must be greater than 0\n");
    #endif

    #if UPWIND_ORDER == 2
        double dx = 3*L_x/M_PI/RESOLUTION;
        double dz = 3*L_z/M_PI/RESOLUTION;
        //double dx = 100;
        //double dz = 100;


    #endif
 
    int nx = L_x/dx + 1.5; // Number of grid points in x-direction
    int nz = L_z/dz + 1.5; // Number of grid points in z-direction 

    printf("nx*dx = %.10f\n", nx*dx/R_SUN);
    printf("nz*dz = %.10f\n", nz*dz/R_SUN);

    struct BackgroundVariables *bg;
    struct ForegroundVariables2D *fg, *fg_previous, *tmp_ptr;

    allocate_background_struct(nz, dz, &bg);
    allocate_foreground_struct_2D(nz, nx, dz, dx, &fg);
    allocate_foreground_struct_2D(nz, nx, dz, dx, &fg_previous);

    solar_s_background_initialization(bg);
    save_background(bg);
    

    //initialize_foreground_struct_zeros(fg_previous);
    initialize_foreground_struct_pertubation(fg_previous);
    //initialize_foreground_struct_zeros(fg_previous);
    // Copy foreground_variables to fg_previous

    // SOME ERROR HERE
    //deep_copy_foreground_2D(fg_previous, fg);

    double dt = 0.1;

    save_foreground(fg_previous, 0, 0.0);
    double t = 0.0;

    int save_nr = 1;
    
    double epsilon = 1e-6; // Small number for comparing doubles
    
    while (t < T)
    {
        one_time_step(bg, fg_previous, fg, dt);
        t += dt;
        
        if (fabs(fmod(t, SAVE_INTERVAL)) < epsilon)
        {
        save_foreground(fg, save_nr, t);
        save_nr++;
        }
        

        // pointer swap
        tmp_ptr = fg_previous;
        fg_previous = fg;
        fg = tmp_ptr;
        
        printf("t = %.1f\n", t);
        

    }
    save_foreground(fg_previous, save_nr, t);
    
    deallocate_background_struct(bg);
    deallocate_foreground_struct_2D(fg);
    deallocate_foreground_struct_2D(fg_previous);
    return 0;
}