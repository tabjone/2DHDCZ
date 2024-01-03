#include "./one_time_step_3D/one_time_step_3D.h"
#include "data_structures/data_structures.h"
#include "global_parameters.h"
#include "functions_3D.h"

int main_3D(int argc, char *argv[], struct MpiInfo *mpi_info)
{
    /*
    This is a GOD function. It should be simplified and split into smaller functions.
    */
    FLOAT_P dt, dt_last, t_since_save;

    int save_nr;

    // Declare the background variables, foreground variables and grid info
    struct BackgroundVariables *bg = NULL;
    struct ForegroundVariables3D *fg = NULL, *fg_previous = NULL, *tmp_ptr = NULL;
    struct GridInfo3D *grid_info = NULL;
    struct PrecalculatedVariables3D *precalc = NULL;

    FLOAT_P t, first_t;
    #if LOAD == 0
        initialize_simulation_3D(&bg, &fg, &fg_previous, &grid_info, &precalc, mpi_info, &save_nr, &t, &first_t);
        
    #elif LOAD == 1
        load_simulation_3D();
        
    #endif // LOAD

    #if SAVE_RHS == 1
        save_rhs(fg_previous, bg, grid_info, mpi_info, precalc, 0);
        save_elliptic_vars(fg_previous, bg, grid_info, mpi_info, precalc, 0);
    #endif // SAVE_RHS
    t_since_save = 0.0;
    dt_last = 0.0;
    while (t < T)
    {
        dt = one_time_step_3D(bg, fg_previous, fg, grid_info, mpi_info, precalc, dt_last, first_t == t);
        t += dt;
        //break;
        
        

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
            save_foreground_3D(fg, grid_info, mpi_info, save_nr, t);
            #if SAVE_RHS == 1
                save_rhs(fg, bg, grid_info, mpi_info, precalc, save_nr);
                save_elliptic_vars(fg, bg, grid_info, mpi_info, precalc, save_nr);
            #endif // SAVE_RHS
            save_nr++;
            t_since_save = 0.0;
        }
        else if (SAVE_ALL == 1)
        {
            save_foreground_3D(fg, grid_info, mpi_info, save_nr, t);
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

    save_foreground_3D(fg_previous, grid_info, mpi_info, save_nr, t);
    #if SAVE_RHS == 1
        save_elliptic_vars(fg_previous, bg, grid_info, mpi_info, precalc, save_nr);
        save_rhs(fg_previous, bg, grid_info, mpi_info, precalc, save_nr);
    #endif // SAVE_RHS

    deallocate_grid_info_struct_3D(grid_info);
    deallocate_background_struct(bg);
    deallocate_precalculated_data_struct_3D(&precalc);
    deallocate_foreground_struct_3D(fg_previous);
    deallocate_foreground_struct_3D(fg);

    return 0;
}