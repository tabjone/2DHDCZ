#include "solver/one_time_step/one_time_step_2D/one_time_step_2D.h"
#include "data_structures/data_structures.h"
#include "functions_2D.h"
#include "global_float_precision.h"
#include "io_operations/io_operations.h"
#include <stdlib.h>

#include "global_constants.h"
#include "global_parameters.h"
#include "global_boundary.h"
#include <math.h>

int main_2D(int argc, char *argv[], struct MpiInfo *mpi_info)
{
    // Initialize random number generator
    srand(0);

    FLOAT_P dt, dt_last, t_since_save;

    int save_nr;

    // Declare the background variables, foreground variables and grid info
    struct BackgroundVariables *bg = NULL;
    struct ForegroundVariables2D *fg = NULL, *fg_previous = NULL, *tmp_ptr = NULL;
    struct GridInfo2D *grid_info = NULL;
    struct PrecalculatedVariables2D *precalc = NULL;

    FLOAT_P t, first_t;
    #if LOAD == 0
        initialize_simulation_2D(&bg, &fg, &fg_previous, &grid_info, &precalc, mpi_info, &save_nr, &t, &first_t);
        
    #elif LOAD == 1
        load_simulation_2D(&bg, &fg, &fg_previous, &grid_info, &precalc, mpi_info, &save_nr, &t, &first_t);
    #endif // LOAD

    t_since_save = 0.0;
    dt_last = 0.0;

    while (t < T)
    {
        dt = one_time_step_2D(bg, fg_previous, fg, grid_info, mpi_info, precalc, dt_last, first_t == t);
        t += dt; 

        if (mpi_info->rank == 0)
            printf("t = %.3f, dt=%.3f\n", t, dt);
        
        t_since_save += dt;
        dt_last = dt;
        
        if (dt_last < 0.005)
        {
            break;
        }
        
        if (t_since_save > SAVE_INTERVAL && SAVE_ALL == 0)
        {
            save_foreground_2D(fg, grid_info, mpi_info, save_nr, t);
            save_nr++;
            t_since_save = 0.0;
        }
        else if (SAVE_ALL == 1)
        {
            save_foreground_2D(fg, grid_info, mpi_info, save_nr, t);
            save_nr++;
        }

        // pointer swap
        tmp_ptr = fg_previous;
        fg_previous = fg;
        fg = tmp_ptr;
    }

    save_foreground_2D(fg_previous, grid_info, mpi_info, save_nr, t);

    deallocate_grid_info_struct_2D(grid_info);
    deallocate_background_struct(bg);
    deallocate_precalculated_data_struct_2D(&precalc);
    deallocate_foreground_struct_2D(fg_previous);
    deallocate_foreground_struct_2D(fg);

    return 0;
}