#include <stdbool.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "iterative_solver_module/iterative_solver_2D/iterative_solver_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"
#include "MPI_module/MPI_module.h"
#include "solver/rhs_functions/rhs_functions_2D/rhs_functions_2D.h"

void solve_elliptic_equation_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc)
{
    /*
    Solves the elliptic equation for the pressure field

    Parameters
    ----------
    bg : struct BackgroundVariables
        Background variables
    fg_prev : struct ForegroundVariables2D
        Foreground variables at previous time step
    fg : struct ForegroundVariables2D
        Foreground variables at current time step
    grid_info : struct GridInfo2D
        Grid parameters
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;

    FLOAT_P **rhs;
    allocate_2D_array(&rhs, nz, ny);

    // Solve inside the grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            rhs[i][j] = rhs_elliptic_eq_2D(bg, fg, grid_info, precalc, i+nz_ghost, j);
        }
    }

    iterative_solver_2D(rhs, fg->p1, fg_prev->p1, nz, nz_ghost, ny, dz, dy, mpi_info);

    deallocate_2D_array(rhs);
}