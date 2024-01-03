#include <stdbool.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "iterative_solver_module/iterative_solver_1D/iterative_solver_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_struct_1D.h"
#include "MPI_module/MPI_module.h"
#include "solver/rhs_functions/rhs_functions_1D/rhs_functions_1D.h"

void solve_elliptic_equation_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg_prev, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables1D *precalc)
{
    /*
    Solves the elliptic equation for the pressure field

    Parameters
    ----------
    bg : struct BackgroundVariables
        Background variables
    fg_prev : struct ForegroundVariables1D
        Foreground variables at previous time step
    fg : struct ForegroundVariables1D
        Foreground variables at current time step
    grid_info : struct GridInfo1D
        Grid parameters
    */

    // Getting grid info
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;

    FLOAT_P **rhs;
    allocate_1D_array(&rhs, nz);

    // Solve inside the grid
    for (int i = 0; i < nz; i++)
    {
        rhs[i][j] = rhs_elliptic_eq_1D(bg, fg, grid_info, precalc, i+nz_ghost);
    }

    iterative_solver_1D(rhs, fg->p1, fg_prev->p1, nz, nz_ghost, dz, mpi_info);

    deallocate_1D_array(rhs);
}