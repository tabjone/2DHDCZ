#include <stdbool.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "iterative_solver_module/iterative_solver_3D/iterative_solver_3D.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/precalculated_data/precalculated_data_3D/precalculated_data_struct_3D.h"
#include "MPI_module/MPI_module.h"
#include "solver/rhs_functions/rhs_functions_3D/rhs_functions_3D.h"

void solve_elliptic_equation_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables3D *precalc)
{
    /*
    Solves the elliptic equation for the pressure field

    Parameters
    ----------
    bg : struct BackgroundVariables
        Background variables
    fg_prev : struct ForegroundVariables3D
        Foreground variables at previous time step
    fg : struct ForegroundVariables3D
        Foreground variables at current time step
    grid_info : struct GridInfo3D
        Grid parameters
    */

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dx = grid_info->dx;

    FLOAT_P ***rhs;
    allocate_3D_array(&rhs, nz, ny, nx);

    // Solve inside the grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                rhs[i][j][k] = rhs_elliptic_eq_3D(bg, fg, grid_info, precalc, i+nz_ghost, j, k);
            }
        }
    }

    iterative_solver_3D(rhs, fg->p1, fg_prev->p1, nz, nz_ghost, ny, nx, dz, dy, dx, mpi_info);

    deallocate_3D_array(rhs);
}