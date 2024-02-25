#include "global_float_precision.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "MPI_module/MPI_module.h"
#include "curve_fitting/curve_fitting.h"

#include "array_utilities/array_memory_management/array_memory_management.h"
#include "solver/boundary/boundary_2D/boundary_2D.h"

void update_vertical_boundary_pressure_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    communicate_2D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost, ny);

    int nz_full = grid_info->nz_full;
    /*
    FLOAT_P *damping_factor;
    allocate_1D_array(&damping_factor, nz_full);
    calculate_damping_2D(damping_factor, bg, grid_info, mpi_info);

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = damping_factor[i] * fg->p1[i][j];
        }
    }

    deallocate_1D_array(damping_factor);
    */
    
    // All non periodic boundary conditions should extrapolate ghost cells
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        //int nz_full = grid_info->nz_full;
        // If not periodic boundary extrapolate ghost cells at top and bottom
        if (!mpi_info->has_neighbor_below)
        {
            extrapolate_2D_array_symmetric_down(fg->p1, nz_ghost, ny);
        }
        if (!mpi_info->has_neighbor_above)
        {
            extrapolate_2D_array_symmetric_up(fg->p1, nz_full, nz_ghost, ny);
        }
    #endif // VERTICAL_BOUNDARY_TYPE
}