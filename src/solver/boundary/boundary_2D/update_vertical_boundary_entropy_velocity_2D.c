#include "global_float_precision.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "MPI_module/MPI_module.h"
#include "curve_fitting/curve_fitting.h"

void update_vertical_boundary_entropy_velocity_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    communicate_2D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->vy, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->vz, mpi_info, nz, nz_ghost, ny);
    
    // All non periodic boundary conditions should extrapolate ghost cells
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        int nz_full = grid_info->nz_full;
        // If not periodic boundary extrapolate ghost cells at top and bottom
        if (!mpi_info->has_neighbor_below)
        {
            extrapolate_2D_array_antisymmetric_down(fg->vy, nz_ghost, ny);
            extrapolate_2D_array_constant_down(fg->vz, nz_ghost, ny);
            extrapolate_2D_array_constant_down(fg->s1, nz_ghost, ny);
        }
        if (!mpi_info->has_neighbor_above)
        {
            extrapolate_2D_array_antisymmetric_up(fg->vy, nz_full, nz_ghost, ny);
            extrapolate_2D_array_constant_up(fg->vz, nz_full, nz_ghost, ny);
            extrapolate_2D_array_constant_up(fg->s1, nz_full, nz_ghost, ny);
        }
    #endif // VERTICAL_BOUNDARY_TYPE
}