#include "global_float_precision.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "MPI_module/MPI_module.h"
#include "curve_fitting/curve_fitting.h"

void update_vertical_boundary_entropy_velocity_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info)
{
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    int nx = grid_info->nx;

    communicate_3D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost, ny, nx);
    communicate_3D_ghost_above_below(fg->vx, mpi_info, nz, nz_ghost, ny, nx);
    communicate_3D_ghost_above_below(fg->vy, mpi_info, nz, nz_ghost, ny, nx);
    communicate_3D_ghost_above_below(fg->vz, mpi_info, nz, nz_ghost, ny, nx);
    
    // All non periodic boundary conditions should extrapolate ghost cells
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        int nz_full = grid_info->nz_full;
        // If not periodic boundary extrapolate ghost cells at top and bottom
        if (!mpi_info->has_neighbor_below)
        {
            extrapolate_3D_array_antisymmetric_down(fg->vy, nz_ghost, ny, nx);
            extrapolate_3D_array_antisymmetric_down(fg->vx, nz_ghost, ny, nx);
            extrapolate_3D_array_constant_down(fg->vz, nz_ghost, ny, nx);
            extrapolate_3D_array_constant_down(fg->s1, nz_ghost, ny, nx);
        }
        if (!mpi_info->has_neighbor_above)
        {
            extrapolate_3D_array_antisymmetric_up(fg->vy, nz_full, nz_ghost, ny, nx);
            extrapolate_3D_array_antisymmetric_up(fg->vx, nz_full, nz_ghost, ny, nx);
            extrapolate_3D_array_constant_up(fg->vz, nz_full, nz_ghost, ny, nx);
            extrapolate_3D_array_constant_up(fg->s1, nz_full, nz_ghost, ny, nx);
        }
    #endif // VERTICAL_BOUNDARY_TYPE
}