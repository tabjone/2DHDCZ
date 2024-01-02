#include "global_float_precision.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "MPI_module/MPI_module.h"
#include "curve_fitting/curve_fitting.h"

void update_vertical_boundary_pressure_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info)
{
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;

    communicate_1D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost);
    
    // All non periodic boundary conditions should extrapolate ghost cells
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        int nz_full = grid_info->nz_full;
        // If not periodic boundary extrapolate ghost cells at top and bottom
        if (!mpi_info->has_neighbor_below)
        {
            extrapolate_1D_array_constant_down(fg->p1, nz_ghost);
        }
        if (!mpi_info->has_neighbor_above)
        {
            extrapolate_1D_array_constant_up(fg->p1, nz_full, nz_ghost);
        }
    #endif // VERTICAL_BOUNDARY_TYPE
}