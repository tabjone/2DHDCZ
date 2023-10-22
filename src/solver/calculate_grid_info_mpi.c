#include "functions.h"
#include <mpi.h>

void calculate_grid_info_mpi(struct MpiInfo *mpi_info, struct GridInfo *grid_info)
{
    // Calculating the size of the grid
    FLOAT_P L_z = (R_END - R_START)*R_SUN;
    FLOAT_P L_y = Y_SIZE*R_SUN;

    // Calculating the size of the grid cells
    FLOAT_P dy = L_y/(NY - 1);
    FLOAT_P dz = L_z/(NZ - 1);

    // Calculating the number of ghost cells
    int nz_ghost;
    if (UPWIND_ORDER >= CENTRAL_ORDER)
    {
        nz_ghost = UPWIND_ORDER;
    }
    else
    {
        nz_ghost = CENTRAL_ORDER;
    }

    // Number of cells to be divided among processes
    int total_cells = NZ;

    // Average number of cells per process
    int avg_cells_per_process = total_cells / mpi_info->size;

    // Extra cells that can't be evenly distributed among processes
    int extra_cells = total_cells % mpi_info->size;

    int my_offset, my_nz;

    // If current rank is less than extra_cells, then it takes one of the extra cells
    if (mpi_info->rank < extra_cells) {
        my_offset = mpi_info->rank * (avg_cells_per_process + 1);
        my_nz = avg_cells_per_process + 1;
    } else {
        my_offset = mpi_info->rank * avg_cells_per_process + extra_cells;
        my_nz = avg_cells_per_process;
    }

    FLOAT_P z_offset = (my_offset) * dz;

    FLOAT_P z0, z1, y0, y1;

    

    z0 = R_SUN * R_START + (my_offset) * dz;
    z1 = R_SUN * R_START + (my_offset + my_nz-1) * dz;
    y0 = 0.0;
    y1 = R_SUN * Y_SIZE;

    int nz_full = my_nz + 2*nz_ghost;

    allocate_grid_info_struct(grid_info, my_nz, nz_ghost, nz_full, NY, dz, dy, z0, z1, z_offset, y0, y1);
}