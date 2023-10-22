#include "functions.h"

void calculate_grid_info(struct GridInfo *grid_info, struct MpiInfo *mpi_info)
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

    // Calculating the number of cells in the full grid
    int nz_full = NZ + 2*nz_ghost;
    FLOAT_P z0 = R_SUN * R_START;
    FLOAT_P z1 = R_SUN * R_END;
    FLOAT_P y0 = 0.0;
    FLOAT_P y1 = R_SUN * Y_SIZE;

    FLOAT_P z_offset = 0.0;

    allocate_grid_info_struct(grid_info, NZ, nz_ghost, nz_full, NY, dz, dy, z0, z1, z_offset, y0, y1);
    mpi_info->size = 1;
}