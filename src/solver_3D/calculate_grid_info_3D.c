#include "functions_3D.h"

void calculate_grid_info_3D(struct GridInfo3D **grid_info, struct MpiInfo *mpi_info)
{

    // Calculating the size of the grid
    FLOAT_P L_z = (R_END - R_START)*R_SUN;
    FLOAT_P L_y = Y_SIZE*R_SUN;
    FLOAT_P L_x = X_SIZE*R_SUN;

    // Calculating the size of the grid cells
    FLOAT_P dx = L_x/(NX - 1);
    FLOAT_P dy = L_y/(NY - 1);
    FLOAT_P dz = L_z/(NZ - 1);

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
    FLOAT_P x0 = 0.0;
    FLOAT_P x1 = R_SUN * X_SIZE;

    FLOAT_P z_offset = 0.0;


    allocate_grid_info_struct_3D(grid_info, NZ, nz_ghost, nz_full, NY, NX, dz, dy, dx, z0, z1, z_offset, y0, y1, x0, x1);
    mpi_info->size = 1;
}