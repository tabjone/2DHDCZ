#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"
#include "MPI_module/mpi_info_struct.h"
#include <hdf5.h>
#include "spacial_derivatives_module/spacial_derivatives_module.h"
#include "array_utilities/array_memory_management/array_memory_management.h"

#include "global_parameters.h"
#include "global_float_precision.h"

void save_derivatives_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, struct MpiInfo *mpi_info, int snap_number)
{
    if (mpi_info->rank == 0)
    {
        printf("Saving derivatives for snapshot %d\n", snap_number);
    }

    // Getting grid info
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;

    // File path
    char file_path[150];

    snprintf(file_path, sizeof(file_path), "%s%s/rhs/rhs%d_%d.h5", SAVE_DIR, RUN_NAME, snap_number, mpi_info->rank);

    hsize_t dims[2] = {nz, ny};
    herr_t status;

    hid_t file, group_variables;
    hid_t dataspace_2d;

    
    // Create new default file 
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
    }

   // Create group for variables
    group_variables = H5Gcreate2(file, "/rhs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_variables < 0) {
        fprintf(stderr, "Failed to create variables group\n");
    }

    dataspace_2d = H5Screate_simple(2, dims, NULL);
    if (dataspace_2d < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }
    
    // Creating pointers to foreground arrays
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;
    FLOAT_P *one_over_rho0 = precalc->one_over_rho0;
    FLOAT_P *rho0 = bg->rho0;

    
    // Create arrays for derivatives
    FLOAT_P **pressure, **advection, **rhs, **viscosity;
    allocate_2D_array(&pressure, nz, ny);
    allocate_2D_array(&advection, nz, ny);
    allocate_2D_array(&rhs, nz, ny);
    allocate_2D_array(&viscosity, nz, ny);

    FLOAT_P dp1_dy, dvy_dy, dvy_dz;
    FLOAT_P dd_vy_ddy, dd_vz_dydz, dd_vy_ddz;

    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            dp1_dy = central_first_derivative_y_2D(p1, i, j, ny, precalc->one_over_2dy);
            dvy_dy = upwind_first_derivative_y_2D(vy, vy, i, j, ny, precalc->one_over_dy, precalc->one_over_2dy);
            dvy_dz = upwind_first_derivative_z_2D(vy, vz, i, j, precalc->one_over_dz, precalc->one_over_2dz);
            dd_vy_ddy = central_second_derivative_y_2D(vy, i, j, ny, precalc->one_over_dydy);
            dd_vz_dydz = central_second_derivative_yz_2D(vz, i, j, ny, precalc->one_over_4dydz);
            dd_vy_ddz = central_second_derivative_z_2D(vy, i, j, precalc->one_over_dzdz);

            pressure[i-nz_ghost][j] = -one_over_rho0[i]*dp1_dy;
            advection[i-nz_ghost][j] = -vy[i][j]*dvy_dy - vz[i][j]*dvy_dz;
            viscosity[i-nz_ghost][j] = precalc->VIS_COEFF_over_rho0[i]*
            (
                4.0/3.0*dd_vy_ddy + 1.0/3.0 * dd_vz_dydz + dd_vy_ddz
            );

            rhs[i-nz_ghost][j] = pressure[i-nz_ghost][j] + advection[i-nz_ghost][j] + viscosity[i-nz_ghost][j];
        }
    }
    
    create_write_dataset(group_variables, "pressure", H5_FLOAT_P, dataspace_2d, pressure[0], "");
    create_write_dataset(group_variables, "advection", H5_FLOAT_P, dataspace_2d, advection[0], "");
    create_write_dataset(group_variables, "rhs", H5_FLOAT_P, dataspace_2d, rhs[0], "");
    create_write_dataset(group_variables, "viscosity", H5_FLOAT_P, dataspace_2d, viscosity[0], "");

    status = H5Gclose(group_variables);
    if (status < 0) {
        fprintf(stderr, "Failed to close group\n");
    }
    status = H5Sclose(dataspace_2d);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
    
    deallocate_2D_array(pressure);
    deallocate_2D_array(advection);
    deallocate_2D_array(rhs);
    deallocate_2D_array(viscosity);
    
}