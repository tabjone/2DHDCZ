#include "io_functions.h"
#include "../spacial_derivatives/spacial_derivatives.h"

void save_elliptic_vars(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc, int save_nr)
{
    /*
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    FLOAT_P **dvy_dy, **dvz_dz, **dd_vy_ddy, **dd_vz_ddz, **dd_vy_dydz, **dd_vz_dydz, **dvy_dy_sqrd, **dvz_dz_sqrd;
    FLOAT_P **dvy_dz, **dvz_dy;
    FLOAT_P **first_term, **second_term;
    FLOAT_P **rhs, **rhs_manual;
    FLOAT_P **third_term;
    
    allocate_2D_array(&dvy_dy, nz, ny);
    allocate_2D_array(&dvz_dz, nz, ny);
    allocate_2D_array(&dd_vy_ddy, nz, ny);
    allocate_2D_array(&dd_vz_ddz, nz, ny);
    allocate_2D_array(&dd_vy_dydz, nz, ny);
    allocate_2D_array(&dd_vz_dydz, nz, ny);
    allocate_2D_array(&dvy_dy_sqrd, nz, ny);
    allocate_2D_array(&dvz_dz_sqrd, nz, ny);
    allocate_2D_array(&first_term, nz, ny);
    allocate_2D_array(&second_term, nz, ny);
    allocate_2D_array(&dvy_dz, nz, ny);
    allocate_2D_array(&dvz_dy, nz, ny);
    allocate_2D_array(&rhs, nz, ny);
    allocate_2D_array(&rhs_manual, nz, ny);
    allocate_2D_array(&third_term, nz, ny);

    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;
    FLOAT_P *grad_rho0 = precalc->grad_rho0;

    FLOAT_P drho1_dz;

    int i_;
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            i_ = i + nz_ghost;
            drho1_dz = central_first_derivative_z(rho1, precalc, i_, j);

            dvy_dy[i][j]      = central_first_derivative_y(vy, precalc, i_, j);
            dvz_dz[i][j]      = central_first_derivative_z(vz, precalc, i_, j);
            dd_vy_ddy[i][j]   = central_second_derivative_y(vy, precalc, i_, j);
            dd_vz_ddz[i][j]   = central_second_derivative_z(vz, precalc, i_, j);
            dd_vy_dydz[i][j]  = central_second_derivative_yz(vy, precalc, i_, j);
            dd_vz_dydz[i][j]  = central_second_derivative_yz(vz, precalc, i_, j);
            dvy_dz[i][j]      = central_first_derivative_z(vy, precalc, i_, j);
            dvz_dy[i][j]      = central_first_derivative_y(vz, precalc, i_, j);
            dvy_dy_sqrd[i][j] = dvy_dy[i][j]*dvy_dy[i][j];
            dvz_dz_sqrd[i][j] = dvz_dz[i][j]*dvz_dz[i][j];

            first_term[i][j]  = -rho0[i_] * ( vy[i_][j]*dd_vy_ddy[i][j] + vz[i_][j]*dd_vz_ddz[i][j] + dvy_dy_sqrd[i][j] + dvz_dz_sqrd[i][j]
                          +2*dvz_dy[i][j]*dvy_dz[i][j] + vy[i_][j]*dd_vz_dydz[i][j] + vz[i_][j]*dd_vy_dydz[i][j] );
            second_term[i][j] = -grad_rho0[i_] * (vy[i_][j]*dvz_dy[i][j] + vz[i_][j]*dvz_dz[i][j]);

            third_term[i][j] = -g[i_] * drho1_dz - precalc->grad_g[i_]*rho1[i_][j];

            rhs_manual[i][j] = first_term[i][j] + second_term[i][j] + third_term[i][j];

            rhs[i][j] = rhs_elliptic_eq_2D(bg, fg, grid_info, precalc, i_, j);
        }
    }

    // File path
    char file_path[150];

    #if MPI_ON == 0
        snprintf(file_path, sizeof(file_path), "%s%s/elliptic_vars/elliptic_vars%d.h5", SAVE_DIR, RUN_NAME, save_nr);
    #elif MPI_ON == 1
        snprintf(file_path, sizeof(file_path), "%s%s/elliptic_vars/elliptic_vars%d_%d.h5", SAVE_DIR, RUN_NAME, save_nr, mpi_info->rank);
    #endif // MPI_ON

    // Create file with no groups
    hsize_t dims[2] = {nz, ny};
    hid_t file, dataspace_2d;
    herr_t status;

    // Create new default file
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
    }

    // Create dataspace for 2D arrays
    dataspace_2d = H5Screate_simple(2, dims, NULL);
    if (dataspace_2d < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }

    // Create datasets
    create_write_dataset(file, "dvy_dy", H5_FLOAT_P, dataspace_2d, dvy_dy[0], "");
    create_write_dataset(file, "dvz_dz", H5_FLOAT_P, dataspace_2d, dvz_dz[0], "");
    create_write_dataset(file, "dvz_dy", H5_FLOAT_P, dataspace_2d, dvz_dy[0], "");
    create_write_dataset(file, "dvy_dz", H5_FLOAT_P, dataspace_2d, dvy_dz[0], "");
    create_write_dataset(file, "dd_vy_ddy", H5_FLOAT_P, dataspace_2d, dd_vy_ddy[0], "");
    create_write_dataset(file, "dd_vz_ddz", H5_FLOAT_P, dataspace_2d, dd_vz_ddz[0], "");
    create_write_dataset(file, "dd_vy_dydz", H5_FLOAT_P, dataspace_2d, dd_vy_dydz[0], "");
    create_write_dataset(file, "dd_vz_dydz", H5_FLOAT_P, dataspace_2d, dd_vz_dydz[0], "");
    create_write_dataset(file, "dvy_dy_sqrd", H5_FLOAT_P, dataspace_2d, dvy_dy_sqrd[0], "");
    create_write_dataset(file, "dvz_dz_sqrd", H5_FLOAT_P, dataspace_2d, dvz_dz_sqrd[0], "");
    create_write_dataset(file, "first_term", H5_FLOAT_P, dataspace_2d, first_term[0], "");
    create_write_dataset(file, "second_term", H5_FLOAT_P, dataspace_2d, second_term[0], "");
    create_write_dataset(file, "third_term", H5_FLOAT_P, dataspace_2d, third_term[0], "");
    create_write_dataset(file, "rhs", H5_FLOAT_P, dataspace_2d, rhs[0], "");
    create_write_dataset(file, "rhs_manual", H5_FLOAT_P, dataspace_2d, rhs_manual[0], "");

    // Close dataspace
    status = H5Sclose(dataspace_2d);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    // Close file
    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }

    // Deallocate arrays
    deallocate_2D_array(dvz_dy);
    deallocate_2D_array(dvy_dz);
    deallocate_2D_array(dvy_dy);
    deallocate_2D_array(dvz_dz);
    deallocate_2D_array(dd_vy_ddy);
    deallocate_2D_array(dd_vz_ddz);
    deallocate_2D_array(dd_vy_dydz);
    deallocate_2D_array(dd_vz_dydz);
    deallocate_2D_array(dvy_dy_sqrd);
    deallocate_2D_array(dvz_dz_sqrd);
    deallocate_2D_array(first_term);
    deallocate_2D_array(second_term);
    deallocate_2D_array(rhs);
    deallocate_2D_array(rhs_manual);
    deallocate_2D_array(third_term);
    */
}