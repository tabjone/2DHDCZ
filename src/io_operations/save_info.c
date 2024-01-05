#include "MPI_module/mpi_info_struct.h"
#include "hdf5.h"
#include <stdio.h>
#include <string.h>
#include "global_float_precision.h"
#include "global_constants.h"
#include "global_parameters.h"
#include "io_operations/general_io/general_io.h"

void save_simulation_info(struct MpiInfo *mpi_info)
{
    if (mpi_info->rank == 0)
    {
        const char* root_header = "This is the root header";
        const char* constants_header = "This is the constants header";
        const char* parameters_header = "This is the parameters header";

        char file_path[150];
        snprintf(file_path, sizeof(file_path), "%s%s/info.h5", SAVE_DIR, RUN_NAME);

        hid_t file;
        hid_t group_constants, group_parameters;
        hid_t dataspace_constants, dataspace_parameters;

        herr_t status;

        // Create new default file 
        file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file < 0) {
            fprintf(stderr, "Failed to create file\n");
        }

        // Adding root header
        add_string_attribute(file, "header", root_header);
        
        // Create group for grid_data
        group_parameters = H5Gcreate2(file, "/global_parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (group_parameters < 0) {
            fprintf(stderr, "Failed to create global_parameters group\n");
        }

        add_string_attribute(group_parameters, "header", parameters_header);

        dataspace_parameters = H5Screate(H5S_SCALAR);
        if (dataspace_parameters < 0) {
            fprintf(stderr, "Failed to create dataspace\n");
        }

        FLOAT_P max_dt = MAX_DT;
        FLOAT_P save_interval = SAVE_INTERVAL;
        FLOAT_P cfl_cut = CFL_CUT;
        int upwind_order = UPWIND_ORDER;
        int central_order = CENTRAL_ORDER;
        int time_order = TIME_ORDER;
        int float_precision = FLOAT_PRECISION;
        int units = UNITS;
        FLOAT_P cz_start = CZ_START;
        FLOAT_P r_start = R_START;
        FLOAT_P r_end = R_END;
        FLOAT_P x_size = X_SIZE;
        FLOAT_P y_size = Y_SIZE;
        int nx = NX;
        int ny = NY;
        int nz = NZ;
        FLOAT_P iter_tol = ITERATIVE_SOLVER_TOLERANCE;
        int iter_max_iter = ITERATIVE_SOLVER_MAX_ITERATIONS;
        FLOAT_P del_ad = DEL_AD;
        FLOAT_P gamma = GAMMA;
        FLOAT_P nabla_ad = NABLA_AD;
        FLOAT_P superad_param = SUPERAD_PARAM;
        int gravity_on = GRAVITY_ON;
        int advection_on = ADVECTION_ON;
        int gas_pressure_on = GAS_PRESSURE_ON;
        int viscosity_on = VISCOSITY_ON;
        int thermal_diffusivity_on = THERMAL_DIFFUSIVITY_ON;
        int bfield_on = BFIELD_ON;

        // Adding constants
        create_write_dataset(group_parameters, "MAX_DT", H5_FLOAT_P, dataspace_parameters, &max_dt, "s");
        create_write_dataset(group_parameters, "SAVE_INTERVAL", H5_FLOAT_P, dataspace_parameters, &save_interval, "s");
        create_write_dataset(group_parameters, "CFL_CUT", H5_FLOAT_P, dataspace_parameters, &cfl_cut, "");
        create_write_dataset(group_parameters, "UPWIND_ORDER", H5T_NATIVE_INT, dataspace_parameters, &upwind_order, "Order of upwind scheme");
        create_write_dataset(group_parameters, "CENTRAL_ORDER", H5T_NATIVE_INT, dataspace_parameters, &central_order, "Order of central scheme");
        create_write_dataset(group_parameters, "TIME_ORDER", H5T_NATIVE_INT, dataspace_parameters, &time_order, "Order of temporal scheme");
        create_write_dataset(group_parameters, "FLOAT_PRECISION", H5T_NATIVE_INT, dataspace_parameters, &float_precision, "0 for float, 1 for double, 2 for long double");
        create_write_dataset(group_parameters, "UNITS", H5T_NATIVE_INT, dataspace_parameters, &units, "0 for cgs, 1 for SI");
        create_write_dataset(group_parameters, "CZ_START", H5_FLOAT_P, dataspace_parameters, &cz_start, "Beginning of CZ in solar radii");
        create_write_dataset(group_parameters, "R_START", H5_FLOAT_P, dataspace_parameters, &r_start, "Beginning of domain in z-direction in solar radii");
        create_write_dataset(group_parameters, "R_END", H5_FLOAT_P, dataspace_parameters, &r_end, "End of domain in z-direction in solar radii");
        create_write_dataset(group_parameters, "X_SIZE", H5_FLOAT_P, dataspace_parameters, &x_size, "Size of domain in x-direction in solar radii");
        create_write_dataset(group_parameters, "Y_SIZE", H5_FLOAT_P, dataspace_parameters, &y_size, "Size of domain in y-direction in solar radii");
        create_write_dataset(group_parameters, "NX", H5T_NATIVE_INT, dataspace_parameters, &nx, "Number of grid points in x-direction");
        create_write_dataset(group_parameters, "NY", H5T_NATIVE_INT, dataspace_parameters, &ny, "Number of grid points in y-direction");
        create_write_dataset(group_parameters, "NZ", H5T_NATIVE_INT, dataspace_parameters, &nz, "Number of grid points in z-direction");
        create_write_dataset(group_parameters, "ITER_TOL", H5_FLOAT_P, dataspace_parameters, &iter_tol, "Gauss-Seidel tolerance");
        create_write_dataset(group_parameters, "ITER_MAX_ITER", H5T_NATIVE_INT, dataspace_parameters, &iter_max_iter, "Gauss-Seidel max iterations");
        create_write_dataset(group_parameters, "DEL_AD", H5_FLOAT_P, dataspace_parameters, &del_ad, "Adiabatic temperature gradient");
        create_write_dataset(group_parameters, "GAMMA", H5_FLOAT_P, dataspace_parameters, &gamma, "Superadiabatic parameter");
        create_write_dataset(group_parameters, "NABLA_AD", H5_FLOAT_P, dataspace_parameters, &nabla_ad, "Adiabatic temperature gradient");
        create_write_dataset(group_parameters, "SUPERAD_PARAM", H5_FLOAT_P, dataspace_parameters, &superad_param, "Superadiabacicity parameter in CZ");
        create_write_dataset(group_parameters, "GRAVITY_ON", H5T_NATIVE_INT, dataspace_parameters, &gravity_on, "0 for gravity off, 1 for gravity on");
        create_write_dataset(group_parameters, "ADVECTION_ON", H5T_NATIVE_INT, dataspace_parameters, &advection_on, "0 for advection off, 1 for advection on");
        create_write_dataset(group_parameters, "GAS_PRESSURE_ON", H5T_NATIVE_INT, dataspace_parameters, &gas_pressure_on, "0 for gas pressure off, 1 for gas pressure on");
        create_write_dataset(group_parameters, "VISCOSITY_ON", H5T_NATIVE_INT, dataspace_parameters, &viscosity_on, "0 for viscosity off, 1 for viscosity on");
        create_write_dataset(group_parameters, "THERMAL_DIFFUSIVITY_ON", H5T_NATIVE_INT, dataspace_parameters, &thermal_diffusivity_on, "0 for thermal diffusivity off, 1 for thermal diffusivity on");
        create_write_dataset(group_parameters, "BFIELD_ON", H5T_NATIVE_INT, dataspace_parameters, &bfield_on, "0 for B-field off, 1 for B-field on");

        status = H5Sclose(dataspace_parameters);
        if (status < 0) {
            fprintf(stderr, "Failed to close dataspace\n");
        }

        status = H5Gclose(group_parameters);
        if (status < 0) {
            fprintf(stderr, "Failed to close group\n");
        }

        // Create group for constants
        group_constants = H5Gcreate2(file, "/constants", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (group_constants < 0) {
            fprintf(stderr, "Failed to create constants group\n");
        }

        add_string_attribute(group_constants, "header", constants_header);

        dataspace_constants = H5Screate(H5S_SCALAR);
        if (dataspace_constants < 0) {
            fprintf(stderr, "Failed to create dataspace\n");
        }

        FLOAT_P eta = ETA;
        FLOAT_P mu = MU;
        FLOAT_P k_b = K_B;
        FLOAT_P m_u = M_U;
        FLOAT_P r_sun = R_SUN;
        FLOAT_P m_sun = M_SUN;
        FLOAT_P viscosity_coeff = VISCOSITY_COEFF;
        FLOAT_P G_ = G;


        // Adding constants
        create_write_dataset(group_constants, "ETA", H5_FLOAT_P, dataspace_constants, &eta, "magnetic diffusivity [cm2 s-1]");
        create_write_dataset(group_constants, "MU", H5_FLOAT_P, dataspace_constants, &mu, "mean molecular weight");
        create_write_dataset(group_constants, "K_B", H5_FLOAT_P, dataspace_constants, &k_b, "boltzmann constant [cm2 g s-2 K-1]");
        create_write_dataset(group_constants, "M_U", H5_FLOAT_P, dataspace_constants, &m_u, "atomic mass constant [g]");
        create_write_dataset(group_constants, "R_SUN", H5_FLOAT_P, dataspace_constants, &r_sun, "solar radius [cm]");
        create_write_dataset(group_constants, "M_SUN", H5_FLOAT_P, dataspace_constants, &m_sun, "solar mass [g]");
        create_write_dataset(group_constants, "VISCOSITY_COEFF", H5_FLOAT_P, dataspace_constants, &viscosity_coeff, "viscosity coefficient [P] (Poise)");
        create_write_dataset(group_constants, "G", H5_FLOAT_P, dataspace_constants, &G_, "gravitational constant [cm3 g-1 s-2]");

        status = H5Sclose(dataspace_constants);
        if (status < 0) {
            fprintf(stderr, "Failed to close dataspace\n");
        }

        status = H5Gclose(group_constants);
        if (status < 0) {
            fprintf(stderr, "Failed to close group\n");
        }

        // Close file
        status = H5Fclose(file);
        if (status < 0) {
            fprintf(stderr, "Failed to close file\n");
        }
    }
}