#include "global_float_precision.h"
#include "io_operations/background_io/background_io.h"
#include "global_parameters.h"
#include "global_constants.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "curve_fitting/interpolation/interpolation.h"
#include <math.h>
#include "solar_s_initialization.h"

void get_initial_values_from_solar_s(FLOAT_P *p0_initial, FLOAT_P *T0_initial, FLOAT_P *rho0_initial, FLOAT_P *m_initial, FLOAT_P r_integration_start)
{
    // Allocating arrays for solar_s data
    int solar_s_size = 2482;
    FLOAT_P *r_over_R_solar_s, *rho_solar_s, *p_solar_s, *T_solar_s;

    allocate_1D_array(&r_over_R_solar_s, solar_s_size);
    allocate_1D_array(&rho_solar_s, solar_s_size);
    allocate_1D_array(&p_solar_s, solar_s_size);
    allocate_1D_array(&T_solar_s, solar_s_size);

    // Reading solar_s data
    read_solar_s_data("additional_files/solar_s.h5", r_over_R_solar_s, rho_solar_s, p_solar_s, T_solar_s, solar_s_size);

    // Since the solar_s datapoints most likely will not match the grid points exactly, we need to interpolate to get the exact values at the grid point where we start the integration

    // Finding the index of the first datapoint that is inside the integration domain
    int integration_start_index = 0;
    while (r_over_R_solar_s[integration_start_index] > r_integration_start/R_SUN)
    {
        integration_start_index++;
    }

    // Getting values before and after integration start
    FLOAT_P r_before = r_over_R_solar_s[integration_start_index]  * R_SUN; // inside domain
    FLOAT_P r_after = r_over_R_solar_s[integration_start_index-1] * R_SUN; // outside domain
    FLOAT_P p0_before = p_solar_s[integration_start_index];
    FLOAT_P p0_after = p_solar_s[integration_start_index-1];
    FLOAT_P T0_before = T_solar_s[integration_start_index];
    FLOAT_P T0_after = T_solar_s[integration_start_index-1];

    // Interpolating to get the exact value at integration_start
    FLOAT_P r_initial = r_integration_start;

    FLOAT_P p0 = interpolate_1D_linear(r_before, r_after, p0_before, p0_after, r_initial);
    FLOAT_P T0 = interpolate_1D_linear(r_before, r_after, T0_before, T0_after, r_initial);
    FLOAT_P rho0 = get_rho_ideal_gas(p0, T0);
    
    // Finding the initial mass
    FLOAT_P C = 4 * M_PI * pow(R_SUN,3);
    FLOAT_P m0 = 0.0;
    FLOAT_P dr;
    // Backward integration of Solar S data (data is in reverse order)
    for (int i = solar_s_size-2; i > integration_start_index-1; i--)
    {
        dr = r_over_R_solar_s[i] - r_over_R_solar_s[i+1];
        m0 += C * 0.5 * (rho_solar_s[i]*pow(r_over_R_solar_s[i], 2) + rho_solar_s[i+1]*pow(r_over_R_solar_s[i+1], 2)) * dr;
    }

    // Small mass between integration_start and r_initial
    dr = r_initial/R_SUN - r_over_R_solar_s[integration_start_index];
    m0 += C * 0.5 * (rho0*pow(r_initial/R_SUN, 2) + rho_solar_s[integration_start_index]*pow(r_over_R_solar_s[integration_start_index], 2)) * dr;

    // Deallocating solar_s arrays
    deallocate_1D_array(r_over_R_solar_s);
    deallocate_1D_array(rho_solar_s);
    deallocate_1D_array(p_solar_s);
    deallocate_1D_array(T_solar_s);

    // Setting initial values
    *p0_initial = p0;
    *T0_initial = T0;
    *rho0_initial = rho0;
    *m_initial = m0;
}