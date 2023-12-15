#include "solar_s_initialization.h"
#include <math.h>

void solar_s_background_initialization(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int grid_nz_full, int grid_nz_ghost, FLOAT_P grid_dz, FLOAT_P grid_z0, FLOAT_P grid_z1, int grid_nz)
{
    // Start by initializing the background radius array
    for (int l = 0; l < grid_nz_full; l++)
    {
        bg->r[l] = grid_z0 + (l-grid_nz_ghost) * grid_dz;
    }

    // Creating radius array for integration
    // If we need more points we pad around this array
    FLOAT_P *r_integration;
    int nz_ghost = grid_nz_ghost;
    int nz_full = NZ + 2*nz_ghost;


    allocate_1D_array(&r_integration, nz_full);
    
    for (int l = 0; l < nz_full; l++)
    {
        r_integration[l] = R_START*R_SUN + (l-nz_ghost) * grid_dz;
    }

    // Loading in the solar S data
    int solar_s_size = 2482;
    FLOAT_P *r_over_R_solar_s, *rho_solar_s, *p_solar_s, *T_solar_s;

    allocate_1D_array(&r_over_R_solar_s, solar_s_size);
    allocate_1D_array(&rho_solar_s, solar_s_size);
    allocate_1D_array(&p_solar_s, solar_s_size);
    allocate_1D_array(&T_solar_s, solar_s_size);

    // Reading solar_s data
    read_solar_s_data("additional_files/solar_s.h5", r_over_R_solar_s, rho_solar_s, p_solar_s, T_solar_s, solar_s_size);

    // Finding the index closest to integration start and interpolating to get the exact value
    int integration_start_index = 0;
    while (r_over_R_solar_s[integration_start_index] > R_END)
    {
        integration_start_index++;
    }
    FLOAT_P integration_start = R_END * R_SUN;

    FLOAT_P r_before = r_over_R_solar_s[integration_start_index]  * R_SUN;
    FLOAT_P r_after = r_over_R_solar_s[integration_start_index-1] * R_SUN;
    FLOAT_P rho0_before = rho_solar_s[integration_start_index];
    FLOAT_P rho0_after = rho_solar_s[integration_start_index-1];
    FLOAT_P p0_before = p_solar_s[integration_start_index];
    FLOAT_P p0_after = p_solar_s[integration_start_index-1];
    FLOAT_P T0_before = T_solar_s[integration_start_index];
    FLOAT_P T0_after = T_solar_s[integration_start_index-1];

    // Interpolating to get the exact value at integration_start
    FLOAT_P r_initial = integration_start;
    FLOAT_P rho0_initial = rho0_before + (rho0_after - rho0_before) / (r_after - r_before) * (integration_start - r_before);
    FLOAT_P p0_initial = p0_before + (p0_after - p0_before) / (r_after - r_before) * (integration_start - r_before);
    FLOAT_P T0_initial = T0_before + (T0_after - T0_before) / (r_after - r_before) * (integration_start - r_before);

    FLOAT_P C = 4 * M_PI * pow(R_SUN,3);
    FLOAT_P m_initial = 0.0;
    FLOAT_P dr;
    // Backward integration of Solar S data (data is in reverse order)
    for (int i = solar_s_size-2; i > integration_start_index-1; i--)
    {
        dr = r_over_R_solar_s[i] - r_over_R_solar_s[i+1];
        m_initial += C * 0.5 * (rho_solar_s[i]*pow(r_over_R_solar_s[i], 2) + rho_solar_s[i+1]*pow(r_over_R_solar_s[i+1], 2)) * dr;
    }

    // Small mass between integration_start and r_initial
    m_initial += C * 0.5 * (rho0_initial*pow(integration_start/R_SUN, 2) + rho_solar_s[integration_start_index]*pow(r_over_R_solar_s[integration_start_index], 2)) * (integration_start/R_SUN - r_over_R_solar_s[integration_start_index]);

    #if UNITS == 1
        // Convert from CGS to SI
        r_initial *= 1e-2; // NOPE, With current setup I need to do this for all values in the array for both bg->r and r_integration
        p0_initial *= 1e3;
        rho0_initial *= 1e-1;
        m_initial *= 1e-3;
    #endif

    // Deallocating solar_s arrays
    deallocate_1D_array(r_over_R_solar_s);
    deallocate_1D_array(rho_solar_s);
    deallocate_1D_array(p_solar_s);
    deallocate_1D_array(T_solar_s);

    // Create arrays for integration variables
    FLOAT_P *p0, *T0, *rho0, *grad_s0, *g, *m;
    allocate_1D_array(&p0, nz_full);
    allocate_1D_array(&T0, nz_full);
    allocate_1D_array(&rho0, nz_full);
    allocate_1D_array(&grad_s0, nz_full);
    allocate_1D_array(&g, nz_full);
    allocate_1D_array(&m, nz_full);

    FLOAT_P r_star = K_B / (MU * M_U);
    FLOAT_P c_p = r_star /(1.0-1.0/GAMMA);

    FLOAT_P dp_dr, dT_dr, ds_dr, dm_dr, k, nabla_star;

    // Setting initial values for integration arrays
    p0[nz_full-1] = p0_initial;
    T0[nz_full-1] = T0_initial;
    rho0[nz_full-1] = rho0_initial;

    rho0[nz_full-1] = M_U * MU / K_B * p0[nz_full-1]/T0[nz_full-1];
    m[nz_full-1] = m_initial;

    k = get_k_value(r_initial);
    nabla_star = NABLA_AD + k;
    dp_dr = - G * m_initial /pow(r_initial,2) * rho0_initial;

    grad_s0[nz_full-1] = c_p * (nabla_star - NABLA_AD) / p0_initial * dp_dr;
    g[nz_full-1] = G * m[nz_full-1] / pow(r_integration[nz_full-1],2);
    
    // Integrating downward
    dr = r_integration[nz_full-2] - r_integration[nz_full-1];

    for (int j = nz_full-2; j >= 0; j--)
    {
        k = get_k_value(r_integration[j+1]);
        nabla_star = NABLA_AD + k;

        dm_dr = 4 * M_PI * pow(r_integration[j+1],2) * rho0[j+1];
        dp_dr = - G * m[j+1] /pow(r_integration[j+1],2) * rho0[j+1];
        dT_dr = nabla_star * T0[j+1]/p0[j+1] * dp_dr;
        ds_dr = c_p * (nabla_star - NABLA_AD) / p0[j+1] * dp_dr;

        grad_s0[j] = ds_dr;
        m[j] = dm_dr * dr + m[j+1];
        p0[j] = dp_dr * dr + p0[j+1];
        T0[j] = dT_dr * dr + T0[j+1];
        rho0[j] = M_U * MU / K_B * p0[j]/T0[j]; //ideal gas law

        g[j] = G * m[j] / pow(r_integration[j],2); 
    }

    #if CONSTANT_BACKGROUND == 1
        // Constant background
        for (int i = 0; i < grid_nz_full; i++)
        {
            rho0[i] = CONSTANT_BACKGROUND_DENSITY;
            p0[i] = CONSTANT_BACKGROUND_PRESSURE;
            T0[i] = CONSTANT_BACKGROUND_TEMPERATURE;
            grad_s0[i] = 0.0;
            g[i] = 0.0;
        }
    #endif // CONSTANT_BACKGROUND

    // Setting the background variables
    int l = 0;
    while (bg->r[0] > r_integration[l])
    {
        l++;
    }

    for (int i = 0; i < grid_nz_full; i++)
    {
        bg->r[i] = r_integration[l+i];
        bg->p0[i] = p0[l+i];
        bg->T0[i] = T0[l+i];
        bg->rho0[i] = rho0[l+i];
        bg->grad_s0[i] = grad_s0[l+i];
        bg->g[i] = g[l+i];
    }

    // Deallocating integration arrays
    deallocate_1D_array(r_integration);
    deallocate_1D_array(p0);
    deallocate_1D_array(T0);
    deallocate_1D_array(rho0);
    deallocate_1D_array(grad_s0);
    deallocate_1D_array(g);
    deallocate_1D_array(m);
}