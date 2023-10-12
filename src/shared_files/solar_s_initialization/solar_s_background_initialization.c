#include "solar_s_initialization.h"

void solar_s_background_initialization(struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    // Creating arrays for solar_s data and allocating memory
    int solar_s_size = 2482;
    FLOAT_P *r_over_R_solar_s, *rho_solar_s, *p_solar_s, *T_solar_s;

    allocate_1D_array(&r_over_R_solar_s, solar_s_size);
    allocate_1D_array(&rho_solar_s, solar_s_size);
    allocate_1D_array(&p_solar_s, solar_s_size);
    allocate_1D_array(&T_solar_s, solar_s_size);

    // Reading solar_s data
    read_solar_s_data("additional_files/solar_s.h5", r_over_R_solar_s, rho_solar_s, p_solar_s, T_solar_s, solar_s_size);

    // Setting CZ start and finding the closest radius to that in the solar_s data
    int cz_start_index = 0;
    while (r_over_R_solar_s[cz_start_index] > CZ_START)
    {
        cz_start_index++;
    }
    cz_start_index--;

    // Using CGS units for Solar S data, then converting to SI if UNITS == 1
    FLOAT_P R_SUN_CGS = 6.957e10; // solar radius in CGS

    // Variables for start of integration
    FLOAT_P r_i = r_over_R_solar_s[cz_start_index] * R_SUN_CGS;
    FLOAT_P p0_i = p_solar_s[cz_start_index];
    FLOAT_P T0_i = T_solar_s[cz_start_index];
    FLOAT_P rho0_i = rho_solar_s[cz_start_index];
    
    // Finding the mass at the start of CZ by cumulative trapezoidal integration
    FLOAT_P m_i = 0.0; 
    FLOAT_P dr;
    FLOAT_P C = 4*M_PI*pow(R_SUN_CGS,3);
    for (int i = solar_s_size-2; i >= cz_start_index; i--) {
        dr =  r_over_R_solar_s[i]-r_over_R_solar_s[i+1];
        m_i += C * 0.5 * (rho_solar_s[i]*pow(r_over_R_solar_s[i], 2) + rho_solar_s[i+1]*pow(r_over_R_solar_s[i+1], 2)) * dr;
    }

    #if UNITS == 1
        // Convert from CGS to SI
        r_i *= 1e-2;
        p0_i *= 1e3;
        rho0_i *= 1e-1;
        m_i *= 1e-3;
    #endif

    // Deallocating solar_s arrays
    deallocate_1D_array(r_over_R_solar_s);
    deallocate_1D_array(rho_solar_s);
    deallocate_1D_array(p_solar_s);
    deallocate_1D_array(T_solar_s);

    // Creating struct for holding upward integration variables
    struct IntegrationVariables i_var_up;
    i_var_up.N = 100; // Initial size of arrays
    i_var_up.N_increment = 2; // For dynamic array resizing

    // Allocating and setting initial values for integration arrays
    allocate_integration_variables(&i_var_up);
    set_initial_integration_values(&i_var_up, r_i, p0_i, T0_i, rho0_i, m_i);

    // Creating struct for holding downward integration variables
    struct IntegrationVariables i_var_down;
    i_var_down.N = 100;
    i_var_down.N_increment = 2;

    // Allocating and setting initial values for integration arrays
    allocate_integration_variables(&i_var_down);
    set_initial_integration_values(&i_var_down, r_i, p0_i, T0_i, rho0_i, m_i);

    // Integration variables for last point
    FLOAT_P dp_dr, nabla_star, k;

    // Upward integration
    int i = 0;
    while (i_var_up.r[i] < R_SUN * R_END)
    {   
        integrate_one_step(&i_var_up, i, false);
        i++;
    }

    // Calculating entropy gradient in last point
    k = get_k_value(i_var_up.r[i]);
    nabla_star = NABLA_AD + k;

    dp_dr = - G * i_var_up.m[i] /pow(i_var_up.r[i],2) * i_var_up.rho0[i];
    i_var_up.grad_s0[i] = (nabla_star-NABLA_AD) / (i_var_up.rho0[i] * i_var_up.T0[i]) * dp_dr;

    // Downward integration
    int j = 0;
    while (i_var_down.r[j] > R_SUN * R_START)
    {   
        integrate_one_step(&i_var_down, j, true); 
        j++;
    }

    // Calculating entropy gradient in last point
    k = get_k_value(i_var_down.r[j]);
    nabla_star = NABLA_AD + k;

    dp_dr = - G * i_var_down.m[j] /pow(i_var_down.r[j],2) * i_var_down.rho0[j];

    i_var_down.grad_s0[j] = (nabla_star-NABLA_AD) / (i_var_down.rho0[j] * i_var_down.T0[j]) * dp_dr;

    // Total number of points in upward and downward integration
    int nz = i+j;

    // Combining the upward and downward integration arrays
    FLOAT_P *r = malloc(sizeof(FLOAT_P)*(nz));
    FLOAT_P *p0 = malloc(sizeof(FLOAT_P)*(nz));
    FLOAT_P *T0 = malloc(sizeof(FLOAT_P)*(nz));
    FLOAT_P *rho0 = malloc(sizeof(FLOAT_P)*(nz));
    FLOAT_P *grad_s0 = malloc(sizeof(FLOAT_P)*(nz));
    FLOAT_P *g = malloc(sizeof(FLOAT_P)*(nz));

    for (int jj = j; jj > 0; jj--)
    {
        r[j-jj] = i_var_down.r[jj];
        p0[j-jj] = i_var_down.p0[jj];
        T0[j-jj] = i_var_down.T0[jj];
        grad_s0[j-jj] = i_var_down.grad_s0[jj];
        rho0[j-jj] = i_var_down.rho0[jj];
        g[j-jj] = G * i_var_down.m[jj] / pow(i_var_down.r[jj],2);
    }

    for (int ii = j; ii < nz; ii++)
    {
        r[ii] = i_var_up.r[ii-j];
        p0[ii] = i_var_up.p0[ii-j];
        T0[ii] = i_var_up.T0[ii-j];
        grad_s0[ii] = i_var_up.grad_s0[ii-j];
        rho0[ii] = i_var_up.rho0[ii-j];
        g[ii] = G * i_var_up.m[ii-j] / pow(i_var_up.r[ii-j],2);
    }

    // Deallocating integration arrays
    deallocate_integration_variables(&i_var_up);
    deallocate_integration_variables(&i_var_down);

    // Initialize background radius array
    for (i = 0; i < grid_info->nz_full-grid_info->nz_ghost; i++)
    {
        bg->r[i+grid_info->nz_ghost] = R_SUN * R_START + i * R_SUN * (R_END - R_START) / (grid_info->nz-1.0);
    }

    // Interpolating the background variables to the grid
    FLOAT_P x0, x1;
    for (i = grid_info->nz_ghost; i < grid_info->nz_full-grid_info->nz_ghost; i++)
    {
        // Handle edge cases
        if (bg->r[i] < r[0])
        {
            bg->p0[i] = p0[0];
            bg->T0[i] = T0[0];
            bg->rho0[i] = rho0[0];
            bg->grad_s0[i] = grad_s0[0];
            bg->g[i] = g[0];
        }
        else if (bg->r[i] > r[nz-1])
        {
            bg->p0[i] = p0[nz-1];
            bg->T0[i] = T0[nz-1];
            bg->rho0[i] = rho0[nz-1];
            bg->grad_s0[i] = grad_s0[nz-1];
            bg->g[i] = g[nz-1];
        }
        else
        {
            int k = 0;
            while (bg->r[i] > r[k])
            {
                k++;
            }
            x0 = r[k-1];
            x1 = r[k];
            bg->p0[i] = interpolate_1D_linear(x0, x1, p0[k-1], p0[k], bg->r[i]);
            bg->T0[i] = interpolate_1D_linear(x0, x1, T0[k-1], T0[k], bg->r[i]);
            bg->rho0[i] = interpolate_1D_linear(x0, x1, rho0[k-1], rho0[k], bg->r[i]);
            bg->grad_s0[i] = interpolate_1D_linear(x0, x1, grad_s0[k-1], grad_s0[k], bg->r[i]);
            bg->g[i] = interpolate_1D_linear(x0, x1, g[k-1], g[k], bg->r[i]);
        }
    }

    // Extrapolating background variables to ghost cells
    extrapolate_background(bg, grid_info);

    #if VERTICAL_BOUNDARY_TYPE == 2
        // Periodic vertical boundary with constant values for background
        for (int i = 0; i < grid_info->nz_full; i++)
        {
            bg->rho0[i] = 1.0e-1;
            bg->grad_s0[i] = 0.0;
            bg->T0[i] = 1.0e6;
            bg->g[i] = 0.0;
            bg->p0[i] = 1.0e13;
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    // Pre-calculate 1/rho0 and eta/(4*pi*rho0*T0)
    for (i = 0; i < grid_info->nz_full; i++)
    {
        bg->one_over_rho0[i] = 1.0/bg->rho0[i];
        bg->eta_over_four_pi_rho0_T0[i] = ETA/(4*M_PI*bg->rho0[i]*bg->T0[i]);

    }

    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    // Pre-calculate grad_g and grad_rho0
    // First for endpoints
    bg->grad_g[nz_ghost] = (bg->g[nz_ghost+1]-bg->g[nz_ghost])/(bg->r[nz_ghost+1]-bg->r[nz_ghost]);
    bg->grad_rho0[nz_ghost] = (bg->rho0[nz_ghost+1]-bg->rho0[nz_ghost])/(bg->r[nz_ghost+1]-bg->r[nz_ghost]);

    bg->grad_g[nz_full-nz_ghost-1] = (bg->g[nz_full-nz_ghost-1]-bg->g[nz_full-nz_ghost-2])/(bg->r[nz_full-nz_ghost-1]-bg->r[nz_full-nz_ghost-2]);
    bg->grad_rho0[nz_full-nz_ghost-1] = (bg->rho0[nz_full-nz_ghost-1]-bg->rho0[nz_full-nz_ghost-2])/(bg->r[nz_full-nz_ghost-1]-bg->r[nz_full-nz_ghost-2]);

    // Then for the rest of the grid using central derivative second order
    for (i = nz_ghost+1; i < nz_full-nz_ghost-1; i++)
    {
        bg->grad_g[i] = (bg->g[i+1]-bg->g[i-1])/(bg->r[i+1]-bg->r[i-1]);
        bg->grad_rho0[i] = (bg->rho0[i+1]-bg->rho0[i-1])/(bg->r[i+1]-bg->r[i-1]);
    }

    // Deallocating full integration arrays
    free(r);
    free(p0);
    free(T0);
    free(rho0);
    free(grad_s0);
    free(g);
}