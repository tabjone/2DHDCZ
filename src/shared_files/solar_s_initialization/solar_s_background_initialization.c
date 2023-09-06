#include "solar_s_initialization.h"
#include "../array_memory_management/array_memory_management.h"
#include "../..//global_parameters.h"
#include "../../global_constants.h"
#include "../structs/structs.h"
#include <math.h>
#include <stdlib.h>
#include "../interpolation/interpolation.h"

void solar_s_background_initialization(struct BackgroundVariables *bg)
{
    // Creating arrays for solar_s data
    int solar_s_size = 2482;
    double *r_over_R_solar_s, *rho_solar_s, *p_solar_s, *T_solar_s;

    allocate_1D_array(&r_over_R_solar_s, solar_s_size);
    allocate_1D_array(&rho_solar_s, solar_s_size);
    allocate_1D_array(&p_solar_s, solar_s_size);
    allocate_1D_array(&T_solar_s, solar_s_size);

    // Reading solar_s data
    read_solar_s_data("../additional_files/solar_s.h5", r_over_R_solar_s, rho_solar_s, p_solar_s, T_solar_s, solar_s_size);

    // Setting CZ start at 0.7 solar radii and finding the closest radius to that in the solar_s data
    int cz_start_index = 0;
    while (r_over_R_solar_s[cz_start_index] > 0.7)
    {
        cz_start_index++;
    }
    cz_start_index--;

    // Variables for start of integration
    double r_i = r_over_R_solar_s[cz_start_index] * R_SUN;
    double p0_i = p_solar_s[cz_start_index];
    double T0_i = T_solar_s[cz_start_index];
    double rho0_i = rho_solar_s[cz_start_index];
    
    // Finding the mass at the start of CZ by cumulative trapezoidal integration
    double m_i = 0.0; 
    double dr;
    double C = 4*M_PI*pow(R_SUN,3);
    for (int i = solar_s_size-2; i >= cz_start_index; i--) {
        dr =  r_over_R_solar_s[i]-r_over_R_solar_s[i+1];
        m_i += C * 0.5 * (rho_solar_s[i]*pow(r_over_R_solar_s[i], 2) + rho_solar_s[i+1]*pow(r_over_R_solar_s[i+1], 2)) * dr;
    }

    // Deallocating solar_s arrays
    deallocate_1D_array(r_over_R_solar_s);
    deallocate_1D_array(rho_solar_s);
    deallocate_1D_array(p_solar_s);
    deallocate_1D_array(T_solar_s);

    // Now starting integration
    
    int N = 100; // Number of integration steps, this will dynamically increase if needed

    // Creating arrays for upward integration from bottom of CZ
    double *r_up = malloc(sizeof(double)*N);
    double *p_up = malloc(sizeof(double)*N);
    double *m_up = malloc(sizeof(double)*N);
    double *T_up = malloc(sizeof(double)*N);
    double *rho_up = malloc(sizeof(double)*N);
    double *s_up = malloc(sizeof(double)*N);
    double *grad_s0_up = malloc(sizeof(double)*N);

    // Creating arrays for downward integration from bottom of CZ
    double *r_down = malloc(sizeof(double)*N);
    double *p_down = malloc(sizeof(double)*N);
    double *m_down = malloc(sizeof(double)*N);
    double *T_down = malloc(sizeof(double)*N);
    double *rho_down = malloc(sizeof(double)*N);
    double *s_down = malloc(sizeof(double)*N);
    double *grad_s0_down = malloc(sizeof(double)*N);

    // Setting the specific heat capacity at constant pressure
    double c_p = p0_i /((rho0_i * T0_i)*(1.0-1.0/GAMMA));

    // Setting the first element of the arrays to the values at the start of CZ
    r_up[0] = r_i;
    p_up[0] = p0_i;
    m_up[0] = m_i;
    T_up[0] = T0_i;
    rho_up[0] = rho0_i;
    s_up[0] = c_p;

    r_down[0] = r_i;
    p_down[0] = p0_i;
    m_down[0] = m_i;
    T_down[0] = T0_i;
    rho_down[0] = rho0_i;
    s_down[0] = c_p;

    // Integration variables
    double dr1, dr2, dr3, dr4;
    double dm_dr, dp_dr, dT_dr, ds_dr;
    double nabla_star;

    double k;

    // Upward integration
    int i = 0;
    while (r_up[i] < R_SUN * R_END)
    {   
        if (r_up[i]>=r_i)
        {
            k = K;
        }
        else
        {
            k = 0.0;
        }

        nabla_star = NABLA_AD + k;

        dm_dr = 4 * M_PI * pow(r_up[i],2) *rho_up[i];
        dp_dr = - G * m_up[i] /pow(r_up[i],2) * rho_up[i];
        dT_dr = nabla_star * T_up[i]/p_up[i] * dp_dr;
        ds_dr = (nabla_star-NABLA_AD) / (rho_up[i] * T_up[i]) * dp_dr;

        if (i % N == 0 && i != 0)
        {
            r_up = realloc(r_up, sizeof(double)*(i*N));
            p_up = realloc(p_up, sizeof(double)*(i*N));
            m_up = realloc(m_up, sizeof(double)*(i*N));
            T_up = realloc(T_up, sizeof(double)*(i*N));
            rho_up = realloc(rho_up, sizeof(double)*(i*N));
            s_up = realloc(s_up, sizeof(double)*(i*N));
            grad_s0_up = realloc(grad_s0_up, sizeof(double)*(i*N));
        }

        dr1 = fabs(p * m_up[i]/dm_dr);
        dr2 = fabs(p * p_up[i]/dp_dr);
        dr3 = fabs(p * T_up[i]/dT_dr);
        dr4 = fabs(p * s_up[i]/ds_dr);

        if (dr1 < dr2 && dr1 < dr3 && dr1 < dr4) {
            dr = dr1;
        } else if (dr2 < dr1 && dr2 < dr3 && dr2 < dr4) {
            dr = dr2;
        } else if (dr3 < dr1 && dr3 < dr2 && dr3 < dr4) {
            dr = dr3;
        } else {
            dr = dr4;
        }
        grad_s0_up[i] = ds_dr;

        r_up[i+1] = r_up[i] + dr;
        m_up[i+1] = m_up[i] + dm_dr * dr;
        p_up[i+1] = p_up[i] + dp_dr * dr;
        T_up[i+1] = T_up[i] + dT_dr * dr;
        s_up[i+1] = s_up[i] + ds_dr * dr;
        rho_up[i+1] = M_U * MU / K_B * p_up[i+1]/T_up[i+1]; //ideal gas law
         
        i++;
    }

    // Calculating last point of entropy gradient
    dp_dr = - G * m_up[i] /pow(r_up[i],2) * rho_up[i];
    grad_s0_up[i] = (nabla_star-NABLA_AD) / (rho_up[i] * T_up[i]) * dp_dr;

    // Downward integration
    N = 100;
    int j = 0;
    while (r_down[j] > R_SUN * R_START)
    {   
        if (r_down[j]>=r_i)
        {
            k = K;
        }
        else
        {
            k = 0.0;
        }

        nabla_star = NABLA_AD + k;

        dm_dr = 4 * M_PI * pow(r_down[j],2) *rho_down[j];
        dp_dr = - G * m_down[j] /pow(r_down[j],2) * rho_down[j];
        dT_dr = nabla_star * T_down[j]/p_down[j] * dp_dr;
        ds_dr = (nabla_star-NABLA_AD) / (rho_down[j] * T_down[j]) * dp_dr;

        if (j % N == 0 && j != 0)
        {
            r_down = realloc(r_down, sizeof(double)*(j*N));
            p_down = realloc(p_down, sizeof(double)*(j*N));
            m_down = realloc(m_down, sizeof(double)*(j*N));
            T_down = realloc(T_down, sizeof(double)*(j*N));
            rho_down = realloc(rho_down, sizeof(double)*(j*N));
            s_down = realloc(s_down, sizeof(double)*(j*N));
            grad_s0_down = realloc(grad_s0_down, sizeof(double)*(j*N));
        }

        dr1 = fabs(p * m_down[j]/dm_dr);
        dr2 = fabs(p * p_down[j]/dp_dr);
        dr3 = fabs(p * T_down[j]/dT_dr);
        dr4 = fabs(p * s_down[j]/ds_dr);

        if (dr1 < dr2 && dr1 < dr3 && dr1 < dr4) {
            dr = -dr1;
        } else if (dr2 < dr1 && dr2 < dr3 && dr2 < dr4) {
            dr = -dr2;
        } else if (dr3 < dr1 && dr3 < dr2 && dr3 < dr4) {
            dr = -dr3;
        } else {
            dr = -dr4;
        }
        grad_s0_down[j] = ds_dr;

        r_down[j+1] = r_down[j] + dr;
        m_down[j+1] = m_down[j] + dm_dr * dr;
        p_down[j+1] = p_down[j] + dp_dr * dr;
        T_down[j+1] = T_down[j] + dT_dr * dr;
        s_down[j+1] = s_down[j] + ds_dr * dr;
        rho_down[j+1] = M_U * MU / K_B * p_down[j+1]/T_down[j+1]; //ideal gas law
         
        j++;
    }

    // Calculating last point of entropy gradient
    dp_dr = - G * m_down[j] /pow(r_down[j],2) * rho_down[j];
    grad_s0_down[j] = (nabla_star-NABLA_AD) / (rho_down[j] * T_down[j]) * dp_dr;

    int nz_upwdown = i+j;

    // Combining the upward and downward integration arrays
    double *r = malloc(sizeof(double)*(nz_upwdown));
    double *p0 = malloc(sizeof(double)*(nz_upwdown));
    double *T0 = malloc(sizeof(double)*(nz_upwdown));
    double *rho0 = malloc(sizeof(double)*(nz_upwdown));
    double *grad_s0 = malloc(sizeof(double)*(nz_upwdown));
    double *g = malloc(sizeof(double)*(nz_upwdown));

    for (int jj = j; jj > 0; jj--)
    {
        r[j-jj] = r_down[jj];
        p0[j-jj] = p_down[jj];
        T0[j-jj] = T_down[jj];
        grad_s0[j-jj] = grad_s0_down[jj];
        rho0[j-jj] = rho_down[jj];
        g[j-jj] = -G * m_down[jj] / pow(r_down[jj],2);
    }

    for (int ii = j; ii < i+j; ii++)
    {
        r[ii] = r_up[ii-j];
        p0[ii] = p_up[ii-j];
        T0[ii] = T_up[ii-j];
        grad_s0[ii] = grad_s0_up[ii-j];
        rho0[ii] = rho_up[ii-j];
        g[ii] = -G * m_up[ii-j] / pow(r_up[ii-j],2);
    }

    // Deallocating integration arrays
    free(r_up);
    free(p_up);
    free(m_up);
    free(T_up);
    free(s_up);
    free(rho_up);
    free(grad_s0_up);

    free(r_down);
    free(p_down);
    free(m_down);
    free(T_down);
    free(s_down);
    free(rho_down);
    free(grad_s0_down);

    // Initialize background radius array
    for (i = 0; i < bg->nz; i++)
    {
        bg->r[i] = R_SUN * R_START + i * R_SUN * (R_END - R_START) / (bg->nz-1.0);
    }

    // Interpolating the background variables to the grid
    for (i = 0; i < bg->nz; i++)
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
        else if (bg->r[i] > r[nz_upwdown-1])
        {
            bg->p0[i] = p0[nz_upwdown-1];
            bg->T0[i] = T0[nz_upwdown-1];
            bg->rho0[i] = rho0[nz_upwdown-1];
            bg->grad_s0[i] = grad_s0[nz_upwdown-1];
            bg->g[i] = g[nz_upwdown-1];
        }
        else
        {
            int k = 0;
            while (bg->r[i] > r[k])
            {
                k++;
            }
            double x0 = r[k-1];
            double x1 = r[k];
            bg->p0[i] = interpolate_1D_linear(x0, x1, p0[k-1], p0[k], bg->r[i]);
            bg->T0[i] = interpolate_1D_linear(x0, x1, T0[k-1], T0[k], bg->r[i]);
            bg->rho0[i] = interpolate_1D_linear(x0, x1, rho0[k-1], rho0[k], bg->r[i]);
            bg->grad_s0[i] = interpolate_1D_linear(x0, x1, grad_s0[k-1], grad_s0[k], bg->r[i]);
            bg->g[i] = interpolate_1D_linear(x0, x1, g[k-1], g[k], bg->r[i]);
        }
    }

    // Deallocating integration arrays
    free(r);
    free(p0);
    free(T0);
    free(rho0);
    free(grad_s0);
    free(g);
}