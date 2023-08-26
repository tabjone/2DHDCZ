#include "solar_s_initialization.h"
#include "../array_memory_management/array_memory_management.h"
#include "../..//global_parameters.h"
#include "../../global_constants.h"
#include "../../hd_solver/2D/structs/structs.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void solar_s_background_initialization(struct BackgroundVariables *bg)
{
    float p = 0.1; // number for determining the step size

    // Creating arrays for solar_s data
    int solar_s_size = 2482;
    double *r_over_R_solar_s, *c_s_solar_s, *rho_solar_s, *p_solar_s, *Gamma_1_solar_s, *T_solar_s;

    allocate_1D_array(&r_over_R_solar_s, solar_s_size);
    allocate_1D_array(&c_s_solar_s, solar_s_size);
    allocate_1D_array(&rho_solar_s, solar_s_size);
    allocate_1D_array(&p_solar_s, solar_s_size);
    allocate_1D_array(&Gamma_1_solar_s, solar_s_size);
    allocate_1D_array(&T_solar_s, solar_s_size);

    // Reading solar_s data
    read_solar_s_data("../additional_files/solar_s.h5", r_over_R_solar_s, c_s_solar_s, rho_solar_s, p_solar_s, Gamma_1_solar_s, T_solar_s, solar_s_size);

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
    deallocate_1D_array(c_s_solar_s);
    deallocate_1D_array(rho_solar_s);
    deallocate_1D_array(p_solar_s);
    deallocate_1D_array(Gamma_1_solar_s);
    deallocate_1D_array(T_solar_s);

    // Now starting integration
    
    int N = 100; // Number of integration steps, this will dynamically increase if needed

    // Creating arrays for upward integration from bottom of CZ
    double *r_up = malloc(sizeof(double)*N);
    double *p_up = malloc(sizeof(double)*N);
    double *m_up = malloc(sizeof(double)*N);
    double *T_up = malloc(sizeof(double)*N);
    double *rho_up = malloc(sizeof(double)*N);

    // Creating arrays for downward integration from bottom of CZ
    double *r_down = malloc(sizeof(double)*N);
    double *p_down = malloc(sizeof(double)*N);
    double *m_down = malloc(sizeof(double)*N);
    double *T_down = malloc(sizeof(double)*N);
    double *rho_down = malloc(sizeof(double)*N);

    // Setting the first element of the arrays to the values at the start of CZ
    r_up[0] = r_i;
    p_up[0] = p0_i;
    m_up[0] = m_i;
    T_up[0] = T0_i;
    rho_up[0] = rho0_i;

    r_down[0] = r_i;
    p_down[0] = p0_i;
    m_down[0] = m_i;
    T_down[0] = T0_i;
    rho_down[0] = rho0_i;

    // Integration variables
    double dr1, dr2, dr3;
    double dm_dr, dp_dr, dT_dr;
    double nabla_star;

    double k;

    // Upward integration
    int i = 0;
    while (r_up[i] < R_SUN * R_END)
    {   
        if (r_up[i]>=r_i)
        {
            k = 0.01;
        }
        else
        {
            k = 0.0;
        }

        nabla_star = NABLA_AD + k;

        dm_dr = 4 * M_PI * pow(r_up[i],2) *rho_up[i];
        dp_dr = - G * m_up[i] /pow(r_up[i],2) * rho_up[i];
        dT_dr = nabla_star * T_up[i]/p_up[i] * dp_dr;

        if (i % N == 0 && i != 0)
        {
            r_up = realloc(r_up, sizeof(double)*(i*N));
            p_up = realloc(p_up, sizeof(double)*(i*N));
            m_up = realloc(m_up, sizeof(double)*(i*N));
            T_up = realloc(T_up, sizeof(double)*(i*N));
            rho_up = realloc(rho_up, sizeof(double)*(i*N));
        }

        dr1 = fabs(p * m_up[i]/dm_dr);
        dr2 = fabs(p * p_up[i]/dp_dr);
        dr3 = fabs(p * T_up[i]/dT_dr);

        if (dr1 < dr2 && dr1 < dr3) {
            dr = dr1;
        } else if (dr2 < dr1 && dr2 < dr3) {
            dr = dr2;
        } else {
            dr = dr3;
        }

        r_up[i+1] = r_up[i] + dr;
        m_up[i+1] = m_up[i] + dm_dr * dr;
        p_up[i+1] = p_up[i] + dp_dr * dr;
        T_up[i+1] = T_up[i] + dT_dr * dr;
        rho_up[i+1] = M_U * MU / K_B * p_up[i+1]/T_up[i+1]; //ideal gas law
         
        i++;
    }

    // Downward integration
    N = 100;
    int j = 0;
    while (r_down[j] > R_SUN * R_START)
    {   
        if (r_down[j]>=r_i)
        {
            k = 0.01;
        }
        else
        {
            k = 0.0;
        }

        nabla_star = NABLA_AD + k;

        dm_dr = 4 * M_PI * pow(r_down[j],2) *rho_down[j];
        dp_dr = - G * m_down[j] /pow(r_down[j],2) * rho_down[j];
        dT_dr = nabla_star * T_down[j]/p_down[j] * dp_dr;

        if (j % N == 0 && j != 0)
        {
            r_down = realloc(r_down, sizeof(double)*(j*N));
            p_down = realloc(p_down, sizeof(double)*(j*N));
            m_down = realloc(m_down, sizeof(double)*(j*N));
            T_down = realloc(T_down, sizeof(double)*(j*N));
            rho_down = realloc(rho_down, sizeof(double)*(j*N));
        }

        dr1 = fabs(p * m_down[j]/dm_dr);
        dr2 = fabs(p * p_down[j]/dp_dr);
        dr3 = fabs(p * T_down[j]/dT_dr);

        if (dr1 < dr2 && dr1 < dr3) {
            dr = -dr1;
        } else if (dr2 < dr1 && dr2 < dr3) {
            dr = -dr2;
        } else {
            dr = -dr3;
        }

        r_down[j+1] = r_down[j] + dr;
        m_down[j+1] = m_down[j] + dm_dr * dr;
        p_down[j+1] = p_down[j] + dp_dr * dr;
        T_down[j+1] = T_down[j] + dT_dr * dr;
        rho_down[j+1] = M_U * MU / K_B * p_down[j+1]/T_down[j+1]; //ideal gas law
         
        j++;
    }

    allocate_background_struct(i+j, bg);
    
    for (int jj = j; jj > 0; jj--)
    {
        bg->r[j-jj] = r_down[jj];
        bg->p0[j-jj] = p_down[jj];
        bg->T0[j-jj] = T_down[jj];
        bg->rho0[j-jj] = rho_down[jj];
    }

    printf("i: %d\n", i);
    for (int ii = j; ii < i+j; ii++)
    {
        bg->r[ii] = r_up[ii-j];
        bg->p0[ii] = p_up[ii-j];
        bg->T0[ii] = T_up[ii-j];
        bg->rho0[ii] = rho_up[ii-j];
    }

    bg->nz = i+j;

    





    free(r_up);
    free(p_up);
    free(m_up);
    free(T_up);
    free(rho_up);

    free(r_down);
    free(p_down);
    free(m_down);
    free(T_down);
    free(rho_down);

    


    //while (r[i])


    //printf("sum: %f\n", m_i/M_SUN);

    //printf("%d\n" ,solar_s_size-cz_start_index);
    //dobule m_i = 4.0 * PI * rho0_i * pow(r_cz, 3.0) / 3.0;


    //printf("r_cz: %f\n", r_cz);

    //printf("cz_start_index: %d\n", cz_start_index);
    //Then we print the radius
    //printf("r_over_R_solar_s[cz_start_index]: %f\n", r_over_R_solar_s[cz_start_index]);

    //printf("Hello world!\n");

    //calculate_pressure_scale_height(r_over_R_solar_s, p_solar_s, H_solar_s, solar_s_size);
    //calculate_superadiabcicity_parameter(p_solar_s, T_solar_s, superad_param_solar_s, del_ad, solar_s_size);
    //calculate_entropy_gradient(p_solar_s, rho_solar_s, T_solar_s, Gamma_1_solar_s, H_solar_s, superad_param_solar_s, grad_s_solar_s, solar_s_size);
    //calculate_gravitational_acceleration(r_over_R_solar_s, rho_solar_s, g_solar_s, solar_s_size);



   //deallocate_1D_array(H_solar_s);
    //deallocate_1D_array(superad_param_solar_s);
    //deallocate_1D_array(grad_s_solar_s);
    //deallocate_1D_array(g_solar_s);
}