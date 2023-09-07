#include "solar_s_initialization.h"

void set_initial_integration_values(struct IntegrationVariables *i_var, double r_i, double p0_i, double T0_i, double rho0_i, double m_i)
    {
    // Setting the specific heat capacity at constant pressure for ideal gas
    double c_p = p0_i /((rho0_i * T0_i)*(1.0-1.0/GAMMA));

    // Setting the initial values for the integration variables
    i_var->r[0] = r_i;
    i_var->p0[0] = p0_i;
    i_var->T0[0] = T0_i;
    i_var->rho0[0] = rho0_i;
    i_var->s0[0] = c_p;
    i_var->m[0] = m_i;
    }