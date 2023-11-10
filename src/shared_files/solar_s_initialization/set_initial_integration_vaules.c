#include "solar_s_initialization.h"

void set_initial_integration_values(struct IntegrationVariables *i_var, FLOAT_P r_i, FLOAT_P p0_i, FLOAT_P T0_i, FLOAT_P rho0_i, FLOAT_P m_i)
{
    // Setting the specific heat capacity at constant pressure for ideal gas
    FLOAT_P r_star = K_B / (MU * M_U);
    FLOAT_P c_p = r_star /(1.0-1.0/GAMMA);

    // Setting the initial values for the integration variables
    i_var->r[0] = r_i;
    i_var->p0[0] = p0_i;
    i_var->T0[0] = T0_i;
    i_var->rho0[0] = rho0_i;
    i_var->s0[0] = c_p;
    i_var->m[0] = m_i;

    i_var->rho0[0] = M_U * MU / K_B * p0_i/T0_i;
}