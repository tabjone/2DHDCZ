#include "solar_s_initialization.h"

void deallocate_integration_variables(struct IntegrationVariables *i_var)
{
    free(i_var->r);
    free(i_var->p0);
    free(i_var->T0);
    free(i_var->rho0);
    free(i_var->s0);
    free(i_var->grad_s0);
    free(i_var->m);
}