#include "solar_s_initialization.h"

void allocate_integration_variables(struct IntegrationVariables *i_var)
{
    i_var->r = malloc(sizeof(FLOAT_P)*i_var->N);
    i_var->rho0 = malloc(sizeof(FLOAT_P)*i_var->N);
    i_var->p0 = malloc(sizeof(FLOAT_P)*i_var->N);
    i_var->m = malloc(sizeof(FLOAT_P)*i_var->N);
    i_var->T0 = malloc(sizeof(FLOAT_P)*i_var->N);
    i_var->s0 = malloc(sizeof(FLOAT_P)*i_var->N);
    i_var->grad_s0 = malloc(sizeof(FLOAT_P)*i_var->N);
}