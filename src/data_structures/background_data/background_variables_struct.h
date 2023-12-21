#ifndef BACKGROUND_VARIABLES_STRUCT_H_
#define BACKGROUND_VARIABLES_STRUCT_H_

#include "global_float_precision.h"

struct BackgroundVariables
{
    FLOAT_P *r;
    FLOAT_P *p0;
    FLOAT_P *T0;
    FLOAT_P *rho0;
    FLOAT_P *grad_s0;
    FLOAT_P *g;
};

#endif // BACKGROUND_VARIABLES_STRUCT_H_