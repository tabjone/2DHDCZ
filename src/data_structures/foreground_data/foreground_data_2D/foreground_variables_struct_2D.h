#ifndef FOREGROUND_VARIABLES_STRUCT_2D_H_
#define FOREGROUND_VARIABLES_STRUCT_2D_H_

#include "global_float_precision.h"

struct ForegroundVariables2D
{
    FLOAT_P **p1;
    FLOAT_P **rho1;
    FLOAT_P **T1;
    FLOAT_P **s1;
    FLOAT_P **vy;
    FLOAT_P **vz;
};

#endif // FOREGROUND_VARIABLES_STRUCT_2D_H_