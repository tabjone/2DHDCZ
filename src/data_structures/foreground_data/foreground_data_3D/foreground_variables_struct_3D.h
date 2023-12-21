#ifndef FOREGROUND_VARIABLES_STRUCT_3D_H
#define FOREGROUND_VARIABLES_STRUCT_3D_H

#include "global_float_precision.h"

struct ForegroundVariables3D
{
    FLOAT_P ***p1;
    FLOAT_P ***rho1;
    FLOAT_P ***T1;
    FLOAT_P ***s1;
    FLOAT_P ***vx;
    FLOAT_P ***vy;
    FLOAT_P ***vz;
};

#endif // FOREGROUND_VARIABLES_STRUCT_3D_H