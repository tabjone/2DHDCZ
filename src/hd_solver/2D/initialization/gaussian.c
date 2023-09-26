#include "initialization.h"

FLOAT_P gaussian(FLOAT_P x, FLOAT_P y, FLOAT_P x0, FLOAT_P y0, FLOAT_P sigma_x, FLOAT_P sigma_y, FLOAT_P A) 
{
    return A * exp(-(x - x0) * (x - x0) / (2 * sigma_x * sigma_x) 
                 -(y - y0) * (y - y0) / (2 * sigma_y * sigma_y));
}