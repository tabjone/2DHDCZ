#include "initialization.h"

double gaussian(double x, double y, double x0, double y0, double sigma_x, double sigma_y, double A) 
{
    return A * exp(-(x - x0) * (x - x0) / (2 * sigma_x * sigma_x) 
                 -(y - y0) * (y - y0) / (2 * sigma_y * sigma_y));
}