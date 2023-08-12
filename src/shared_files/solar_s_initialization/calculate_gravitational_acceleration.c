#include <math.h>
#include <stdio.h>

void calculate_gravitational_acceleration(double *r_over_R, double *rho, double *g, int nz)
{
    /*
    Calculates the gravitational acceleration using Newtons formula for entire grid

    Parameters
    ----------
    r_over_R : double array
        Pointer to radius array normalized by the solar radius
    rho : double array
        Pointer to density array of background state
    g : double array
        Pointer to gravitational acceleration array of background state
    nz : int
        Number of grid points in z-direction
    */

    // Radius of the sun in cgs units
    double R_sun = 6.957e10;
    double M_sun = 1.989e33; // Solar mass in cgs units
    double G = 6.674e-8; // Gravitational constant in cgs units

    // Solar S goes from 0 to 1, so nz-1 is r=0, then there is no gravity
    g[nz-1] = 0.0;

    double integral = 0.0;
    double dz;
    int i;
    for (int j = 1; j < nz; j++)
    {
        i = nz - 1 - j; // We want to start from the bottom
        dz = r_over_R[i] - r_over_R[i+1];
        integral += dz * rho[i] * r_over_R[i]*r_over_R[i];
        g[i] = 4 * M_PI * G * R_sun/ pow(r_over_R[i], 2) * integral;
    }
}