#include "hdf5.h"
#include "../../../shared_files/shared_files.h"

double rhs_dvz_dt_2D_hd(double **rho1, double **p1, double **vx, double **vz, double *rho0, double *g, int i, int j, double dx, double dz)
{
    /*
    Calculates the right hand side of the equation dvz/dt = ... for the 2D HD solver.

    Parameters
    ----------
    rho1 : double**
        Pointer to 2D density pertubation array
    p1 : double**
        Pointer to 2D pressure pertubation array
    vx : double**
        Pointer to 2D x-velocity array
    vz : double**
        Pointer to 2D z-velocity array
    rho0 : double*
        Pointer to 1D background density array
    g : double*
        Pointer to 1D gravitational acceleration array
    i : int
        The index of the z-direction
    j : int
        The index of the x-direction
    dx : double
        The grid spacing in the x-direction
    dz : double
        The grid spacing in the z-direction

    Returns
    -------
    double
        The right hand side of the equation dvz/dt = ...
    */

    double first_term, second_term, third_term;

    // Force due to gas pressure potential, dp1/dz
    first_term = - 1/rho0[i] * central_first_derivative_second_order(p1[i-1][j], p1[i+1][j], dz);

    // Force due to advection (I think), vx*dvz/dx + vz*dvz/dz
    // Forward derivative for negative velocities, backward derivative for positive velocities
    second_term = 0.0;
    if (vx[i][j] >= 0.0)
    {
        second_term -= vx[i][j]*backward_first_derivative_first_order(vz[i][j], vz[i][j-1], dx);
    }
    else
    {
        second_term -= vx[i][j]*forward_first_derivative_first_order(vz[i][j], vz[i][j+1], dx);
    }

    if (vz[i][j] >= 0)
    {
        second_term -= vz[i][j]*backward_first_derivative_first_order(vz[i][j], vz[i-1][j], dz);
    }
    else
    {
        second_term -= vz[i][j]*forward_first_derivative_first_order(vz[i][j], vz[i+1][j], dz);
    }

    // Force due to gravity
    third_term = - rho1[i][j]/rho0[i] * g[i];

    return first_term + second_term + third_term;
}