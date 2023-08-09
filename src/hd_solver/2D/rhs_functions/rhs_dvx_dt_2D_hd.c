#include "../../../shared_files/shared_functions.h"

double rhs_dvx_dt_2D_hd(double **rho1, double **p1, double **vx, double **vz, double *rho0, double *g, int i, int j, double dx, double dz)
{
    /*
    Calculates the right hand side of the equation dvx/dt = ... for the 2D HD solver.

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
        The right hand side of the equation dvx/dt = ...
    */
    
    double first_term, second_term;

    // Force due to gas pressure potential, dp1/dx
    first_term = - 1/rho0[i] * central_derivative_second_order(p1[i][j-1], p1[i][j+1], dx);
    
    // Force due to advection (I think), vx*dvx/dx + vz*dvx/dz
    // Forward derivative for negative velocities, backward derivative for positive velocities
    second_term = 0.0;
    if (vx[i][j] >= 0.0)
    {
        second_term -= vx[i][j]*backward_first_derivative_first_order(vx[i][j], vx[i][j-1], dx);
    }
    else
    {
        second_term -= vx[i][j]*forward_first_derivative_first_order(vx[i][j], vx[i][j+1], dx);
    }

    if (vz[i][j] >= 0)
    {
        second_term -= vz[i][j]*backward_first_derivative_first_order(vx[i][j], vx[i-1][j], dz);
    }
    else
    {
        second_term -= vz[i][j]*forward_first_derivative_first_order(vx[i][j], vx[i+1][j], dz);
    }

    return first_term + second_term;
}