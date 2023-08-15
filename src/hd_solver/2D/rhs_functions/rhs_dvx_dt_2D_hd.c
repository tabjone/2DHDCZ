#include "hdf5.h"
#include "../../../shared_files/shared_files.h"
#include "../boundaries/boundary_functions_2D_hd.h"

double rhs_dvx_dt_2D_hd(double **p1, double **vx, double **vz, double *rho0, int i, int j, double dx, double dz, int nx)
{
    /*
    Calculates the right hand side of the equation dvx/dt = ... for the 2D HD solver.

    Parameters
    ----------
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

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, nx);
    int j_minus2 = periodic_boundary(j-2, nx);
    int j_plus = periodic_boundary(j+1, nx);
    int j_plus2 = periodic_boundary(j+2, nx);

    // Force due to gas pressure potential, dp1/dx
    first_term = - 1/rho0[i] * central_first_derivative_second_order(p1[i][j_minus], p1[i][j_plus], dx);
    
    // Force due to advection (I think), vx*dvx/dx + vz*dvx/dz
    // Forward derivative for negative velocities, backward derivative for positive velocities
    second_term = 0.0;
    if (vx[i][j] >= 0.0)
    {
        second_term -= vx[i][j]*backward_first_derivative_second_order(vx[i][j], vx[i][j_minus], vx[i][j_minus2], dx);
    }
    else
    {
        second_term -= vx[i][j]*forward_first_derivative_second_order(vx[i][j], vx[i][j_plus], vx[i][j_plus2], dx);
    }

    if (vz[i][j] >= 0)
    {
        second_term -= vz[i][j]*backward_first_derivative_second_order(vx[i][j], vx[i-1][j], vx[i-2][j], dz);
    }
    else
    {
        second_term -= vz[i][j]*forward_first_derivative_second_order(vx[i][j], vx[i+1][j], vx[i+2][j], dz);
    }

    return first_term + second_term;
}