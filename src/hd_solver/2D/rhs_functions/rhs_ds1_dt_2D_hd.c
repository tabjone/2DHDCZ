#include "hdf5.h"
#include "../../../shared_files/shared_files.h"
#include "../boundaries/boundary_functions_2D_hd.h"

double rhs_ds1_dt_2D_hd(double **s1, double *grad_s0, double **vx, double **vz, double *T0, double *rho0, int i, int j, double dx, double dz, int nx)
{
    /*
    Calculates the right hand side of the equation ds1/dt = ... for the 2D HD solver.

    Parameters
    ----------
    s1 : double**
        Pointer to 2D entropy pertubation array
    s0 : double*
        Pointer to 1D background entropy array
    vx : double**
        Pointer to 2D x-velocity array
    vz : double**
        Pointer to 2D z-velocity array
    T0 : double*
        Pointer to 1D background temperature array
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
        The right hand side of the equation ds1/dt = ...
    */

    double bracket_term;
    double central_ds1_dz, central_T0_dz;
    double upwind_ds1_dx, upwind_ds1_dz;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, nx);
    int j_plus = periodic_boundary(j+1, nx); 

    // Central derivatives in the bracket terms
    central_ds1_dz = central_first_derivative_second_order(s1[i-1][j], s1[i+1][j], dz);
    central_T0_dz = central_first_derivative_second_order(T0[i-1], T0[i+1], dz);

    // Upwind derivatives for the second term
    // Forward derivative for negative velocities, backward derivative for positive velocities
    if (vx[i][j] >= 0)
    {
        upwind_ds1_dx = backward_first_derivative_first_order(s1[i][j], s1[i][j_minus], dx);
    }
    else
    {
        upwind_ds1_dx = forward_first_derivative_first_order(s1[i][j], s1[i][j_plus], dx);
    }

    if (vz[i][j] >= 0)
    {
        upwind_ds1_dz = backward_first_derivative_first_order(s1[i][j], s1[i-1][j], dz);
    }
    else
    {
        upwind_ds1_dz = forward_first_derivative_first_order(s1[i][j], s1[i+1][j], dz);
    }

    bracket_term = T0[i]*central_second_derivative_second_order(s1[i][j], s1[i][j_minus], s1[i][j_plus], dx)
                 + T0[i]*central_second_derivative_second_order(s1[i][j], s1[i-1][j], s1[i+1][j], dz)
                 + central_ds1_dz * central_T0_dz;


    return 1.866e6/(rho0[i] * T0[i]) * bracket_term - vx[i][j] * upwind_ds1_dx - vz[i][j] * upwind_ds1_dz
           - vz[i][j] * grad_s0[i];

}