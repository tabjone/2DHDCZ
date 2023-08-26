#include "../../../shared_files/derivatives/spacial_derivatives/spacial_derivative_functions.h"
#include "../boundaries/boundary_functions_2D_hd.h"
#include "../../../global_parameters.h"
#include "../structs/structs.h"

double rhs_dvx_dt_2D_hd(struct BackgroundVariables *background_variables, struct ForegroundVariables *foreground_variables, int i, int j, double dx, double dz, int nx)
{
    double **p1 = foreground_variables->p1;
    double **vx = foreground_variables->vx;
    double **vz = foreground_variables->vz;
    double *rho0 = background_variables->rho0;

    double dp1_dx, dvx_dx, dvx_dz;


    double first_term, second_term;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, nx);
    int j_minus2 = periodic_boundary(j-2, nx);
    int j_plus = periodic_boundary(j+1, nx);
    int j_plus2 = periodic_boundary(j+2, nx);

    #if UPWIND_ORDER == 1
    if (vx[i][j] >= 0)
    {
        dvx_dx = backward_first_derivative_first_order(vx[i][j], vx[i][j_minus], dx);
    }
    else
    {
        dvx_dx = forward_first_derivative_first_order(vx[i][j], vx[i][j_plus], dx);
    }
    if (vz[i][j] >= 0)
    {
        dvx_dz = backward_first_derivative_first_order(vx[i][j], vx[i-1][j], dz);
    }
    else
    {
        dvx_dz = forward_first_derivative_first_order(vx[i][j], vx[i+1][j], dz);
    }
    #elif UPWIND_ORDER == 2
    if (vx[i][j] >= 0.0)
    {
        dvx_dx = backward_first_derivative_second_order(vx[i][j], vx[i][j_minus], vx[i][j_minus2], dx);
    }
    else
    {
        dvx_dx = forward_first_derivative_second_order(vx[i][j], vx[i][j_plus], vx[i][j_plus2], dx);
    }
    if (vz[i][j] >= 0.0)
    {
        dvx_dz = backward_first_derivative_second_order(vx[i][j], vx[i-1][j], vx[i-2][j], dz);
    }
    else
    {
        dvx_dz = forward_first_derivative_second_order(vx[i][j], vx[i+1][j], vx[i+2][j], dz);
    }
    #endif

    #if CENTRAL_ORDER == 2
    dp1_dx = central_first_derivative_second_order(p1[i][j_minus], p1[i][j_plus], dx);
    #endif

    return -1.0/rho0[i] * dp1_dx - vx[i][j]*dvx_dx - vz[i][j]*dvx_dz;
}