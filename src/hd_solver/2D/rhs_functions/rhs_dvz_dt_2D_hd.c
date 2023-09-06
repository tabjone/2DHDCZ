#include "../../../shared_files/derivatives/spacial_derivatives/spacial_derivative_functions.h"
#include "../boundaries/boundary_functions_2D_hd.h"
#include "../../../global_parameters.h"
#include "../../../shared_files/structs/structs.h"

double rhs_dvz_dt_2D_hd(struct BackgroundVariables *background_variables, struct ForegroundVariables2D *foreground_variables, int i, int j)
{
    int nx = foreground_variables->nx;
    
    double **rho1 = foreground_variables->rho1;
    double **p1 = foreground_variables->p1;
    double **vx = foreground_variables->vx;
    double **vz = foreground_variables->vz;
    double *rho0 = background_variables->rho0;
    double *g = background_variables->g;

    double dp1_dz, dvz_dx, dvz_dz;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, nx);
    int j_minus2 = periodic_boundary(j-2, nx);
    int j_plus = periodic_boundary(j+1, nx);
    int j_plus2 = periodic_boundary(j+2, nx);

    #if UPWIND_ORDER == 1
    if (vx[i][j] >= 0)
    {
        dvz_dx = backward_first_derivative_first_order(vz[i][j], vz[i][j_minus], dx);
    }
    else
    {
        dvz_dx = forward_first_derivative_first_order(vz[i][j], vz[i][j_plus], dx);
    }
    if (vz[i][j] >= 0)
    {
        dvz_dz = backward_first_derivative_first_order(vz[i][j], vz[i-1][j], dz);
    }
    else
    {
        dvz_dz = forward_first_derivative_first_order(vz[i][j], vz[i+1][j], dz);
    }
    #elif UPWIND_ORDER == 2
    if (vx[i][j] >= 0.0)
    {
        dvz_dx = backward_first_derivative_second_order(vz[i][j], vz[i][j_minus], vz[i][j_minus2], dx);
    }
    else
    {
        dvz_dx = forward_first_derivative_second_order(vz[i][j], vz[i][j_plus], vz[i][j_plus2], dx);
    }

    if (vz[i][j] >= 0)
    {
        dvz_dz = backward_first_derivative_second_order(vz[i][j], vz[i-1][j], vz[i-2][j], dz);
    }
    else
    {
        dvz_dz = forward_first_derivative_second_order(vz[i][j], vz[i+1][j], vz[i+2][j], dz);
    }

    #endif

    #if CENTRAL_ORDER == 2
    dp1_dz = central_first_derivative_second_order(p1[i-1][j], p1[i+1][j], dz);
    #endif

    return - 1.0/rho0[i] * dp1_dz - vx[i][j]*dvz_dx - vz[i][j]*dvz_dz - rho1[i][j]/rho0[i] * g[i];
}