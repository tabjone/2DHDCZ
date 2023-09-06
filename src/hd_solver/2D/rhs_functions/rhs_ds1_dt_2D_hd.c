#include "../../../shared_files/derivatives/spacial_derivatives/spacial_derivative_functions.h"
#include "../boundaries/boundary_functions_2D_hd.h"
#include "../../../global_parameters.h"
#include "../../../shared_files/structs/structs.h"

double rhs_ds1_dt_2D_hd(struct BackgroundVariables *background_variables, struct ForegroundVariables2D *foreground_variables, int i, int j)
{
    int nx = foreground_variables->nx;
    double central_ds1_dz, central_T0_dz;
    double upwind_ds1_dx, upwind_ds1_dz;

    double *grad_s0 = background_variables->grad_s0;

    double **vx = foreground_variables->vx;
    double **vz = foreground_variables->vz;
    double **s1 = foreground_variables->s1;

    // print the value of vx[0][0] from the struct

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, nx);
    int j_minus2 = periodic_boundary(j-2, nx);
    int j_plus = periodic_boundary(j+1, nx); 
    int j_plus2 = periodic_boundary(j+2, nx);

    #if UPWIND_ORDER == 1
    if (vx >= 0)
    {
        upwind_ds1_dx = backward_first_derivative_first_order(s1[i][j], s1[i][j_minus], dx);
    }
    else
    {
        upwind_ds1_dx = forward_first_derivative_first_order(s1[i][j], s1[i][j_plus], dx);
    }
    if (vz >= 0)
    {
        upwind_ds1_dz = backward_first_derivative_first_order(s1[i][j], s1[i-1][j], dz);
    }
    else
    {
        upwind_ds1_dz = forward_first_derivative_first_order(s1[i][j], s1[i+1][j], dz);
    }
    #elif UPWIND_ORDER == 2
    if (vx[i][j] >= 0)
    {
        upwind_ds1_dx = backward_first_derivative_second_order(s1[i][j], s1[i][j_minus], s1[i][j_minus2], dx);
    }
    else
    {
        upwind_ds1_dx = forward_first_derivative_second_order(s1[i][j], s1[i][j_plus], s1[i][j_plus2], dx);
    }

    if (vz[i][j] >= 0)
    {
        upwind_ds1_dz = backward_first_derivative_second_order(s1[i][j], s1[i-1][j], s1[i-2][j], dz);
    }
    else
    {
        upwind_ds1_dz = forward_first_derivative_second_order(s1[i][j], s1[i+1][j], s1[i+2][j], dz);
    }
    #endif

    return - vx[i][j] * upwind_ds1_dx 
           - vz[i][j] * upwind_ds1_dz
           - vz[i][j] * grad_s0[i];;
}