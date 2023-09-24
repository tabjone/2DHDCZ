#include "solve_diff_eqs.h"

double rhs_ds1_dt_horizontal_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j)
{
    /*
    int nx = fg->nx;
    double upwind_ds1_dx;

    double **vx = fg->vx;
    double **s1 = fg->s1;

    double dx = fg->dx;

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
    #elif UPWIND_ORDER == 2
    if (vx[i][j] >= 0)
    {
        upwind_ds1_dx = backward_first_derivative_second_order(s1[i][j], s1[i][j_minus], s1[i][j_minus2], dx);
    }
    else
    {
        upwind_ds1_dx = forward_first_derivative_second_order(s1[i][j], s1[i][j_plus], s1[i][j_plus2], dx);
    }
    #else
    // Print error message
    printf("Error: UPWIND_ORDER must be 1 or 2\n");
    #endif

    return - vx[i][j] * upwind_ds1_dx;*/
    printf("THIS SHOULD NOT BE USED!\n");
    return 1.0;
}