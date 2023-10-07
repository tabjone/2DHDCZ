#include "solve_elliptic_equation.h"

FLOAT_P rhs_elliptic_eq(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j)
{
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;

    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    int ny = grid_info->ny;

    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);

    // I write this in to avoid bugs, but it is not necessary
    FLOAT_P vy_up_right = vy[i+1][j_plus];
    FLOAT_P vy_up_left = vy[i+1][j_minus];
    FLOAT_P vy_down_right = vy[i-1][j_plus];
    FLOAT_P vy_down_left = vy[i-1][j_minus];

    FLOAT_P vz_up_right = vz[i+1][j_plus];
    FLOAT_P vz_up_left = vz[i+1][j_minus];
    FLOAT_P vz_down_right = vz[i-1][j_plus];
    FLOAT_P vz_down_left = vz[i-1][j_minus];

    #if CENTRAL_ORDER == 2
    // First derivatives
    FLOAT_P drho1g_dz = central_first_derivative_second_order(rho1[i-1][j]*g[i-1], rho1[i+1][j]*g[i+1], dz);
    FLOAT_P drho0_dz = central_first_derivative_second_order(rho0[i-1], rho0[i+1], dz);

    FLOAT_P dvz_dy = central_first_derivative_second_order(vz[i][j_minus], vz[i][j_plus], dy);
    FLOAT_P dvy_dz = central_first_derivative_second_order(vy[i-1][j], vy[i+1][j], dz);

    // Mixed derivatives
    FLOAT_P dd_vy_vz_dydz = central_mixed_derivative_second_order(vy_up_right*vz_up_right, vy_down_right*vz_down_right, vy_up_left*vz_up_left, vy_down_left*vz_down_left, dy, dz);

    // Second derivatives
    FLOAT_P dd_vz_ddz = central_second_derivative_second_order(vz[i][j], vz[i-1][j], vz[i+1][j], dz);
    #else
    // Print error message
    printf("Error: CENTRAL_ORDER must be 2\n");
    #endif

    #if GRAVITY_ON == 1 OR WHATEVER I CALLED IT

    return - drho1g_dz 
           - rho0[i]*(dd_vy_vz_dydz+vz[i][j]*dd_vz_ddz+2*dvy_dz*dvz_dy)-vy[i][j]*dvy_dz*drho0_dz;
}