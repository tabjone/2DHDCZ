#include "solve_elliptic_equation.h"

FLOAT_P rhs_elliptic_eq(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j)
{
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;

    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **vx = fg->vx;
    FLOAT_P **vz = fg->vz;

    int nx = grid_info->nx;

    FLOAT_P dx = grid_info->dx;
    FLOAT_P dz = grid_info->dz;

    int j_minus = periodic_boundary(j-1, nx);
    int j_plus = periodic_boundary(j+1, nx);

    // I write this in to avoid bugs, but it is not necessary
    FLOAT_P vx_up_right = vx[i+1][j_plus];
    FLOAT_P vx_up_left = vx[i+1][j_minus];
    FLOAT_P vx_down_right = vx[i-1][j_plus];
    FLOAT_P vx_down_left = vx[i-1][j_minus];

    FLOAT_P vz_up_right = vz[i+1][j_plus];
    FLOAT_P vz_up_left = vz[i+1][j_minus];
    FLOAT_P vz_down_right = vz[i-1][j_plus];
    FLOAT_P vz_down_left = vz[i-1][j_minus];

    #if CENTRAL_ORDER == 2
    // First derivatives
    FLOAT_P drho1g_dz = central_first_derivative_second_order(rho1[i-1][j]*g[i-1], rho1[i+1][j]*g[i+1], dz);
    FLOAT_P drho0_dz = central_first_derivative_second_order(rho0[i-1], rho0[i+1], dz);

    FLOAT_P dvz_dx = central_first_derivative_second_order(vz[i][j_minus], vz[i][j_plus], dx);
    FLOAT_P dvx_dz = central_first_derivative_second_order(vx[i-1][j], vx[i+1][j], dz);

    // Mixed derivatives
    FLOAT_P dd_vx_vz_dxdz = central_mixed_derivative_second_order(vx_up_right*vz_up_right, vx_down_right*vz_down_right, vx_up_left*vz_up_left, vx_down_left*vz_down_left, dx, dz);

    // Second derivatives
    FLOAT_P dd_vz_ddz = central_second_derivative_second_order(vz[i][j], vz[i-1][j], vz[i+1][j], dz);
    #else
    // Print error message
    printf("Error: CENTRAL_ORDER must be 2\n");
    #endif

    return - drho1g_dz 
           - rho0[i]*(dd_vx_vz_dxdz+vz[i][j]*dd_vz_ddz+2*dvx_dz*dvz_dx)-vx[i][j]*dvx_dz*drho0_dz;
    

}