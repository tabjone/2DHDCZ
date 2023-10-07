#include "solve_elliptic_equation.h"

FLOAT_P rhs_elliptic_eq_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j)
{
    // THIS IS NOT USED!!!!
    // THIS IS THE VERTICAL BOUNDARY
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;
    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **vx = fg->vx;
    
    int nx = grid_info->nx;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dz = grid_info->dz;

    int j_minus = periodic_boundary(j-1, nx);
    int j_plus = periodic_boundary(j+1, nx);

    #if CENTRAL_ORDER == 2
    FLOAT_P central_rho_g_dz = central_first_derivative_second_order(rho1[i-1][j]*g[i-1], rho1[i+1][j]*g[i+1], dz);
    FLOAT_P central_vx_dx = central_first_derivative_second_order(vx[i][j_minus], vx[i][j_plus], dx);
    FLOAT_P central_rho0_vx_dx = central_first_derivative_second_order(rho0[i]*vx[i][j_minus], rho0[i]*vx[i][j_plus], dx);
    FLOAT_P central_second_vx_dx = central_second_derivative_second_order(vx[i][j], vx[i][j_minus], vx[i][j_plus], dx);
    #else
    // Print error message
    printf("Error: CENTRAL_ORDER must be 2\n");
    #endif

    printf("central_rho_g_dz=%f\n", central_rho_g_dz);
    
    return - central_rho_g_dz - rho0[i]*vx[i][j]*central_second_vx_dx - central_vx_dx*central_rho0_vx_dx;
           

}