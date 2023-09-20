#include "solve_elliptic_equation.h"

double rhs_elliptic_eq_horizontal_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j)
{
    double *rho0 = bg->rho0;
    double *g = bg->g;
    double **rho1 = fg->rho1;
    double **vx = fg->vx;
    int nx = fg->nx;

    double dx = fg->dx;
    double dz = fg->dz;

    int j_minus = periodic_boundary(j-1, nx);
    int j_plus = periodic_boundary(j+1, nx);

    #if CENTRAL_ORDER == 2
    double central_rho_g_dz = central_first_derivative_second_order(rho1[i-1][j]*g[i-1], rho1[i+1][j]*g[i+1], dz);
    double central_vx_dx = central_first_derivative_second_order(vx[i][j_minus], vx[i][j_plus], dx);
    double central_rho0_vx_dx = central_first_derivative_second_order(rho0[i]*vx[i][j_minus], rho0[i]*vx[i][j_plus], dx);
    double central_second_vx_dx = central_second_derivative_second_order(vx[i][j], vx[i][j_minus], vx[i][j_plus], dx);
    #else
    // Print error message
    printf("Error: CENTRAL_ORDER must be 2\n");
    #endif
    
    return - central_rho_g_dz - rho0[i]*vx[i][j]*central_second_vx_dx - central_vx_dx*central_rho0_vx_dx;
           

}