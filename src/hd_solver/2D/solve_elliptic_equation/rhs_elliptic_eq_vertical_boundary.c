#include "rhs_functions.h"

double rhs_elliptic_eq_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j)
{
    double **vx = fg->vx;

    double **rho1 = fg->rho1;
    double *rho0 = bg->rho0;
    double *g = bg->g;

    int nx = fg->nx;

    int j_minus = periodic_boundary(j-1, nx);
    int j_plus = periodic_boundary(j+1, nx);

    // THIS IS NOW HARD CODED TO BE SECOND ORDER CENTRAL DIFFERENCES, CHANGE THIS!!!
    // THIS IS NOW HARD CODED TO BE SECOND ORDER CENTRAL DIFFERENCES, CHANGE THIS!!!

    double central_rho_g_dz = central_first_derivative_second_order(rho1[i-1][j]*g[i-1], rho1[i+1][j]*g[i+1], dz);
    double central_second_vx_dx = central_second_derivative_second_order(vx[i][j], vx[i][j_minus], vx[i][j_plus], dx);

    // THIS IS NOW HARD CODED TO BE SECOND ORDER CENTRAL DIFFERENCES, CHANGE THIS!!!
    // THIS IS NOW HARD CODED TO BE SECOND ORDER CENTRAL DIFFERENCES, CHANGE THIS!!!
 
    double central_vx_dx = central_first_derivative_second_order(vx[i][j_minus], vx[i][j_plus], dx);

    double central_rho0_vx_dx = central_first_derivative_second_order(rho0[i]*vx[i][j_minus], rho0[i]*vx[i][j_plus], dx);

    // THIS IS NOW HARD CODED TO BE SECOND ORDER CENTRAL DIFFERENCES, CHANGE THIS!!!

    return - central_rho_g_dz - rho0[i]*vx[i][j]*central_second_vx_dx - central_vx_dx*central_rho0_vx_dx;
           

}