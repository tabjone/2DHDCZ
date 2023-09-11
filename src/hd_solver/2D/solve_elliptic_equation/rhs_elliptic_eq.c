#include "solve_elliptic_equation.h"

double rhs_elliptic_eq(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j)
{
    double *rho0 = bg->rho0;
    double *g = bg->g;

    double **rho1 = fg->rho1;
    double **vx = fg->vx;
    double **vz = fg->vz;

    int nx = fg->nx;

    int j_minus = periodic_boundary(j-1, nx);
    int j_plus = periodic_boundary(j+1, nx);

    #if CENTRAL_ORDER == 2
    // First derivatives
    double central_rho_g_dz = central_first_derivative_second_order(rho1[i-1][j]*g[i-1], rho1[i+1][j]*g[i+1], dz);

    double central_rho0_vx_dx = central_first_derivative_second_order(rho0[i]*vx[i][j_minus], rho0[i]*vx[i][j_plus], dx);
    double central_rho0_vx_dz = central_first_derivative_second_order(rho0[i-1]*vx[i-1][j], rho0[i+1]*vx[i+1][j], dz);

    double central_rho0_vz_dx = central_first_derivative_second_order(rho0[i]*vz[i][j_minus], rho0[i]*vz[i][j_plus], dx);
    double central_rho0_vz_dz = central_first_derivative_second_order(rho0[i-1]*vz[i-1][j], rho0[i+1]*vz[i+1][j], dz);

    double central_vx_dx = central_first_derivative_second_order(vx[i][j_minus], vx[i][j_plus], dx);
    double central_vz_dx = central_first_derivative_second_order(vz[i][j_minus], vz[i][j_plus], dx);

    double central_vx_dz = central_first_derivative_second_order(vx[i-1][j], vx[i+1][j], dz);
    double central_vz_dz = central_first_derivative_second_order(vz[i-1][j], vz[i+1][j], dz);

    // Mixed derivatives
    double central_vz_dxdz = central_mixed_derivative_second_order(vz[i+1][j_plus], vz[i-1][j_plus], vz[i+1][j_minus], vz[i-1][j_minus], dx, dy);

    double central_vx_dxdz = central_mixed_derivative_second_order(vx[i+1][j_plus], vx[i-1][j_plus], vx[i+1][j_minus], vx[i-1][j_minus], dx, dy);
    // Second derivatives
    double central_second_vx_dx = central_second_derivative_second_order(vx[i][j], vx[i][j_minus], vx[i][j_plus], dx);
    double central_second_vz_dz = central_second_derivative_second_order(vz[i][j], vz[i+1][j], vz[i-1][j], dz);
    #else
    // Print error message
    printf("Error: CENTRAL_ORDER must be 2\n");
    #endif

    return - central_rho_g_dz - rho0[i]*vx[i][j]*(central_second_vx_dx + central_vz_dxdz)
           - rho0[i]*vz[i][j]*(central_second_vz_dz + central_vx_dxdz)
           - central_vx_dx*central_rho0_vx_dx - central_vz_dx*central_rho0_vx_dz
           - central_vx_dz*central_rho0_vz_dx - central_vz_dz*central_rho0_vz_dz;

}